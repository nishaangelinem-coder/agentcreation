/* monitor.c - real-time patient-monitor simulator.
 *
 * Threads (fixed priorities, highest first):
 *   T_ALARM  (80)  alarm manager: event-driven (blocking receive) or periodic polling variant
 *   T_WDG    (75)  software watchdog over the acquisition/processing chain, 20 ms period
 *   T_ACQ    (70)  500 Hz acquisition (simulated ADC: synthetic ECG + PPG), absolute-time periodic
 *   T_PROC   (60)  QRS detection, SpO2 estimation, alarm-condition evaluation (one job per sample)
 *   T_UI     (--)  display / operator, time-sharing class, never real-time
 *
 * Alarm latency  L = t_fire - t_rel(k*)  where k* is the sample whose processing first makes the
 * alarm condition true, t_rel its nominal (ideal) release instant and t_fire the instant the alarm
 * manager has driven the annunciation output.                                                   */
#define _GNU_SOURCE
#include "alarm.h"
#include "chan.h"
#include "common.h"
#include "ecg_proc.h"
#include "rt_util.h"
#include "siggen.h"
#include "spo2_proc.h"
#include <fcntl.h>
#include <getopt.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

enum { PRIO_ALARM = 80, PRIO_WDG = 75, PRIO_ACQ = 70, PRIO_PROC = 60 };
enum { EV_SET = 1, EV_CLEAR = 0, EV_ACK = 2 };
typedef struct { uint8_t id, kind; uint64_t k; int64_t t_rel; } alarm_ev_t;
typedef struct { uint8_t id, kind, out_prio; uint64_t k; int64_t t_rel, t_evt; } alarm_log_t;

static struct {
    int fifo, poll, poll_ms, cpu, ui, probe, mlock, lo_persist;
    double duration, stall_at, stall_len, deadline_ms; uint64_t seed; const char *out, *scenario;
} cfg = { .fifo = 1, .poll = 0, .poll_ms = 5, .cpu = 1, .ui = 0, .probe = 0, .mlock = -1, .lo_persist = 0,
          .duration = 30, .stall_at = -1, .stall_len = 0.4, .deadline_ms = 10, .seed = 1, .out = "results/run", .scenario = "clinical" };

static const segment_t SC_CLINICAL[] = {
    { 20, RH_NSR, 75, 97 }, { 10, RH_ASYSTOLE, 0, 97 }, { 20, RH_NSR, 75, 97 }, { 20, RH_NSR, 180, 97 },
    { 20, RH_NSR, 75, 97 }, { 25, RH_NSR, 32, 97 },     { 20, RH_NSR, 75, 97 }, { 25, RH_NSR, 75, 80 },
    { 20, RH_NSR, 75, 97 }, { 10, RH_LEADS_OFF, 75, 97 }, { 20, RH_NSR, 75, 97 } };
static const segment_t SC_NORMAL[] = { { 1e9, RH_NSR, 75, 97 } };

static chan_t ch_samp, ch_alarm;
static siggen_t gen; static ecg_proc_t ecg; static spo2_proc_t spo2; static alarm_mgr_t mgr;
static int64_t t0_ns; static uint64_t n_total;
static atomic_int running = 1;
static atomic_ullong hb_acq, hb_proc;
static atomic_int out_prio;                               /* annunciation output (audio/visual pattern) */
static _Atomic int64_t latched_since[AL_N];               /* >0: latched, condition ended at this time */
static atomic_uint active_mask;
static int fd_annun;

/* polling variant: condition table shared between T_PROC/T_WDG/T_UI and T_ALARM */
static struct { pthread_mutex_t mu; uint8_t pend_set[AL_N], cur[AL_N], dirty[AL_N], pend_ack[AL_N];
                uint64_t set_k[AL_N]; int64_t set_trel[AL_N]; } pt;

/* display snapshot, shared between T_PROC (RT) and T_UI (non-RT) under a PI mutex */
#define DISP_N 2000
static struct { pthread_mutex_t mu; float ecg[DISP_N]; int w; double hr, spo2; int hr_ok, spo2_ok; uint64_t k; } snap;

/* pre-allocated logs (no allocation or I/O in the real-time phase) */
static int32_t *lg_jit, *lg_pipe, *lg_cacq, *lg_cproc;
static alarm_log_t *lg_al; static atomic_int n_al; static int cap_al;
static int64_t wc_alarm_ns, sum_alarm_ns; static long n_alarm_jobs;
static int64_t wc_wdg_ns, wc_ui_cs_ns;                   /* observed execution-time maxima (for the response-time bound) */
static int64_t stall_start_ns;

static void post_alarm(alarm_id_t id, int kind, uint64_t k, int64_t t_rel) {
    if (!cfg.poll) { alarm_ev_t ev = { (uint8_t)id, (uint8_t)kind, k, t_rel }; chan_push(&ch_alarm, &ev); return; }
    pthread_mutex_lock(&pt.mu);
    if (kind == EV_ACK) pt.pend_ack[id] = 1;
    else {
        pt.cur[id] = (uint8_t)kind; pt.dirty[id] = 1;
        if (kind == EV_SET && !pt.pend_set[id]) { pt.pend_set[id] = 1; pt.set_k[id] = k; pt.set_trel[id] = t_rel; }
    }
    pthread_mutex_unlock(&pt.mu);
}

static void log_alarm(uint8_t id, uint8_t kind, uint64_t k, int64_t t_rel, int64_t t_evt) {
    int i = atomic_fetch_add(&n_al, 1);
    if (i < cap_al) lg_al[i] = (alarm_log_t){ id, kind, (uint8_t)mgr.out, k, t_rel, t_evt };
}

/* Executed by T_ALARM only. */
static void alarm_handle(uint8_t id, uint8_t kind, uint64_t k, int64_t t_rel) {
    if (kind == EV_ACK) { if (alarm_mgr_ack(&mgr, id)) { atomic_store(&latched_since[id], 0); log_alarm(id, 3, k, 0, now_ns()); } }
    else {
        int fired = alarm_mgr_event(&mgr, (alarm_id_t)id, kind == EV_SET, k, t_rel, 0);
        atomic_store(&out_prio, mgr.out);
        if (fired) {
            unsigned char pat = (unsigned char)mgr.out;
            if (write(fd_annun, &pat, 1) < 0) { }           /* drive the annunciator (simulated device) */
            int64_t t = now_ns(); mgr.s[id].t_fire = t;
            log_alarm(id, 1, k, t_rel, t);
        } else if (kind == EV_CLEAR) {
            int64_t t = now_ns(); log_alarm(id, 0, k, t_rel, t);
            atomic_store(&latched_since[id], mgr.s[id].latched ? t : 0);
        }
        if (kind == EV_SET) atomic_store(&latched_since[id], 0);
    }
    unsigned m = 0; for (int i = 0; i < AL_N; i++) if (mgr.s[i].active) m |= 1u << i;
    atomic_store(&active_mask, m); atomic_store(&out_prio, mgr.out);
}

static void *th_alarm(void *arg) {
    (void)arg;
    if (!cfg.poll) {
        alarm_ev_t ev;
        while (chan_pop(&ch_alarm, &ev) == 0) {
            int64_t c0 = thread_cpu_ns();
            alarm_handle(ev.id, ev.kind, ev.k, ev.t_rel);
            int64_t c = thread_cpu_ns() - c0; if (c > wc_alarm_ns) wc_alarm_ns = c; sum_alarm_ns += c; n_alarm_jobs++;
        }
    } else {
        int64_t next = t0_ns, per = (int64_t)cfg.poll_ms * 1000000LL;
        while (atomic_load(&running)) {
            sleep_until_ns(next); next += per;
            int64_t c0 = thread_cpu_ns();
            uint8_t ps[AL_N], cur[AL_N], dirty[AL_N], ack[AL_N]; uint64_t sk[AL_N]; int64_t st[AL_N];
            pthread_mutex_lock(&pt.mu);
            memcpy(ps, pt.pend_set, AL_N); memcpy(cur, pt.cur, AL_N); memcpy(dirty, pt.dirty, AL_N); memcpy(ack, pt.pend_ack, AL_N);
            memcpy(sk, pt.set_k, sizeof sk); memcpy(st, pt.set_trel, sizeof st);
            memset(pt.pend_set, 0, AL_N); memset(pt.dirty, 0, AL_N); memset(pt.pend_ack, 0, AL_N);
            pthread_mutex_unlock(&pt.mu);
            for (int i = 0; i < AL_N; i++) {
                if (ps[i]) alarm_handle((uint8_t)i, EV_SET, sk[i], st[i]);
                if (dirty[i] && !cur[i]) alarm_handle((uint8_t)i, EV_CLEAR, 0, 0);
                if (ack[i]) alarm_handle((uint8_t)i, EV_ACK, 0, 0);
            }
            int64_t c = thread_cpu_ns() - c0; if (c > wc_alarm_ns) wc_alarm_ns = c; sum_alarm_ns += c; n_alarm_jobs++;
        }
    }
    return NULL;
}

static void *th_acq(void *arg) {
    (void)arg;
    int64_t next = t0_ns;
    for (uint64_t k = 0; k < n_total; k++, next += PERIOD_NS) {
        sleep_until_ns(next);
        int64_t t = now_ns(), c0 = thread_cpu_ns();
        sample_t s; siggen_step(&gen, k, &s);               /* stands in for the ADC/SPI read */
        s.t_rel_ns = next; s.t_acq_ns = t;
        if (chan_push(&ch_samp, &s) != 0) post_alarm(AL_DATA_LOSS, EV_SET, k, next);
        atomic_fetch_add(&hb_acq, 1);
        lg_jit[k] = (int32_t)(t - next); lg_cacq[k] = (int32_t)(thread_cpu_ns() - c0);
    }
    chan_close(&ch_samp);
    return NULL;
}

static void *th_proc(void *arg) {
    (void)arg;
    sample_t s; uint8_t prev[AL_N] = { 0 };
    int stalled = 0; uint64_t stall_k = cfg.stall_at >= 0 ? (uint64_t)(cfg.stall_at * FS_ECG) : UINT64_MAX;
    while (chan_pop(&ch_samp, &s) == 0) {
        int64_t c0 = thread_cpu_ns();
        if (s.k >= stall_k && !stalled) {                    /* fault injection: runaway loop in T_PROC */
            stalled = 1; stall_start_ns = now_ns();
            int64_t until = stall_start_ns + (int64_t)(cfg.stall_len * 1e9);
            while (now_ns() < until) { }
            c0 = thread_cpu_ns();
        }
        ecg_proc_step(&ecg, s.ecg_mv);
        if (s.has_ppg) spo2_proc_step(&spo2, s.red, s.ir);
        uint8_t c[AL_N] = { 0 };
        c[AL_ASYSTOLE] = ecg.cond_asys; c[AL_TACHY] = ecg.cond_tachy; c[AL_BRADY] = ecg.cond_brady;
        c[AL_DESAT] = spo2.cond_desat; c[AL_SPO2_LOW] = spo2.cond_low; c[AL_LEADS_OFF] = ecg.leads_off; c[AL_LOW_PERF] = spo2.cond_perf;
        for (int i = 0; i <= AL_LOW_PERF; i++)
            if (c[i] != prev[i]) { post_alarm((alarm_id_t)i, c[i] ? EV_SET : EV_CLEAR, s.k, s.t_rel_ns); prev[i] = c[i]; }
        pthread_mutex_lock(&snap.mu);
        snap.ecg[snap.w] = s.ecg_mv; snap.w = (snap.w + 1) % DISP_N; snap.k = s.k;
        snap.hr = ecg.hr; snap.hr_ok = ecg.hr_valid; snap.spo2 = spo2.spo2; snap.spo2_ok = spo2.valid;
        pthread_mutex_unlock(&snap.mu);
        atomic_fetch_add(&hb_proc, 1);
        lg_cproc[s.k] = (int32_t)(thread_cpu_ns() - c0);
        int64_t d = now_ns() - s.t_rel_ns; lg_pipe[s.k] = d > INT32_MAX ? INT32_MAX : (int32_t)d;
    }
    return NULL;
}

static void *th_wdg(void *arg) {
    (void)arg;
    int64_t next = t0_ns + 100000000LL; unsigned long long la = 0, lp = 0; int miss = 0, raised = 0;
    while (atomic_load(&running)) {
        sleep_until_ns(next); next += 20000000LL;
        int64_t c0 = thread_cpu_ns();
        unsigned long long a = atomic_load(&hb_acq), p = atomic_load(&hb_proc);
        if (a == n_total) break;
        miss = (a == la || p == lp) ? miss + 1 : 0; la = a; lp = p;
        if (miss >= 5 && !raised) { raised = 1; post_alarm(AL_SYS_STALL, EV_SET, a, stall_start_ns); }
        if (miss == 0 && raised)  { raised = 0; post_alarm(AL_SYS_STALL, EV_CLEAR, a, 0); }
        int64_t c = thread_cpu_ns() - c0; if (c > wc_wdg_ns) wc_wdg_ns = c;
    }
    return NULL;
}

static void render(char *buf, size_t cap, const float *e, int w, double hr, int hr_ok, double sp, int sp_ok, uint64_t k) {
    enum { COLS = 100, ROWS = 9 }; char grid[ROWS][COLS + 1];
    for (int r = 0; r < ROWS; r++) { memset(grid[r], ' ', COLS); grid[r][COLS] = 0; }
    for (int c = 0; c < COLS; c++) {
        float lo = 1e9f, hi = -1e9f;
        for (int j = 0; j < DISP_N / COLS; j++) { float v = e[(w + c * (DISP_N / COLS) + j) % DISP_N]; if (v < lo) lo = v; if (v > hi) hi = v; }
        int r0 = (int)((1.4f - hi) / 2.0f * (ROWS - 1) + 0.5f), r1 = (int)((1.4f - lo) / 2.0f * (ROWS - 1) + 0.5f);
        if (r0 < 0) r0 = 0; if (r1 > ROWS - 1) r1 = ROWS - 1; if (r0 > ROWS - 1) r0 = ROWS - 1; if (r1 < 0) r1 = 0;
        for (int r = r0; r <= r1; r++) grid[r][c] = (r1 - r0 >= 2) ? '|' : '-';
    }
    int p = atomic_load(&out_prio); unsigned m = atomic_load(&active_mask);
    static const char *col[] = { "\x1b[42;30m", "\x1b[46;30m", "\x1b[43;30m", "\x1b[41;97;5m" };
    static const char *pn[] = { " NO ALARM ", " LOW PRIORITY ", " MEDIUM PRIORITY ", " *** HIGH PRIORITY *** " };
    int n = snprintf(buf, cap, "\x1b[H\x1b[2J RT PATIENT MONITOR SIMULATOR   t=%7.1f s   %s%s\x1b[0m\n\x1b[32m", (double)k / FS_ECG, col[p], pn[p]);
    for (int r = 0; r < ROWS; r++) n += snprintf(buf + n, cap - (size_t)n, " %s\n", grid[r]);
    n += snprintf(buf + n, cap - (size_t)n, "\x1b[0m\n  \x1b[1;32mHR  ");
    n += hr_ok ? snprintf(buf + n, cap - (size_t)n, "%3.0f bpm", hr) : snprintf(buf + n, cap - (size_t)n, "--- bpm");
    n += snprintf(buf + n, cap - (size_t)n, "\x1b[0m      \x1b[1;36mSpO2  ");
    n += sp_ok ? snprintf(buf + n, cap - (size_t)n, "%3.0f %%", sp) : snprintf(buf + n, cap - (size_t)n, "--- %%");
    n += snprintf(buf + n, cap - (size_t)n, "\x1b[0m\n\n  Active alarms: ");
    for (int i = 0; i < AL_N; i++) if (m & (1u << i)) n += snprintf(buf + n, cap - (size_t)n, "%s%s\x1b[0m ", col[alarm_prio[i]], alarm_name[i]);
    snprintf(buf + n, cap - (size_t)n, "\n");
}

static void *th_ui(void *arg) {
    (void)arg;
    static float e[DISP_N]; static char buf[8192];
    while (atomic_load(&running)) {
        usleep(100000);
        pthread_mutex_lock(&snap.mu);
        int64_t cs0 = thread_cpu_ns();
        memcpy(e, snap.ecg, sizeof e); int w = snap.w; double hr = snap.hr, sp = snap.spo2; int ho = snap.hr_ok, so = snap.spo2_ok; uint64_t k = snap.k;
        int64_t cs = thread_cpu_ns() - cs0; if (cs > wc_ui_cs_ns) wc_ui_cs_ns = cs;
        pthread_mutex_unlock(&snap.mu);
        render(buf, sizeof buf, e, w, hr, ho, sp, so, k);
        if (cfg.ui) { fputs(buf, stdout); fflush(stdout); }
        /* scripted operator: acknowledge a latched alarm 3 s after its condition has ended */
        int64_t t = now_ns();
        for (int i = 0; i < AL_N; i++) { int64_t ls = atomic_load(&latched_since[i]); if (ls > 0 && t - ls > 3000000000LL) post_alarm((alarm_id_t)i, EV_ACK, 0, 0); }
    }
    return NULL;
}

static void usage(void) {
    fprintf(stderr, "monitor [--policy fifo|other] [--alarm event|poll] [--poll-ms N] [--scenario clinical|probe|normal]\n"
                    "        [--duration S] [--seed N] [--cpu N|-1] [--mlock 0|1] [--inject-stall T] [--deadline-ms D] [--ui] [--out PREFIX]\n");
    exit(2);
}

int main(int argc, char **argv) {
    static const struct option lo[] = { { "policy", 1, 0, 'p' }, { "alarm", 1, 0, 'a' }, { "poll-ms", 1, 0, 'P' }, { "scenario", 1, 0, 's' },
        { "duration", 1, 0, 'd' }, { "seed", 1, 0, 'S' }, { "cpu", 1, 0, 'c' }, { "mlock", 1, 0, 'm' }, { "inject-stall", 1, 0, 'i' },
        { "deadline-ms", 1, 0, 'D' }, { "ui", 0, 0, 'u' }, { "out", 1, 0, 'o' }, { 0, 0, 0, 0 } };
    int o;
    while ((o = getopt_long(argc, argv, "", lo, NULL)) != -1) switch (o) {
        case 'p': cfg.fifo = !strcmp(optarg, "fifo"); break;
        case 'a': cfg.poll = !strcmp(optarg, "poll"); break;
        case 'P': cfg.poll_ms = atoi(optarg); break;
        case 's': cfg.scenario = optarg; break;
        case 'd': cfg.duration = atof(optarg); break;
        case 'S': cfg.seed = strtoull(optarg, NULL, 10); break;
        case 'c': cfg.cpu = atoi(optarg); break;
        case 'm': cfg.mlock = atoi(optarg); break;
        case 'i': cfg.stall_at = atof(optarg); break;
        case 'D': cfg.deadline_ms = atof(optarg); break;
        case 'u': cfg.ui = 1; break;
        case 'o': cfg.out = optarg; break;
        default: usage();
    }
    if (cfg.mlock < 0) cfg.mlock = cfg.fifo;
    cfg.probe = !strcmp(cfg.scenario, "probe");
    const segment_t *sc = SC_NORMAL; int nsc = 1;
    if (!strcmp(cfg.scenario, "clinical")) {
        sc = SC_CLINICAL; nsc = (int)(sizeof SC_CLINICAL / sizeof SC_CLINICAL[0]);
        double tot = 0; for (int i = 0; i < nsc; i++) tot += sc[i].dur_s; cfg.duration = tot;
    }
    n_total = (uint64_t)(cfg.duration * FS_ECG);

    siggen_init(&gen, sc, nsc, cfg.seed, 0.02, 1.0);
    if (cfg.probe) { siggen_enable_probe(&gen, 4, 8, 24); cfg.lo_persist = 3; }
    ecg_proc_init(&ecg, cfg.lo_persist); spo2_proc_init(&spo2); alarm_mgr_init(&mgr);
    chan_init(&ch_samp, sizeof(sample_t), 512, cfg.fifo); chan_init(&ch_alarm, sizeof(alarm_ev_t), 256, cfg.fifo);
    rt_mutex_init(&pt.mu, cfg.fifo); rt_mutex_init(&snap.mu, cfg.fifo);
    lg_jit = calloc(n_total, 4); lg_pipe = calloc(n_total, 4); lg_cacq = calloc(n_total, 4); lg_cproc = calloc(n_total, 4);
    cap_al = (int)(n_total / 4) + 1024; lg_al = calloc((size_t)cap_al, sizeof *lg_al);
    fd_annun = open("/dev/null", O_WRONLY);
    if (!lg_jit || !lg_pipe || !lg_cacq || !lg_cproc || !lg_al || fd_annun < 0) { perror("init"); return 1; }
    if (cfg.mlock && rt_lock_memory() != 0) fprintf(stderr, "warning: memory not locked\n");

    t0_ns = now_ns() + 200000000LL;
    pthread_t ta, tw, tq, tp, tu;
    if (rt_thread_create(&ta, cfg.fifo, PRIO_ALARM, cfg.cpu, th_alarm, NULL)) return 1;
    if (rt_thread_create(&tw, cfg.fifo, PRIO_WDG, cfg.cpu, th_wdg, NULL)) return 1;
    if (rt_thread_create(&tp, cfg.fifo, PRIO_PROC, cfg.cpu, th_proc, NULL)) return 1;
    if (rt_thread_create(&tq, cfg.fifo, PRIO_ACQ, cfg.cpu, th_acq, NULL)) return 1;
    if (rt_thread_create(&tu, 0, 0, cfg.cpu, th_ui, NULL)) return 1;

    pthread_join(tq, NULL); pthread_join(tp, NULL);
    usleep(50000);
    atomic_store(&running, 0); chan_close(&ch_alarm);
    pthread_join(ta, NULL); pthread_join(tw, NULL); pthread_join(tu, NULL);
    if (cfg.ui) printf("\x1b[0m\n");

    /* ---- non-real-time phase: dump logs ---- */
    char fn[512]; FILE *f;
    snprintf(fn, sizeof fn, "%s_samples.csv", cfg.out); f = fopen(fn, "w"); if (!f) { perror(fn); return 1; }
    fprintf(f, "k,jitter_ns,pipeline_ns,c_acq_ns,c_proc_ns\n");
    for (uint64_t k = 0; k < n_total; k++) fprintf(f, "%llu,%d,%d,%d,%d\n", (unsigned long long)k, lg_jit[k], lg_pipe[k], lg_cacq[k], lg_cproc[k]);
    fclose(f);
    int na = atomic_load(&n_al); if (na > cap_al) na = cap_al;
    snprintf(fn, sizeof fn, "%s_alarms.csv", cfg.out); f = fopen(fn, "w"); if (!f) { perror(fn); return 1; }
    fprintf(f, "alarm,event,k,t_s,latency_us,out_prio\n");
    static const char *kn[] = { "clear", "fire", "?", "ack" };
    long nfire = 0, nmiss = 0; double lmax = 0;
    for (int i = 0; i < na; i++) {
        alarm_log_t *a = &lg_al[i];
        double lat = (a->kind == 1 && a->t_rel) ? (double)(a->t_evt - a->t_rel) / 1e3 : -1;
        if (a->kind == 1 && a->id != AL_SYS_STALL) { nfire++; if (lat > lmax) lmax = lat; if (lat > cfg.deadline_ms * 1e3) nmiss++; }
        fprintf(f, "%s,%s,%llu,%.6f,%.1f,%d\n", alarm_name[a->id], kn[a->kind], (unsigned long long)a->k, (double)(a->t_evt - t0_ns) / 1e9, lat, a->out_prio);
    }
    fclose(f);
    snprintf(fn, sizeof fn, "%s_truth.csv", cfg.out); f = fopen(fn, "w"); if (!f) { perror(fn); return 1; }
    fprintf(f, "segment,rhythm,hr,spo2,k_start,k_last_r\n");
    for (int i = 0; i < nsc && i < 64; i++) fprintf(f, "%d,%d,%.0f,%.0f,%llu,%llu\n", i, sc[i].rhythm, sc[i].hr, sc[i].spo2,
                                                   (unsigned long long)gen.seg_start_k[i], (unsigned long long)gen.seg_last_r_k[i]);
    fclose(f);
    snprintf(fn, sizeof fn, "%s_summary.json", cfg.out); f = fopen(fn, "w"); if (!f) { perror(fn); return 1; }
    fprintf(f, "{\"policy\":\"%s\",\"alarm_path\":\"%s\",\"poll_ms\":%d,\"scenario\":\"%s\",\"duration_s\":%.1f,\"seed\":%llu,\"cpu\":%d,"
               "\"mlock\":%d,\"deadline_ms\":%.1f,\"samples\":%llu,\"sample_drops\":%llu,\"sample_q_high_water\":%zu,\"alarm_ev_drops\":%llu,"
               "\"alarm_fires\":%ld,\"deadline_misses\":%ld,\"latency_max_us\":%.1f,\"wcet_alarm_ns\":%lld,\"mean_alarm_ns\":%.0f,\"stall_at_s\":%.1f,"
               "\"wcet_wdg_ns\":%lld,\"wcet_ui_cs_ns\":%lld}\n",
            cfg.fifo ? "fifo" : "other", cfg.poll ? "poll" : "event", cfg.poll_ms, cfg.scenario, cfg.duration, (unsigned long long)cfg.seed, cfg.cpu,
            cfg.mlock, cfg.deadline_ms, (unsigned long long)n_total, (unsigned long long)ch_samp.drops, ch_samp.high_water,
            (unsigned long long)ch_alarm.drops, nfire, nmiss, lmax, (long long)wc_alarm_ns, n_alarm_jobs ? (double)sum_alarm_ns / n_alarm_jobs : 0.0, cfg.stall_at,
            (long long)wc_wdg_ns, (long long)wc_ui_cs_ns);
    fclose(f);
    if (!cfg.ui) printf("%s/%s %s: fires=%ld max_latency=%.1f us misses(>%.0f ms)=%ld drops=%llu\n", cfg.fifo ? "fifo" : "other",
                        cfg.poll ? "poll" : "event", cfg.scenario, nfire, lmax, cfg.deadline_ms, nmiss, (unsigned long long)ch_samp.drops);
    return 0;
}
