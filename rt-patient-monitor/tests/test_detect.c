/* test_detect.c - deterministic off-target verification of detection algorithms and alarm logic.
 * Runs the signal source through the same processing code as the real-time build, at full speed.
 * Emits one CSV row per test case: id,description,metric,value,limit,verdict                     */
#include "../src/alarm.h"
#include "../src/ecg_proc.h"
#include "../src/siggen.h"
#include "../src/spo2_proc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static FILE *out; static int n_fail;
static double g_match_from_s = 3.0;   /* beats before this time are excluded from Se/PPV (threshold learning) */
static void report(const char *id, const char *desc, const char *metric, double v, const char *op, double lim) {
    int ok = !strcmp(op, "<=") ? v <= lim : !strcmp(op, ">=") ? v >= lim : v == lim;
    if (!ok) n_fail++;
    fprintf(out, "%s,\"%s\",%s,%.4f,%s %.4g,%s\n", id, desc, metric, v, op, lim, ok ? "PASS" : "FAIL");
    printf("%-10s %-4s %-62s %s = %.4f (%s %.4g)\n", id, ok ? "PASS" : "FAIL", desc, metric, v, op, lim);
}

typedef struct {
    long tp, fn, fp;                 /* QRS detection vs ground truth          */
    double hr_err_max, spo2_err_max; /* steady-state numeric errors            */
    int64_t first_set[AL_N];         /* first sample index each condition true */
    long n_set[AL_N];
    double mwi_peak;
} res_t;

static void run(const segment_t *seg, int nseg, uint64_t seed, double noise, double amp, int lo_persist,
                double settle_s, res_t *r, siggen_t *gout) {
    static siggen_t g; static ecg_proc_t e; static spo2_proc_t s;
    siggen_init(&g, seg, nseg, seed, noise, amp); ecg_proc_init(&e, lo_persist); spo2_proc_init(&s);
    memset(r, 0, sizeof *r); for (int i = 0; i < AL_N; i++) r->first_set[i] = -1;
    double total = 0; for (int i = 0; i < nseg; i++) total += seg[i].dur_s;
    int64_t N = (int64_t)(total * FS_ECG);
    enum { MAXB = 20000 }; static int64_t truth[MAXB], det[MAXB]; int nt = 0, nd = 0;
    uint8_t prev[AL_N] = {0};
    int64_t seg_k0 = 0; int cur = 0;
    for (int64_t k = 0; k < N; k++) {
        sample_t sm; siggen_step(&g, (uint64_t)k, &sm);
        if (g.cur != cur) { cur = g.cur; seg_k0 = k; }
        ecg_proc_step(&e, sm.ecg_mv);
        if (sm.has_ppg) spo2_proc_step(&s, sm.red, sm.ir);
        if (e.mwi > r->mwi_peak) r->mwi_peak = e.mwi;
        if (sm.r_peak && nt < MAXB) truth[nt++] = k;
        if (e.qrs_flag && nd < MAXB) det[nd++] = e.qrs_idx;
        uint8_t c[AL_N] = {0};
        c[AL_ASYSTOLE] = e.cond_asys; c[AL_TACHY] = e.cond_tachy; c[AL_BRADY] = e.cond_brady;
        c[AL_DESAT] = s.cond_desat; c[AL_SPO2_LOW] = s.cond_low; c[AL_LEADS_OFF] = e.leads_off; c[AL_LOW_PERF] = s.cond_perf;
        for (int i = 0; i < AL_N; i++) { if (c[i] && !prev[i]) { r->n_set[i]++; if (r->first_set[i] < 0) r->first_set[i] = k; } prev[i] = c[i]; }
        if ((double)(k - seg_k0) / FS_ECG > settle_s && sm.rhythm == RH_NSR) {
            if (e.hr_valid) { double d = fabs(e.hr - sm.true_hr); if (d > r->hr_err_max) r->hr_err_max = d; }
            if (s.valid)    { double d = fabs(s.spo2 - sm.true_spo2); if (d > r->spo2_err_max) r->spo2_err_max = d; }
        }
    }
    /* beat matching: a detection within [0, 150 ms] after a true R peak (filter group delay) */
    int j = 0; int64_t k_ok = (int64_t)(g_match_from_s * FS_ECG);
    static char used[MAXB]; memset(used, 0, sizeof used);
    for (int i = 0; i < nt; i++) {
        if (truth[i] < k_ok) continue;
        while (j < nd && det[j] < truth[i] - FS_ECG / 20) j++;
        if (j < nd && det[j] <= truth[i] + FS_ECG * 15 / 100) { r->tp++; used[j] = 1; j++; } else r->fn++;
    }
    for (int i = 0; i < nd; i++) if (!used[i] && det[i] >= k_ok + FS_ECG * 15 / 100) r->fp++;
    if (gout) *gout = g;
}

static double maxd(double a, double b) { return a > b ? a : b; }

int main(int argc, char **argv) {
    out = fopen(argc > 1 ? argv[1] : "results/verification.csv", "w");
    if (!out) { perror("open"); return 2; }
    fprintf(out, "test_id,description,metric,value,limit,verdict\n");
    res_t r; siggen_t g;

    if (argc > 2 && !strcmp(argv[2], "calib")) {
        segment_t s[] = { { 30, RH_NSR, 75, 97 } };
        for (double a = 0.1; a <= 2.01; a += 0.1) { run(s, 1, 1, 0.02, a / 1.1, 0, 5, &r, NULL); printf("R=%.2f mV mwi_peak=%.3e tp=%ld fn=%ld\n", a, r.mwi_peak, r.tp, r.fn); }
        segment_t z[] = { { 30, RH_ASYSTOLE, 0, 97 } };
        for (double n = 0.02; n <= 0.11; n += 0.02) { run(z, 1, 1, n, 1, 0, 5, &r, NULL); printf("noise=%.2f mV mwi_peak=%.3e\n", n, r.mwi_peak); }
        return 0;
    }

    /* ---- TC-ECG-01/02: QRS detection and heart-rate accuracy over rate and amplitude ---- */
    { long tp = 0, fn = 0, fp = 0; double hre = 0, hre_rel = 0;
      const double hrs[] = { 40, 50, 60, 80, 100, 120, 150, 180 }, amps[] = { 0.5, 1.0, 2.0, 4.0 };
      for (unsigned a = 0; a < 4; a++) for (unsigned h = 0; h < 8; h++) for (uint64_t sd = 1; sd <= 3; sd++) {
          segment_t s[] = { { 120, RH_NSR, hrs[h], 97 } };
          run(s, 1, sd * 100 + h * 10 + a, 0.02, amps[a], 0, 15, &r, NULL);
          tp += r.tp; fn += r.fn; fp += r.fp;
          /* generator beat-to-beat variability (2 %) bounds the instantaneous 4-beat error */
          hre = maxd(hre, r.hr_err_max); hre_rel = maxd(hre_rel, r.hr_err_max / hrs[h] * 100.0);
      }
      report("TC-ECG-01", "QRS sensitivity, 40-180 bpm, R 0.55-4.4 mV (96 records, 3.2 h)", "Se_pct", 100.0 * tp / (tp + fn), ">=", 99.5);
      report("TC-ECG-01b", "QRS positive predictivity, same records", "PPV_pct", 100.0 * tp / (tp + fp), ">=", 99.5);
      report("TC-ECG-02", "Heart-rate error (4-beat mean) vs programmed rate, worst case", "err_pct", hre_rel, "<=", 10.0);
      printf("           beats=%ld fn=%ld fp=%ld max_abs_hr_err=%.2f bpm\n", tp + fn, fn, fp, hre);
    }
    /* ---- TC-ASY-01: asystole time-to-alarm from last QRS ---- */
    { double tmin = 1e9, tmax = 0; int missed = 0;
      const double amps[] = { 0.5, 1.0, 2.0 }, hrs[] = { 50, 75, 120 };
      for (uint64_t sd = 1; sd <= 20; sd++) for (int a = 0; a < 3; a++) for (int h = 0; h < 3; h++) {
          segment_t s[] = { { 20 + (double)sd * 0.137, RH_NSR, hrs[h], 97 }, { 12, RH_ASYSTOLE, 0, 97 }, { 10, RH_NSR, 75, 97 } };
          run(s, 3, sd * 7 + a, 0.02, amps[a], 0, 15, &r, &g);
          if (r.first_set[AL_ASYSTOLE] < 0) { missed++; continue; }
          double d = (double)(r.first_set[AL_ASYSTOLE] - (int64_t)g.seg_last_r_k[1]) / FS_ECG;
          if (d < tmin) tmin = d; if (d > tmax) tmax = d;
      }
      report("TC-ASY-01", "Asystole: missed alarms in 180 episodes", "missed", missed, "==", 0);
      report("TC-ASY-01b", "Asystole: max time-to-alarm from last R wave", "t_max_s", tmax, "<=", 4.5);
      report("TC-ASY-01c", "Asystole: min time-to-alarm (no premature alarm)", "t_min_s", tmin, ">=", 4.0);
    }
    /* ---- TC-ASY-02: asystole must not be masked by noise adaptation ---- */
    { int missed = 0; double tmax = 0;
      for (uint64_t sd = 1; sd <= 20; sd++) {
          segment_t s[] = { { 20, RH_NSR, 75, 97 }, { 30, RH_ASYSTOLE, 0, 97 } };
          run(s, 2, sd, 0.05, 1.0, 0, 15, &r, &g);       /* 0.05 mV rms noise, ~0.3 mV p-p */
          if (r.first_set[AL_ASYSTOLE] < 0) missed++;
          else tmax = maxd(tmax, (double)(r.first_set[AL_ASYSTOLE] - (int64_t)g.seg_last_r_k[1]) / FS_ECG);
          if (r.n_set[AL_ASYSTOLE] != 1) missed++;        /* must stay asserted for the whole episode */
      }
      report("TC-ASY-02", "Asystole with 0.05 mV rms noise: missed/interrupted alarms in 20 episodes", "missed", missed, "==", 0);
      /* ECG reconnected while the patient is already in asystole: thresholds learnt on noise only */
      int missed2 = 0; double t2 = 0;
      for (uint64_t sd = 1; sd <= 20; sd++) {
          segment_t s[] = { { 15, RH_NSR, 75, 97 }, { 5, RH_LEADS_OFF, 75, 97 }, { 30, RH_ASYSTOLE, 0, 97 } };
          run(s, 3, sd, 0.03, 1.0, 0, 15, &r, &g);
          if (r.first_set[AL_ASYSTOLE] < 0) missed2++;
          else t2 = maxd(t2, (double)(r.first_set[AL_ASYSTOLE] - (int64_t)g.seg_start_k[2]) / FS_ECG);
      }
      report("TC-ASY-03", "Asystole present at lead reconnection: missed alarms in 20 episodes", "missed", missed2, "==", 0);
      report("TC-ASY-03b", "Asystole present at lead reconnection: max time-to-alarm", "t_max_s", t2, "<=", 10.0);
    }
    /* ---- TC-TAC-01 / TC-BRA-01: rate alarms ---- */
    { double tmax = 0; int missed = 0; const double tg[] = { 160, 180, 220 };
      for (uint64_t sd = 1; sd <= 20; sd++) for (int h = 0; h < 3; h++) {
          segment_t s[] = { { 20 + (double)sd * 0.11, RH_NSR, 75, 97 }, { 20, RH_NSR, tg[h], 97 } };
          run(s, 2, sd, 0.02, 1.0, 0, 15, &r, &g);
          if (r.first_set[AL_TACHY] < 0) missed++; else tmax = maxd(tmax, (double)(r.first_set[AL_TACHY] - (int64_t)g.seg_start_k[1]) / FS_ECG);
      }
      report("TC-TAC-01", "Tachycardia (75 to 160/180/220 bpm): missed alarms in 60 episodes", "missed", missed, "==", 0);
      report("TC-TAC-01b", "Tachycardia: max time-to-alarm from rate change", "t_max_s", tmax, "<=", 10.0);
      tmax = 0; missed = 0; const double bg[] = { 30, 35 };
      for (uint64_t sd = 1; sd <= 20; sd++) for (int h = 0; h < 2; h++) {
          segment_t s[] = { { 20 + (double)sd * 0.11, RH_NSR, 75, 97 }, { 30, RH_NSR, bg[h], 97 } };
          run(s, 2, sd, 0.02, 1.0, 0, 15, &r, &g);
          if (r.first_set[AL_BRADY] < 0) missed++; else tmax = maxd(tmax, (double)(r.first_set[AL_BRADY] - (int64_t)g.seg_start_k[1]) / FS_ECG);
      }
      report("TC-BRA-01", "Bradycardia (75 to 30/35 bpm): missed alarms in 40 episodes", "missed", missed, "==", 0);
      report("TC-BRA-01b", "Bradycardia: max time-to-alarm from rate change", "t_max_s", tmax, "<=", 10.0);
    }
    /* ---- TC-SPO-01/02: SpO2 accuracy and desaturation alarms ---- */
    { double emax = 0;
      for (double sp = 70; sp <= 100; sp += 5) for (double hr = 50; hr <= 150; hr += 50) {
          segment_t s[] = { { 60, RH_NSR, hr, sp } };
          run(s, 1, (uint64_t)(sp + hr), 0.02, 1.0, 0, 10, &r, NULL); emax = maxd(emax, r.spo2_err_max);
      }
      report("TC-SPO-01", "SpO2 estimate error, 70-100 %, 50-150 bpm", "err_max_pct", emax, "<=", 2.0);
      double t1 = 0, t2 = 0; int missed = 0, wrong = 0;
      for (uint64_t sd = 1; sd <= 20; sd++) {
          segment_t s[] = { { 20 + (double)sd * 0.13, RH_NSR, 75, 97 }, { 30, RH_NSR, 75, 80 } };
          run(s, 2, sd, 0.02, 1.0, 0, 15, &r, &g);
          if (r.first_set[AL_DESAT] < 0) missed++; else t1 = maxd(t1, (double)(r.first_set[AL_DESAT] - (int64_t)g.seg_start_k[1]) / FS_ECG);
          segment_t q[] = { { 20 + (double)sd * 0.13, RH_NSR, 75, 97 }, { 30, RH_NSR, 75, 88 } };
          run(q, 2, sd, 0.02, 1.0, 0, 15, &r, &g);
          if (r.first_set[AL_SPO2_LOW] < 0) missed++; else t2 = maxd(t2, (double)(r.first_set[AL_SPO2_LOW] - (int64_t)g.seg_start_k[1]) / FS_ECG);
          if (r.n_set[AL_DESAT]) wrong++;                 /* 88 % must not raise the <85 % alarm */
      }
      report("TC-SPO-02", "Desaturation alarms: missed or wrong-level alarms in 40 episodes", "errors", missed + wrong, "==", 0);
      report("TC-SPO-02b", "DESAT (97 to 80 %): max time-to-alarm", "t_max_s", t1, "<=", 10.0);
      report("TC-SPO-02c", "SPO2_LOW (97 to 88 %): max time-to-alarm", "t_max_s", t2, "<=", 15.0);
    }
    /* ---- TC-LO-01: leads-off ---- */
    { double dmax = 0; int bad = 0;
      for (uint64_t sd = 1; sd <= 20; sd++) {
          segment_t s[] = { { 20, RH_NSR, 75, 97 }, { 15, RH_LEADS_OFF, 75, 97 }, { 30, RH_NSR, 75, 97 } };
          g_match_from_s = 20 + 15 + 3.0;       /* judge recovery after reconnection + re-learning */
          run(s, 3, sd, 0.02, 1.0, 0, 15, &r, &g); g_match_from_s = 3.0;
          if (r.first_set[AL_LEADS_OFF] < 0) { bad++; continue; }
          dmax = maxd(dmax, fabs((double)(r.first_set[AL_LEADS_OFF] - (int64_t)g.seg_start_k[1]) / FS_ECG - 0.5));
          for (int i = 0; i < AL_N; i++) if (i != AL_LEADS_OFF && r.n_set[i]) bad++;   /* no secondary false alarms */
          if (r.fn > 1 || r.fp > 0) bad++;                                                       /* detector recovers after reconnection */
      }
      report("TC-LO-01", "Leads-off: missed alarms, secondary false alarms or failed recovery (20 episodes)", "errors", bad, "==", 0);
      report("TC-LO-01b", "Leads-off: deviation from 0.5 s persistence time", "dev_s", dmax, "<=", 0.004);
    }
    /* ---- TC-PERF-01: loss of pulse on PPG during asystole ---- */
    { int missed = 0; double t = 0;
      for (uint64_t sd = 1; sd <= 10; sd++) {
          segment_t s[] = { { 20, RH_NSR, 75, 97 }, { 20, RH_ASYSTOLE, 0, 97 } };
          run(s, 2, sd, 0.02, 1.0, 0, 15, &r, &g);
          if (r.first_set[AL_LOW_PERF] < 0) missed++; else t = maxd(t, (double)(r.first_set[AL_LOW_PERF] - (int64_t)g.seg_start_k[1]) / FS_ECG);
          if (r.n_set[AL_DESAT] || r.n_set[AL_SPO2_LOW]) missed++;   /* no spurious SpO2 values without a pulse */
      }
      report("TC-PERF-01", "Low-perfusion alarm on pulse loss: errors in 10 episodes", "errors", missed, "==", 0);
      report("TC-PERF-01b", "Low-perfusion alarm: max time-to-alarm", "t_max_s", t, "<=", 10.0);
    }
    /* ---- TC-FA-01: false alarms on normal rhythm ---- */
    { long fa = 0; double hours = 0;
      for (uint64_t sd = 1; sd <= 24; sd++) {
          double hr = 50 + (double)((sd * 37) % 90);            /* 50..139 bpm */
          segment_t s[] = { { 600, RH_NSR, hr, 93 + (double)(sd % 7) } };
          run(s, 1, sd, 0.03, 0.5 + 0.25 * (double)(sd % 5), 0, 15, &r, NULL);
          for (int i = 0; i < AL_N; i++) fa += r.n_set[i];
          hours += 600.0 / 3600.0;
      }
      report("TC-FA-01", "False alarms on 4 h of normal rhythm (50-139 bpm, 93-99 %)", "false_alarms", (double)fa, "==", 0);
      (void)hours;
    }
    /* ---- TC-ALM-01..04: alarm-system logic ---- */
    { alarm_mgr_t m; int bad = 0;
      alarm_mgr_init(&m);
      if (alarm_mgr_output(&m) != PR_NONE) bad++;
      if (!alarm_mgr_event(&m, AL_LOW_PERF, 1, 1, 0, 10) || m.out != PR_LOW) bad++;
      if (!alarm_mgr_event(&m, AL_SPO2_LOW, 1, 2, 0, 20) || m.out != PR_MED) bad++;
      if (!alarm_mgr_event(&m, AL_ASYSTOLE, 1, 3, 0, 30) || m.out != PR_HIGH) bad++;
      if (alarm_mgr_event(&m, AL_ASYSTOLE, 1, 4, 0, 40)) bad++;                 /* no duplicate annunciation */
      report("TC-ALM-01", "Priority arbitration: highest active priority drives the output", "errors", bad, "==", 0);
      bad = 0;
      alarm_mgr_event(&m, AL_ASYSTOLE, 0, 5, 0, 50);
      if (!m.s[AL_ASYSTOLE].active || m.out != PR_HIGH) bad++;                  /* latched after condition ends */
      alarm_mgr_event(&m, AL_SPO2_LOW, 0, 6, 0, 60);
      if (m.s[AL_SPO2_LOW].active) bad++;                                       /* non-latching alarm clears */
      report("TC-ALM-02", "Latching: high-priority physiological alarm persists, non-latching clears", "errors", bad, "==", 0);
      bad = 0;
      alarm_mgr_event(&m, AL_TACHY, 1, 7, 0, 70);
      if (alarm_mgr_ack(&m, AL_TACHY)) bad++;                                   /* cannot dismiss an ongoing condition */
      if (!m.s[AL_TACHY].active) bad++;
      if (!alarm_mgr_ack(&m, AL_ASYSTOLE) || m.s[AL_ASYSTOLE].active) bad++;    /* ended condition is released */
      if (m.out != PR_HIGH) bad++;                                              /* TACHY still active */
      alarm_mgr_event(&m, AL_TACHY, 0, 8, 0, 80); alarm_mgr_ack(&m, AL_TACHY);
      if (m.out != PR_LOW) bad++;
      report("TC-ALM-03", "Acknowledge: releases ended alarms only; ongoing conditions cannot be dismissed", "errors", bad, "==", 0);
      bad = 0;
      if (!alarm_mgr_event(&m, AL_ASYSTOLE, 1, 9, 0, 90)) bad++;                /* re-fires after release */
      if (m.s[AL_ASYSTOLE].n_fired != 2 || m.s[AL_ASYSTOLE].t_fire != 90) bad++;
      report("TC-ALM-04", "Re-annunciation of a recurring condition after release", "errors", bad, "==", 0);
    }
    fclose(out);
    printf("\n%s (%d failing)\n", n_fail ? "VERIFICATION FAILED" : "ALL VERIFICATION TESTS PASSED", n_fail);
    return n_fail ? 1 : 0;
}
