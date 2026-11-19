#include "siggen.h"
#include <math.h>
#include <string.h>

#define TWO_PI 6.283185307179586
#define T_INF  1e18

static inline uint64_t rng_next(uint64_t *s) {            /* xorshift64* */
    uint64_t x = *s; x ^= x >> 12; x ^= x << 25; x ^= x >> 27; *s = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static inline double rng_uni(uint64_t *s) { return (double)(rng_next(s) >> 11) * (1.0 / 9007199254740992.0); }
static inline double rng_gauss(uint64_t *s) {              /* bounded-time approx. N(0,1) */
    double a = 0; for (int i = 0; i < 12; i++) a += rng_uni(s); return a - 6.0;
}

/* PQRST complex as a sum of Gaussians in time relative to the R peak (McSharry-style morphology). */
static double beat(double dt, double s) {
    static const double a[5]  = { 0.12, -0.10, 1.10, -0.22, 0.30 };
    static const double mu[5] = { -0.180, -0.030, 0.0, 0.030, 0.260 };
    static const double sd[5] = { 0.025, 0.008, 0.010, 0.009, 0.045 };
    if (dt < -0.5 || dt > 0.6) return 0.0;
    double v = 0;
    for (int i = 0; i < 5; i++) {
        double m = (i == 0 || i == 4) ? mu[i] * s : mu[i];
        double w = (i == 0 || i == 4) ? sd[i] * s : sd[i];
        double d = (dt - m) / w;
        if (d > -6 && d < 6) v += a[i] * exp(-0.5 * d * d);
    }
    return v;
}
static double circ(double a, double b) { double d = fabs(a - b); return d > 0.5 ? 1.0 - d : d; }
static double pulse_shape(double ph) {
    double d1 = circ(ph, 0.25) / 0.07, d2 = circ(ph, 0.55) / 0.10;
    return exp(-0.5 * d1 * d1) + 0.35 * exp(-0.5 * d2 * d2);
}

void siggen_init(siggen_t *g, const segment_t *seg, int nseg, uint64_t seed, double noise_sd, double amp) {
    memset(g, 0, sizeof *g);
    g->seg = seg; g->nseg = nseg; g->cur = 0; g->seg_end_t = seg[0].dur_s;
    g->rng = seed * 0x9E3779B97F4A7C15ULL + 0x1234567ULL; if (!g->rng) g->rng = 1;
    g->noise_sd = noise_sd; g->amp = amp;
    g->t_prev_r = -T_INF; g->rr_cur = 60.0 / (seg[0].hr > 0 ? seg[0].hr : 75.0);
    g->t_next_r = (seg[0].rhythm == RH_ASYSTOLE) ? T_INF : 0.3 * rng_uni(&g->rng) + 0.1;
    g->perf = (seg[0].rhythm == RH_ASYSTOLE) ? 0.0 : 1.0;
    g->probe_next = UINT64_MAX;
}
void siggen_enable_probe(siggen_t *g, int len, int gmin, int gmax) {
    g->probe_len = len; g->probe_gap_min = gmin; g->probe_gap_max = gmax; g->probe_next = 5 * FS_ECG;
}

void siggen_step(siggen_t *g, uint64_t k, sample_t *o) {
    double t = (double)k / FS_ECG;
    while (t >= g->seg_end_t && g->cur < g->nseg - 1) {
        g->cur++; g->seg_end_t += g->seg[g->cur].dur_s;
        if (g->cur < 64) { g->seg_start_k[g->cur] = k; g->seg_last_r_k[g->cur] = g->k_last_r; }
        /* entering asystole: a beat whose P wave has already started is completed (no truncated complex) */
        if (g->seg[g->cur].rhythm == RH_ASYSTOLE) { if (g->t_next_r - t > 0.3) g->t_next_r = T_INF; }
        else if (g->t_next_r >= T_INF) g->t_next_r = t + 0.3;
    }
    const segment_t *sg = &g->seg[g->cur];
    int beating = (sg->rhythm != RH_ASYSTOLE);

    o->r_peak = 0;
    if (t >= g->t_next_r) {
        o->r_peak = 1; g->k_last_r = k; g->t_prev_r = g->t_next_r;
        double rr = 60.0 / sg->hr * (1.0 + 0.02 * rng_gauss(&g->rng));
        if (rr < 0.2) rr = 0.2;
        g->rr_cur = rr; g->t_next_r = g->t_prev_r + rr;
        if (sg->rhythm == RH_ASYSTOLE) { g->t_next_r = T_INF; if (g->cur < 64) g->seg_last_r_k[g->cur] = k; }
    }
    double s = sqrt(g->rr_cur / 0.8); if (s < 0.6) s = 0.6; if (s > 1.2) s = 1.2;
    double ecg = g->amp * (beat(t - g->t_prev_r, s) + (g->t_next_r < T_INF ? beat(t - g->t_next_r, s) : 0.0));
    ecg += 0.05 * sin(TWO_PI * 0.25 * t) + 0.02 * sin(TWO_PI * 50.0 * t) + g->noise_sd * rng_gauss(&g->rng);

    o->probe = 0;
    if (sg->rhythm == RH_LEADS_OFF) ecg = ECG_RAIL_MV;
    if (g->probe_len > 0) {
        if (g->probe_rem > 0) { g->probe_rem--; ecg = ECG_RAIL_MV; o->probe = 1; }
        else if (k >= g->probe_next) {
            ecg = ECG_RAIL_MV; o->probe = 1; g->probe_rem = g->probe_len - 1;
            int gap = g->probe_gap_min + (int)(rng_uni(&g->rng) * (g->probe_gap_max - g->probe_gap_min + 1));
            g->probe_next = k + g->probe_len + gap;
        }
    }
    o->k = k; o->ecg_mv = (float)ecg; o->rhythm = (uint8_t)sg->rhythm;
    o->true_hr = beating ? (float)sg->hr : 0.0f; o->true_spo2 = (float)sg->spo2;

    o->has_ppg = (k % PPG_DIV == 0);
    if (o->has_ppg) {
        const double dt = 1.0 / FS_PPG;
        g->ppg_phase += dt / g->rr_cur; if (g->ppg_phase >= 1.0) g->ppg_phase -= 1.0;
        g->perf += ((beating ? 1.0 : 0.0) - g->perf) * (dt / 0.5);
        double R = (110.0 - sg->spo2) / 25.0, p = 0.02 * g->perf * pulse_shape(g->ppg_phase);
        o->ir  = (float)(1.0 - p     + 0.0004 * rng_gauss(&g->rng));
        o->red = (float)(1.0 - R * p + 0.0004 * rng_gauss(&g->rng));
    } else { o->ir = o->red = 0; }
}
