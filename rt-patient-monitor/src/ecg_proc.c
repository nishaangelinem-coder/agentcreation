#include "ecg_proc.h"
#include <math.h>
#include <string.h>

#define REFRACT   (FS_ECG / 5)        /* 200 ms */
#define TWAVE_WIN (FS_ECG * 36 / 100) /* 360 ms */
#define CONFIRM   (FS_ECG * 8 / 100)  /* 80 ms without a higher value confirms a peak */
#define LEARN_LEN (2 * FS_ECG)
#define SETTLE    (FS_ECG / 2)
#define SAT_MV    4.8
/* Absolute floor of the detection threshold on the integrator output. It corresponds to a QRS of
 * roughly 0.2 mV (integrator peak for a 1 mV QRS is ~1.9e-3, scaling with amplitude squared), so
 * that noise-only input can never be adapted into "beats" (asystole must not be masked). */
#define THR_FLOOR 7.5e-5

static void bq_design(biquad_t *f, double fc, int highpass) {
    double w0 = 2.0 * M_PI * fc / FS_ECG, c = cos(w0), al = sin(w0) / (2.0 * M_SQRT1_2), a0 = 1.0 + al;
    if (highpass) { f->b0 = (1 + c) / 2 / a0; f->b1 = -(1 + c) / a0; f->b2 = f->b0; }
    else          { f->b0 = (1 - c) / 2 / a0; f->b1 =  (1 - c) / a0; f->b2 = f->b0; }
    f->a1 = -2 * c / a0; f->a2 = (1 - al) / a0; f->z1 = f->z2 = 0;
}
static inline double bq(biquad_t *f, double x) {
    double y = f->b0 * x + f->z1;
    f->z1 = f->b1 * x - f->a1 * y + f->z2;
    f->z2 = f->b2 * x - f->a2 * y;
    return y;
}

static void start_learning(ecg_proc_t *e, int64_t from) {
    e->learning = 1; e->learn_start = from; e->learn_end = from + SETTLE + LEARN_LEN;
    e->learn_max = 0; e->learn_sum = 0; e->learn_cnt = 0;
    e->rr_n = 0; e->rr_i = 0; e->hr_valid = 0; e->have_qrs = 0;
    e->cand_v = 0; e->cand_is_peak = 0; e->sb_v = 0;
    e->cond_tachy = e->cond_brady = 0;   /* rate conditions are meaningless without valid ECG */
}

void ecg_proc_init(ecg_proc_t *e, int lo_persist) {
    memset(e, 0, sizeof *e);
    bq_design(&e->hp, 5.0, 1); bq_design(&e->lp, 15.0, 0);
    e->lo_persist = lo_persist > 0 ? lo_persist : FS_ECG / 2;
    start_learning(e, 0);
}

static void accept_qrs(ecg_proc_t *e, double v, int64_t i, double w) {
    if (v > 4.0 * e->spki) v = 4.0 * e->spki;              /* bound the influence of one outlier */
    e->spki = w * v + (1.0 - w) * e->spki;
    if (e->have_qrs) {
        double rr = (double)(i - e->last_qrs) / FS_ECG;
        if (rr > 3.0) { e->rr_n = 0; e->rr_i = 0; e->hr_valid = 0; }   /* pause: restart rate estimate */
        else {
            e->rr[e->rr_i] = rr; e->rr_i = (e->rr_i + 1) & 3; if (e->rr_n < 4) e->rr_n++;
            if (e->rr_n == 4) {
                e->hr = 240.0 / (e->rr[0] + e->rr[1] + e->rr[2] + e->rr[3]); e->hr_valid = 1;
                if (e->hr > HR_TACHY_SET) e->cond_tachy = 1; else if (e->hr < HR_TACHY_CLR) e->cond_tachy = 0;
                if (e->hr < HR_BRADY_SET) e->cond_brady = 1; else if (e->hr > HR_BRADY_CLR) e->cond_brady = 0;
            }
        }
    }
    e->have_qrs = 1; e->last_qrs = i; e->last_qrs_v = v; e->ref_i = i;
    e->sb_v = 0; e->qrs_flag = 1; e->qrs_idx = i;
}

static void classify(ecg_proc_t *e, double v, int64_t i) {
    double thr = e->npki + 0.25 * (e->spki - e->npki); if (thr < THR_FLOOR) thr = THR_FLOOR;
    if (e->have_qrs && i - e->last_qrs < REFRACT) return;
    int twave = e->have_qrs && (i - e->last_qrs < TWAVE_WIN) && (v < 0.5 * e->last_qrs_v);
    if (v > thr && !twave) accept_qrs(e, v, i, 0.125);
    else {
        e->npki = 0.125 * v + 0.875 * e->npki;
        /* search-back candidates exclude the T-wave window: otherwise the T wave of the last beat
         * before a pause is accepted as a beat and delays the asystole alarm by ~0.5 s */
        if (v > e->sb_v && (!e->have_qrs || i - e->last_qrs >= TWAVE_WIN)) { e->sb_v = v; e->sb_i = i; }
    }
}

void ecg_proc_step(ecg_proc_t *e, double x) {
    int64_t n = e->n;
    e->qrs_flag = 0;

    /* 1. front-end saturation / leads-off supervision; hold last good value to avoid filter transients */
    int sat = fabs(x) >= SAT_MV;
    if (sat) { e->sat_run++; e->inhibit_until = n + FS_ECG * 3 / 10; x = e->last_good; }
    else     { e->sat_run = 0; e->last_good = x; }
    e->leads_off = (uint8_t)(e->sat_run >= e->lo_persist);

    /* 2. band-pass 5-15 Hz, 5-point derivative, squaring, moving-window integration */
    double y = bq(&e->lp, bq(&e->hp, x));
    e->d[4] = e->d[3]; e->d[3] = e->d[2]; e->d[2] = e->d[1]; e->d[1] = e->d[0]; e->d[0] = y;
    double dv = (2 * e->d[0] + e->d[1] - e->d[3] - 2 * e->d[4]) * 0.125;
    double sq = dv * dv;
    e->mwi_sum += sq - e->mwi_buf[e->mi]; e->mwi_buf[e->mi] = sq; e->mi = (e->mi + 1) % MWI_N;
    if (e->mi == 0) { double s = 0; for (int j = 0; j < MWI_N; j++) s += e->mwi_buf[j]; e->mwi_sum = s; }
    double m = e->mwi_sum / MWI_N; e->mwi = m;

    /* 3. signal invalid: restart threshold learning after the blanking interval */
    if (sat || n < e->inhibit_until) {
        if (!e->learning || e->learn_start < e->inhibit_until) start_learning(e, e->inhibit_until);
    }
    if (e->learning) {
        if (n >= e->learn_start + SETTLE) { if (m > e->learn_max) e->learn_max = m; e->learn_sum += m; e->learn_cnt++; }
        if (n >= e->learn_end) {
            e->spki = 0.5 * e->learn_max; e->npki = e->learn_cnt ? 0.5 * e->learn_sum / e->learn_cnt : 0;
            e->learning = 0; e->ref_i = n; e->cand_v = m; e->cand_is_peak = 0; e->cand_i = n;
        }
    } else {
        /* 4. local-maximum search with adaptive dual thresholds */
        if (m > e->cand_v) { e->cand_v = m; e->cand_i = n; e->cand_is_peak = 1; }
        else if (n - e->cand_i >= CONFIRM) {
            if (e->cand_is_peak) classify(e, e->cand_v, e->cand_i);
            e->cand_v = m; e->cand_i = n; e->cand_is_peak = 0;
        }
        /* 5. search-back for a missed beat */
        if (e->have_qrs && e->rr_n >= 2 && e->sb_v > 0) {
            double rrm = 0; for (int j = 0; j < e->rr_n; j++) rrm += e->rr[j]; rrm /= e->rr_n;
            double thr2 = 0.5 * (e->npki + 0.25 * (e->spki - e->npki)); if (thr2 < THR_FLOOR) thr2 = THR_FLOOR;
            if ((double)(n - e->last_qrs) > 1.66 * rrm * FS_ECG && e->sb_v > thr2 && e->sb_i - e->last_qrs >= REFRACT)
                accept_qrs(e, e->sb_v, e->sb_i, 0.25);
        }
    }

    /* 6. asystole: no QRS for ASYS_SEC while the ECG is valid */
    e->cond_asys = (uint8_t)(!e->learning && !e->leads_off && (n - e->ref_i) >= (int64_t)(ASYS_SEC * FS_ECG));
    if (e->hr_valid && n - e->last_qrs > 3 * FS_ECG) e->hr_valid = 0;
    e->n = n + 1;
}
