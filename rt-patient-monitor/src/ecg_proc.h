/* ecg_proc.h - Pan-Tompkins-style QRS detector (500 Hz) and rhythm alarm conditions. */
#ifndef ECG_PROC_H
#define ECG_PROC_H
#include "common.h"
#define MWI_N 75                 /* 150 ms integration window */
typedef struct { double b0, b1, b2, a1, a2, z1, z2; } biquad_t;
typedef struct {
    biquad_t hp, lp; double d[5]; double mwi_buf[MWI_N]; int mi; double mwi_sum;
    int64_t n;
    double cand_v; int64_t cand_i; int cand_is_peak;
    double spki, npki; int learning; int64_t learn_start, learn_end; double learn_max, learn_sum; int64_t learn_cnt;
    int64_t last_qrs, ref_i; double last_qrs_v; int have_qrs;
    double rr[4]; int rr_n, rr_i; double hr; int hr_valid;
    double sb_v; int64_t sb_i;
    int sat_run; int64_t inhibit_until; double last_good;
    int lo_persist;
    /* outputs */
    int qrs_flag; int64_t qrs_idx;
    uint8_t leads_off, cond_asys, cond_tachy, cond_brady;
    double mwi;                    /* last integrator output (for calibration/tests) */
} ecg_proc_t;

#define ASYS_SEC      4.0
#define HR_TACHY_SET  150.0
#define HR_TACHY_CLR  140.0
#define HR_BRADY_SET  40.0
#define HR_BRADY_CLR  45.0

void ecg_proc_init(ecg_proc_t *e, int lo_persist_samples);
void ecg_proc_step(ecg_proc_t *e, double x_mv);
#endif
