/* spo2_proc.h - ratio-of-ratios SpO2 estimator (100 Hz) with desaturation / perfusion conditions. */
#ifndef SPO2_PROC_H
#define SPO2_PROC_H
#include "common.h"
#define SPO2_W (3 * FS_PPG)
typedef struct {
    double xr[SPO2_W], xi[SPO2_W], ar[SPO2_W], ai[SPO2_W];
    double sxr, sxi, sar, sai; int idx; int64_t n;
    double spo2, pi; int valid;
    int cnt_desat, cnt_low, cnt_perf;
    uint8_t cond_desat, cond_low, cond_perf;
} spo2_proc_t;
#define SPO2_DESAT_LIM 85.0
#define SPO2_DESAT_SEC 5
#define SPO2_LOW_LIM   90.0
#define SPO2_LOW_SEC   10
#define PERF_MIN       0.0015
#define PERF_SEC       3
void spo2_proc_init(spo2_proc_t *s);
void spo2_proc_step(spo2_proc_t *s, double red, double ir);
#endif
