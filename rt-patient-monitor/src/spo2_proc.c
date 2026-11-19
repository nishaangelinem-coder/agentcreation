#include "spo2_proc.h"
#include <math.h>
#include <string.h>

void spo2_proc_init(spo2_proc_t *s) { memset(s, 0, sizeof *s); }

void spo2_proc_step(spo2_proc_t *s, double red, double ir) {
    int i = s->idx;
    /* DC: moving average over 3 s; AC: residual, RMS over 3 s */
    s->sxr += red - s->xr[i]; s->xr[i] = red;
    s->sxi += ir  - s->xi[i]; s->xi[i] = ir;
    int filled = s->n + 1 >= SPO2_W;
    double dcr = s->sxr / SPO2_W, dci = s->sxi / SPO2_W;
    double acr = filled ? red - dcr : 0, aci = filled ? ir - dci : 0;
    s->sar += acr * acr - s->ar[i]; s->ar[i] = acr * acr;
    s->sai += aci * aci - s->ai[i]; s->ai[i] = aci * aci;
    s->idx = (i + 1) % SPO2_W;
    if (s->idx == 0) {                       /* periodic exact recomputation: no drift of running sums */
        double a = 0, b = 0, c = 0, d = 0;
        for (int j = 0; j < SPO2_W; j++) { a += s->xr[j]; b += s->xi[j]; c += s->ar[j]; d += s->ai[j]; }
        s->sxr = a; s->sxi = b; s->sar = c; s->sai = d;
    }
    s->n++;
    int warm = s->n >= 2 * SPO2_W;
    if (warm && dcr > 0 && dci > 0) {
        double rr = sqrt(fmax(s->sar, 0) / SPO2_W) / dcr, ri = sqrt(fmax(s->sai, 0) / SPO2_W) / dci;
        s->pi = ri;
        if (ri >= PERF_MIN) {
            double v = 110.0 - 25.0 * (rr / ri);
            s->spo2 = v > 100 ? 100 : (v < 0 ? 0 : v); s->valid = 1; s->cnt_perf = 0;
        } else { s->valid = 0; s->cnt_perf++; }
    } else { s->valid = 0; s->cnt_perf = 0; }

    s->cnt_desat = (s->valid && s->spo2 < SPO2_DESAT_LIM) ? s->cnt_desat + 1 : 0;
    s->cnt_low   = (s->valid && s->spo2 < SPO2_LOW_LIM)   ? s->cnt_low + 1   : 0;
    s->cond_desat = (uint8_t)(s->cnt_desat >= SPO2_DESAT_SEC * FS_PPG);
    s->cond_low   = (uint8_t)(s->cnt_low   >= SPO2_LOW_SEC * FS_PPG);
    s->cond_perf  = (uint8_t)(s->cnt_perf  >= PERF_SEC * FS_PPG);
}
