/* dump_wave.c - off-target trace of the clinical scenario for figures: k,ecg,mwi,qrs,hr,spo2,conditions */
#include "../src/ecg_proc.h"
#include "../src/siggen.h"
#include "../src/spo2_proc.h"
#include <stdio.h>
int main(void) {
    static const segment_t sc[] = { { 20, RH_NSR, 75, 97 }, { 10, RH_ASYSTOLE, 0, 97 }, { 20, RH_NSR, 75, 97 }, { 20, RH_NSR, 180, 97 },
        { 20, RH_NSR, 75, 97 }, { 25, RH_NSR, 32, 97 }, { 20, RH_NSR, 75, 97 }, { 25, RH_NSR, 75, 80 },
        { 20, RH_NSR, 75, 97 }, { 10, RH_LEADS_OFF, 75, 97 }, { 20, RH_NSR, 75, 97 } };
    static siggen_t g; static ecg_proc_t e; static spo2_proc_t s;
    siggen_init(&g, sc, 11, 11, 0.02, 1.0); ecg_proc_init(&e, 0); spo2_proc_init(&s);
    printf("k,ecg,mwi,qrs,hr,hr_ok,ir,spo2,spo2_ok,asys,tachy,brady,desat,low,lo,perf\n");
    for (long k = 0; k < 210 * FS_ECG; k++) {
        sample_t sm; siggen_step(&g, (unsigned long)k, &sm); ecg_proc_step(&e, sm.ecg_mv);
        if (sm.has_ppg) spo2_proc_step(&s, sm.red, sm.ir);
        printf("%ld,%.4f,%.3e,%d,%.1f,%d,%.5f,%.2f,%d,%d,%d,%d,%d,%d,%d,%d\n", k, sm.ecg_mv, e.mwi, e.qrs_flag, e.hr, e.hr_valid,
               sm.has_ppg ? sm.ir : -1.0, s.spo2, s.valid, e.cond_asys, e.cond_tachy, e.cond_brady, s.cond_desat, s.cond_low, e.leads_off, s.cond_perf);
    }
    return 0;
}
