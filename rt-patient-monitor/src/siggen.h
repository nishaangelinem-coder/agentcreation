/* siggen.h - deterministic synthetic ECG + dual-wavelength PPG source with scripted scenarios. */
#ifndef SIGGEN_H
#define SIGGEN_H
#include "common.h"
typedef struct { double dur_s; rhythm_t rhythm; double hr; double spo2; } segment_t;

typedef struct {
    const segment_t *seg; int nseg; int cur; double seg_end_t;
    double t_prev_r, t_next_r, rr_cur;
    double ppg_phase, perf;
    uint64_t rng; double noise_sd, amp;
    /* latency-probe overlay: short front-end saturation bursts at random intervals */
    int probe_len, probe_gap_min, probe_gap_max, probe_rem; uint64_t probe_next;
    /* ground truth bookkeeping */
    uint64_t k_last_r; uint64_t seg_start_k[64]; uint64_t seg_last_r_k[64];
} siggen_t;

void siggen_init(siggen_t *g, const segment_t *seg, int nseg, uint64_t seed, double noise_sd, double amp);
void siggen_enable_probe(siggen_t *g, int len, int gap_min, int gap_max);
void siggen_step(siggen_t *g, uint64_t k, sample_t *out);
#endif
