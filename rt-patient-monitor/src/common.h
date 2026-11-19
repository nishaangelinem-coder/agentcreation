/* common.h - shared types for the real-time patient-monitor simulator.
 * Portable POSIX.1-2008 C11; intended targets: QNX Neutrino (SDP 7.1/8.0) and Linux. */
#ifndef COMMON_H
#define COMMON_H
#include <stdint.h>

#define FS_ECG      500          /* ECG sampling rate, Hz            */
#define PPG_DIV     5            /* PPG sampled every 5th tick: 100 Hz */
#define FS_PPG      (FS_ECG / PPG_DIV)
#define PERIOD_NS   (1000000000LL / FS_ECG)
#define ECG_RAIL_MV 5.0f         /* simulated front-end saturation    */

typedef enum { RH_NSR = 0, RH_ASYSTOLE = 1, RH_LEADS_OFF = 2 } rhythm_t;

typedef struct {
    uint64_t k;          /* sample index                                   */
    int64_t  t_rel_ns;   /* nominal release time of this sample            */
    int64_t  t_acq_ns;   /* actual time the acquisition job started        */
    float    ecg_mv;
    float    red, ir;    /* PPG channels (valid if has_ppg)                */
    uint8_t  has_ppg;
    uint8_t  r_peak;     /* ground truth: R peak at this sample            */
    uint8_t  rhythm;     /* ground truth rhythm                            */
    uint8_t  probe;      /* ground truth: probe (saturation) sample        */
    float    true_hr, true_spo2;
} sample_t;

/* Alarm conditions. Priorities follow the IEC 60601-1-8 HIGH/MEDIUM/LOW scheme. */
typedef enum {
    AL_ASYSTOLE = 0, AL_TACHY, AL_BRADY, AL_DESAT, AL_SPO2_LOW,
    AL_LEADS_OFF, AL_LOW_PERF, AL_SYS_STALL, AL_DATA_LOSS, AL_N
} alarm_id_t;

typedef enum { PR_NONE = 0, PR_LOW = 1, PR_MED = 2, PR_HIGH = 3 } alarm_prio_t;

#endif
