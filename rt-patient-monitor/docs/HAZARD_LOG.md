# Hazard log and traceability report

*Generated 2026-11-19 by `scripts/gen_hazard_log.py` from executed verification results. Do not edit by hand.*

> **Scope and status.** This document is *styled after* the risk-management file required by ISO 14971 and the
> traceability chain of IEC 62304 (hazardous situation -> software item -> software cause -> risk control measure ->
> verification). It belongs to a research/teaching **simulator**, not a medical device. It has not been reviewed by a
> clinical or regulatory expert, severity/probability ratings are illustrative engineering judgements with no field
> data behind them, and nothing here constitutes evidence of conformity to any standard.

## 1. Risk acceptability scheme

Severity S: 1 negligible, 2 minor, 3 serious, 4 critical, 5 catastrophic. Probability of harm P: 1 improbable ... 5 frequent.
Following IEC 62304, the probability that a software fault *occurs* is taken as 1; P rates the likelihood that the
resulting hazardous situation leads to harm, before and after the risk control measures.
Risk index R = S x P: **>= 12 unacceptable**, **6-10 tolerable** (reduce as far as possible), **<= 5 acceptable**.

Software safety classification (IEC 62304, illustrative): the alarm path (acquisition, processing, alarm manager,
watchdog) is treated as **Class C**; the display/operator thread as Class B, segregated by scheduling class and priority.

## 2. Hazard log

Verification summary: **34 PASS, 2 FAIL** out of 36 executed checks.

### H-01 - Life-threatening arrhythmia (asystole) not annunciated

| | |
|---|---|
| Software cause | QRS detector adapts its threshold to noise and counts noise as beats; search-back accepts a T wave as a beat |
| Hazardous situation | Patient in cardiac standstill, monitor shows a heart rate and stays silent |
| Harm | Death or irreversible injury from delayed resuscitation |
| Software item | `ecg_proc` |
| Initial risk | S5 x P3 = 15 (Unacceptable) |
| Risk control measures | Absolute threshold floor (0.2 mV-equivalent); T-wave window excluded from search-back; asystole timer on validated beats only |
| Control type | Inherent safety by design |
| Requirements | SRS-ECG-03, SRS-ALM-01 |
| Verification | TC-ASY-01 (PASS), TC-ASY-01b (PASS), TC-ASY-01c (PASS), TC-ASY-02 (PASS), TC-ASY-03 (PASS), TC-ASY-03b (PASS), TC-RT-04 (PASS), TC-RT-04b (PASS), TC-RT-04c (PASS) |
| Residual risk | S5 x P1 = 5 (Acceptable) |
| Status | **Verified** |
| Residual-risk note | Ventricular fibrillation is not classified by this prototype (out of scope); coarse VF may be counted as beats. |

### H-02 - Alarm annunciated too late

| | |
|---|---|
| Software cause | Alarm path delayed by lower-criticality software (display, logging, other processes) through CPU contention |
| Hazardous situation | Alarm condition detected internally but output delayed by an unbounded time |
| Harm | Delayed clinical response |
| Software item | `monitor (scheduling)` |
| Initial risk | S4 x P3 = 12 (Unacceptable) |
| Risk control measures | Fixed-priority preemptive scheduling with alarm manager highest; event-driven (non-polled) alarm path; 10 ms hard deadline verified under overload |
| Control type | Inherent safety by design |
| Requirements | SRS-RT-01, SRS-RT-02 |
| Verification | TC-INSP-01 (PASS), TC-RT-01 (PASS), TC-RT-01b (FAIL) |
| Residual risk | S4 x P1 = 4 (Acceptable) |
| Status | **OPEN: control not verified** |
| Residual-risk note | Bound demonstrated empirically on the evaluation platform only; a certified RTOS/board latency bound is required for a product. |

### H-03 - Alarm annunciated too late

| | |
|---|---|
| Software cause | Priority inversion: display thread holds a mutex needed by the processing thread and is pre-empted by unrelated load |
| Hazardous situation | Processing thread blocked for an unbounded time |
| Harm | Delayed clinical response |
| Software item | `chan / monitor` |
| Initial risk | S4 x P3 = 12 (Unacceptable) |
| Risk control measures | Priority-inheritance protocol on every shared mutex; critical sections limited to fixed-size copies |
| Control type | Inherent safety by design |
| Requirements | SRS-RT-02 |
| Verification | TC-INSP-01 (PASS), TC-RT-01 (PASS), TC-RT-01b (FAIL) |
| Residual risk | S4 x P1 = 4 (Acceptable) |
| Status | **OPEN: control not verified** |

### H-04 - Alarm annunciated too late

| | |
|---|---|
| Software cause | Page fault, heap allocation or file I/O inside a real-time thread |
| Hazardous situation | Real-time thread suspended on a kernel service with unbounded duration |
| Harm | Delayed clinical response |
| Software item | `monitor / rt_util` |
| Initial risk | S4 x P3 = 12 (Unacceptable) |
| Risk control measures | All buffers pre-allocated before the real-time phase; memory locked (mlockall); logs written only after the run |
| Control type | Inherent safety by design |
| Requirements | SRS-RT-03 |
| Verification | TC-INSP-01 (PASS) |
| Residual risk | S4 x P1 = 4 (Acceptable) |
| Status | **Verified** |

### H-05 - Physiological data silently lost or stale

| | |
|---|---|
| Software cause | Processing thread hangs (runaway loop) or acquisition stops |
| Hazardous situation | Waveform frozen, no physiological alarm can be generated |
| Harm | Missed alarm; death or serious injury |
| Software item | `monitor (watchdog)` |
| Initial risk | S5 x P2 = 10 (Tolerable) |
| Risk control measures | Independent watchdog thread above the supervised threads; SYS_STALL HIGH-priority technical alarm within 200 ms |
| Control type | Protective measure |
| Requirements | SRS-WDG-01 |
| Verification | TC-WDG-01 (PASS), TC-WDG-01b (PASS) |
| Residual risk | S5 x P1 = 5 (Acceptable) |
| Status | **Verified** |
| Residual-risk note | A stalled alarm-manager thread or a kernel hang is not covered: an external hardware watchdog is required. |

### H-06 - Physiological data silently lost

| | |
|---|---|
| Software cause | Sample queue overflow under overload |
| Hazardous situation | Gaps in the ECG seen by the detector; beats missed or mis-timed |
| Harm | False or missed rate alarms |
| Software item | `chan` |
| Initial risk | S3 x P3 = 9 (Tolerable) |
| Risk control measures | Bounded queue with drop counter; DATA_LOSS technical alarm on any drop; queue sized for 1 s of backlog |
| Control type | Protective measure |
| Requirements | SRS-ACQ-02 |
| Verification | TC-RT-03 (PASS) |
| Residual risk | S3 x P1 = 3 (Acceptable) |
| Status | **Verified** |

### H-07 - False asystole alarm after lead reconnection

| | |
|---|---|
| Software cause | Saturation transient drives the adaptive threshold far above true QRS amplitude |
| Hazardous situation | Repeated false HIGH-priority alarms |
| Harm | Alarm fatigue leading to a missed true alarm |
| Software item | `ecg_proc` |
| Initial risk | S3 x P4 = 12 (Unacceptable) |
| Risk control measures | Hold-last-value during saturation; blanking and threshold re-learning after reconnection; bounded single-beat influence on the threshold |
| Control type | Inherent safety by design |
| Requirements | SRS-ALM-06, SRS-ALM-10 |
| Verification | TC-ASY-03 (PASS), TC-ASY-03b (PASS), TC-FA-01 (PASS), TC-LO-01 (PASS), TC-LO-01b (PASS), TC-RT-04 (PASS), TC-RT-04b (PASS), TC-RT-04c (PASS) |
| Residual risk | S3 x P2 = 6 (Tolerable) |
| Status | **Verified** |

### H-08 - Asystole masked by leads-off condition

| | |
|---|---|
| Software cause | ECG rhythm alarms are inhibited while electrodes are disconnected |
| Hazardous situation | Patient deteriorates while ECG is invalid |
| Harm | Delayed clinical response |
| Software item | `ecg_proc / spo2_proc` |
| Initial risk | S5 x P2 = 10 (Tolerable) |
| Risk control measures | LEADS_OFF technical alarm within 0.5 s; independent pulse supervision through PPG perfusion (LOW_PERF); asystole detected within 10 s of reconnection even if thresholds were learnt on noise |
| Control type | Protective measure + information for safety |
| Requirements | SRS-ALM-06, SRS-ALM-07 |
| Verification | TC-ASY-03 (PASS), TC-ASY-03b (PASS), TC-LO-01 (PASS), TC-LO-01b (PASS), TC-PERF-01 (PASS), TC-PERF-01b (PASS), TC-RT-04 (PASS), TC-RT-04b (PASS), TC-RT-04c (PASS) |
| Residual risk | S5 x P1 = 5 (Acceptable) |
| Status | **Verified** |
| Residual-risk note | LEADS_OFF is MEDIUM and LOW_PERF is LOW priority; escalation when both are present is a recommended extension. |

### H-09 - Wrong SpO2 value displayed

| | |
|---|---|
| Software cause | Ratio-of-ratios evaluated on noise when the pulse is absent or perfusion is very low (R tends to 1, i.e. 85 %) |
| Hazardous situation | Plausible but meaningless SpO2 numeric shown |
| Harm | Hypoxaemia not recognised or unnecessary intervention |
| Software item | `spo2_proc` |
| Initial risk | S4 x P3 = 12 (Unacceptable) |
| Risk control measures | Perfusion-index validity gate: numeric blanked and SpO2 limit alarms suppressed, LOW_PERF raised |
| Control type | Inherent safety by design |
| Requirements | SRS-ALM-07, SRS-SPO-01 |
| Verification | TC-PERF-01 (PASS), TC-PERF-01b (PASS), TC-SPO-01 (PASS) |
| Residual risk | S4 x P1 = 4 (Acceptable) |
| Status | **Verified** |
| Residual-risk note | Empirical calibration curve and motion artefact are not modelled. |

### H-10 - Hypoxaemia not annunciated

| | |
|---|---|
| Software cause | Desaturation detector fails or is delayed |
| Hazardous situation | Sustained SpO2 below limit without alarm |
| Harm | Hypoxic injury |
| Software item | `spo2_proc` |
| Initial risk | S4 x P2 = 8 (Tolerable) |
| Risk control measures | Two-level SpO2 limits with bounded persistence times (5 s / 10 s) |
| Control type | Inherent safety by design |
| Requirements | SRS-ALM-04, SRS-ALM-05 |
| Verification | TC-RT-04 (PASS), TC-RT-04b (PASS), TC-RT-04c (PASS), TC-SPO-02 (PASS), TC-SPO-02b (PASS), TC-SPO-02c (PASS) |
| Residual risk | S4 x P1 = 4 (Acceptable) |
| Status | **Verified** |

### H-11 - High-priority alarm hidden by lower-priority alarm

| | |
|---|---|
| Software cause | Annunciation follows the most recent instead of the most severe alarm |
| Hazardous situation | HIGH-priority condition presented with a MEDIUM/LOW signal |
| Harm | Delayed clinical response |
| Software item | `alarm` |
| Initial risk | S4 x P2 = 8 (Tolerable) |
| Risk control measures | Priority arbitration over all active alarms on every state change |
| Control type | Inherent safety by design |
| Requirements | SRS-ALM-08 |
| Verification | TC-ALM-01 (PASS) |
| Residual risk | S4 x P1 = 4 (Acceptable) |
| Status | **Verified** |

### H-12 - Transient life-threatening event unnoticed

| | |
|---|---|
| Software cause | Alarm clears automatically when the condition ends before staff attend |
| Hazardous situation | Self-terminating asystole or tachycardia leaves no trace |
| Harm | Missed diagnosis; recurrence without preparation |
| Software item | `alarm` |
| Initial risk | S4 x P3 = 12 (Unacceptable) |
| Risk control measures | Latching of HIGH-priority physiological alarms until operator acknowledgement; ongoing conditions cannot be dismissed |
| Control type | Inherent safety by design |
| Requirements | SRS-ALM-09 |
| Verification | TC-ALM-02 (PASS), TC-ALM-03 (PASS), TC-ALM-04 (PASS) |
| Residual risk | S4 x P1 = 4 (Acceptable) |
| Status | **Verified** |
| Residual-risk note | Audio-pause and alarm-limit configuration are not implemented. |

### H-13 - Excessive false alarms

| | |
|---|---|
| Software cause | Detector false positives or negatives on normal rhythm; rate estimate unstable |
| Hazardous situation | Frequent non-actionable alarms |
| Harm | Alarm fatigue leading to a missed true alarm |
| Software item | `ecg_proc / spo2_proc` |
| Initial risk | S3 x P4 = 12 (Unacceptable) |
| Risk control measures | Detector performance requirements; 4-beat averaging with hysteresis on rate limits; persistence times on SpO2 limits |
| Control type | Inherent safety by design |
| Requirements | SRS-ECG-01, SRS-ECG-02, SRS-ALM-10 |
| Verification | TC-ECG-01 (PASS), TC-ECG-01b (PASS), TC-ECG-02 (PASS), TC-FA-01 (PASS) |
| Residual risk | S3 x P2 = 6 (Tolerable) |
| Status | **Verified** |
| Residual-risk note | Verified on synthetic signals only; clinical false-alarm rate is unknown. |

### H-14 - Rate alarms mistimed

| | |
|---|---|
| Software cause | Sampling jitter or drift distorts RR intervals and persistence timers |
| Hazardous situation | Heart-rate and timing-based conditions evaluated on a wrong time base |
| Harm | False or late alarms |
| Software item | `monitor (acquisition)` |
| Initial risk | S3 x P2 = 6 (Tolerable) |
| Risk control measures | Absolute-time periodic release; timers counted in sample indices, not wall-clock reads |
| Control type | Inherent safety by design |
| Requirements | SRS-ACQ-01, SRS-RT-04 |
| Verification | TC-RT-02 (FAIL) |
| Residual risk | S3 x P1 = 3 (Acceptable) |
| Status | **OPEN: control not verified** |

## 3. Software requirements and verification status

| ID | Requirement | Tests | Status |
|---|---|---|---|
| SRS-ACQ-01 | ECG shall be acquired at 500 Hz and PPG at 100 Hz by a periodic task released on an absolute time base (no cumulative drift). | TC-RT-02 | **FAIL** |
| SRS-ACQ-02 | No acquired sample shall be lost between acquisition and processing; any loss shall raise the HIGH-priority technical alarm DATA_LOSS. | TC-RT-03 | **PASS** |
| SRS-ECG-01 | QRS detection sensitivity and positive predictivity shall each be >= 99.5 % for R-wave amplitudes 0.5-4.4 mV and heart rates 40-180 bpm on the reference signal set. | TC-ECG-01, TC-ECG-01b | **PASS** |
| SRS-ECG-02 | Displayed heart rate shall be within +/-10 % of the true rate. | TC-ECG-02 | **PASS** |
| SRS-ECG-03 | The QRS detection threshold shall have an absolute floor so that signals below 0.15 mV (noise, isoelectric line) are never counted as beats. | TC-ASY-02, TC-ASY-03, TC-ASY-03b | **PASS** |
| SRS-ALM-01 | ASYSTOLE (HIGH) shall be raised when no QRS is detected for 4.0 s; time to alarm from the last QRS shall be >= 4.0 s and <= 4.5 s. | TC-ASY-01, TC-ASY-01b, TC-ASY-01c, TC-RT-04, TC-RT-04b, TC-RT-04c | **PASS** |
| SRS-ALM-02 | TACHY (HIGH) shall be raised when the 4-beat mean heart rate exceeds 150 bpm, within 10 s of the rate change. | TC-RT-04, TC-RT-04b, TC-RT-04c, TC-TAC-01, TC-TAC-01b | **PASS** |
| SRS-ALM-03 | BRADY (HIGH) shall be raised when the 4-beat mean heart rate falls below 40 bpm, within 10 s of the rate change. | TC-BRA-01, TC-BRA-01b, TC-RT-04, TC-RT-04b, TC-RT-04c | **PASS** |
| SRS-ALM-04 | DESAT (HIGH) shall be raised when SpO2 < 85 % persists for 5 s, within 10 s of a step desaturation. | TC-RT-04, TC-RT-04b, TC-RT-04c, TC-SPO-02, TC-SPO-02b, TC-SPO-02c | **PASS** |
| SRS-ALM-05 | SPO2_LOW (MEDIUM) shall be raised when SpO2 < 90 % persists for 10 s, within 15 s of a step change. | TC-SPO-02, TC-SPO-02b, TC-SPO-02c | **PASS** |
| SRS-ALM-06 | LEADS_OFF (MEDIUM, technical) shall be raised after 0.5 s of front-end saturation; ECG rhythm alarms shall be inhibited while the ECG is invalid and the detector shall re-learn its thresholds after reconnection. | TC-ASY-03, TC-ASY-03b, TC-LO-01, TC-LO-01b, TC-RT-04, TC-RT-04b, TC-RT-04c | **PASS** |
| SRS-ALM-07 | SpO2 shall be declared invalid (no numeric, no SpO2 limit alarms) when the perfusion index is below the minimum, and LOW_PERF (LOW, technical) shall be raised within 10 s. | TC-PERF-01, TC-PERF-01b | **PASS** |
| SRS-ALM-08 | The annunciation output shall always reflect the highest priority among the active alarms (HIGH > MEDIUM > LOW). | TC-ALM-01 | **PASS** |
| SRS-ALM-09 | HIGH-priority physiological alarms shall latch until acknowledged; acknowledgement shall not dismiss an alarm whose condition is still present; a recurring condition shall be re-annunciated. | TC-ALM-02, TC-ALM-03, TC-ALM-04 | **PASS** |
| SRS-ALM-10 | No alarm shall be raised on normal sinus rhythm (50-139 bpm, SpO2 93-99 %) with nominal noise. | TC-FA-01 | **PASS** |
| SRS-SPO-01 | SpO2 shall be estimated within +/-2 % (absolute) over 70-100 %. | TC-SPO-01 | **PASS** |
| SRS-RT-01 | Alarm-system latency: the annunciation output shall be driven no later than 10 ms after the nominal release instant of the sample that makes the alarm condition true (hard deadline), including under CPU, memory and I/O overload from non-critical software. | TC-RT-01, TC-RT-01b | **FAIL** |
| SRS-RT-02 | Safety-related threads shall run under fixed-priority preemptive scheduling with the alarm manager at the highest priority; every mutex shared with a lower-priority thread shall use priority inheritance. | TC-INSP-01, TC-RT-01, TC-RT-01b | **FAIL** |
| SRS-RT-03 | The real-time phase shall perform no dynamic memory allocation and no file-system I/O, and process memory shall be locked in RAM. | TC-INSP-01 | **PASS** |
| SRS-RT-04 | Acquisition release jitter shall not exceed 1 ms. | TC-RT-02 | **FAIL** |
| SRS-WDG-01 | A stall of the acquisition or processing thread longer than 100 ms shall raise the HIGH-priority technical alarm SYS_STALL within 200 ms of the start of the stall, independently of the stalled thread. | TC-WDG-01, TC-WDG-01b | **PASS** |

## 4. Verification results (as executed)

| Test | Description | Metric | Value | Limit | Verdict |
|---|---|---|---|---|---|
| TC-ECG-01 | QRS sensitivity, 40-180 bpm, R 0.55-4.4 mV (96 records, 3.2 h) | Se_pct | 99.874 | >= 99.5 | **PASS** |
| TC-ECG-01b | QRS positive predictivity, same records | PPV_pct | 100 | >= 99.5 | **PASS** |
| TC-ECG-02 | Heart-rate error (4-beat mean) vs programmed rate, worst case | err_pct | 4.6572 | <= 10 | **PASS** |
| TC-ASY-01 | Asystole: missed alarms in 180 episodes | missed | 0 | == 0 | **PASS** |
| TC-ASY-01b | Asystole: max time-to-alarm from last R wave | t_max_s | 4.132 | <= 4.5 | **PASS** |
| TC-ASY-01c | Asystole: min time-to-alarm (no premature alarm) | t_min_s | 4.088 | >= 4 | **PASS** |
| TC-ASY-02 | Asystole with 0.05 mV rms noise: missed/interrupted alarms in 20 episodes | missed | 0 | == 0 | **PASS** |
| TC-ASY-03 | Asystole present at lead reconnection: missed alarms in 20 episodes | missed | 0 | == 0 | **PASS** |
| TC-ASY-03b | Asystole present at lead reconnection: max time-to-alarm | t_max_s | 6.798 | <= 10 | **PASS** |
| TC-TAC-01 | Tachycardia (75 to 160/180/220 bpm): missed alarms in 60 episodes | missed | 0 | == 0 | **PASS** |
| TC-TAC-01b | Tachycardia: max time-to-alarm from rate change | t_max_s | 2.402 | <= 10 | **PASS** |
| TC-BRA-01 | Bradycardia (75 to 30/35 bpm): missed alarms in 40 episodes | missed | 0 | == 0 | **PASS** |
| TC-BRA-01b | Bradycardia: max time-to-alarm from rate change | t_max_s | 7.64 | <= 10 | **PASS** |
| TC-SPO-01 | SpO2 estimate error, 70-100 %, 50-150 bpm | err_max_pct | 0.6839 | <= 2 | **PASS** |
| TC-SPO-02 | Desaturation alarms: missed or wrong-level alarms in 40 episodes | errors | 0 | == 0 | **PASS** |
| TC-SPO-02b | DESAT (97 to 80 %): max time-to-alarm | t_max_s | 6.98 | <= 10 | **PASS** |
| TC-SPO-02c | SPO2_LOW (97 to 88 %): max time-to-alarm | t_max_s | 12.22 | <= 15 | **PASS** |
| TC-LO-01 | Leads-off: missed alarms, secondary false alarms or failed recovery (20 episodes) | errors | 0 | == 0 | **PASS** |
| TC-LO-01b | Leads-off: deviation from 0.5 s persistence time | dev_s | 0.002 | <= 0.004 | **PASS** |
| TC-PERF-01 | Low-perfusion alarm on pulse loss: errors in 10 episodes | errors | 0 | == 0 | **PASS** |
| TC-PERF-01b | Low-perfusion alarm: max time-to-alarm | t_max_s | 8.27 | <= 10 | **PASS** |
| TC-FA-01 | False alarms on 4 h of normal rhythm (50-139 bpm, 93-99 %) | false_alarms | 0 | == 0 | **PASS** |
| TC-ALM-01 | Priority arbitration: highest active priority drives the output | errors | 0 | == 0 | **PASS** |
| TC-ALM-02 | Latching: high-priority physiological alarm persists, non-latching clears | errors | 0 | == 0 | **PASS** |
| TC-ALM-03 | Acknowledge: releases ended alarms only; ongoing conditions cannot be dismissed | errors | 0 | == 0 | **PASS** |
| TC-ALM-04 | Re-annunciation of a recurring condition after release | errors | 0 | == 0 | **PASS** |
| TC-RT-01 | Alarm latency, proposed configuration under overload (7161 alarms, 3 runs): maximum | L_max_ms | 3.7299 | <= 10 | **PASS** |
| TC-RT-01b | Pipeline latency of every sample (release to end of processing), proposed configuration under overload (148500 samples): samples above 10 ms | late_samples | 1 | == 0 | **FAIL** |
| TC-RT-02 | Acquisition release jitter, proposed configuration under overload (148500 samples): maximum | J_max_ms | 8.7194 | <= 1 | **FAIL** |
| TC-RT-03 | Samples or alarm events lost, proposed configuration, all probe runs | lost | 0 | == 0 | **PASS** |
| TC-WDG-01 | Runaway processing thread injected under overload (5 runs): undetected stalls | undetected | 0 | == 0 | **PASS** |
| TC-WDG-01b | Runaway processing thread: maximum time from stall start to SYS_STALL annunciation | t_max_ms | 118.001 | <= 200 | **PASS** |
| TC-RT-04 | Clinical scenario in real time, proposed configuration (3 runs): runs with missing, extra or repeated alarms | bad_runs | 0 | == 0 | **PASS** |
| TC-RT-04b | Clinical scenario in real time: smallest margin to the time-to-alarm limits | margin_s | 0.006 | >= 0 | **PASS** |
| TC-RT-04c | Clinical scenario in real time: maximum alarm-system latency | L_max_ms | 0.1644 | <= 10 | **PASS** |
| TC-INSP-01 | Automated inspection: heap, stdio or sleep calls inside real-time thread bodies; PI on all shared mutexes; memory locking | violations | 0 | == 0 | **PASS** |

## 5. Traceability matrix

See `docs/traceability.csv` (hazard, software item, requirement, test, verdict).
