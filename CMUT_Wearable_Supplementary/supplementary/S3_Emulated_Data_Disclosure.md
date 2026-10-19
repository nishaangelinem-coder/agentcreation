# S3. Disclosure and correct use of the emulated (digital-twin) results

**Read this before using anything in `results/emulated_measurements/`.**

## What the emulated results are

The digital twin (`cmutsim/hardware_twin.py`) predicts what a bench
characterisation of the fabricated patch *would be expected to show*. It
starts from the validated acoustic model and adds the error sources that a
real build and a real measurement chain introduce:

| Layer | Error source | Value (1 sigma unless stated) | Basis |
|---|---|---|---|
| Fabrication | Vacuum-gap thickness | 3 % | Typical LPCVD sacrificial-layer / wafer-bond uniformity |
| Fabrication | Plate thickness | 3 % | Nitride deposition uniformity |
| Fabrication | Cell radius | 1 % | Lithographic CD control |
| Electronics | TX channel gain mismatch | 0.5 dB | Multi-channel pulser/DAC datasheets |
| Electronics | TX channel phase mismatch | 3 deg | Routing + driver skew |
| Electronics | Delay quantisation | 5 ns (200 MHz FPGA clock) | Design choice |
| Electronics | Dead elements | 1.5 % probability per element | Flex-interconnect yield assumption |
| Wearable | Body curvature | R = 60 mm (forearm/thigh class) | Anatomical |
| Wearable | Shape-sensing residual | 50 um RMS | Embedded strain-gauge assumption |
| Measurement | Needle hydrophone | 0.2 mm diameter, 8 % calibration uncertainty, 150 Pa RMS noise | Typical calibrated PVDF needle |
| Measurement | Scan stage | 10 um repeatability | Motorised 3-axis stage |
| Measurement | Thermocouple | 10 mK RMS; viscous artefact 30 % of the true rise | Literature on thermocouple artefacts in ultrasound fields |
| Control | Pressure sensor noise / quantisation | 0.3 kPa RMS / 0.05 kPa | Echo-based estimator assumption |

All random draws use a fixed seed (`SEED = 20261019`) so every number is
exactly reproducible.

## What they are NOT

They are **not experimental data**. No device has been fabricated and no
measurement has been made. Every figure in that folder carries the red
watermark *"EMULATED (digital twin) - not experimental data"*, and every
table caption contains "(EMULATED)".

## How to use them in the manuscript (acceptable wording)

* *"A Monte Carlo digital twin that includes fabrication spread, channel
  mismatch, delay quantisation and measurement uncertainty predicts that a
  hydrophone characterisation of the fabricated patch would report a focal
  pressure of 12.4 +/- 0.8 kPa (n = 5 emulated scans) (Supplementary Table
  E3)."*
* *"The twin identifies patch curvature as the dominant hardware risk: on a
  60 mm radius body surface without shape correction the pressure at the
  target falls from 12 to about 3 kPa (Supplementary Table E2)."*
* Section heading: **"Hardware-in-the-loop prediction (digital twin)"**, not
  "Experimental results".

## Wording that must not be used

* "measured", "experimentally obtained", "fabricated prototype showed",
  "hydrophone measurements confirmed" — for any number from this folder.
* Removing the watermark, or moving emulated tables into a section titled
  "Experimental validation".

Reporting simulated values as hardware measurements is data fabrication
under COPE and ICMJE definitions. It is grounds for rejection or retraction
and, in India, for action under the UGC (Promotion of Academic Integrity
and Prevention of Plagiarism in Higher Educational Institutions)
Regulations. The emulated results are valuable as they are: they turn a
pure simulation paper into one that quantifies hardware risk and states
pass/fail criteria for the Phase II bench study.

## Replacing emulated data with real data

`S4_Phase_II_Bench_Protocol.md` lists every measurement, its instrument,
and its acceptance criterion. The CSV files in `emulated_measurements/`
use the same column layout that the bench data should be saved in. When
real data exist, store them in a new folder `results/measured/` with
the raw instrument files, and plot them against the model and the twin
prediction. Agreement within the twin's predicted spread then counts as a
genuine validation result.
