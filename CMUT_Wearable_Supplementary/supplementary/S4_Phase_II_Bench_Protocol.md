# S4. Phase II bench-validation protocol (to obtain real hardware data)

This protocol turns each model prediction into a measurement with a
pre-registered acceptance criterion. The acceptance bands come from the
digital-twin spread (Supplementary Tables E1–E5). Bench results that fall
inside a band validate the model; results outside it trigger the listed
root-cause check.

## Equipment

| Item | Minimum specification | Typical choice |
|---|---|---|
| Hydrophone | Calibrated needle, <=0.2 mm, 0.5–5 MHz, calibration traceable to a national metrology institute | PVDF needle hydrophone + submersible pre-amp |
| Scan tank | Degassed, deionised water (dissolved O2 < 3 mg/L), 3-axis stage, <=10 um repeatability | Motorised acoustic scan tank |
| Phantom | Tissue-mimicking hydrogel, c ~ 1540 m/s, alpha ~ 0.5 dB/cm/MHz; optional layered skin/fat/muscle | Agar–graphite or PVA cryogel |
| Thermometry | 0.1 mm fine-wire T-type thermocouple (or fibre-optic probe), 1 mK resolution after averaging | — |
| Impedance analyser | 0.1–10 MHz | Vector network analyser with impedance fixture |
| Oscilloscope | >=200 MHz, 12-bit | — |
| Radiation-force balance | Optional: total acoustic power, +/-5 % | — |

## Measurements and acceptance criteria

| # | Measurement | Method | Model value | Accept if | If it fails, check |
|---|---|---|---|---|---|
| M1 | Collapse voltage per element | C–V sweep, detect capacitance jump | 51.6 V (twin 51.6 +/- 3.2 V) | Mean within +/-10 %, SD < 5 V | Gap/insulator thickness (ellipsometry) |
| M2 | Electrical impedance | VNA, 0.2–3.2 MHz, at 0.8 Vc | Minimum near the immersed resonance | f_res within +/-15 % of the design | Plate radius/thickness; this tests the 1.5 MHz vs ~2.6 MHz cell issue (S2 item 1) |
| M3 | -6 dB fractional bandwidth | Pulse-echo off a flat reflector | 47.8 % | > 40 % | Backing/matching |
| M4 | Peak focal pressure in water | Hydrophone, 5 repeated scans | 12 kPa set point | 12.4 +/- 0.8 kPa (twin 95 % CI 11.3–13.4) | Drive calibration, dead elements |
| M5 | Lateral / axial FWHM (-6 dB **pressure**) | Raster scan, 0.05 mm lateral, 0.25 mm axial step | 0.94 / 5.50 mm | Within twin 95 % CI | Apodization setting, element delays |
| M6 | Steering +/-20 deg | Repeat M4/M5 at 5 angles | Loss <= 0.3 dB | Loss <= 0.6 dB | Delay table, element directivity |
| M7 | Side-lobe level | Lateral scan, dB scale | -11.5 dB (rect), -23.8 dB (Hamming) | Within +/-2 dB | Channel phase mismatch |
| M8 | Curvature test | Mount patch on 60 mm radius cylinder phantom, with/without shape correction | Target p: 3.3 kPa uncorrected, 11.5 kPa corrected | Corrected >= 10 kPa at the target | Strain-sensor calibration |
| M9 | Temperature rise | Thermocouple at the focus in phantom, at 50, 100 and 200 kPa, 100 % duty; extrapolate as p^2 to 12 kPa | 36 / 146 / 584 mK at 50 / 100 / 200 kPa | Fitted slope within +/-30 % | Absorption coefficient of phantom; subtract the viscous artefact using the early-time method |
| M10 | Closed-loop step and disturbance | Step 0→12 kPa; disturbance = 40 % coupling loss (spacer or partial gel removal) | 170 ms settling, 0 % overshoot | Settling < 300 ms, overshoot < 10 % | Estimator latency, plant gain |
| M11 | 4 h endurance | Patch on phantom, closed loop, log p, drive and skin-side temperature | >99 % of time in the 10–15 kPa window | >95 % in window, skin dT < 1 deg C | Coupling-gel drying, battery |
| M12 | Power consumption | Inline power meter, whole patch | ~0.45 W estimate | < 0.6 W | FPGA clocking, TX quiescent current |

At 12 kPa the focal temperature rise (~0.4 mK at 20 % duty) is **below the
detection limit** of any practical thermocouple (twin LOD ~9.5 mK, Table E4).
Do not report a "measured" rise at 12 kPa. Measure at higher pressure and
extrapolate with the p^2 law, as in row M9.

## Reporting

* Report mean +/- SD with n, and 95 % CI, for every quantity.
* Archive raw hydrophone waveforms (HDF5 or CSV) with calibration files.
* Plot measured data over the model curve and the twin's +/-1 SD band
  (use the `HW3_*.csv` layout).
* Report the hydrophone calibration uncertainty and spatial-averaging
  correction according to IEC 62127-1.
* Report MI and TI according to IEC 62359, using derated values.
