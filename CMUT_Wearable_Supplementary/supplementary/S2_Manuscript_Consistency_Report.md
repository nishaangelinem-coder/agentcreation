# S2. Manuscript consistency report (pre-submission audit)

This is an internal audit for the authors, not a file to upload as-is. It
lists every point where the independent reproduction (`run_all.py`)
disagrees with the manuscript, or where a Q1 reviewer is likely to object.
Items are ordered by severity.

**Overall:** most headline numbers reproduce well:

* Peak pressure, Isppa, Ispta and MI agree to within 2 %.
* Steering loss is −0.28 dB here vs −0.33 dB in the manuscript.
* Apodization side lobes agree to within 0.4 dB for all four windows.
* Focal ΔT is 0.42 mK here vs 0.37–0.40 mK.
* The PI loop gives 0 % overshoot and about 0.04 Pa steady-state error.

The main claims of feasibility stand. The items below need correcting or
clarifying before submission.

---

## A. Must fix (a reviewer will catch these)

### A1. The CMUT cell does not resonate at 1.5 MHz
With a = 30 µm and t = 1.5 µm (Si₃N₄, E = 250 GPa), the clamped-plate model
gives:
* vacuum resonance 7.2 MHz;
* immersed resonance 3.07 MHz at 0 V and 2.60 MHz at 0.8 V_c, with
  fluid-loading factor β = 4.3.

Placing the biased, immersed resonance at 1.5 MHz needs **a ≈ 37.7 µm**
(then V_c ≈ 33 V), or a thinner plate. A 1.5 MHz drive is still inside the
low-frequency skirt of a broadband CMUT, so the design is not invalid. But
the text "We size the cell for a transmit resonance near the array
operating frequency" is not supported by the stated dimensions.
**Action:** either change a/t in Table I or state that the array runs
below the cell resonance and give the sensitivity penalty. See Fig S1(c)
and Table S1.

### A2. Focal-spot widths are −3 dB (intensity) widths, labelled as FWHM
The reproduction gives:

| | −6 dB (pressure FWHM) | −3 dB (intensity FWHM) | Manuscript |
|---|---|---|---|
| Lateral | 0.94 mm | **0.69 mm** | 0.69 mm |
| Axial | 5.50 mm | **3.95 mm** | 3.74 mm |

The manuscript's 0.69 mm matches the −3 dB width exactly. The paraxial
estimate (≈1.2 λF/D = 0.77 mm) is not valid at F/0.63. An independent
point-source sum gives 0.92–0.96 mm (−6 dB).
**Action:** write "intensity FWHM (−3 dB)" in the text, Table III,
Table IV and Table XII, or report the pressure FWHM of 0.94 mm. The
< 2 mm targeting requirement is still met either way.

The same convention explains Table VI. The −3 dB main-lobe widths here are
Hann 1.02 mm (manuscript 1.03 mm) and Hamming 0.97 mm (manuscript 0.86 mm).

### A3. Nominal apodization is inconsistent
Table III says the nominal apodization is **Hamming**, but the 0.69 mm
width and the 12 kPa / 4.45 mW·cm⁻² figures correspond to **rectangular**
weighting. Hamming widens the −3 dB spot to 0.97 mm and costs **4.9 dB of
focal gain**, which the drive must supply. **Action:** state one nominal
window and quote its numbers consistently. Mention the gain penalty.

### A4. "Thermal index" is numerically equal to ΔT
Table VII lists TI = 0.0001 / 0.0004 / 0.0009 / 0.0019, which equals the ΔT
column. TI is a power ratio (IEC 62359), not a temperature. The IEC
formulae give:

| Duty cycle | TIS | TIB |
|---|---|---|
| 20 % | 6.6 × 10⁻⁴ | 5.7 × 10⁻³ |
| 100 % | 3.3 × 10⁻³ | 2.9 × 10⁻² |

These are still far below 1 (Table S7). **Action:** replace the TI
column with TIS/TIB, or drop it.

### A5. Grating-lobe table (Table V)
* At 0.50 λ pitch there is **no grating lobe** in visible space for ±20°
  steering (sin θ_g = sin 20° − 2 < −1). The −19.8 dB value must be a side
  lobe. The reproduction gives a −16.0 dB peak side lobe on the 20 mm arc.
* Table V shows a *lower* lobe at 0.75 λ (−22.2 dB) than at 0.50 λ, which
  is non-monotonic and will attract questions. The reproduction gives
  −13.0 dB at 0.75 λ (grating lobe at −82°), −12.0 dB at 1.0 λ (at −41°)
  and −5.0 dB at 1.5 λ (at −19°). The 1.0 λ and 1.5 λ values agree well
  with the manuscript.
* Table I says 0.49 λ but Table V says 0.50 λ = 0.51 mm. Use one value.

**Action:** relabel the column "peak side/grating lobe" and recompute.

### A6. Steering widths at ±10°
Table IV shows the focal width shrinking from 0.69 mm (0°) to 0.60 mm
(±10°), then returning to 0.69 mm (±20°). Physically the width grows
roughly as 1/cos θ. The reproduction gives 0.693 → 0.696 → 0.706 mm (−3 dB).
**Action:** recompute Table IV.

### A7. Closed-loop figure shows negative pressure
Fig 18 dips to about −40 kPa after the disturbance. A pressure amplitude
cannot be negative, and a coupling loss lowers the gain rather than
reversing it. With the PI gains as stated and a 40 % coupling loss, the
reproduction gives a dip of only 1.2 kPa and recovery in 181 ms (Fig S11).
**Action:**
* define the disturbance;
* clip the plant output at 0;
* add anti-windup;
* state the plant parameters (K, τ, θ), which are not given anywhere.

### A8. Figure/table numbering order
* Fig 22 appears before Fig 21.
* Table II appears after Table VI.
* The Fig 17 title has a typo: "20%% duty".
* Fig 5(a) shows a 12 × 12 2-D matrix, but the array is a 64-element
  linear array.

---

## B. Should fix (strengthens the paper)

### B1. Robustness of the PI loop
With ±50 % uncertainty in plant gain and lag, plus a 2–10 ms estimator
delay, the median settling time is 173 ms (95th percentile 244 ms). The
worst-case overshoot, however, reaches **56 %** (Table S10, Fig S11b). That
happens at high gain combined with a long delay. For a skin-coupled
wearable, where coupling varies, consider:
* gain scheduling from the echo amplitude; or
* lowering K_p to about 2.5 (trade-off: slower settling).

### B2. Patch curvature is the dominant hardware risk (new)
The digital twin (Table E2, Fig E2) shows the effect of a 60 mm radius
body surface with no shape correction:
* the focus moves **5.3 mm toward the skin**;
* the pressure at the 20 mm target falls from **12 to 3.3 kPa**, outside
  the 10–15 kPa window.

Shape-sensed delay correction (50 µm residual) restores 11.4 ± 0.2 kPa.
This belongs in the Discussion and in the system architecture (Fig 1). It
is a strong, specific contribution for a wearable-array paper.

### B3. Thermal verification at 12 kPa is not measurable
At 20 % duty the focal ΔT is about 0.4 mK. Even at 100 % duty it is about
2 mK, below the ~9.5 mK detection limit of a phantom thermocouple
(Table E4). Phase II must therefore measure at 50–200 kPa and extrapolate
with the p² law (protocol S4, row M9). Phrase the Phase II plan this way.

### B4. ">4 h on a wearable power budget" has no numbers
The surface pressure needed is only 2.8 kPa, which takes ≈0.51 V AC at a
41 V bias (5.4 kPa/V). The TX dynamic power is therefore negligible
(< 1 mW). The budget is dominated by the FPGA and channel quiescent power,
estimated at ≈0.45 W, so a 4 h session needs ≈1.8 Wh, or ≈600 mAh at
3.7 V (Table S12). Add this to support the claim.

### B5. Radiation-force displacement
The 3-D Kelvin-solution estimate gives on-axis static displacements of:

| Layer | Displacement |
|---|---|
| Skin | 0.002 nm |
| Fat | 0.041 nm |
| Muscle | 0.011 nm |
| Tumour | 0.005 nm |

Peak strains are 10⁻¹⁰–10⁻⁹. The manuscript gives 0.042 nm in the tumour.
Both agree on the conclusion that the stimulus is sub-nanometre and
non-destructive, but the method behind Table IX is not described.
**Action:** state the displacement model.

Note that the *time-averaged* body force at 12 kPa is only about
0.14 N/m³. Reviewers will ask whether such a small force can engage
mechanotransduction. The Limitations section already raises this. Consider
quantifying it against published LIPUS mechanobiology doses, which are
typically 30 mW/cm² SATA, i.e. about 30× the Ispta here.

### B6. Assumptions that must be stated
These quantities are absent from the manuscript and were assumed here:
* element elevation height (6 mm);
* cell pitch (66 µm);
* plate material constants;
* insulation layer (0.15 µm SiO₂);
* number of cells per element (540);
* element C₀ (68 pF);
* the FOPDT plant.

Add them to Table I or to a supplementary table. Without them the results
cannot be reproduced.

---

## C. References to verify before submission

* **[1] Ciucci et al., Science 392 (6796), eads9412, 2026** and
  **[2] Paltzer et al., same issue.** The central motivation rests on [1].
  Check that both entries (authors, article number, DOI) are correct
  against the publisher's site.
* [4] Puts et al. and [5] Liu et al.: check the year/volume against the DOI.
* For acoustic output reporting, add IEC 62127-1 (hydrophone measurement)
  and the AIUM/NEMA Output Display Standard.

---

## D. What the reproduction confirms

| Claim | Status |
|---|---|
| 12 kPa focus at 20 mm achievable | Confirmed (peak at 19.8 mm, homogeneous; 19.8–20.0 mm layered) |
| Steering ±20° with < 0.33 dB loss | Confirmed (−0.28 dB) |
| Targeting < 2 mm | Confirmed (0.94 mm −6 dB; focal error < 0.25 mm) |
| ΔT < 1 °C | Confirmed (0.42 mK, 20 % duty; 1.0 mK with fat properties, worst case) |
| MI < 0.5 | Confirmed (0.0098) |
| PI settling ≈150 ms, no overshoot | Confirmed for the nominal plant (170 ms, 0 %); not robust to +50 % gain combined with a long delay |
| Layered tissue compensable | Confirmed (+0.44 dB drive; focal shift ≤ 0.2 mm) |
