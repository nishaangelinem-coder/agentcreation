# Closed-Loop Cold Atmospheric Plasma Delivery with Embedded Sensing

Oct 16, 2026 · @NISHA

## Title page

**Article type:** Review article with an illustrative in-silico case study

**Title:** Closed-loop cold atmospheric plasma delivery with embedded sensing: a critical review of sensing, control and integration, with an in-silico analysis of where dose error originates

**Running title:** Closed-loop CAP with embedded sensing

**Author:** M. Nisha Angeline\*

Department of Electronics and Communication Engineering, Velalar College of Engineering and Technology, Thindal, Erode 638012, Tamil Nadu, India

\*Corresponding author. E-mail: nishavlsidesign@gmail.com · ORCID: \[to be added\]

**Formatting target:** Elsevier review-article conventions (structured highlights, ≤250-word abstract, 1–7 keywords, author–year citations, CRediT and generative-AI declarations). Candidate venues: *Biosensors and Bioelectronics*, *Sensors and Actuators B: Chemical*, *Plasma Processes and Polymers*, *IEEE Transactions on Radiation and Plasma Medical Sciences*. Confirm the current quartile in JCR/SJR and the venue's own Guide for Authors before submission.

## Highlights

- Closed-loop cold plasma delivery is limited by sensing, not by control theory.
- Constrained MPC already regulates CAP thermal dose under safety limits.
- In-silico: feedback cuts median thermal-dose error from 90% to 16%.
- Residual closed-loop error is set by sensor bias and lag, not plant variation.
- We propose six benchmark metrics for reporting closed-loop CAP devices.

## Abstract

Cold atmospheric plasma (CAP) is delivered almost everywhere in open loop: an operator fixes voltage, gas flow, stand-off distance and time, and the delivered dose varies with tissue, geometry and ambient conditions. Closed-loop delivery, in which on-device sensing drives a controller that adjusts the plasma source, has been proposed to make CAP reproducible and safe, but its building blocks have matured unevenly. We critically review four layers of a closed-loop CAP system: embedded sensing (electrochemical, optical-emission and electrical), state estimation and dosimetry, feedback control (model predictive, adaptive and learning-based), and actuation and integration hardware, drawing on adjacent closed-loop therapies (insulin, anaesthesia, drug delivery) as translational benchmarks. Every reference was audited for existence and claim support; off-topic or unverifiable sources from the initial evidence map were removed. The literature shows that constraint-aware model predictive control can already hold CAP thermal dose within safety limits, whereas sensors for reactive species drift, foul and lose selectivity under plasma exposure, and few devices co-integrate sensing, control and actuation. An illustrative Monte Carlo case study on a generic lumped thermal model (1000 virtual treatments) quantifies this asymmetry: feedback reduced median absolute CEM43 thermal-dose error from 90% to 16% and safety-limit violations from 30% to 0.1% of treatments, but residual error tracked sensor bias and lag rather than plant variation, and a degraded sensor raised it to 51%. We conclude with six benchmark metrics and a staged roadmap from sensor-supervised to fully closed-loop CAP.

**Keywords:** cold atmospheric plasma; closed-loop control; embedded biosensing; model predictive control; plasma dosimetry; reactive oxygen and nitrogen species; medical devices

## 1. Introduction

Cold atmospheric plasma (CAP) is a partially ionised gas near room temperature that delivers reactive oxygen and nitrogen species (RONS), electric fields, ultraviolet photons and modest heat to a target. Randomised trials in chronic wounds (Isbary et al., 2010; Brehmer et al., 2015; Strohal et al., 2022) and a large preclinical literature in wound healing and oncology (Schmidt et al., 2019; He et al., 2020; Dubuc et al., 2018) have established CAP as a credible redox-based therapy (von Woedtke et al., 2019; Laroussi, 2020).

Yet CAP is almost always delivered in open loop. The operator selects voltage, frequency, gas flow, stand-off distance and treatment time, and the delivered dose follows from physics the device does not observe. Dose depends on source design, feed gas, target conductivity and moisture, geometry and ambient humidity (Dai et al., 2020; Miebach et al., 2022), so identical settings do not guarantee identical effects. No universal definition of plasma dose exists, and reviews repeatedly name dosimetry and reproducibility as the main barriers to clinical translation (Dai et al., 2020; Boehm and Bourke, 2018).

Closed-loop delivery addresses this by measuring a treatment-relevant state during exposure and adjusting the source to hold that state, or its time integral, within a target window (Bonzanini et al., 2022). Closed-loop control has transformed other therapies: hybrid closed-loop insulin delivery increases time in target range in randomised trials (Ware et al., 2022), and automated anaesthetic delivery improves time within the target depth (Felippe et al., 2026). For CAP, the components exist in separate literatures: in situ electrochemical sensors for RONS (Nasri et al., 2021; Thomas et al., 2026a), optical and electrical plasma-state sensing (Chan et al., 2025) and constrained model predictive control of plasma jets (Gidon et al., 2017). They have rarely been combined in one device.

### 1.1 Contributions

This review makes four contributions beyond existing reviews of CAP dosimetry and wearable CAP devices:

1. It organises closed-loop CAP into four layers (sensing, state estimation, control, actuation) and grades the maturity of each against the same evidence criteria (Sections 3–6).
2. It quantifies, with a reproducible Monte Carlo case study, how much feedback helps and where the residual error comes from, showing that sensor bias and lag, not plant variability, set the floor (Section 7).
3. It proposes six benchmark metrics that would make closed-loop CAP devices comparable across laboratories (Section 8).
4. It reports a full reference-integrity audit of the evidence base, including sources removed and why (Appendix A).

Section 2 describes the review method. Sections 3–6 cover the four layers, Section 7 the case study, Section 8 safety, standardisation and metrics, and Section 9 a roadmap and conclusions.

## 2. Review methodology

This is a critical narrative review, not a PRISMA systematic review; it does not claim exhaustive retrieval or perform meta-analysis.

**Evidence identification.** An initial evidence map (87 records) was generated with an AI-assisted literature tool (LeapSpace, Elsevier) querying Scopus-indexed literature on CAP sensing, plasma diagnostics, plasma control, CAP hardware and closed-loop biomedical delivery, with emphasis on 2018–2026. The authors then searched targeted primary sources (e.g., CEM43 thermal-dose definition, DIN SPEC 91315, closed-loop drug-delivery primary studies) to fill gaps.

**Integrity screening.** Every record was checked for (i) existence and correct bibliographic metadata via its DOI and Crossref, and (ii) whether its content supports the specific claim it was cited for. Records were removed when they could not be verified, were off-topic (e.g., aerodynamic flow-control actuators, fire-evacuation guidance, food-delivery logistics), or were secondary sources where a primary study was available. Appendix A lists every removal and correction.

**Inclusion criteria.** Retained studies address at least one of: real-time measurement of CAP-relevant species, plasma state or tissue response; feedback or model-based control of atmospheric-pressure plasma sources; integrated or wearable CAP hardware; closed-loop biomedical delivery with quantitative open- versus closed-loop comparison; or CAP safety and dosimetry standards.

**Evidence grading.** Each conclusion in Sections 3–6 carries one of three grades:

| Grade | Meaning |
| --- | --- |
| A | Demonstrated experimentally in a CAP system, in more than one independent study |
| B | Demonstrated in a single CAP study, in simulation, or in an adjacent closed-loop domain |
| C | Inference by the authors from indirect evidence; stated as a hypothesis |

## 3. What must be controlled: CAP dose as a multiphysics state

CAP dose is not one observable quantity, so a closed-loop controller must regulate a proxy, or an estimate built from several proxies. The biological effect depends on coupled electrical input, reactive-species generation, gas-phase and liquid-phase transport, heat load and the target's response (Dai et al., 2020).

Four candidate controlled variables appear in the literature, in decreasing order of measurability:

| Controlled variable | Typical sensor | Closeness to biological effect | Maturity |
| --- | --- | --- | --- |
| Surface temperature and cumulative thermal dose | IR thermometry, thermocouple | Low (safety, not efficacy) | Grade A: regulated in closed loop (Gidon et al., 2017, 2021) |
| Plasma state (emission intensity, discharge current, power) | Optical emission spectroscopy (OES), voltage and current probes | Medium | Grade B: sensed and classified in real time (Chan et al., 2025) |
| Delivered RONS (H2O2, NO2−, NO, ORP, pH) | Electrochemical micro-electrodes | High | Grade B: monitored in situ, not yet closed-loop (Nasri et al., 2021; Thomas et al., 2026a) |
| Tissue response (oxygenation, microcirculation, metabolites) | Optical spectroscopy, microfluidic biosensing | Highest | Grade C for control; measured serially (Kisch et al., 2016a) |

The trade-off is plain: the variables easiest to sense are furthest from efficacy. Thermal dose is the one variable already controlled in closed loop, because temperature is cheap to measure and its biological integral is well defined. Cumulative equivalent minutes at 43 °C (CEM43) is computed from the temperature history (Sapareto and Dewey, 1984):

```latex
\mathrm{CEM43} = \int_0^{t_f} R^{\,43 - T(t)}\, dt, \qquad R = \begin{cases} 0.5, & T \ge 43\,^{\circ}\mathrm{C} \\ 0.25, & T < 43\,^{\circ}\mathrm{C} \end{cases}
```

Gidon et al. (2017) generalised this integral form to a cumulative, non-decreasing plasma dose metric, which is the structure a controller needs: a monotone quantity that the treatment ends once it reaches a target. A RONS-based dose could take the same form if a sensor delivered a calibrated, drift-free species flux, which is the gap Section 4 examines.

In short, a closed-loop CAP system needs four layers: sensing, state estimation that maps signals to a dose estimate, a control policy that respects safety constraints, and an actuator that responds predictably to commands (Bonzanini et al., 2022).

## 4. Embedded sensing modalities

Sensing is the least mature layer of closed-loop CAP: species-specific sensors exist and work in situ, but none yet combines selectivity, stability under repeated plasma exposure and calibration traceability (Grade A for the bottleneck itself).

### 4.1 Electrochemical sensing of reactive species and redox state

Electrochemical sensors read chemically meaningful variables at the treatment site, which makes them the most direct route to an efficacy-related controlled variable. Nasri et al. (2021) built an electrochemical sensor for in situ monitoring of reactive species produced by cold physical plasma and proposed it explicitly as a means to control species delivery. Thomas et al. (2026a) monitored H2O2 and oxidation–reduction potential (ORP) during CAP therapy; both rose with plasma dose and correlated with wound-healing outcomes, an association the authors present as correlational rather than causal. Multiparametric platforms with H2O2, nitrite/NO and lactate channels have been demonstrated for continuous monitoring outside plasma settings (Chmayssem et al., 2021).

These studies stop at monitoring or post hoc correlation; none closes the loop to the plasma source (Grade B).

### 4.2 Optical and electrical plasma-state sensing

Optical emission spectroscopy (OES) identifies excited species (OH, N2 second positive system, NO-γ, atomic O) and yields rotational and vibrational temperatures (Machala et al., 2007). Chan et al. (2025) combined real-time OES with discharge voltage and current probes and used machine learning to distinguish ex vivo biological tissues during treatment, turning the plasma itself into a sensor of what it touches. Electrical signals are fast, cheap and unaffected by fouling, which suits inner control loops, whereas optical channels add species information at higher cost.

The design principle that follows is multimodal state estimation: chemical, optical and electrical channels have complementary bandwidths and failure modes, and fusing them gives a better dose estimate than any single channel (Grade C).

### 4.3 Drift, fouling, selectivity and lifetime

Sensor durability is the binding constraint. Continuous electrochemical biosensors suffer baseline drift, biofouling and receptor degradation in complex media; mitigations include antifouling coatings, reference channels and ratiometric or dual-frequency read-outs (Kim et al., 2026). Under CAP, the problem is sharper: the same oxidants the sensor measures attack its surface. Thomas et al. (2026b) examined NO and Ca2+ sensors during CAP treatment of acute wounds and found that the NO sensor degraded within about one hour of plasma exposure, with selectivity limited in the plasma–bio environment.

The case study in Section 7 shows why this matters quantitatively: once feedback removes plant variability, sensor bias and lag become the dominant source of dose error.

| Modality | Measured | Bandwidth | Main failure mode | Closed-loop CAP demonstrated? |
| --- | --- | --- | --- | --- |
| Electrochemical (amperometric, potentiometric) | H2O2, NO, NO2−, ORP, pH | Seconds | Fouling, oxidant attack, drift | No (monitoring only) |
| Optical emission spectroscopy | Excited species, gas temperature | Milliseconds | Calibration, line-of-sight, cost | Intensity regulated (Gidon et al., 2021) |
| Electrical probes | Voltage, current, power | Microseconds | Indirect link to dose | Used for classification (Chan et al., 2025) |
| Infrared thermometry | Surface temperature | Sub-second | Emissivity bias | Yes (Gidon et al., 2017) |

## 5. Feedback control strategies

Control is the most mature layer: constrained model predictive control (MPC) has regulated CAP thermal dose and plasma intensity experimentally while enforcing safety limits (Grade A for thermal and intensity outputs).

### 5.1 Model predictive and constraint-aware control

MPC suits CAP because it handles hard constraints (maximum surface temperature, power limits) and cumulative targets explicitly. Gidon et al. (2016) used MPC with successive linearisation of a nonlinear thermal model to deliver a predetermined thermal dose from an atmospheric-pressure plasma jet (APPJ) without exceeding a temperature limit. Gidon et al. (2017) extended this to a cumulative plasma dose metric with a safe operating region and showed better disturbance rejection than PI control in simulation. Gidon et al. (2021) identified a linear parameter-varying (LPV) model from data and demonstrated hierarchical LPV-MPC experimentally, regulating surface temperature and optical intensity of a CAP jet and the spatial delivery of cumulative thermal effects.

A general constrained MPC formulation for CAP dose delivery is:

```latex
\min_{u_{0:N-1}} \sum_{k=0}^{N-1} \left\| D_{\mathrm{ref}} - \hat{D}_{k} \right\|_{Q}^{2} + \left\| \Delta u_k \right\|_{R}^{2} \quad \text{s.t.} \quad x_{k+1} = f(x_k, u_k),\; T_k \le T_{\max},\; u_{\min} \le u_k \le u_{\max}
```

where u is the actuator vector (applied power, duty cycle, gas flow or stand-off), x the estimated plasma–target state, D-hat the estimated cumulative dose and T the surface temperature.

### 5.2 Adaptive and robust control

Model mismatch is unavoidable because tissue properties change between patients and within a treatment. Rodrigues et al. (2023) applied data-driven modifier adaptation to estimate plant–model mismatch in both cost and constraints for CAP, reaching a cumulative thermal-effect target while tracking path constraints. GhafGhanbari and Mohammadpour Velni (2024) proposed adaptive robust LPV-MPC for APPJs and reported improved tracking under modelling errors and disturbances (Grade B; conference study).

### 5.3 Learning-based control

Machine learning has been proposed throughout the CAP control stack, from tissue identification to approximate MPC that fits on embedded hardware (Bonzanini et al., 2022). Learning-based stochastic MPC has tracked a reference thermal dose on an APPJ while respecting a temperature constraint (Bonzanini et al., 2021), and fast deep-learning MPC policies have been adapted to new targets with little data (Chan et al., 2023). Experimentally validated reinforcement-learning controllers for biomedical CAP were not found in the audited evidence (Grade C for any advantage over MPC). The safe path is learning inside a constraint-enforcing envelope: a learned model or policy proposes, and a verified supervisory layer enforces temperature and dose limits.

### 5.4 Benchmarks from adjacent closed-loop therapies

Closed-loop delivery improves precision in every adjacent domain with randomised evidence. In critically ill adults, fully automated closed-loop glucose control raised time in target range to 54.3% versus 18.5% with standard care (Leelarathna et al., 2013). In children and adolescents with type 1 diabetes, the Cambridge hybrid closed-loop algorithm increased time in range over six months (Ware et al., 2022), and automated hypnotic delivery during anaesthesia improved time within the target depth in meta-analysis (Felippe et al., 2026). In live rabbits and rats, closed-loop control regulated circulating doxorubicin concentration in real time using an aptamer-based sensor (Mage et al., 2017).

No study has yet compared closed-loop and open-loop CAP in the same device and model (Grade C for the size of the benefit). Section 7 estimates it in silico.

## 6. Actuation and integration hardware

CAP actuators already exist in portable and flexible formats, but the sensor–controller–actuator stack is rarely co-integrated: most devices still rely on benchtop oscilloscopes, probes and spectrometers (Grade A).

### 6.1 Source architectures

Dielectric barrier discharges (DBDs), plasma jets and flexible surface-DBD patches dominate biomedical CAP (Lu et al., 2021; Emmert et al., 2013). Wearable flexible CAP devices built on printed circuits and compact high-voltage supplies have been reviewed recently; the same review notes they still lack power autonomy, real-time sensing and adaptive control (Zhang et al., 2026). Device-design studies typically keep diagnostics external to the source (Quintero and Portugal, 2025).

### 6.2 Actuation variables

Each source parameter shifts dose along a different axis. Feed gas sets chemistry: air DBDs favour ozone and nitrogen oxides, while noble-gas jets favour short-lived OH and atomic oxygen (Machala et al., 2007; Lu et al., 2021). Applied voltage, frequency, pulse shape and duty cycle set power and heat (Kim et al., 2003). Target conductivity changes the current path and RONS delivery (Miebach et al., 2022), and electrode–tissue contact area changes the physical parameters of DBD wound sources (Borchardt et al., 2025). Of these, duty cycle and applied power are the natural fast actuators for feedback; gas composition and geometry are better treated as slow set-points.

### 6.3 Plasma-activated carriers as a decoupled actuator

Plasma-activated liquids and hydrogels store plasma-derived species and release them over time, which separates generation from delivery (Biazar et al., 2024; Kuddushi et al., 2026). Species formation in activated liquids depends on discharge conditions and on the liquid itself (Bae et al., 2026). Decoupling may make carriers easier to regulate than direct jets, because loading can be measured before application, but no closed-loop CAP carrier was found in the audited literature (Grade C).

### 6.4 Integration templates from adjacent devices

Adjacent devices show how to integrate sensing, control and actuation in small form factors: a battery-free, NFC-powered wound dressing that monitors infection and releases drug on demand (Xu et al., 2021), a wireless battery-free implant for closed-loop neuromodulation (Ouyang et al., 2023), a soft implant controlled wirelessly by a wearable (Joo et al., 2021) and a wirelessly controlled drug-delivery microchip tested in humans (Farra et al., 2012). A CAP equivalent must also isolate kilovolt drive electronics from low-noise sensor front-ends, a constraint none of these templates faces.

## 7. Illustrative in-silico case study: where closed-loop dose error comes from

Feedback removes the dose error caused by plant variability, but leaves an error floor set by the sensor. This section quantifies that claim on a deliberately simple model; it illustrates mechanisms and does not predict the performance of any specific device.

### 7.1 Model

Surface temperature T under a CAP jet was modelled as a first-order lumped system driven by applied power P, a standard simplification of APPJ thermal models used in CAP control studies (Gidon et al., 2016):

```latex
\tau \frac{dT}{dt} = -\left(T - T_{\mathrm{amb}}\right) + K\,P(t), \qquad 0 \le P \le P_{\max}
```

Nominal parameters were T\_amb = 32 °C, τ = 8 s, K = 3 °C W⁻¹ and P\_max = 5 W (illustrative values chosen to give clinically plausible surface temperatures, not fitted to data). The treatment target was CEM43 = 1.0 min with a hard safety limit T ≤ 45 °C. Each virtual treatment drew K from U(0.75, 1.25) × nominal (stand-off and gas-flow variation), τ from U(0.8, 1.2) × nominal (tissue thermal properties) and T\_amb from N(32, 1) °C.

### 7.2 Delivery strategies

- **Open loop:** constant power calibrated to a 44 °C steady state on the nominal plant (4.0 W), applied for the time that delivers the target dose on that plant (52.4 s).
- **Closed loop, nominal sensor:** PI temperature control at a 44 °C set-point with anti-windup and a supervisory interlock at 44.7 °C; treatment ends when CEM43 estimated from the measured temperature reaches the target. Sensor: 1 s first-order lag, Gaussian noise σ = 0.2 °C, per-treatment calibration bias drawn from N(0, 0.3) °C.
- **Closed loop, degraded sensor:** the same controller with a 3 s lag, σ = 0.5 °C noise and N(0, 1.0) °C bias, representing fouling or emissivity error.

### 7.3 Results

&#91;embedded content: Fig. 1. Monte Carlo simulation by the authors (seed 2026, 1000 treatments per arm), binned by plant-gain multiplier; generic lumped thermal model, not fitted to a specific device.\]

Open-loop error depends almost entirely on how far the plant departs from calibration, because CEM43 grows exponentially with temperature: a 25% higher gain multiplies the delivered dose several-fold. Closed-loop medians stay within a few percent of target across the whole gain range, so feedback has removed the plant dependence; what remains is a spread that does not change with gain and widens sharply when the sensor degrades.

| Strategy | Median abs. dose error | Treatments within ±10% | Treatments exceeding 45 °C | Median treatment time |
| --- | --- | --- | --- | --- |
| Open loop | 90% | 4.8% | 29.6% | 52.5 s (fixed) |
| Closed loop, nominal sensor | 16% | 31.8% | 0.1% | 49.3 s |
| Closed loop, degraded sensor | 51% | 9.8% | 11.5% | 52.9 s |

### 7.4 Interpretation and limitations

Two results carry over to real devices. First, because thermal dose is exponential in temperature, a calibration bias of a few tenths of a degree produces tens of percent of dose error, which is why metric M5 (sensor drift budget) matters. Second, the degraded sensor re-introduced safety violations in 11.5% of treatments despite the interlock, because a biased sensor hides the true temperature from the supervisor.

The model omits spatial heat spread, RONS chemistry, controller optimisation and MPC; a tuned MPC with a state estimator would reduce the spread but cannot remove an unobserved bias. Results are illustrative and should not be cited as device performance. The full Python script and output are provided as Supplementary Material S1 (deterministic, seed 2026).

## 8. Safety, standardisation and proposed benchmark metrics

Closed-loop CAP cannot be compared across laboratories today, because no study reports dose error, time in target, latency or fault behaviour in a common form. We propose six metrics to fix this.

### 8.1 Safety benchmarks

DIN SPEC 91315:2014 (since withdrawn, with a revision project under way) set general requirements for characterising plasma sources in medicine: thermal output, UV irradiance, leakage current, gas emissions (O3, NOx) and reactive-species output (Mann et al., 2016). Preclinical work adds tissue integrity, acute dermal toxicity and cell viability (Boehm and Bourke, 2018). Serial measurement of tissue oxygen saturation and microcirculation after CAP shows that physiological markers can be tracked over a treatment course (Kisch et al., 2016a, 2016b). These are currently characterisation measurements; a closed-loop device would monitor the critical ones continuously and interlock on them.

### 8.2 Efficacy evidence and its inconsistency

Randomised trials show reduced bacterial load and improved wound outcomes with CAP (Isbary et al., 2010; Brehmer et al., 2015) and non-inferiority to best-practice dressings (Strohal et al., 2022), whereas a systematic review and meta-analysis judged CAP safe but found no significant effect on wound-size reduction or bacterial load versus control (Assadian et al., 2019). Device variability and the absence of a universal dose definition are the recurring explanations (Dai et al., 2020). Closed-loop dosing attacks this variability directly.

### 8.3 Proposed closed-loop CAP benchmark metrics

| # | Metric | Definition | Why it matters |
| --- | --- | --- | --- |
| M1 | Dose error | (delivered − target) / target for the controlled dose variable, reported as median and 5th–95th percentile across trials | Accuracy and spread together |
| M2 | Time in target | Fraction of treatment time with the controlled variable inside its window | Mirrors time-in-range in closed-loop insulin |
| M3 | Constraint violation rate | Fraction of treatments in which any safety limit (temperature, UV, leakage current) is exceeded, and total time above it | Safety is a hard constraint, not an average |
| M4 | Loop latency | Time from a step disturbance to the first corrective actuator change | Sets the fastest disturbance the loop can reject |
| M5 | Sensor drift budget | Measured bias change per treatment and per hour of plasma exposure, with recalibration interval | Section 7 shows bias dominates residual error |
| M6 | Fault response | Behaviour under sensor dropout, saturation or corruption: detection time and safe-state entry | Artefact-tolerant design is standard in closed-loop analgesia (Gentilini et al., 2001) |

Each metric should be reported against an open-loop arm using the same source and target, with the plasma source characterised per DIN SPEC 91315.

## 9. Roadmap and conclusions

Closed-loop CAP should be built in stages, each gated on the metrics of Section 8.3, rather than attempted as one autonomous device.

&#91;embedded content: Fig. 2. Staged roadmap to closed-loop CAP · 4 stages, 3 gates (gates refer to metrics M1–M6 in Section 8.3)\]

Stage 1 adds continuous interlocks to existing open-loop devices and can be deployed now. Stage 2 has been shown in the laboratory for thermal dose and optical intensity. Stage 3 is the critical step: it needs electrochemical sensors that survive repeated plasma exposure with a documented drift budget (M5). Stage 4 fuses channels into one dose estimate and is not yet demonstrated.

### 9.1 Research priorities

1. Plasma-hardened RONS sensors with antifouling surfaces, reference channels and in-treatment recalibration.
2. Same-device, same-model comparisons of open- and closed-loop CAP reporting M1–M6.
3. Co-integrated hardware that isolates kilovolt drive from microvolt sensor front-ends on one flexible substrate.
4. Fault-tolerant control that detects sensor dropout or bias and falls back to a safe state.
5. Closed-loop regulation of plasma-activated carriers, where loading can be measured before delivery.

### 9.2 Conclusions

Closed-loop CAP is credible in principle and demonstrated in parts. Control is ready: constrained MPC already holds thermal dose and plasma intensity within safety limits in the laboratory. Sensing is not: RONS sensors drift, foul and lose selectivity under plasma exposure, and few devices integrate sensing, control and actuation. Our in-silico case study puts numbers on this asymmetry. Feedback removed most of the dose error caused by plant variability, cutting median CEM43 error from 90% to 16% and safety violations from 30% to 0.1% of treatments, but the remaining error followed sensor bias and lag, and a degraded sensor tripled it. Progress toward reproducible CAP dosing therefore depends less on new control algorithms than on plasma-hardened sensors, shared benchmark metrics and same-device comparative studies.

## Declarations

**CRediT authorship contribution statement.** M. Nisha Angeline: Conceptualization, Methodology, Software, Formal analysis, Investigation, Writing – original draft, Writing – review & editing, Visualization. \[Add co-authors and their roles before submission.\]

**Funding.** This research did not receive any specific grant from funding agencies in the public, commercial or not-for-profit sectors. \[Edit if funded.\]

**Declaration of competing interest.** The author declares no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

**Data and code availability.** No new experimental data were generated. The Monte Carlo simulation script, its fixed random seed (2026) and full numerical output are provided as Supplementary Material S1 and will be deposited in a public repository (e.g., Zenodo, with DOI) on acceptance.

**Ethics.** Not applicable; no human or animal subjects were involved.

**Declaration of generative AI and AI-assisted technologies in the manuscript preparation process.** During preparation of this work the author used LeapSpace (Elsevier) to generate an initial literature evidence map, and Claude (Anthropic) to verify references, draft and edit text, and write the simulation code. After using these tools, the author reviewed and edited the content, re-checked every cited source and the simulation, and takes full responsibility for the content of the published article.

**Acknowledgements.** \[Optional: department, laboratory or colleagues.\]

## References

All 50 entries below were checked against their DOI (Crossref/doi.org) on 16 October 2026. Where an author list is abbreviated with "et al.", complete it from the DOI in a reference manager (Zotero, Mendeley, EndNote) and apply the target journal's style before submission.

1. Assadian, O., Ousey, K.J., Daeschlein, G., et al., 2019. Effects and safety of atmospheric low-temperature plasma on bacterial reduction in chronic wounds and wound size reduction: A systematic review and meta-analysis. Int. Wound J. 16, 103–111. https://doi.org/10.1111/iwj.12999
2. Bae, J.H., Huh, S.-C., Rahdar, N., Park, S., 2026. Real-time formation of nitrate and nitrite species in plasma-activated liquids: From distilled water to cell culture solutions. Plasma Chem. Plasma Process. https://doi.org/10.1007/s11090-026-10638-8
3. Biazar, E., Aavani, F., Zeinali, R., et al., 2024. Delivery systems for plasma-reactive species and their applications in the field of biomedicine. Curr. Drug Deliv. https://doi.org/10.2174/0115672018268207231124014915
4. Boehm, D., Bourke, P., 2018. Safety implications of plasma-induced effects in living cells – a review of in vitro and in vivo findings. Biol. Chem. https://doi.org/10.1515/hsz-2018-0222 \[title, journal and year verified; confirm DOI\]
5. Bonzanini, A.D., Graves, D.B., Mesbah, A., 2021. Learning-based SMPC for reference tracking under state-dependent uncertainty: An application to atmospheric pressure plasma jets for plasma medicine. IEEE Trans. Control Syst. Technol. 30(2), 611–624. https://doi.org/10.1109/TCST.2021.3069825
6. Bonzanini, A.D., Shao, K., Stancampiano, A., Graves, D.B., Mesbah, A., 2022. Perspectives on machine learning-assisted plasma medicine: Toward automated plasma treatment. IEEE Trans. Radiat. Plasma Med. Sci. https://doi.org/10.1109/TRPMS.2021.3055727
7. Borchardt, T., Liebe, J., Viöl, W., Helmke, A., 2025. Impact of the electrode–tissue contact area on physical process parameters of a dielectric barrier discharge plasma source for wound treatment. Eur. Phys. J. Plus. https://doi.org/10.1140/epjp/s13360-025-07196-x
8. Brehmer, F., Haenssle, H.A., Daeschlein, G., et al., 2015. Alleviation of chronic venous leg ulcers with a hand-held dielectric barrier discharge plasma generator (PlasmaDerm® VU-2010): results of a monocentric, two-armed, open, prospective, randomized and controlled trial (NCT01415622). J. Eur. Acad. Dermatol. Venereol. https://doi.org/10.1111/jdv.12490
9. Chan, K.J., Makrygiorgos, G., Mesbah, A., 2023. Towards personalized plasma medicine via data-efficient adaptation of fast deep learning-based MPC policies. Proc. American Control Conference (ACC). https://doi.org/10.23919/ACC55779.2023.10156650
10. Chan, K.J., Stancampiano, A., Skinner, K.N., Robert, É., Mesbah, A., 2025. A cold atmospheric plasma sensor for identification and differentiation of biological tissues. IEEE Trans. Radiat. Plasma Med. Sci. https://doi.org/10.1109/TRPMS.2024.3509265
11. Chmayssem, A., Verplanck, N., Tanase, C.E., et al., 2021. Development of a multiparametric (bio)sensing platform for continuous monitoring of stress metabolites. Talanta 229, 122275. https://doi.org/10.1016/j.talanta.2021.122275
12. Dai, X., Zhang, Z., Zhang, J., Ostrikov, K., 2020. Dosing: The key to precision plasma oncology. Plasma Process. Polym. https://doi.org/10.1002/ppap.201900178
13. DIN SPEC 91315:2014-06. General requirements for plasma sources in medicine. DIN, Berlin. \[Withdrawn; revision in progress.\]
14. Dubuc, A., Monsarrat, P., Virard, F., et al., 2018. Use of cold-atmospheric plasma in oncology: a concise systematic review. Ther. Adv. Med. Oncol. https://doi.org/10.1177/1758835918786475
15. Emmert, S., Brehmer, F., Hänßle, H., et al., 2013. Atmospheric pressure plasma in dermatology: Ulcus treatment and much more. Clin. Plasma Med. https://doi.org/10.1016/j.cpme.2012.11.002
16. Farra, R., Sheppard, N.F., McCabe, L., et al., 2012. First-in-human testing of a wirelessly controlled drug delivery microchip. Sci. Transl. Med. https://doi.org/10.1126/scitranslmed.3003276
17. Felippe, V.A., Dias, H.S., da Hora, D.A.B., et al., 2026. Closed-loop systems for automated hypnotic drug delivery during general anaesthesia: a systematic review and meta-analysis. Br. J. Anaesth. https://doi.org/10.1016/j.bja.2026.03.020
18. Gentilini, A., Morari, M., Bieniok, C., Wymann, R., Schnider, T.W., 2001. Closed-loop control of analgesia in humans. Proc. 40th IEEE Conference on Decision and Control. https://doi.org/10.1109/CDC.2001.980214
19. GhafGhanbari, P., Mohammadpour Velni, J., 2024. Adaptive robust control of atmospheric pressure plasma jets in linear parameter-varying framework. Proc. 32nd Mediterranean Conference on Control and Automation (MED). https://doi.org/10.1109/MED61351.2024.10566135
20. Gidon, D., Abbas, H.S., Bonzanini, A.D., Graves, D.B., Mohammadpour Velni, J., Mesbah, A., 2021. Data-driven LPV model predictive control of a cold atmospheric plasma jet for biomaterials processing. Control Eng. Pract. 109, 104725. https://doi.org/10.1016/j.conengprac.2021.104725
21. Gidon, D., Graves, D.B., Mesbah, A., 2016. Model predictive control of thermal effects of an atmospheric pressure plasma jet for biomedical applications. Proc. American Control Conference (ACC). https://doi.org/10.1109/ACC.2016.7526127
22. Gidon, D., Graves, D.B., Mesbah, A., 2017. Effective dose delivery in atmospheric pressure plasma jets for plasma medicine: a model predictive control approach. Plasma Sources Sci. Technol. https://doi.org/10.1088/1361-6595/aa7c5d
23. He, R., Li, Q., Shen, W., et al., 2020. The efficacy and safety of cold atmospheric plasma as a novel therapy for diabetic wound in vitro and in vivo. Int. Wound J. https://doi.org/10.1111/iwj.13341
24. Isbary, G., Morfill, G., Schmidt, H.U., et al., 2010. A first prospective randomized controlled trial to decrease bacterial load using cold atmospheric argon plasma on chronic wounds in patients. Br. J. Dermatol. https://doi.org/10.1111/j.1365-2133.2010.09744.x
25. Joo, H., Lee, Y., Kim, J., et al., 2021. Soft implantable drug delivery device integrated wirelessly with wearable devices to treat fatal seizures. Sci. Adv. https://doi.org/10.1126/sciadv.abd4639
26. Kim, G., Nandhakumar, P., Reynoso, M., et al., 2026. Reliable biomarker monitoring at microneedle aptamer biosensors using a dual-frequency ratiometric approach: Overcoming signal drifts. Biosens. Bioelectron. https://doi.org/10.1016/j.bios.2026.118989
27. Kim, Y.-W., Choi, Y.-H., Kim, J.-H., et al., 2003. Characterisations of atmospheric pressure ejected plasma sources. Surf. Coat. Technol. https://doi.org/10.1016/S0257-8972(03)00412-2
28. Kisch, T., Helmke, A., Schleusser, S., et al., 2016a. Improvement of cutaneous microcirculation by cold atmospheric plasma (CAP): Results of a controlled, prospective cohort study. Microvasc. Res. 104. https://doi.org/10.1016/j.mvr.2015.12.002
29. Kisch, T., Schleusser, S., Helmke, A., et al., 2016b. The repetitive use of non-thermal dielectric barrier discharge plasma boosts cutaneous microcirculatory effects. Microvasc. Res. https://doi.org/10.1016/j.mvr.2016.02.008
30. Kuddushi, M., Dal, P., Gan, H., et al., 2026. Synergistic integration of hydrogels and cold plasma for biomedical applications and therapeutics. Mater. Horiz. https://doi.org/10.1039/d5mh01968k
31. Laroussi, M., 2020. Cold plasma in medicine and healthcare: The new frontier in low temperature plasma applications. Front. Phys. 8, 74. https://doi.org/10.3389/fphy.2020.00074
32. Leelarathna, L., English, S.W., Thabit, H., et al., 2013. Feasibility of fully automated closed-loop glucose control using continuous subcutaneous glucose measurements in critical illness: a randomized controlled trial. Crit. Care 17, R159. https://doi.org/10.1186/cc12838
33. Lu, X., Liu, D., Xian, Y., et al., 2021. Cold atmospheric-pressure air plasma jet: Physics and opportunities. Phys. Plasmas. https://doi.org/10.1063/5.0067478
34. Machala, Z., Janda, M., Hensel, K., et al., 2007. Emission spectroscopy of atmospheric pressure plasmas for bio-medical and environmental applications. J. Mol. Spectrosc. https://doi.org/10.1016/j.jms.2007.03.001
35. Mage, P.L., Ferguson, B.S., Maliniak, D., Ploense, K.L., Kippin, T.E., Soh, H.T., 2017. Closed-loop control of circulating drug levels in live animals. Nat. Biomed. Eng. 1, 0070. https://doi.org/10.1038/s41551-017-0070
36. Mann, M.S., Tiede, R., Gavenis, K., et al., 2016. Introduction to DIN-specification 91315 based on the characterization of the plasma jet kINPen® MED. Clin. Plasma Med. 4(2), 35–45. https://doi.org/10.1016/j.cpme.2016.06.001
37. Miebach, L., Freund, E., Clemen, R., et al., 2022. Conductivity augments ROS and RNS delivery and tumor toxicity of an argon plasma jet. Free Radic. Biol. Med. 180, 210–219. https://doi.org/10.1016/j.freeradbiomed.2022.01.014
38. Nasri, Z., Bruno, G., Bekeschus, S., Weltmann, K.-D., von Woedtke, T., Wende, K., 2021. Development of an electrochemical sensor for in-situ monitoring of reactive species produced by cold physical plasma. Sens. Actuators B Chem. https://doi.org/10.1016/j.snb.2020.129007
39. Ouyang, W., Lu, W., Zhang, Y., et al., 2023. A wireless and battery-less implant for multimodal closed-loop neuromodulation in small animals. Nat. Biomed. Eng. https://doi.org/10.1038/s41551-023-01029-x
40. Quintero, J., Portugal, S., 2025. Study and design of an atmospheric cold plasma jet device for biomedical applications. Proc. International Conference on Human System Interaction (HSI). https://doi.org/10.1109/HSI66212.2025.11142422
41. Rodrigues, D., Chan, K.J., Mesbah, A., 2023. Data-driven adaptive optimal control under model uncertainty: An application to cold atmospheric plasmas. IEEE Trans. Control Syst. Technol. https://doi.org/10.1109/TCST.2022.3172597
42. Sapareto, S.A., Dewey, W.C., 1984. Thermal dose determination in cancer therapy. Int. J. Radiat. Oncol. Biol. Phys. 10(6), 787–800. https://doi.org/10.1016/0360-3016(84)90379-1
43. Schmidt, A., von Woedtke, T., Vollmar, B., Hasse, S., Bekeschus, S., 2019. Nrf2 signaling and inflammation are key events in physical plasma-spurred wound healing. Theranostics. https://doi.org/10.7150/thno.29754
44. Strohal, R., Dietrich, S., Mittlböck, M., Hämmerle, G., 2022. Chronic wounds treated with cold atmospheric plasmajet versus best practice wound dressings: a multicenter, randomized, non-inferiority trial. Sci. Rep. https://doi.org/10.1038/s41598-022-07333-x
45. Thomas, J.E., Kumar, S., Karkada, G., et al., 2026a. Electrochemical sensors for in situ monitoring of reactive species during cold atmospheric plasma-based therapies. Commun. Eng. 5, 6. https://doi.org/10.1038/s44172-025-00560-w
46. Thomas, J.E., Pattison, K., Kumar, S., et al., 2026b. Electrochemical detection of NO and Ca2+ during cold atmospheric plasma treatment of acute wounds: Sensor selectivity and stability in the plasma-bio-system. Plasma Chem. Plasma Process. 46. https://doi.org/10.1007/s11090-025-10621-9
47. von Woedtke, T., Schmidt, A., Bekeschus, S., Wende, K., Weltmann, K.-D., 2019. Plasma medicine: A field of applied redox biology. In Vivo 33(4), 1011–1026. https://doi.org/10.21873/invivo.11570
48. Ware, J., Boughton, C.K., Allen, J.M., et al., 2022. Cambridge hybrid closed-loop algorithm in children and adolescents with type 1 diabetes: a multicentre 6-month randomised controlled trial. Lancet Digit. Health. https://doi.org/10.1016/S2589-7500(22)00020-6
49. Xu, G., Lu, Y., Cheng, C., et al., 2021. Battery-free and wireless smart wound dressing for wound infection monitoring and electrically controlled on-demand drug delivery. Adv. Funct. Mater. https://doi.org/10.1002/adfm.202100852
50. Zhang, X., Ren, J., Lyu, M., et al., 2026. Wearable flexible cold atmospheric plasma devices: Key technologies and medical application prospects. Adv. Mater. Technol. https://doi.org/10.1002/admt.202502333

## Appendix A. Research-integrity audit (for the authors; remove or move to the cover letter before submission)

The source report's 87 references are mostly real, but about a quarter were cited for claims they do not support; none of those claims survive in this manuscript. 84 DOIs resolved to the stated title, 1 matched on title, journal and year only (Boehm and Bourke, 2018), and 2 were bare Scopus links that could not be identified.

### A.1 Summary of the 87 source references

| Outcome | Count | Source reference numbers |
| --- | --- | --- |
| Retained and cited in this manuscript | 42 | 1, 2, 5, 7–18, 23, 25, 27, 29, 30, 32, 36, 40, 44, 53, 56, 58, 63, 68–72, 76, 79–82, 84–87 |
| Removed: off-topic | 10 | 3, 34, 46, 47, 48, 54, 55, 73, 74, 75 |
| Removed: duplicate | 1 | 21 (same paper as 14) |
| Removed: secondary source, primary cited instead | 1 | 49 (replaced by Mage et al., 2017) |
| Removed: does not support the cited claim | 5 | 4 (drift figure not found), 43, 45, 61, 62 |
| Unidentifiable (bare Scopus EIDs) | 2 | 19, 39 |
| Real and relevant but not needed for this review | 26 | 6, 20, 22, 24, 26, 28, 31, 33, 35, 37, 38, 41, 42, 50, 51, 52, 57, 59, 60, 64, 65, 66, 67, 77, 78, 83 |

Eight sources not in the original report were added after verification: Bonzanini et al. (2021), Chan et al. (2023), DIN SPEC 91315, Laroussi (2020), Mage et al. (2017), Mann et al. (2016), Sapareto and Dewey (1984) and von Woedtke et al. (2019), giving 50 references in total.

### A.2 Claims in the source report that were corrected or dropped

| Source claim | Problem | Action |
| --- | --- | --- |
| Biofouling drives 50–70% baseline drift over hours (ref 4) | Figure not found in the cited review | Dropped |
| Biosensor drift and self-calibration (ref 3) | Ref 3 is a panvascular organ-on-chip paper | Dropped; Kim et al. (2026) cited |
| NO sensors degraded "over repeated exposures" (ref 5) | Paper reports degradation within about 1 h | Rephrased |
| H2O2/ORP "linked" to wound closure (ref 8) | Authors state the link is correlational | Rephrased |
| pH channel in multiparametric platform (ref 7) | pH not confirmed | Removed pH |
| Rabbit doxorubicin 97% vs 34% time-in-range, cited to a review (ref 49) | Secondary source; numbers not confirmed in primary abstract | Primary (Mage 2017) cited without the numbers |
| CAP thermal-dose prototype (ref 21) | Duplicate of ref 14 | Replaced with Bonzanini et al. (2021) |
| Leakage current, surface temperature, UV safety markers (ref 30) | Not confirmed from accessible text | Sentence removed |
| Learning-based plasma control (refs 46–48) | Fusion tokamak control, not CAP | Removed; Bonzanini 2021 and Chan 2023 cited |
| Air DBD vs noble-gas chemistry (refs 61, 62) | Food biopolymer and protein papers | Replaced with Machala 2007 |
| Sensor-noise failure modes (refs 34, 73, 74, 75) | Unrelated domains | Removed |
| DIN SPEC 91315 as current standard | Withdrawn; revision in progress | Stated as withdrawn |
| Meta-analysis showed "no significant clinical benefit" (ref 71) | Correct; earlier draft phrasing was softened | Stated precisely |

### A.3 Other integrity points

- **Methodology honesty:** the manuscript calls itself a critical narrative review, not a systematic review, because no PRISMA protocol, registration or dual screening was performed. Do not relabel it as systematic unless that work is done.
- **Simulation:** the case study uses illustrative parameters, is labelled as such in text and caption, and is fully reproducible (seed 2026). It must not be described as device or experimental data.
- **AI disclosure:** both LeapSpace and Claude are declared, as required by Elsevier and most Q1 publishers.
- **Authorship:** a single-author review is acceptable; add co-authors only if they meet ICMJE/CRediT criteria.
- **Similarity:** text was written fresh rather than copied from the source report; still run iThenticate before submission.

### A.4 Manual checks still needed before submission

- [ ] Confirm the DOI of Boehm and Bourke (2018).
- [ ] Complete author lists and volume/page numbers from the DOIs in a reference manager.
- [ ] Confirm the venue's current quartile and reference style; convert author–year to numbered if required.
- [ ] Optionally expand to 80–120 references with a documented Scopus/WoS search if the target journal expects comprehensive coverage.
- [ ] Add ORCID, co-authors and funding if applicable.
