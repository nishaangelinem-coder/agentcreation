# Results

## A virtual instrument with exact ground truth

The simulator reproduces the qualitative structure of real data (Fig. 2). Simulated spectra show complementary b/y ladders with enhanced cleavage N-terminal to Pro, neutral-loss satellites and noise peaks interleaved with fragments (Fig. 2a). Training and test peptides have matched length distributions (Fig. 2b) but zero sequence overlap. LC-MS1 maps display isotope envelopes and elution profiles (Fig. 2c), and case–control difference maps place phenotype signal on both library and non-library ions (Fig. 2d). Site effects shift the total-ion chromatogram (Fig. 2e). Thirty percent of ions are not identifiable by design (Fig. 2f), and the abundance distribution is truncated by a limit of detection (Fig. 2g).

In the MS1 cohorts, ID-free alignment matched detected features to the correct consensus feature with ≥ 99.99 % accuracy. The fraction of missing values was 9.7 ± 2.7 % for identification-free quantification and 25.5 ± 1.2 % for DDA identification-first quantification. Table 1 summarises the corpora and Table 2 the model configurations.

![**Fig. 1. From identification-first to spectrum-native proteomics.** (a) The classical workflow loses unidentified spectra, variants, PTMs and non-canonical peptides at the search step, and creates missing values through stochastic DDA identification. (b) Spectrum-native models share an encoder that may be pretrained on unlabeled corpora and support database-free sequencing, identification-free phenotype prediction and representation learning.](../figures/fig1_paradigm.png){width=100%}

![**Fig. 2. The virtual instrument.** (a) Simulated MS2 spectrum of LVNEDTPAGLK with annotated b (blue) and y (red) ions; grey, noise and neutral-loss peaks. (b) Peptide-length distributions of training (proteome T) and test (proteome P) spectra. (c) Rendered LC-MS1 map of one run (64 RT × 256 m/z bins). (d) Mean case − control difference map; circles, library biomarkers; squares, unidentifiable biomarkers. (e) Total-ion chromatograms of the four acquisition sites. (f) MS1 feature landscape coloured by identifiability. (g) Abundance distribution and limit of detection (LOD).](../figures/fig2_virtual_instrument.png){width=100%}

{{TABLE1}}

**Table 2. Model and training configurations.**

| Item | SpecNovo | FeatureFormer |
|---|---|---|
| Input | ≤ 80 peaks (m/z, √I) + precursor (M, z) | 300 aligned MS1 features (value or missing) |
| Embedding | sinusoidal m/z (λ 10⁻³–10⁴ Da) + complementary m/z + intensity | feature ID + 0.1·sin(m/z, RT) + value × w_j (FT-style) |
| Encoder | 3 × pre-norm transformer, d = 128, 4 heads, FF 512 | 2 × pre-norm transformer, d = 64, 4 heads, FF 128 |
| Decoder / head | 3-layer autoregressive decoder + prefix/residual mass + ion-evidence MLP (8→32→1) | [CLS] logit + zero-initialised additive token evidence |
| Parameters | 1.44 M | 0.14 M |
| Objective | cross-entropy, label smoothing 0.05 | SSL: masked-intensity MSE (30 %); fine-tune: cross-entropy, smoothing 0.1 |
| Optimiser | AdamW, lr 5 × 10⁻⁴ (3 × 10⁻⁴ for evidence stage), cosine, batch 128 | AdamW, one-cycle, lr 10⁻³ (head) / 2 × 10⁻⁴ (body), batch 32, 150 steps |
| Inference | beam K = 5, precursor-mass constraint 0.05 Da | single forward pass |
| Hardware | 3 CPU threads; final model ≈ 1.45 h wall-clock (mass-aware + evidence stages) | 1 CPU thread per seed |

![**Fig. 3. Model architectures.** (a) SpecNovo encodes peaks and their complementary masses, decodes residues conditioned on prefix and residual mass, and adds a differentiable ion-evidence logit; decoding is mass-constrained. (b) FeatureFormer tokenises identification-free MS1 features, is pretrained by masked-intensity modelling on unlabeled runs and fine-tuned for phenotype; gradient × input yields feature-level attributions.](../figures/fig3_architecture.png){width=100%}

## Database-free sequencing recovers peptides that database search cannot

On canonical peptides, target–decoy database search was essentially perfect: recall 0.999 at 1 % FDR, with 0.03 % false identifications (Table 3, Fig. 4c). This is the regime for which identification-first proteomics was designed. Introducing a single amino-acid variant, or drawing peptides from a proteome absent from the database, reduced database-search recall to 0.000 (Table 3). Target–decoy competition correctly rejected the near-miss target matches, so no false identifications passed 1 % FDR; but every variant and novel peptide was lost.

Both de novo methods were, as expected, insensitive to database membership. SpecNovo recall was 0.416, 0.389 and 0.387, and spectrum-graph recall 0.475, 0.475 and 0.468, on the canonical, SAAV and novel sets respectively.

At full coverage the classical spectrum-graph algorithm had higher peptide recall and amino-acid precision (0.834 vs 0.687 on inDB). SpecNovo, however, produced better-ranked confidence scores, with precision–coverage AUC of 0.663 vs 0.630 (Fig. 4b). Its precision was higher whenever predictions were filtered by confidence:

| Coverage (inDB) | SpecNovo precision | Spectrum-graph precision |
|---|---|---|
| 10 % | 0.874 | 0.801 |
| 25 % | 0.787 | 0.716 |
| 50 % | 0.660 | 0.612 |

The same held for SAAV spectra, where precision at 10 % coverage was 0.905 vs 0.743. In practice, this is the regime that matters for database-free discovery, where only high-confidence de novo calls are followed up.

{{TABLE3}}

![**Fig. 4. De novo sequencing benchmark.** (a) Validation token accuracy (teacher forcing) as a function of training spectra for the base transformer, the mass-aware variant and the ion-evidence variant (initialised from the mass-aware checkpoint). (b) Peptide precision–coverage curves on inDB. (c) Peptide recall of database search, spectrum-graph DP and SpecNovo on canonical, SAAV and novel-proteome spectra; hatched, database false identifications at 1 % FDR. (d, e) Peptide recall versus fragment observation probability and chemical-noise peaks. (f) Ablation: amino-acid precision and peptide recall on inDB; MA, mass-aware; EV, ion evidence; +Nk, additional training spectra.](../figures/fig4_denovo.png){width=100%}

## Inductive biases, not capacity, determine data efficiency

The ablation (Table 5, Fig. 4a,f) is the most instructive result of this study. A standard Casanovo-style encoder–decoder reached a teacher-forced validation token accuracy of 0.578 after 448k training spectra. Under mass-constrained beam search, however, it achieved only 0.8 % peptide recall and 0.146 amino-acid precision. Inspection of greedy decodes showed that the model had learned the sequence prior (C-terminal K/R, residue frequencies) but very little of the spectrum: the decodes were plausible tryptic peptides that did not explain the observed peaks.

Complementary-ion encoding and prefix/residual-mass conditioning accelerated early learning (token accuracy 0.386 vs 0.303 at 38k spectra) but saturated at the same level (0.582). Adding the differentiable ion-evidence term broke this plateau within 1,200 further optimisation steps. Token accuracy rose to 0.674, 0.849 and 0.864 after 51k, 102k and 154k additional spectra, and peptide recall rose to 4.3 %, 35.9 % and 41.6 %. The total budget (410k spectra, ≈ 1.45 h of wall-clock time on three CPU threads) is two orders of magnitude below that used to train published de novo transformers (e.g., 30 M spectra for Casanovo [9]). The ablation therefore quantifies how much a physics-informed prior substitutes for data. It also explains why repository-scale corpora are a precondition for purely data-driven "proteomics GPTs".

{{TABLE5}}

## Robustness: where classical algorithms still win

Fragment completeness was the dominant factor (Table 4, Fig. 4d). With nearly complete ladders (*p*~obs~ = 0.95), SpecNovo outperformed the spectrum graph: peptide recall 0.621 vs 0.542, PR-AUC 0.840 vs 0.693. Below *p*~obs~ ≈ 0.9 the ranking reversed, and at *p*~obs~ = 0.5 SpecNovo recalled 9.2 % of peptides against 18.2 % for the spectrum graph. SpecNovo was also more sensitive to chemical noise (Fig. 4e). From λ = 10 to 200 noise peaks its recall fell from 0.438 to 0.248, while the spectrum graph fell only from 0.471 to 0.392.

Database search was robust to both perturbations (≥ 0.996), confirming that the hard part of de novo sequencing is the missing information, not the spectrum quality per se. These results mark the current limits of small spectrum-native models. More training data that covers sparse and noisy spectra, or hybrid models that combine learned scoring with exhaustive graph search, are the obvious next steps.

{{TABLE4}}

## Learned representations and peak-level attribution

We tested whether the SpecNovo encoder learns peptide identity without being trained to do so. We generated 12 independent noisy replicate spectra for each of 40 test peptides (*p*~obs~ = 0.75, λ = 60) and compared mean-pooled encoder embeddings with 0.5 Da binned spectra (Fig. 5a–c, Table 8).

The encoder embedding grouped replicates more compactly: cosine silhouette 0.273 vs 0.222. Its top-1 replicate retrieval (0.908) was slightly below that of binned spectra (0.996). In simulation, replicate spectra share exact fragment m/z values, which favours binned matching; on real data with calibration drift, learned embeddings such as GLEAMS [19] are expected to have the advantage.

Gradient × input attribution showed that SpecNovo bases its predictions on genuine fragment ions (Fig. 5d,e). Annotated b/y-family ions carried 42.2 % of the total intensity but received 80.0 % of the attribution mass, and attribution exceeded the intensity share in 299 of 300 spectra. The model therefore actively down-weights noise peaks, an interpretable behaviour that classical intensity-weighted scores do not provide.

![**Fig. 5. Representations and interpretability of SpecNovo.** (a, b) t-SNE of 480 replicate spectra (40 peptides × 12) using encoder embeddings (a) or binned spectra (b); six peptides highlighted. (c) Top-1 replicate retrieval and cosine silhouette. (d) Example spectrum: intensity (top; fragment ions blue) and gradient × input saliency (bottom; fragment ions orange). (e) Saliency share versus intensity share on annotated fragment ions for 300 inDB spectra; dashed line, identity.](../figures/fig5_representation.png){width=100%}

{{TABLE8}}

## Spectrum-to-phenotype: removing the identification bottleneck

The same theme dominated the MS1 experiments, at the cohort level (Table 6, Fig. 6). With half of the biomarkers unidentifiable, the identification-first pipeline reached AUROC 0.558 with 24 labelled runs and 0.672 with 384 (Fig. 6a). Two losses compound here: unidentified biomarkers, and 25 % missing values from stochastic DDA identification, which minimum-value imputation turns into noise.

Match-between-runs (MBR) recovered the missing-value loss but not the unidentified biomarkers (0.932 at n = 192). The identification-free linear model, which uses all aligned features, was best at every label budget: 0.820 at n = 24, 0.975 at n = 192 and 0.981 at n = 384. The ProFormer-style FeatureFormer improved substantially on identification-first statistics (0.767 vs 0.648 at n = 192; 0.802 vs 0.672 at n = 384). At these sample sizes it did not reach the identification-free linear model, because the simulated phenotype signal is additive and linear models are near-optimal for it.

Self-supervised masked-intensity pretraining on 800 unlabeled runs did not improve in-distribution accuracy (0.714 vs 0.767 at n = 192). It did improve transfer to an unseen acquisition site: leave-one-site-out AUROC rose from 0.651 ± 0.087 to 0.728 ± 0.054, with lower variance (Fig. 6b). Identification-free features also transferred best overall (0.986 ± 0.012), because ID-free RT-shift alignment absorbs site-specific chromatographic drift.

{{TABLE6}}

![**Fig. 6. Spectrum-to-phenotype benchmark.** (a) Test AUROC versus number of labelled runs (mean ± s.d., three seeds; 50 % of biomarkers unidentifiable). (b) Leave-one-site-out AUROC on an unseen acquisition site. (c) AUROC (n = 192) versus fraction of biomarkers not identifiable by MS2 search. (d) FeatureFormer gradient × input attribution versus true effect size; library (circles) and unidentifiable (squares) biomarkers. (e) ROC curves on the unseen site (seed 0; colours as in a).](../figures/fig6_phenotype.png){width=100%}

## The dark-proteome sweep quantifies information loss

Varying the fraction of unidentifiable biomarkers from 0 % to 75 % separated the pipelines cleanly (Table 7, Fig. 6c). Identification-free linear models were invariant to identifiability (0.973 → 0.968). Identification-first AUROC declined monotonically (0.808 → 0.575), approaching chance once three quarters of the biomarkers were dark. MBR matched identification-free performance while most biomarkers were identifiable (0.973 at 0 %) but collapsed at 75 % (0.705 ± 0.055). FeatureFormer+SSL was likewise invariant to identifiability (0.714–0.759), confirming that it uses unidentified signal, but it remained limited by sample size.

Feature attribution recovered true biomarkers above chance (AUROC 0.622), including unidentifiable ones (Fig. 6d). However, a Welch t-test on identification-free features ranked all 20 biomarkers above all background features (AUROC 1.000). For biomarker *discovery* in additive settings, simple statistics on identification-free features therefore remain the reference against which deep attribution must be judged.

{{TABLE7}}

## Validation on real PRIDE and MassIVE-KB data

**Real spectra are sparser and more y-dominated than the default simulator** (Table 9, Fig. 7c). On 126 MassIVE-KB HCD spectra the median b-ion observation rate was 0.29 (IQR 0.14–0.50) against 0.78 in the default simulator, while y-ion observation was 0.67 (0.50–0.86). Summed y/b intensity was 2.7 against 1.7, and b/y-family ions explained 32 % of the total ion current against 52 %.

The proline effect (1.98 vs 1.73) and the number of unexplained peaks (38.5 vs 41) were close to the simulator. The fragment mass error was tighter in the real data (s.d. 4.4 ppm against 8.9 ppm). PRIDE BSA ion-trap CID spectra were different again: about 145 peaks per spectrum, 59 % of the ion current explained, y/b ratio 4.1 and fragment errors of several hundred ppm (Table 9). These measurements identify the simulator's main miscalibration as the overly symmetric b/y fragmentation.

**Table 9. Calibration of the virtual instrument against real annotated spectra (median [IQR] per PSM).**

| Statistic | Simulator (default) | MassIVE-KB HCD (real) | PRIDE BSA CID (real) |
|---|---|---|---|
| PSMs analysed | 500 | 126 | 64 |
| b-ion observation rate | 0.78 [0.67–0.88] | 0.29 [0.14–0.50] | 0.69 [0.50–0.79] |
| y-ion observation rate | 0.78 [0.67–0.88] | 0.67 [0.50–0.86] | 0.83 [0.67–0.89] |
| Summed y/b intensity | 1.66 [1.12–2.27] | 2.66 [1.59–4.78] | 4.08 [2.51–6.44] |
| Pro N-terminal cleavage enhancement | 1.73 [0.91–2.77] | 1.98 [1.25–3.69] | 2.09 [0.98–3.29] |
| Fraction of TIC explained | 0.52 [0.44–0.60] | 0.32 [0.24–0.43] | 0.59 [0.45–0.68] |
| Peaks per spectrum | 62 [55–68] | 53 [38–70] | 145 [114–192] |
| Unexplained peaks | 41 [37–45] | 38.5 [26–55] | 116.5 [90–162] |
| Fragment mass error (s.d.) | 8.9 ppm | 4.4 ppm | 323 ppm (ion trap) |

**Identification on real PRIDE runs reproduces the identification bottleneck.** The binomial-score search separated targets from decoys (Fig. 7b) and accepted 64 PSMs (20 unique BSA peptides) at 1 % FDR. That is 2.0 % of the 3,135 MS2 spectra; under the conservative (*d* + 1)/*t* estimator, 88 PSMs passed at 5 % FDR.

On the MS1 level, feature detection found 25,613 de-isotoped features across the six files (Table 11, Fig. 7e). Of these, 2,355 (9.2 %) were sampled for MS2 and only 48 (0.19 %) carried a 1 % FDR identification. After ID-free alignment of the three samples, 1,342 consensus features were quantified in all three samples, against a single feature with complete identification-first profiles (Fig. 7f).

These numbers overstate the bottleneck in absolute terms. Our minimal search engine identifies far fewer spectra than production engines (MS-GF+, Comet, Andromeda) with rescoring, and a simple centroid-linking feature finder also counts some noise traces. The direction and the order of magnitude, however, agree with the simulation: most MS1 signal is never identified, and identification-first matrices are dominated by missing values.

**Simulation-to-real transfer of de novo sequencing** (Table 10, Fig. 7d). Applied zero-shot to 119 real MassIVE-KB HCD spectra, the simulation-trained SpecNovo recovered 7.6 % of peptides (95 % CI 3.4–12.6 %) with amino-acid precision 0.261. The spectrum-graph baseline, which encodes only fragmentation physics, reached 21.8 % (15.1–29.4 %), a substantial sim-to-real gap for the learned model.

Fine-tuning for 800 steps on the corpus re-simulated with the measured b/y statistics doubled SpecNovo's peptide recall to 15.1 % (9.2–21.8 %; paired difference +7.6 points, 95 % CI +2.5 to +13.4). Amino-acid precision rose to 0.408, and precision–coverage AUC rose to 0.352 (0.215–0.474), on a par with the spectrum graph (0.346, 0.224–0.466). Calibration did not cost simulated performance: peptide recall on the simulated inDB set rose from 0.416 to 0.525, partly from the additional training steps.

On low-resolution ion-trap CID spectra (62 PSMs ≤ 14 residues, 19 peptides), neither the zero-shot nor the HCD-calibrated SpecNovo recovered complete peptides (amino-acid precision 0.105 and 0.144). The spectrum graph, with a 0.5 Da tolerance, recovered 8.1 %. Both de novo approaches therefore remain instrument-specific, and a learned model must be trained on, or calibrated to, the fragmentation mode and mass accuracy it will be applied to.

**Table 10. De novo sequencing of real spectra (95 % bootstrap CI for MassIVE-KB).**

| Data | Method | AA precision | AA recall | Peptide recall | PR-AUC |
|---|---|---|---|---|---|
| MassIVE-KB HCD (n = 119) | SpectrumGraph-DP | 0.593 [0.534–0.651] | 0.511 | 0.218 [0.151–0.294] | 0.346 [0.224–0.466] |
| | SpecNovo, zero-shot (simulation only) | 0.261 [0.211–0.313] | 0.266 | 0.076 [0.034–0.126] | 0.199 [0.086–0.313] |
| | SpecNovo, calibrated simulation | 0.408 [0.350–0.471] | 0.405 | 0.151 [0.092–0.218] | 0.352 [0.215–0.474] |
| PRIDE BSA ion-trap CID (n = 62) | SpectrumGraph-DP (0.5 Da) | 0.292 | 0.334 | 0.081 | 0.113 |
| | SpecNovo, zero-shot | 0.105 | 0.116 | 0.000 | 0.000 |
| | SpecNovo, calibrated (HCD) | 0.144 | 0.154 | 0.000 | 0.000 |

**Table 11. The identification bottleneck on real PRIDE-derived LC-MS runs (three BSA samples × two RT fractions).**

| Quantity | Value |
|---|---|
| MS2 spectra / MS1 scans | 3,135 / 1,676 |
| Target + decoy database entries | 1,725,026 |
| PSMs at 1 % FDR (unique peptides) | 64 (20); 2.0 % of MS2 spectra |
| De-isotoped MS1 features (six files) | 25,613 |
| Features sampled for MS2 | 2,355 (9.2 %) |
| Features with a 1 % FDR identification | 48 (0.19 %) |
| Consensus features after ID-free alignment | 18,636 |
| Quantified in all three samples: ID-free / ID-first | 1,342 / 1 |
| Missing values: ID-free (all features) / ID-first (identified features) | 54.2 % / 55.6 % |
| Pearson r of log-intensities between samples (ID-free, complete features) | 0.44–0.64 |

![**Fig. 7. Validation on real PRIDE and MassIVE-KB data.** (a) Real annotated MassIVE-KB HCD spectrum (b ions blue, y ions red). (b) Target and decoy score distributions for the PRIDE BSA ion-trap search (log scale). (c) Fragmentation statistics of the default simulator against real HCD and CID spectra (median and IQR; ratios divided by 4 for display). (d) De novo sequencing of 119 real HCD spectra: amino-acid precision (light) and peptide recall (dark) for the spectrum-graph baseline, zero-shot SpecNovo and SpecNovo fine-tuned on calibrated simulation. (e) Identification funnel on real MS1 data: detected features, features sampled for MS2 and features with a 1 % FDR identification. (f) Number of features quantified in all three samples with identification-free and identification-first workflows.](../figures/fig7_real_validation.png){width=100%}

# Discussion

The move from identification-first to spectrum-native proteomics is often framed as a contest between classical algorithms and deep learning. The controlled measurements here suggest a more useful framing: **identification is an information bottleneck, and the main gains come from removing that bottleneck.** The choice of model class matters less.

On the MS1 side, the bottleneck alone accounts for a 0.33 AUROC difference at n = 192 (0.648 vs 0.975), and it grows with the size of the dark proteome. This mirrors the ProFormer observation that identification-free MS1 classification is both faster and more accurate than interpretation-based workflows [17]. It also explains why match-between-runs is so valuable in practice [2]: MBR recovers the missing values, but not the biomarkers that were never identified.

On the MS2 side, database search remains unbeatable for canonical peptides but recovers nothing beyond its database. De novo sequencing provides access to variant and novel sequences at a fixed, database-independent rate, which matters for immunopeptidomics, metaproteomics, antibody sequencing and proteogenomic variant detection [9,10].

The second lesson concerns **data efficiency**. A generic sequence-to-sequence transformer trained on ≈ 0.45 M spectra learned a sequence language model but not spectrum interpretation, and was useless under mass-constrained decoding. A small, physics-informed ion-evidence term, analogous to the ion-location features of DeepNovo [15] and to pairwise mass-difference attention [13], made the same architecture competitive within a few hundred additional optimisation steps.

This result is directly relevant to the vision of proteomics GPTs trained on MassNet-style corpora. Corpus scale and inductive bias are substitutes, and at the scales available to most laboratories, inductive bias is cheaper. Conversely, the gap between our CPU-scale model and published models trained on tens of millions of spectra [9,10] is consistent with data scaling, not architectural limits, being the main driver of state-of-the-art de novo accuracy.

Third, **self-supervised pretraining is not a free lunch at small scale.** Masked-intensity pretraining on 800 unlabeled runs did not improve in-distribution phenotype prediction, but it consistently improved cross-site transfer and reduced variance. This pattern is expected if pretraining learns site-invariant co-abundance structure rather than phenotype-specific signal. It suggests that the payoff of corpus-scale pretraining will appear first in robustness and transfer across instruments, laboratories and cohorts.

**Limitations.** The central comparisons are drawn from simulation. The real-data validation is deliberately small: 3,135 PRIDE-derived MS2 spectra and 119 MassIVE-KB spectra, obtained from public mirrors because the PRIDE archive itself could not be reached from the compute environment. It used a minimal search engine and feature finder, and it could not test phenotype prediction, which requires labelled cohorts of tens to hundreds of runs. The real data nonetheless show that the default virtual instrument overestimates b-ion completeness and that learned de novo models suffer a substantial sim-to-real gap, which calibration halves but does not close. The virtual instrument also omits several features of real data:

- fragment-intensity patterns learned from real data (e.g., Prosit or AlphaPeptDeep [20,21]);
- chimeric and DIA spectra;
- isotope and charge-state ambiguity in precursor selection;
- non-tryptic peptides;
- PTMs other than Met oxidation;
- non-additive biology.

The phenotype signal was deliberately additive, which favours linear models. Training budgets were limited to a few CPU-hours, and the SSL corpus contained 800 runs, far from repository scale. The spectrum-graph baseline is a simplified re-implementation, not a production de novo engine.

These choices are intentional: they make every quantity measurable and every experiment reproducible on a laptop. They also mean that the absolute numbers should not be compared with real-data benchmarks. Future work should:

- extend the calibration loop demonstrated here to full MassIVE-KB/PRIDE reprocessing results (e.g., via quantms [31]), including ion-trap CID;
- validate spectrum-to-phenotype models on labelled PRIDE cohorts annotated with SDRF [32];
- add DIA multiplexing;
- benchmark published models (Casanovo, InstaNovo, ProFormer) inside the same controlled environment;
- scale SSL corpora by orders of magnitude.

# Conclusions

We have provided a reproducible, ground-truth benchmark of the transition from identification-first to spectrum-native proteomics, with measurements of what each paradigm gains and loses:

- **Sequencing.** Database search is exact but blind beyond its database. Database-free transformer sequencing becomes viable only with physics-informed inductive biases at modest data scales, offers the best-calibrated confidence, and still trails classical graph search on sparse or noisy spectra.
- **Phenotype.** Removing identification is the dominant effect, and transformer classifiers and self-supervised pretraining add robustness rather than accuracy at cohort sizes of a few hundred.
- **Real data.** On PRIDE-derived runs only 0.19 % of MS1 features and 2.0 % of MS2 spectra were identified by a minimal engine, confirming the bottleneck. On real HCD spectra, calibrating the simulator to measured fragmentation statistics doubled the learned model's peptide recall (7.6 % → 15.1 %) and matched the spectrum-graph baseline in precision–coverage AUC.

These measurements support continued investment in large, shared, machine-learning-ready spectral corpora, and in hybrid models that unite learned representations with the physics of peptide fragmentation.

# Data and code availability

All source code (simulator, models, baselines, real-data validation, evaluation and figure scripts), result files (JSON/pickle), trained model checkpoints and figures are available in the project repository (`spectrum_native_proteomics/`). Every experiment is reproducible on a CPU with Python ≥ 3.11, PyTorch ≥ 2.0, NumPy, SciPy, scikit-learn, pyteomics and matplotlib, using the commands listed in the repository README. Real data: PRIDE-derived BSA runs from github.com/bigbio/quantms-test-datasets (`testdata/lfq_ci/BSA`, with SDRF) and MassIVE-KB annotated spectra from github.com/Noble-Lab/casanovo (`sample_data/`).

# Funding

Not applicable.

# Conflict of interest

The author declares no competing interests.

# Author contributions

M.N.A.: conceptualisation, methodology, software, investigation, formal analysis, visualisation, writing — original draft, review and editing.
