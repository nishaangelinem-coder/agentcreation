---
title: "From identification-first to spectrum-native proteomics: a controlled benchmark, with real-data validation, of transformer models that learn directly from MS1 and MS2 signal"
author:
  - "M. Nisha Angeline^1,\\*^"
date: ""
---

^1^ Department of Electronics and Communication Engineering, Velalar College of Engineering and Technology, Thindal, Erode, Tamil Nadu, India.
^\*^ Corresponding author: nishavlsidesign@gmail.com

**Article type:** Research article (Problem solving protocol / benchmark)
**Suggested journals (Q1):** *Briefings in Bioinformatics*; *Computers in Biology and Medicine*; *Journal of Proteome Research*

# Abstract

**Motivation.** Mass-spectrometry proteomics is moving from an *identification-first* workflow — search MS2 spectra against a database, then do statistics on identified peptides — towards *spectrum-native* models that learn directly from MS1/MS2 signal. Examples include transformer de novo sequencing and ProFormer-style phenotype classifiers. Real data cannot tell us how much information identification discards, because the ground truth of real spectra is itself defined by identification.

**Results.** We built an open, CPU-trainable *virtual instrument* with known ground truth. It simulates b/y fragmentation (proline/aspartate effects, neutral losses, doubly charged ions, stochastic fragment dropout, ppm error, chemical noise) and multi-site LC-MS1 acquisition. In some cohorts a controlled fraction of biomarkers cannot be identified by MS2 search. On 250,000 training spectra we trained SpecNovo, an encoder–decoder transformer with three inductive biases: complementary-ion encoding, prefix-mass decoding and a differentiable ion-evidence term. The ablation showed that the ion-evidence term is decisive. It raised teacher-forced token accuracy from 0.582 to 0.864 and peptide recall from 0.005 to 0.416 at a similar training budget.

Database search was near-perfect on canonical peptides (recall 0.999) but recovered none of the variant or novel-proteome peptides. SpecNovo was database-independent (peptide recall 0.416, 0.389 and 0.387 on canonical, single-amino-acid-variant and novel-proteome spectra). It had the best precision–coverage AUC (0.663 vs 0.630 for a spectrum-graph baseline) and the highest precision at 25 % coverage (0.787 vs 0.716). A classical spectrum-graph algorithm remained more robust to incomplete and noisy spectra.

For phenotype prediction, removing the identification bottleneck was the largest single effect. With 50 % of biomarkers unidentifiable and 192 labelled runs, identification-free features gave AUROC 0.975 ± 0.015 versus 0.648 ± 0.013 for identification-first quantification. On an unseen acquisition site the figures were 0.986 vs 0.627. A ProFormer-style FeatureFormer beat identification-first statistics (0.767 vs 0.648) but not identification-free linear models at these cohort sizes. Masked-intensity pretraining improved site transfer (0.728 vs 0.651) but not in-distribution accuracy.

We validated the benchmark on real data: PRIDE-derived LC-MS/MS runs (3,135 MS2 spectra) and 119 annotated MassIVE-KB spectra. Real HCD spectra were more y-dominated than the default simulator (b-ion observation 0.29 vs 0.78). On real runs only 0.19 % of 25,613 MS1 features were identified at 1 % FDR, and 1,342 features had complete identification-free profiles against one with identification-first quantification. Applied zero-shot to real spectra, SpecNovo recalled 7.6 % of peptides. After fine-tuning on simulation calibrated to the measured fragmentation statistics it recalled 15.1 % (95 % CI 9.2–21.8 %), with precision–coverage AUC equal to the spectrum-graph baseline (0.352 vs 0.346).

**Conclusions.** Spectrum-native learning recovers information that identification-first pipelines lose, especially for variant, novel and dark-proteome signal. At realistic data scales, though, its benefits depend strongly on physics-informed inductive biases and on corpus size. Real-data validation confirms the identification bottleneck and shows that simulator calibration is necessary for learned models to transfer. The benchmark, code and all measurements are released to support "proteomics GPT" development under controlled ground truth.

**Keywords:** de novo peptide sequencing; transformer; spectrum-to-phenotype; self-supervised learning; mass spectrometry; benchmark; explainable AI

# Key points

- An open virtual instrument with exact ground truth quantifies, for the first time in a single framework, what identification-first proteomics loses relative to spectrum-native learning on both MS2 (sequencing) and MS1 (phenotype) tasks.
- Physics-informed inductive biases — above all a differentiable fragment-ion-evidence term — turn an under-trained sequence-to-sequence transformer (0.5 % peptide recall) into a database-free sequencer (41.6 %) with the best-calibrated confidence among tested methods.
- Identification-free MS1 features raise phenotype AUROC from 0.65 to 0.98 when half of the biomarkers are unidentifiable and remain stable from 0 % to 75 % "dark" biomarkers, whereas identification-first performance falls from 0.81 to 0.58.
- Real PRIDE and MassIVE-KB data confirm the bottleneck (0.19 % of MS1 features identified). They also show that calibrating the simulator to measured b/y statistics doubles real-spectrum de novo recall of the learned model.
- At cohort sizes of 24–384 samples, linear models on identification-free features outperform ProFormer-style transformers; self-supervised pretraining chiefly improves cross-site robustness — a quantitative argument for corpus-scale (MassNet-style) pretraining.
