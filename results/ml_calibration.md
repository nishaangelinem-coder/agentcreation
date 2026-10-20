# Stage 2 reliability note

- Model: RandomForest (300 trees, balanced) on 2048-bit Morgan (r=2) fingerprints (sklearn baseline; DeepChem GNN is an optional upgrade).
- Validation: deterministic Bemis-Murcko **scaffold split** (80/20); test scaffolds are unseen in training.
- Headline metric is **AUPRC**, because Tox21 positive classes are rare (test base rates below).
- Probabilities are raw RandomForest vote fractions and are **not isotonically calibrated**; treat them as a relative ranking signal, not absolute activation probabilities.
- Leachate scores are produced by models refit on all labelled Tox21 rows for each endpoint.

| endpoint | n_test | n_pos_test | base_rate | ROC-AUC | AUPRC |
|---|---|---|---|---|---|
| NR-AhR | 1220 | 174 | 0.1426 | 0.8399 | 0.5449 |
| SR-ARE | 943 | 226 | 0.2397 | 0.7154 | 0.4658 |
| SR-p53 | 1258 | 138 | 0.1097 | 0.7591 | 0.3473 |
