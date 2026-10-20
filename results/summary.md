# Pipeline summary — Microplastics & Neuroinflammation

## Predicted polymer ranking (high -> low neuroinflammatory potential)

**PU > PVC > PE > PS > PP**

| polymer   |   index |   receptor_activation_u | dominant_receptor   | dominant_leachate             |   NR-AhR_top3 |   SR-ARE_top3 |   SR-p53_top3 | used_docking   |
|:----------|--------:|------------------------:|:--------------------|:------------------------------|--------------:|--------------:|--------------:|:---------------|
| PU        |  0.976  |                  0.8506 | NR-AhR              | 2,4-Toluenediamine (TDA)      |        0.9967 |        0.9822 |        0.2878 | False          |
| PVC       |  0.7415 |                  0.3345 | SR-ARE              | Bisphenol A (BPA)             |        0.1198 |        0.75   |        0.2478 | False          |
| PE        |  0.6981 |                  0.2968 | SR-ARE              | 4-tert-Octylphenol            |        0.2353 |        0.4322 |        0.2478 | False          |
| PS        |  0.6432 |                  0.2577 | SR-ARE              | Hexabromocyclododecane (HBCD) |        0.1038 |        0.56   |        0.1889 | False          |
| PP        |  0.5706 |                  0.2159 | NR-AhR              | 2,4-di-tert-butylphenol       |        0.2022 |        0.2878 |        0.1422 | False          |

## Stage-2 scaffold-split metrics (imbalance-aware)

| endpoint   |   n_train |   n_test |   n_pos_test |   base_rate_test |   roc_auc |   auprc |
|:-----------|----------:|---------:|-------------:|-----------------:|----------:|--------:|
| NR-AhR     |      5322 |     1220 |          174 |           0.1426 |    0.8399 |  0.5449 |
| SR-ARE     |      4882 |      943 |          226 |           0.2397 |    0.7154 |  0.4658 |
| SR-p53     |      5509 |     1258 |          138 |           0.1097 |    0.7591 |  0.3473 |

## Validation vs experiment

experimental_ranking.csv present but not yet filled in.

## Single most important caveat

**Docking/ML affinity is binding STRENGTH, not activation DIRECTION.** The index assumes the engaged receptors are driven toward a pro-inflammatory state; agonist/antagonist direction must be confirmed per leachate (here the Tox21 label supplies direction for the AhR/stress axes). All virtual-cell kinetic parameters are illustrative, not fitted. Stage-3 receptor docking (TLR4-MD2, NLRP3) could not be run because structure databases were unreachable in the execution sandbox (see docking_status.md); the index therefore rests on the Tox21 ML signal. The ranking is a falsifiable PREDICTION for wet-lab testing, not established fact.
