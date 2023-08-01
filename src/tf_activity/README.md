README
================
Jordan H. Whitlock
2023-08-01

# GTEx SETBP1 decoupleR inputs:

## TF regulation prior: CollecTRI (human)

-   The prior network for
    [decoupleR](https://saezlab.github.io/decoupleR/) was constructed
    using [CollecTRI](https://github.com/saezlab/CollecTRI) for human
    on 230525. The file is located in this github repository as
    `human_prior_tri.csv`

## Expression prior: GTEx expression data (tpm)

-   31 different tissues from human GTEx data, normalized by tpm.
-   input array was already constructed in PANDA network construction
    and is titled `GTEx_exp_files_array.txt`