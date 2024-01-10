README
================
2023-08-02

# The landscape of *SETBP1* gene expression and transcription factor activity across human tissues

## Authors

**Jordan H. Whitlock Elizabeth J. Wilk, Timothy C. Howton, Amanda D.
Clark, Brittany N. Lasseigne**

The University of Alabama at Birmingham (UAB), Heersink School of
Medicine, Department of Cell, Developmental and Integrative Biology
(CDIB)

## Citation
> Whitlock, J. H., Wilk, E. J., Howton, T. C., Clark, A. D., & Lasseigne, B. N. (2024). The landscape of SETBP1 gene expression and transcription factor activity across human tissues. PLOS ONE, 19(1), e0296328. https://doi.org/10.1371/journal.pone.0296328


[The Lasseigne Lab](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="90" height="90">

## Purpose

**The purpose of this research is to investigate the tissue-specific
expression and TF activity landscape of human *SETBP1* in GTEx tissues**
This repository contains the code and accompanying data for our analysis
of the tissue-specific expression and regulation of transcription factor
(TF) *SETBP1* across 31 non-diseased human tissues part of the
Genotype-tissue expression project (GTEx). This project is hypothesis
generating, therefore emphasizing the role of different contexts, such
as tissues and the role they may play in disease when a variant is
introduced.

As part of this project, we also developed an interactive **GTEx TF Activity Web Application** that can be accessed [here](https://lasseignelab.shinyapps.io/gtex_tf_activity/). Our web application enables researchers to investigate the activity of their favorite TFs in order to generate hypotheses about their role in a diseased setting. 

## Supplemental Data Availability:

- **CollecTRI Data:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8222799.svg)](https://doi.org/10.5281/zenodo.8222799)

## Overview
![GTC_Overview_Fig_1 (5)](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/assets/62023125/4110d641-b2b9-4965-a781-af166bc133b4)

## Scripts

Here we provide a framework to investigate the following:

-   [Tissue-specific expression of *SETBP1* and its known
    targets](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/Tissue_Expression)

<!-- -->

    ## src/Tissue_Expression
    ## +-- 01_setbp1_combine_targets.R
    ## +-- 02_setbp1_expression.R
    ## +-- 03_median_expression_plot.R
    ## \-- 04_median_expression_heatmap.R

-   [SETBP1 TF activity
    analysis](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/tf_activity)

<!-- -->

    ## src/tf_activity
    ## +-- 01_decoupleR_analysis.R
    ## +-- 01_decoupleR_array_job.sh
    ## +-- 02_TF_activity_GTEx.Rmd
    ## +-- README.Rmd
    ## \-- README.md

## Dependencies and Resources

This analysis was carried out in Docker using R version 4.1.3. TF
activity inference using a multivariate linear model [(decoupleR)](https://saezlab.github.io/decoupleR/) was run using Docker on UAB’s high-performance computing Cluster.
Bash scripts, including resources used, are included in this repository.
The containers have been made publicly available on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8428932.svg)](https://doi.org/10.5281/zenodo.8428932)

## Additional DOIs
- **Repository:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8225613.svg)](https://doi.org/10.5281/zenodo.8225613)
- **Shiny Application:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8225317.svg)](https://doi.org/10.5281/zenodo.8225317)

## Funding
This work was supported in part by the UAB Lasseigne Lab funds, UAB Pilot Center for Precision Animal Modeling (C-PAM)(1U54OD030167), and JW UAB Predoctoral Training Grant in Cell, Molecular, and Developmental Biology (CMDB T32)(5T32GM008111-35).

## Acknowledgements
The authors thank the Lasseigne Lab members Vishal Oza, Tabea Soelter, Emma Jones, and Victoria Flanary for their feedback throughout this study. In addition, we thank Vishal Oza for his previously published Jaccard Similarity analysis code we adapted and used for this project. We also thank the UAB Biological Data Science group (RRID:SCR_021766) for providing a script for helping to run containers on the UAB high-performance cluster (https://github.com/U-BDS/training_guides/blob/main/run_rstudio_singularity.sh).

## License

[![License](https://img.shields.io/badge/LICENSE-MIT_License-yellow)](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/blob/main/LICENSE)

This repository is licensed under the MIT License; see LICENSE
documentation within this repository for more details.
