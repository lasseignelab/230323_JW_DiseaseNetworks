README
================
2023-08-02

# The Landscape of *SETBP1* tissue-specific gene expression and regulation across human tissues

## Authors

**Jordan H. Whitlock Elizabeth J. Wilk, Timothy C. Howton, Amanda D.
Clark, Brittany N. Lasseigne**

The University of Alabama at Birmingham (UAB), Heersink School of
Medicine, Department of Cell, Developmental and Integrative Biology
(CDIB)

[The Lasseigne Lab](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="90" height="90">

## Purpose

**The purpose of this research is to investigate the tissue-specific
expression and regulatory landscape of human *SETBP1* in GTEx tissues**
This repository contains the code and accompanying data for our analysis
of the tissue-specific expression and regulation of transcription factor
(TF) *SETBP1* across 31 non-diseased human tissues part of the
Genotype-tissue expression project (GTEx). This project is hypothesis
generating, therefore emphasizing the role of different contexts, such
as tissues, and the role they may play in disease when a variant is
introduced.

## Overview
![Overview_Fig (2)](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/assets/62023125/fde343b2-ec31-4dcf-9e75-d435335dfcc1)

## Scripts

Here we provide a framework to investigate the following:

-   [Tissue-specific expresion of *SETBP1* and its known
    targets](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/Tissue_Expression)

<!-- -->

    ## src/Tissue_Expression
    ## +-- 01_setbp1_combine_targets.R
    ## +-- 02_setbp1_expression.R
    ## +-- 03_median_expression_plot.R
    ## \-- 04_median_expression_heatmap.R

-   [Construction of tissue-specific gene regulatory
    networks](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/GTEx_PANDA)

<!-- -->

    ## src/GTEx_PANDA
    ## +-- 01_array_construction.R
    ## +-- 02_PANDA.R
    ## \-- 02_PANDA_array.sh

-   Information on Network inputs and how to produce them can be found
    [here](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/panda_input_construction).

<!-- -->

    ## src/panda_input_construction
    ## +-- 01_Human_TFmotif.Rmd
    ## +-- 02_Human_TF_motif_enrichment.Rmd
    ## \-- 03_Human_PANDA_ppi.Rmd

-   [Tissue-specific community
    detection](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/community_detection)

<!-- -->

    ## src/community_detection
    ## +-- 01_condor_file_inputs.Rmd
    ## +-- 02_condor_array_job.sh
    ## +-- 02_condor_networks_array.R
    ## +-- 03_Setbp1_Community_Identification.Rmd
    ## \-- 04_condor_network_analysis.Rmd

-   [SETBP1 TF activity
    analysis](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/tf_activity)

<!-- -->

    ## src/tf_activity
    ## +-- 01_decoupleR_analysis.R
    ## +-- 01_decoupleR_array_job.sh
    ## +-- 02_TF_activity_GTEx.Rmd
    ## +-- README.Rmd
    ## \-- README.md

-   [TF
    targeting](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/tf_targeting)

<!-- -->

    ## src/tf_targeting
    ## \-- 01_TF_targeting.Rmd

-   [Tissue-specific edge
    identification](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/tissue_specific_edges)

<!-- -->

    ## src/tissue_specific_edges
    ## \-- 01_tissue_specific_edge_identification.Rmd

*However, this code can be adapted and applied to other genes or TFs.*

## Dependencies and Resources

This analysis was carried out in Docker using R version 4.1.3. TF
activity inference using [decoupleR](), [PANDA]() tissue-specific
regulatory network generation, and tissue-specific edge identification
were run using the Docker on UAB’s high performance computing Cluster.
Bash scripts including resources used are included in this repository.
The containers have been made publicly available on Zenodo:

\[put Zenodo buttons here\]

## Funding

List project funding sources.

## Acknowledgements

List project acknowledgements.

## License

[![License](https://img.shields.io/badge/LICENSE-MIT_License-yellow)](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/blob/main/LICENSE)

This repository is licensed under the MIT License, see LICENSE
documentation within this repository for more details.
