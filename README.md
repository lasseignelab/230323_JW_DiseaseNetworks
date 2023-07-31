README
================
2023-07-31

# The Landscape of *SETBP1* tissue-specific gene expression and regulation

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

## Scripts

Here we provide a framework to investigate the following:

-   [Tissue-specific expresion of *SETBP1* and its known
    targets](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/Tissue_Expression)

<!-- -->

    ## src/Tissue_Expression
    ## +-- 01_setbp1_combine_targets.R
    ## +-- 02_setbp1_expression.R
    ## \-- 03_median_expression_plot.R

-   [Construction of tissue-specific gene regulatory
    networks](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/GTEx_PANDA)

<!-- -->

    ## src/GTEx_PANDA
    ## +-- 01_array_construction.R
    ## +-- 02_PANDA_array.sh
    ## +-- PANDA.R
    ## +-- gPANDA_21053445_0.err
    ## +-- gPANDA_21053445_0.out
    ## +-- gPANDA_21053445_1.err
    ## +-- gPANDA_21053445_1.out
    ## +-- gPANDA_21053445_10.err
    ## +-- gPANDA_21053445_10.out
    ## +-- gPANDA_21053445_11.err
    ## +-- gPANDA_21053445_11.out
    ## +-- gPANDA_21053445_12.err
    ## +-- gPANDA_21053445_12.out
    ## +-- gPANDA_21053445_13.err
    ## +-- gPANDA_21053445_13.out
    ## +-- gPANDA_21053445_14.err
    ## +-- gPANDA_21053445_14.out
    ## +-- gPANDA_21053445_15.err
    ## +-- gPANDA_21053445_15.out
    ## +-- gPANDA_21053445_16.err
    ## +-- gPANDA_21053445_16.out
    ## +-- gPANDA_21053445_17.err
    ## +-- gPANDA_21053445_17.out
    ## +-- gPANDA_21053445_18.err
    ## +-- gPANDA_21053445_18.out
    ## +-- gPANDA_21053445_19.err
    ## +-- gPANDA_21053445_19.out
    ## +-- gPANDA_21053445_2.err
    ## +-- gPANDA_21053445_2.out
    ## +-- gPANDA_21053445_20.err
    ## +-- gPANDA_21053445_20.out
    ## +-- gPANDA_21053445_21.err
    ## +-- gPANDA_21053445_21.out
    ## +-- gPANDA_21053445_22.err
    ## +-- gPANDA_21053445_22.out
    ## +-- gPANDA_21053445_23.err
    ## +-- gPANDA_21053445_23.out
    ## +-- gPANDA_21053445_24.err
    ## +-- gPANDA_21053445_24.out
    ## +-- gPANDA_21053445_25.err
    ## +-- gPANDA_21053445_25.out
    ## +-- gPANDA_21053445_26.err
    ## +-- gPANDA_21053445_26.out
    ## +-- gPANDA_21053445_27.err
    ## +-- gPANDA_21053445_27.out
    ## +-- gPANDA_21053445_28.err
    ## +-- gPANDA_21053445_28.out
    ## +-- gPANDA_21053445_29.err
    ## +-- gPANDA_21053445_29.out
    ## +-- gPANDA_21053445_3.err
    ## +-- gPANDA_21053445_3.out
    ## +-- gPANDA_21053445_30.err
    ## +-- gPANDA_21053445_30.out
    ## +-- gPANDA_21053445_4.err
    ## +-- gPANDA_21053445_4.out
    ## +-- gPANDA_21053445_5.err
    ## +-- gPANDA_21053445_5.out
    ## +-- gPANDA_21053445_6.err
    ## +-- gPANDA_21053445_6.out
    ## +-- gPANDA_21053445_7.err
    ## +-- gPANDA_21053445_7.out
    ## +-- gPANDA_21053445_8.err
    ## +-- gPANDA_21053445_8.out
    ## +-- gPANDA_21053445_9.err
    ## \-- gPANDA_21053445_9.out

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
    ## +-- 04_condor_network_analysis.Rmd
    ## +-- condor_21066636_0.condor.err
    ## +-- condor_21066636_0.condor.out
    ## +-- condor_21066636_1.condor.err
    ## +-- condor_21066636_1.condor.out
    ## +-- condor_21066636_10.condor.err
    ## +-- condor_21066636_10.condor.out
    ## +-- condor_21066636_11.condor.err
    ## +-- condor_21066636_11.condor.out
    ## +-- condor_21066636_12.condor.err
    ## +-- condor_21066636_12.condor.out
    ## +-- condor_21066636_13.condor.err
    ## +-- condor_21066636_13.condor.out
    ## +-- condor_21066636_14.condor.err
    ## +-- condor_21066636_14.condor.out
    ## +-- condor_21066636_15.condor.err
    ## +-- condor_21066636_15.condor.out
    ## +-- condor_21066636_16.condor.err
    ## +-- condor_21066636_16.condor.out
    ## +-- condor_21066636_17.condor.err
    ## +-- condor_21066636_17.condor.out
    ## +-- condor_21066636_18.condor.err
    ## +-- condor_21066636_18.condor.out
    ## +-- condor_21066636_19.condor.err
    ## +-- condor_21066636_19.condor.out
    ## +-- condor_21066636_2.condor.err
    ## +-- condor_21066636_2.condor.out
    ## +-- condor_21066636_20.condor.err
    ## +-- condor_21066636_20.condor.out
    ## +-- condor_21066636_21.condor.err
    ## +-- condor_21066636_21.condor.out
    ## +-- condor_21066636_22.condor.err
    ## +-- condor_21066636_22.condor.out
    ## +-- condor_21066636_23.condor.err
    ## +-- condor_21066636_23.condor.out
    ## +-- condor_21066636_24.condor.err
    ## +-- condor_21066636_24.condor.out
    ## +-- condor_21066636_25.condor.err
    ## +-- condor_21066636_25.condor.out
    ## +-- condor_21066636_26.condor.err
    ## +-- condor_21066636_26.condor.out
    ## +-- condor_21066636_27.condor.err
    ## +-- condor_21066636_27.condor.out
    ## +-- condor_21066636_28.condor.err
    ## +-- condor_21066636_28.condor.out
    ## +-- condor_21066636_29.condor.err
    ## +-- condor_21066636_29.condor.out
    ## +-- condor_21066636_3.condor.err
    ## +-- condor_21066636_3.condor.out
    ## +-- condor_21066636_30.condor.err
    ## +-- condor_21066636_30.condor.out
    ## +-- condor_21066636_31.condor.err
    ## +-- condor_21066636_31.condor.out
    ## +-- condor_21066636_4.condor.err
    ## +-- condor_21066636_4.condor.out
    ## +-- condor_21066636_5.condor.err
    ## +-- condor_21066636_5.condor.out
    ## +-- condor_21066636_6.condor.err
    ## +-- condor_21066636_6.condor.out
    ## +-- condor_21066636_7.condor.err
    ## +-- condor_21066636_7.condor.out
    ## +-- condor_21066636_8.condor.err
    ## +-- condor_21066636_8.condor.out
    ## +-- condor_21066636_9.condor.err
    ## \-- condor_21066636_9.condor.out

-   [SETBP1 TF activity
    analysis](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/tf_activity)

<!-- -->

    ## src/tf_activity
    ## +-- 01_GTExSetbp1_decoupleR_inputs.Rmd
    ## +-- 02_decoupleR_analysis.R
    ## +-- 02_decoupleR_array_job.sh
    ## +-- 03_TF_activity_GTEx.Rmd
    ## +-- gtex_decoupleR_21064875_0.err
    ## +-- gtex_decoupleR_21064875_0.out
    ## +-- gtex_decoupleR_21064875_1.err
    ## +-- gtex_decoupleR_21064875_1.out
    ## +-- gtex_decoupleR_21064875_10.err
    ## +-- gtex_decoupleR_21064875_10.out
    ## +-- gtex_decoupleR_21064875_11.err
    ## +-- gtex_decoupleR_21064875_11.out
    ## +-- gtex_decoupleR_21064875_12.err
    ## +-- gtex_decoupleR_21064875_12.out
    ## +-- gtex_decoupleR_21064875_13.err
    ## +-- gtex_decoupleR_21064875_13.out
    ## +-- gtex_decoupleR_21064875_14.err
    ## +-- gtex_decoupleR_21064875_14.out
    ## +-- gtex_decoupleR_21064875_15.err
    ## +-- gtex_decoupleR_21064875_15.out
    ## +-- gtex_decoupleR_21064875_16.err
    ## +-- gtex_decoupleR_21064875_16.out
    ## +-- gtex_decoupleR_21064875_17.err
    ## +-- gtex_decoupleR_21064875_17.out
    ## +-- gtex_decoupleR_21064875_18.err
    ## +-- gtex_decoupleR_21064875_18.out
    ## +-- gtex_decoupleR_21064875_19.err
    ## +-- gtex_decoupleR_21064875_19.out
    ## +-- gtex_decoupleR_21064875_2.err
    ## +-- gtex_decoupleR_21064875_2.out
    ## +-- gtex_decoupleR_21064875_20.err
    ## +-- gtex_decoupleR_21064875_20.out
    ## +-- gtex_decoupleR_21064875_21.err
    ## +-- gtex_decoupleR_21064875_21.out
    ## +-- gtex_decoupleR_21064875_22.err
    ## +-- gtex_decoupleR_21064875_22.out
    ## +-- gtex_decoupleR_21064875_23.err
    ## +-- gtex_decoupleR_21064875_23.out
    ## +-- gtex_decoupleR_21064875_24.err
    ## +-- gtex_decoupleR_21064875_24.out
    ## +-- gtex_decoupleR_21064875_25.err
    ## +-- gtex_decoupleR_21064875_25.out
    ## +-- gtex_decoupleR_21064875_26.err
    ## +-- gtex_decoupleR_21064875_26.out
    ## +-- gtex_decoupleR_21064875_27.err
    ## +-- gtex_decoupleR_21064875_27.out
    ## +-- gtex_decoupleR_21064875_28.err
    ## +-- gtex_decoupleR_21064875_28.out
    ## +-- gtex_decoupleR_21064875_29.err
    ## +-- gtex_decoupleR_21064875_29.out
    ## +-- gtex_decoupleR_21064875_3.err
    ## +-- gtex_decoupleR_21064875_3.out
    ## +-- gtex_decoupleR_21064875_30.err
    ## +-- gtex_decoupleR_21064875_30.out
    ## +-- gtex_decoupleR_21064875_4.err
    ## +-- gtex_decoupleR_21064875_4.out
    ## +-- gtex_decoupleR_21064875_5.err
    ## +-- gtex_decoupleR_21064875_5.out
    ## +-- gtex_decoupleR_21064875_6.err
    ## +-- gtex_decoupleR_21064875_6.out
    ## +-- gtex_decoupleR_21064875_7.err
    ## +-- gtex_decoupleR_21064875_7.out
    ## +-- gtex_decoupleR_21064875_8.err
    ## +-- gtex_decoupleR_21064875_8.out
    ## +-- gtex_decoupleR_21064875_9.err
    ## \-- gtex_decoupleR_21064875_9.out

-   [TF
    targeting](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/tf_targeting)

<!-- -->

    ## src/tf_targeting
    ## \-- 01_TF_targeting.Rmd

-   [Tissue-specific edge
    identification](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/tissue_specific_edges)

<!-- -->

    ## src/tissue_specific_edges
    ## +-- 01_tissue_specific_edges.R
    ## +-- 01_tissue_specific_edges.sh
    ## +-- 02_tissue_specific_edges.Rmd
    ## +-- tse.err
    ## \-- tse.out

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
