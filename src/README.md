README
================
Jordan Whitlock
2023-05-31

# Network Construction and Analysis

## Purpose:

All scripts in this directory are used to build necessary inputs
required for network construction as well as downstream analyses.

### Reproducibility:

**Network Construction** \* PANDA jobs were run using UAB Cheaha
Supercomputer and it’s SLURM scheduler. Jobs were run as an array,
within a
[Docker](https://hub.docker.com/repository/docker/jordanwhitlock/jw_diseasenetworks/general)
container converted Singularity in order to execute. \* Human
Protein-Protein Interaction (PPI) and Transcription Factor Motif
(TF-Motif) inputs were previously generated in [Whitlock et al](). The
github repository for constructing both inputs and associated code is
found
[here](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/tree/main)
for both the [Human
TF-Motif](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/PANDA_input_construction/02_Human_TFmotif.Rmd)
and [Human
PPI](https://github.com/lasseignelab/230227_JW_Setbp1Manuscript/blob/main/src/network_scripts/PANDA_input_construction/04_HumanMouse_ppi.Rmd)
\* Expression data was obtained from GTEx using Recount3 (accessed:
220806-220807) and includes all tissue types except blood and study_na

### Scripts:

#### PANDA Expression Inputs and Network Construction:

This directory contains all scripts needed to construct
[PANDA](https://netzoo.github.io/zooanimals/panda/) regulatory networks.
In addition, the *.err* and *.out* files for each array job are included
here to provide detailed information on the jobs. Prior to building
PANDA networks, the user must first construct the expression inputs
needed and then build the *.txt*, which is a list of expression input
file paths needed for the array job.

    ## GTEx_PANDA
    ## +-- 01_array_construction.R
    ## +-- 02_PANDA_array.sh
    ## \-- PANDA.R