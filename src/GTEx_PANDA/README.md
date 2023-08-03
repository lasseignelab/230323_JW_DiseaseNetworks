README
================
Jordan Whitlock
2023-08-03

# Network Construction and Analysis

## Purpose:

All scripts in this directory are used to build necessary inputs
required for network construction as well as downstream analyses.

### Reproducibility:

**Network Construction** \* PANDA jobs were run using UAB Cheaha
Supercomputer and it’s SLURM scheduler. Jobs were run as an array,
within a
[Docker](https://hub.docker.com/repository/docker/jordanwhitlock/jw_diseasenetworks/general)
container converted Singularity in order to execute.

-   Expression data was obtained from GTEx using Recount3 (accessed:
    220806-220807) and includes all tissue types except blood and
    study_na
-   Human PPI and TF Motif input construction is detailed
    [here](https://github.com/lasseignelab/230323_JW_DiseaseNetworks/tree/main/src/panda_input_construction).

### Scripts:

#### PANDA Expression Inputs and Network Construction:

This directory contains all scripts needed to construct
[PANDA](https://netzoo.github.io/zooanimals/panda/) regulatory networks.
Prior to building PANDA networks, the user must first construct the
expression inputs needed and then build the *.txt*, which is a list of
expression input file paths needed for the array job.

    ## .
    ## +-- 01_array_construction.R
    ## +-- 02_PANDA.R
    ## +-- 02_PANDA_array.sh
    ## +-- README.Rmd
    ## \-- README.md