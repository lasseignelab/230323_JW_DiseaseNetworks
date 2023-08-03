README
================
Jordan Whitlock
2023-05-31

# Network Construction

## Purpose:

All scripts in this directory are used to build necessary inputs
required for network construction as well as downstream analyses.

### Reproducibility:

**Network Construction** \* PANDA jobs were run using UAB Cheaha
Supercomputer and it’s SLURM scheduler. Jobs were run as an array,
within a
[Docker](https://hub.docker.com/repository/docker/jordanwhitlock/jw_diseasenetworks/general)
container converted Singularity in order to execute.

* Expression data was obtained from GTEx using Recount3 (accessed:
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
