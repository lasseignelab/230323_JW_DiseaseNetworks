#!/bin/bash
## run the Rscript 01_tissue_specific_edges.R and schedule this job to SLURM with
## `sbatch 01_tissue_specific_edges.sh`

#SBATCH --job-name=tse
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jbarham3@uab.edu #CHANGE THE EMAIL TO YOUR OWN
#
#Number of tasks needed for this job, the partition, and parameters
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=255000
#SBATCH --nodes=1
#SBATCH --time=150:00:00
#SBATCH --partition=largemem-long #Note: when running brain, can use largemem. Partition info here: https://docs.rc.uab.edu/cheaha/hardware/#partitions
#
#SBATCH --output=tse.out
#SBATCH --error=tse.err


########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
module load R
module load Singularity/3.5.2-GCC-5.4.0-2.26

wd=/data/user/jbarham3/230323_JW_DiseaseNetworks/ #change this to match your path

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER='jbarham3' #change this to your user

cd ${wd}

singularity exec --cleanenv --no-home -B ${wd} ${wd}/bin/downstream_analysis_docker/setbp1_manuscript_1.0.6.sif Rscript --vanilla ${wd}src/tissue_specific_edges/01_tissue_specific_edges.R