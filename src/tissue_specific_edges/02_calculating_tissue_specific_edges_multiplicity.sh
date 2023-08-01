#!/bin/bash
## run the Rscript 02_calculating_tissue_specific_edges_multiplicity.R and schedule this job to SLURM with
## `sbatch 02_calculating_tissue_specific_edges_multiplicity.sh`

#SBATCH --job-name=tse_mult
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jbarham3@uab.edu #CHANGE THE EMAIL TO YOUR OWN
#
#Number of tasks needed for this job, the partition, and parameters
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=800GB
#SBATCH --nodes=1
#SBATCH --time=50:00:00
#SBATCH --partition=largemem #Note: when running brain, can use largemem. Partition info here: https://docs.rc.uab.edu/cheaha/hardware/#partitions
#
#SBATCH --output=tse_mult.out
#SBATCH --error=tse_mult.err


########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
module load R
module load Singularity/3.5.2-GCC-5.4.0-2.26

wd=/data/user/jbarham3/230323_JW_DiseaseNetworks/ #change this to match your path

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER='jbarham3' #change this to your user

cd ${wd}

singularity exec --cleanenv --no-home -B ${wd} ${wd}/bin/downstream_analysis_docker/setbp1_manuscript_1.0.6.sif Rscript --vanilla ${wd}src/tissue_specific_edges/02_calculating_tissue_specific_edges_multiplicity.R