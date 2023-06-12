#!/bin/bash
## run the Rscript PANDA.R and schedule this job to SLURM with
## `sbatch 02_PANDA_array.sh`

#SBATCH --job-name=gPANDA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jwhitlock@uab.edu #CHANGE THE EMAIL TO YOUR OWN
#
#Number of tasks needed for this job, the partition, and parameters
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=255000
#SBATCH --nodes=1
#SBATCH --time=50:00:00
#SBATCH --partition=largemem #partition info here: https://docs.rc.uab.edu/cheaha/hardware/#partitions
#
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-30 #31 tissue types from GTEx networks

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

#load modules
module load Singularity/3.5.2-GCC-5.4.0-2.26

#variables
wd="/data/user/jbarham3/230323_JW_DiseaseNetworks"
src="/data/user/jbarham3/230323_JW_DiseaseNetworks/src" #be sure that your subdirectories are structured the same

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER='jbarham3'

#code to execute docker and script for analysis
cd ${wd}

#array file of cell type specific expression inputs for PANDA
SAMPLE_LIST="${wd}/results/array_inputs/GTEx_exp_files_array.txt" #note: make sure path and file name are correct
SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`) # parantheses instruct bash to create a shell array of strings from SAMPLE_LIST
INPUT=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}` #extracts a single input from the array and prints (using echo) it into INPUT variable, each single input is then assigned an array number by $SLURM_TASK_ID

# NOTE user must have already pulled from docker and built .sif file with singularity below (jordanwhitlock/setbp1_manuscript_panda_1.0.1)
singularity exec --cleanenv --containall -B ${wd} -B /data/project/lasseigne_lab/DATASET_dir/TCGA_GTEx_CCLE_HPA_220807/gtex/ /data/user/jbarham3/230323_JW_DiseaseNetworks/bin/PANDA_construction_docker/jw_diseasenetworks_1.0.0.sif Rscript --vanilla ${src}/GTEx_PANDA/PANDA.R ${INPUT} # here vanilla ensures only the script is run and environment is kept clean