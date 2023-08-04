#!/bin/bash
## run the Rscript 01_decoupleR_analysis.R and schedule this job to SLURM with
## `sbatch 01_decoupleR_array_job.sh`

#SBATCH --job-name=gtex_decoupleR 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER@uab.edu #CHANGE THE EMAIL TO YOUR OWN
#
#Number of tasks needed for this job, the partition, and parameters
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=255000
#SBATCH --nodes=1
#SBATCH --time=50:00:00 #Note: '50:00:00 
#SBATCH --partition=largemem 
#
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-31 

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
module load R
module load Singularity/3.5.2-GCC-5.4.0-2.26

wd="$USER_DATA/230323_JW_DiseaseNetworks" #filepath to repo
src="$USER_DATA/230323_JW_DiseaseNetworks/src" #be sure that your subdirectories are structured the same

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER='$USER'

cd ${wd}

#input expression .Rdata file
SAMPLE_LIST="${wd}/results/array_inputs/GTEx_exp_files_array.txt" #note: make sure path and file name are correct
SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`) # parantheses instruct bash to create a shell array of strings from SAMPLE_LIST
GTEX_FILE=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}` #extracts a single input from the array and prints (using echo) it into INPUT variable, each single input is then assigned an array number by $SLURM_TASK_ID


#input prior network file
PRIOR_NET="/data/project/lasseigne_lab/JordanWhitlock/230323_JW_DiseaseNetworks/data/human_prior_tri.csv"
echo "Prior Network from ${PRIOR_NET}"

#set minimum interaction number
MIN_N=5

# NOTE user must have already pulled from docker and built .sif file with singularity below (jordanwhitlock/setbp1_manuscript_panda_1.0.1)
#singularity exec --cleanenv --no-home -B ${wd} -B /data/project/lasseigne_lab/DATASET_dir/TCGA_GTEx_CCLE_HPA_220807/gtex/ ${wd}/bin/PANDA_construction_docker/jw_diseasenetworks_1.0.0.sif Rscript --vanilla ${wd}/src/tf_activity/01_decoupleR_analysis.R ${GTEX_FILE} ${PRIOR_NET} ${MIN_N} ${wd}
singularity exec --cleanenv --containall \
  -B ${wd} \
  -B /data/project/lasseigne_lab/ \
  ${wd}/bin/PANDA_construction_docker/jw_diseasenetworks_1.0.0.sif \
  Rscript --vanilla ${wd}/src/tf_activity/02_decoupleR_analysis.R ${GTEX_FILE} ${PRIOR_NET} ${MIN_N} ${wd}
