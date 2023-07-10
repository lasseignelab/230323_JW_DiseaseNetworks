#!/bin/bash
## run the Rscript alpaca_networks.R and schedule this job to SLURM with
## `sbatch 02_condor_array_job.sh`

#SBATCH --job-name=condor
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jbarham3@uab.edu #CHANGE THE EMAIL TO YOUR OWN
#
#Number of tasks needed for this job, the partition, and parameters
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=255000
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=express #partition info here: https://docs.rc.uab.edu/cheaha/hardware/#partitions
#
#SBATCH --output=%x_%A_%a.condor.out
#SBATCH --error=%x_%A_%a.condor.err
#SBATCH --array=1-31 # 31 tissue files

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
module load R
module load Singularity/3.5.2-GCC-5.4.0-2.26

wd=/data/user/jbarham3/230323_JW_DiseaseNetworks

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER='jbarham3'


cd "${wd}"

#array file of PANDA network inputs for CONDOR
SAMPLE_LIST="${wd}/results/array_inputs/condor_array.txt" #note: make sure path and file name are correct
SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`) # parantheses instruct bash to create a shell array of strings from SAMPLE_LIST
INPUT=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}` #extracts a single input from the array and prints (using echo) it into INPUT variable, each single input is then assigned an array number by $SLURM_TASK_ID


# Set the path to the folder containing the Rdata files
folder_path="${wd}/results/PANDA"
echo "Folder path: $folder_path"

singularity exec --cleanenv --no-home -B "$wd" "$wd/bin/downstream_analysis_docker/setbp1_manuscript_1.0.6.sif" Rscript --vanilla "$wd/src/community_detection/02_condor_networks_array.R" "$INPUT"
