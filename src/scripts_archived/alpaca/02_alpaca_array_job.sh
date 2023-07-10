#!/bin/bash
## run the Rscript alpaca_networks.R and schedule this job to SLURM with
## `sbatch 02_alpaca_array_job.sh`

#SBATCH --job-name=alpaca
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
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-465 #31 tissue files paired with each other to make 465 pairings of files

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
module load R
module load Singularity/3.5.2-GCC-5.4.0-2.26

wd=/data/user/jbarham3/230323_JW_DiseaseNetworks

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER='jbarham3'


cd "${wd}"

# Set the path to the file containing the pairs of filenames
file_pairs="${wd}/src/community_detection/file_pairs.txt"

# Read the line corresponding to the SLURM_ARRAY_TASK_ID from file_pairs.txt
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$file_pairs")

# Split the line into two variables, file1 and file2
IFS=' ' read -r -a files <<< "$line"
file1="${files[0]}"
file2="${files[1]}"

# Set the path to the folder containing the Rdata files
folder_path="${wd}/results/PANDA"
echo "Folder path: $folder_path"

singularity exec --cleanenv --no-home -B "$wd" "$wd/bin/downstream_analysis_docker/setbp1_manuscript_1.0.6.sif" Rscript --vanilla "$wd/src/community_detection/alpaca_networks_array.R" "$file1" "$file2"
