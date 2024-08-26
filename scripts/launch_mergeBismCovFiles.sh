#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=mergBismCov
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 2
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/output_mergeBismCovFiles.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/error_mergeBismCovFiles.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

covDir1="/home/users/gsantamaria/projects/car_t/data/epigen/bismark_res/2024/new/"
covDir2="/home/users/gsantamaria/projects/car_t/data/wgbs/2024/"
outDir="/home/users/gsantamaria/projects/car_t/data/wgbs/2024_comb/"

# Iterate over the files and run the script for merging them

# Maximum number of parallel processes
max_jobs=11

# Counter to keep track of background processes
job_count=0

covDir1_files=($(ls "$covDir1"))

for file in "${covDir1_files[@]}"; do
  cov1="${covDir1}${file}"
  cov2="${covDir2}${file}"
  outName="${outDir}${file}"
  Rscript mergeBismCovFiles.R --cov1 "$cov1" --cov2 "$cov2" --outName "$outName" &
  # Increment the job counter
  ((job_count++))
  # If we have reached the max number of jobs, wait for them to complete
  if [[ $job_count -ge $max_jobs ]]; then
    wait
    job_count=0
  fi
done

wait

echo "All bismark cov files merged and saved in ${outDir}."

conda deactivate