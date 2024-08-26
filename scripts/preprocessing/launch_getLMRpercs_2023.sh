#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=getLMRpercs
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 8
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/preprocessing/output_launch_getLMRpercs_2023.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/preprocessing/error_launch_getLMRpercs_2023.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

bismFileDir="/home/users/gsantamaria/projects/car_t/data/wgbs/"
lmrFile="/home/users/gsantamaria/projects/car_t/data/d59_d61_d62_M0.hmr"
minCov=5
outTag="2023"
outDir="/home/users/gsantamaria/projects/car_t/results/preprocessing/LMRs/2023/"
num_jobs=10

counter=0

pids=()

for file in "$bismFileDir"/*; do
  if [[ -f "$file" ]]; then
    Rscript getLMRpercs.R "$file" --lmrFile "$lmrFile" --minCov "$minCov" --outDir "$outDir" &
    pids+=($!)

    ((counter++))

    if ((counter >= num_jobs)); then
      # Wait for all background jobs in the array to complete
      for pid in "${pids[@]}"; do
        wait "$pid"
      done
      # Clear the array
      pids=()
      # Reset the counter
      counter=0
    fi

  fi
done

for pid in "${pids[@]}"; do
  wait "$pid"
done

Rscript mergeLMRfiles.R $outDir --outTag $outTag --outDir $outDir

conda deactivate