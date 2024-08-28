#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=liftYngbld
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 1
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/youngblood/output_liftYoungblood.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/youngblood/error_liftYoungblood.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1
# Variables for the pipeline
########################################################################################################################
toLiftDir="/home/users/gsantamaria/projects/car_t/data/youngblood/"
chain="/home/users/gsantamaria/projects/car_t/data/hg19ToHg38.over.chain"
outTag="hg38"
outDir="/home/users/gsantamaria/projects/car_t/data/youngblood_hg38Lift/"
num_jobs=10


# Run the lifting
########################################################################################################################
counter=0

pids=()

for file in "$toLiftDir"/*; do
  if [[ -f "$file" ]] && ! file "$file" | grep -q "tar archive"; then
    Rscript liftMethFiles.R "$file" --chain "$chain" --reorder --compress --outTag "$outTag" --outDir "$outDir" &
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
