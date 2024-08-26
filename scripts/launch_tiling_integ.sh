#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=tilingInteg
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 4
#SBATCH --time=02-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/output_tiling_integ.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/error_tiling_integ.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

#seuratFile="/home/users/gsantamaria/projects/aging/results/integration/blood_int_cca.rds"
#resolution=2
#nPCs=100

Rscript tiling_integ.R


conda deactivate