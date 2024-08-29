#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=getPercMatComb
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 1
#SBATCH --time=00-05:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/preprocessing/tiles/output_getPercMatMan_comb.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/preprocessing/tiles/error_getPercMatMan_comb.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

tileFile="/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_comb/car_t_comb_100bp_tiled.rds"
outName="/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_comb/comb_methPerc_man.csv"

Rscript getPercMatManual.R "$tileFile" --outName "$outName"