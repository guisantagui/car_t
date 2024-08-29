#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=tiling2024
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300GB
#SBATCH -c 4
#SBATCH --time=01-00:30:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/preprocessing/tiles/output_methyl_tiling_2024.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/preprocessing/tiles/error_methyl_tiling_2024.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

urlsFile="/home/users/gsantamaria/projects/car_t/data/PrjCellEx2024_COV_links.txt"
wgbsDir="/home/users/gsantamaria/projects/car_t/data/epigen/bismark_res/2024/new/"
tileSize=100
minCov=5
minCpGs_inTile=5
outName="car_t_2024_new"
outDir="/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_2024/"

#Rscript dwnld_files_from_tab.R $urlsFile --expName $outName --outDir $wgbsDir

Rscript methyl_tiling_exe.R --wgbsDir $wgbsDir --tileSize $tileSize --minCov $minCov --minCpGs_inTile $minCpGs_inTile --outName $outName --outDir $outDir

conda deactivate