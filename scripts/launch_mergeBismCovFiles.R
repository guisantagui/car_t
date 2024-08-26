#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=methTiling
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 4
#SBATCH --time=01-00:30:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/output_methyl_tiling.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/error_methyl_tiling.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

urlsFile1="/home/users/gsantamaria/projects/car_t/data/PrjCellEx2024_COV_links.txt"
expName1="car_t_2023"
urlsFile2="/home/users/gsantamaria/projects/car_t/data/Tcells2023WGBSMeta_Links_collab.xlsx"
wgbsDir="/home/users/gsantamaria/projects/car_t/data/wgbs/2024/"
tileSize=100
minCov=5
minCpGs_inTile=3
outName="car_t_2024"
outDir="/home/users/gsantamaria/projects/car_t/results/preprocessing"

# Download 2023 files and create 2023 samp_info file
Rscript dwnld_files_from_tab.R $urlsFile --expName $outName --outDir $wgbsDir

# Download 2023 files and create 2023 samp_info file
Rscript dwnld_files_from_tab.R $urlsFile --expName $outName --outDir $wgbsDir

Rscript methyl_tiling_exe.R --wgbsDir $wgbsDir --tileSize $tileSize --minCov $minCov --minCpGs_inTile $minCpGs_inTile --outName $outName --outDir $outDir

conda deactivate