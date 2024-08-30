#!/bin/bash

#conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################
file1="/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_2023/car_t_2023_MK_100bp_tiled.rds"
file2="/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_2024/car_t_2024_new_MK_100bp_tiled.rds"
outName="/home/users/gsantamaria/projects/car_t/results/preprocessing/tile_comb/car_t_comb_100bp_tiled.rds"

# Analysis
########################################################################################################################
Rscript mergeTileFiles.R --file1 "$file1" --file2 "$file2" --outName "$outName"

# conda deactivate