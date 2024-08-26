#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=OPLS_RFE
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH --time=01-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/output_OPLS_RFE.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/error_OPLS_RFE.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal


conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

input="/home/users/gsantamaria/projects/car_t/results/preprocessing/percMeth_4OPLS.csv"
sampInfo="/home/users/gsantamaria/projects/car_t/results/preprocessing/sample_info.csv"
respVar='exh_score'
outDir="/home/users/gsantamaria/projects/car_t/results/OPLS_RFE/"
#splitTrainTest='0.66'
splitTrainTest='no'
filtSampsByDisPhen='no'

Rscript OPLS_RFE.R $input --sampMetadata $sampInfo --respVar $respVar --outDir $outDir --splitTrainTest $splitTrainTest --filtSampsByDisPhen $filtSampsByDisPhen

conda deactivate