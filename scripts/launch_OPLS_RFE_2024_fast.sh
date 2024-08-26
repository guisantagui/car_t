#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=OPLS_RFE_fst
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH -c 1
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/output_OPLS_RFE_fast.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/error_OPLS_RFE_fast.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal


conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

input="/home/users/gsantamaria/projects/car_t/results/preprocessing/car_t_2024_100bp_percMeth.csv"
sampInfo="/home/users/gsantamaria/projects/car_t/results/preprocessing/sampInfo_2024.csv"
respVar='exh_score'
worseFeatsProp=0.002
varThrshld=100
outDir="/home/users/gsantamaria/projects/car_t/results/OPLS_RFE_2024_fast/"
#splitTrainTest='0.66'
splitTrainTest='no'
filtSampsByDisPhen='no'

Rscript OPLS_RFE.R $input --sampMetadata $sampInfo --respVar $respVar --worseFeatsProp $worseFeatsProp --varThrshld $varThrshld --outDir $outDir --splitTrainTest $splitTrainTest --filtSampsByDisPhen $filtSampsByDisPhen

conda deactivate