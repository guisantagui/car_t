#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=getYngbld
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4GB
#SBATCH -c 1
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/getData/output_getYoungblood.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/getData/error_getYoungblood.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

cd /home/users/gsantamaria/projects/car_t/data/youngblood/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188325/suppl/GSE188325_RAW.tar