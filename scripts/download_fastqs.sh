#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=dl_fastqs
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4GB
#SBATCH -c 1
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/car_t/scripts/output_download_fastqs.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/car_t/scripts/error_download_fastqs.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate python-3.8

# Define variables
SFTP_USER="andcip91_stanford"
SFTP_SERVER="sftp.genewiz.com"
REMOTE_DIR="./30-970588504/rerun/samples/"
SAMPLE="CD4-CD19-21-67" # This is to be edited for each sample
LOCAL_DIR="/home/users/gsantamaria/projects/car_t/data/epigen/fastqs/2024/"
SFTP_PASSWORD="DDDPvF6AIhvmb5fmgWvi"

# Create an lftp script dynamically
cat <<EOF > sftp_download.lftp
set sftp:auto-confirm yes
open -u $SFTP_USER,$SFTP_PASSWORD sftp://$SFTP_SERVER
lcd $LOCAL_DIR
cd $REMOTE_DIR
mget $SAMPLE
bye
EOF

# Run the lftp script
lftp -f sftp_download.lftp

# Clean up
rm sftp_download.lftp

conda deactivate