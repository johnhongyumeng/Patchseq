#!/bin/bash
#SBATCH --job-name=sst_nwb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=yk2513@nyu.edu
#SBATCH --mail-type=END
#SBATCH --array=1-600

sed -e 's/.$//' 0304_rename_nwb_sst_downloads_a.csv > 0304_linux_rename_nwb_sst_downloads_a.csv

File=$(awk "NR==$SLURM_ARRAY_TASK_ID" 0304_linux_rename_nwb_sst_downloads_a.csv)

mkdir SST_nwb
cd SST_nwb

wget -nc -O $File # download the files
echo $File

# tar -cf sst_nwb.tar.gz *.nwb ### compress all nwb files