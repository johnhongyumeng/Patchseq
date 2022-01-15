#!/bin/bash
#SBATCH --job-name=sst_nwb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --mail-user=yk2513@nyu.edu
#SBATCH --mail-type=END
#SBATCH --array=1-172

File=$(awk "NR==$SLURM_ARRAY_TASK_ID" nwb_sst_L23_downloads.csv)

mkdir SST_nwb
cd SST_nwb
wget $File # download the files

tar -cf sst_nwb.tar.gz *.nwb%0D ### compress all nwb files