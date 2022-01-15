#!/bin/bash
#SBATCH --job-name=pv_nwb
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --mail-user=yk2513@nyu.edu
#SBATCH --mail-type=END
#SBATCH --array=1-145

File=$(awk "NR==$SLURM_ARRAY_TASK_ID" nwb_pv_L23_downloads.csv)

mkdir PV_nwb
cd PV_nwb

wget $File # download the files
echo $File
tar -cf pv_nwb.tar.gz *.nwb%0D ### compress all nwb files