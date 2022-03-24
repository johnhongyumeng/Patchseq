#!/bin/bash
#SBATCH --job-name=0309_pv_AIHW_loop
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --mail-user=yk2513@nyu.edu
#SBATCH --mail-type=END

module purge
module load python/intel/3.8.6

python ./0309_AllenAIHW_Loop_PV.py