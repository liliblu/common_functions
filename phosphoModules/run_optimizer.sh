#!/bin/bash
#SBATCH --partition=fn_short
#SBATCH --job-name=phosphoModules
#SBATCH --output phosModules-%J.log
#SBATCH --time 02:00:00
#SBATCH --mem 50G

std=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" stds.txt)

module unload python
source activate phos
python run_optimizer.py --std $std --nmembers 2 4 6 8 10 12 14 16
