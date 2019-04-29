#!/bin/bash
#SBATCH --partition fn_short
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 48G
#SBATCH --mem 48G
#SBATCH --time 0-10:00:00
#SBATCH --job-name hdbscan
#SBATCH --output hdbscan_parameter_search-%J.log

module unload python
cd /gpfs/data/ruggleslab/home/lmb529/phos_corr/
std=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" parameter_space.txt)
nmembers=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" parameter_space.txt)

source activate phos
python hdbscan_parameter_optimization.py --std ${std} --nmembers ${nmembers}
