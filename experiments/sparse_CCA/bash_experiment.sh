#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=experiments/sparse_CCA/logs/array_%A_%a.out
#SBATCH --error=experiments/sparse_CCA/logs/array_%A_%a.err
#SBATCH --array=1-50
#SBATCH --time=35:00:00
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
#module load libgmp
module load R/4.2.0
module load python

#source activate "r-reticulate"

result_file="new_exp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.csv"
echo "result file is ${result_file}"
cd $SCRATCH/$USER/group-CCA/
Rscript experiments/sparse_CCA/experiment_sparse_CCA.R $SLURM_ARRAY_TASK_ID $result_file $1 $2 $3 $4 $5
# $1 : N
# $2 : r
# $3 : r_pcas
# $4 : criterion (prediction/ correlation) for CV
# $5 : normalized diagonal (0/1)
