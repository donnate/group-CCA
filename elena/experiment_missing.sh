#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=experiments/sparse_CCA/logs/array_%A_%a.out
#SBATCH --error=experiments/sparse_CCA/logs/array_%A_%a.err
#SBATCH --array=1-2
#SBATCH --time=6:00:00
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load R/4.2.0


result_file="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "result file is ${result_file}"
cd $SCRATCH/$USER/group-CCA/
Rscript elena/missing/1.simulation.R $SLURM_ARRAY_TASK_ID $result_file