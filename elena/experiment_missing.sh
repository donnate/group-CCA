#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=elena/logs/array_%A_%a.out
#SBATCH --error=elena/logs/array_%A_%a.err
#SBATCH --array=1-20
#SBATCH --time=35:00:00
#SBATCH --partition=caslake
#SBATCH --mem=20G
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
Rscript elena/missing/1.simulation.R $SLURM_ARRAY_TASK_ID $result_file $1
