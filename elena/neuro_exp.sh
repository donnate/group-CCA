#!/bin/bash

#SBATCH --job-name=neuro_array
#SBATCH --output=elena/group_sparse/logs/neuro_array_%A_%a.out
#SBATCH --error=elena/group_sparse/logs/neuro_array_%A_%a.err
#SBATCH --array=1-16
#SBATCH --time=24:00:00
#SBATCH --partition=caslake
#SBATCH --mem=15G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load R/4.2.0

echo "result file is ${result_file}"
cd $SCRATCH/$USER/group-CCA/
Rscript neuroscience_elena.R $1 $SLURM_ARRAY_TASK_ID
