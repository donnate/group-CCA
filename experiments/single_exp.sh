#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --error=logs/array_%A_%a.err
#SBATCH --array=1-5
#SBATCH --time=35:00:00
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load R
echo $1
echo $2
METH=$2
result_file="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "result file is ${result_file}"
cd $SCRATCH/$USER/gnumap/experiments


if [[ "$METH" == "others" ]]; then
  Rscript experiments/run_experiment.R $1 "others" 
elif [[ "$METH" == "original-ChaoCCA" ]]; then
  Rscript experiments/run_experiment.R $1 $METH $3 $4 $5 20 
elif [[ "$METH" == "genChaoCCA" ]]; then
  Rscript experiments/run_experiment.R $1 $METH $3 $4 $5 20 $6
elif [[ "$METH" == "genCCA" ]]; then
  Rscript experiments/run_experiment.R $1 $METH $3 $4 $5 100 
fi;
