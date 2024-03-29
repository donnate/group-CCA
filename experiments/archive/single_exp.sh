#!/bin/bash
#SBATCH --job-name=single_experiment
#SBATCH --output=experiments/logs/single_experiment.out
#SBATCH --error=experiments/logs/single_experiment.err
#SBATCH --time=2:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=1G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load R/4.2.0
echo $1
echo $2
METH=$2
result_file="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "result file is ${result_file}"
cd $SCRATCH/group-CCA

if [[ "$METH" == "others" ]]; then
  Rscript experiments/run_experiment.R $1 "others" 
elif [[ "$METH" == "original-ChaoCCA" ]]; then
  Rscript experiments/run_experiment.R $1 $METH $3 $4 $5 50 
elif [[ "$METH" == "genChaoCCA" ]]; then
  Rscript experiments/run_experiment.R $1 $METH $3 $4 $5 50 $6
elif [[ "$METH" == "genCCA" ]]; then
  Rscript experiments/run_experiment.R $1 $METH $3 $4 $5 100 
fi;
