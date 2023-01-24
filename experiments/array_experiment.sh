#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --error=logs/array_%A_%a.err
#SBATCH --array=1-5
#SBATCH --time=35:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load R/4.2.0
echo $1
echo $2
id_experiment="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
name_experiment="$3-$1-$2-$4-${id_experiment}-$6"
echo "result file is ${result_file}"
cd $SCRATCH/$USER/group-CCA


# Run one experiment  to create the dataset
Rscript experiments/main_experiment.R $1 $2 $3 $4 $5 $6

#### Launch the analysis with the different methods
sh experiments/single_exp.sh $name_experiment others 

#### Launch the loop for genCCA

for lambda1 in 0.001 0.01 0.1 1 10
do 
  for lambda2 in 0.001 0.01 0.1 1 10
  do 
    for lambda3 in 0.0 0.001 0.01
    do 
       sh experiments/single_exp.sh $name_experiment genCCA $lambda1 $lambda2 $lambda3
       sh experiments/single_exp.sh $name_experiment genChaoCCA $lambda1 $lambda2 $lambda3 smooth
       sh experiments/single_exp.sh $name_experiment genChaoCCA $lambda1 $lambda2 $lambda3 gen
       sh experiments/single_exp.sh $name_experiment original-ChaoCCA $lambda1 $lambda2 $lambda3 
    done
  done
done
