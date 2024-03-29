#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=experiments/logs/array_%A_%a.out
#SBATCH --error=experiments/logs/array_%A_%a.err
#SBATCH --array=1-20
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
name_experiment="$3-$1-$2-$4-${id_experiment}"
echo "name experiment is ${name_experiment}"
cd $SCRATCH/group-CCA


# Run one experiment  to create the dataset
Rscript experiments/main_experiment.R $1 $2 $3 $4 ${id_experiment} ${SLURM_ARRAY_TASK_ID}
wait
#### Launch the analysis with the different methods
sbatch -e "experiments/logs/others_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" experiments/single_exp.sh $name_experiment others 
#### Launch the loop for genCCA

for lambda1 in "0.001" "0.005" "0.01" "0.05" "0.1" "0.5" "1" "10"
do 
  for lambda2 in "0.001" "0.01" "0.1" "1" "10"
  do 
    for lambda3 in 0
    do 
       sbatch -e "experiments/logs/genCCA_${lambda1}_${lambda2}_${lambda3}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" experiments/single_exp.sh $name_experiment genCCA $lambda1 $lambda2 $lambda3
       sbatch -e "experiments/logs/genChao_smooth_${lambda1}_${lambda2}_${lambda3}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" experiments/single_exp.sh $name_experiment genChaoCCA $lambda1 $lambda2 $lambda3 smooth
       sbatch -e "experiments/logs/genChao_sparse_${lambda1}_${lambda2}_${lambda3}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" experiments/single_exp.sh $name_experiment genChaoCCA $lambda1 $lambda2 $lambda3 GEN
       sbatch -e "experiments/logs/Chao_${lambda1}_${lambda2}_${lambda3}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" experiments/single_exp.sh $name_experiment original-ChaoCCA $lambda1 $lambda2 $lambda3 
    done
  done
done
