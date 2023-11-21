#!/bin/bash

# Define the values for the variables
strengths="high medium low"
n_values="100 200 300 500 800 1000"

# Loop over the combinations and launch the sbatch command
for sig_strength in $strengths; do
  for n in $n_values; do
        sbatch elena/simu_group.sh "$n" "$theta"
  done
done