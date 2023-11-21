#!/bin/bash

# Define the values for the variables
#theta_strengths="high medium low"
#n_values="100 200 500 1000 2000 10000"
#theta_strengths="low"
n_values="1000 2000 10000"
theta_strengths="high medium"
for theta in $theta_strengths; do
  for n in $n_values; do
      sbatch elena/rrr_experiment.sh "$n" "$theta"
    #echo "$theta"
    done
done

