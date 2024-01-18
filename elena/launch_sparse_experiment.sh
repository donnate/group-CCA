#!/bin/bash

# Define the values for the variables
theta_strengths="high medium low"
n_values="100 200 500 1000 2000"
p_values="80 100 500 800 1000"
r_values="3"

for theta in $theta_strengths; do
  for n in $n_values; do
    for p in $p_values; do
      for r in $r_values; do
      sbatch elena/sparse_experiment.sh "$n" "$theta" "$p" "$r"
    #echo "$theta"
      done
    done
  done
done
