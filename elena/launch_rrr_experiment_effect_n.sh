#!/bin/bash

# Define the values for the variables
theta_strengths="high medium low"
n_values="100 300 500 1000 2000 10000"
#p_values="80 100 300 500 800"
p_values="500"
q_values="10 50 80"
r_values="2"

for theta in $theta_strengths; do
  for n in $n_values; do
    for p in $p_values; do
      for r in $r_values; do
	  for q in $q_values; do
      sbatch elena/rrr_experiment.sh "$n" "$theta" "$p" "$r" "$q"
    #echo "$theta"
    done
      done
    done
  done
done
~      
