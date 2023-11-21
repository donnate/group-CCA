#!/bin/bash

# Define the values for the variables
#!/bin/bash

# Define the values for the variables
theta_strengths="high medium low"
n_values="200 500 800 1000 2000"
r_values="1 2 3 5 10"

for theta in $theta_strengths; do
  for n in $n_values; do
    for p in $p_values; do
      for r in $r_values; do
      sbatch elena/simu_group.sh "$n" "$theta" "$r
    #echo "$theta"
      done
    done
  done
done
