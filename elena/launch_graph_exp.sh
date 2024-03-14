#!/bin/bash

# Define the values for the variables
#!/bin/bash

# Define the values for the variables
theta_strengths="high medium low"
#in_values="500"
r_values="2"
#p_values="100 300 500 700"
p_values="10 15 20 30 40 50"

for theta in $theta_strengths; do
  for p in $p_values; do
      for r in $r_values; do
      sbatch elena/simu_graph.sh "$theta" "$r" "$p"
    done
  done
done
