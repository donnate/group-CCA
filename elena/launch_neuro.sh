#!/bin/bash

# Define the values for the variables
#!/bin/bash

# Define the values for the variables
test_folds="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
lambda_values="1e-5 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1"

for lambda in $lambda_values; do
      sbatch elena/neuro_exp.sh "$lambda"
done
