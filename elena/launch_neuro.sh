#!/bin/bash

# Define the values for the variables
#!/bin/bash

# Define the values for the variables
test_folds="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16"
lambda_values="0.01 0.1 1 10"

for lambda in $lambda_values; do
  for test_fold in $test_folds; do
      sbatch elena/neuro_exp.sh "$lambda" "$test_fold"
  done
done