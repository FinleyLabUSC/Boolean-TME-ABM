#!/bin/bash

# cmake artifacts
rm -r CMakeFiles/
rm cmake_install.cmake
rm CMakeCache.txt
rm Makefile

cmake .

# Clean and make the project
make clean
make

# Define the base parameters
param2=6
param3=1
param4=1
param5=1

# Array of multipliers for the last three parameters
multipliers=(1 5 10)

# Counter for param1
param1=1

# Set the number of iterations (how many times to increment param1)
iterations=10  # You can set this to any number you want

# Loop for the specified number of iterations
for ((i=1; i<=iterations; i++)); do
  # Iterate over each multiplier for the last three parameters
  for multiplier in "${multipliers[@]}"; do
    # Modify the last three parameters based on the multiplier
    new_param4=$((param4 * multiplier))
    new_param5=$((param5 * multiplier))
    new_param6=$((param3 * multiplier))

    # Run the model with the current combination of parameters
    echo "Running model with param1=$param1, param2=$param2, param4=$new_param4, param5=$new_param5, param6=$new_param6"
    ./runModel test_dir "$param1" "$param2" "$new_param4" "$new_param5" "$new_param6"
  done

  # Increment param1 after each iteration
  ((param1++))
done
