#!/bin/bash

echo "Starting test run..."
make clean
make
./runModel dev_test 0 6 1 1
