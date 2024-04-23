#!/bin/bash

make clean
echo "ATTEMPTING BUILD"
make
echo "STARTING TEST RUN..."
./runModel dev_test 0 6 1 1
