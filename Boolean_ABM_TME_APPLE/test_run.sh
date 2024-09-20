#!/bin/bash

make clean
echo "ATTEMPTING BUILD"
make
echo "STARTING TEST RUNS..."

for i in {52..100}; do
  echo "STARTING TEST RUN $i"
  ./runModel dev_test $i 6 1 1
done
