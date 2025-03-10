#!/bin/bash

# cmake artifacts
rm -r CMakeFiles/
rm cmake_install.cmake
rm CMakeCache.txt
rm Makefile

cmake .


make clean
make

./runModel dev_test 0 6 1 1 1

