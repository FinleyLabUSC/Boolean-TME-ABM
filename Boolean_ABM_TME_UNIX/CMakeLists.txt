cmake_minimum_required(VERSION 3.0)
project(Tumor_Immune_CBM)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-fopenmp -std=c++17 -O3")

include_directories(Model/inc)
file(GLOB SOURCES "Model/src/*.cpp")

add_executable(runModel ${SOURCES})
