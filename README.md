# Boolean-TME-ABM


# Running a Program with CMake and Make

This repository contains a C++ program that builds using CMake and Make. 


## Prerequisites

Before running the program, ensure you have the following installed on your system:

- CMake (https://cmake.org/download/)
- Make (Usually available by default on Unix-like systems)
- C++ compiler (such as g++ for Linux or macOS, or Visual C++ for Windows)

## Build and Run Instructions

Follow these steps to build and run the program:

1. Clone this repository to your local machine:
2. If you are using Apple silicon please work with the model contained in _APPLE or _UNIX if you are using an alternative system. 


```bash
$ git clone https://github.com/FinleyLabUSC/Boolean-TME-ABM.git
$ cd <_APPLE or _UNIX_>
$ cmake .
$ make 

```
 
3. This will build an executable runModel. To run, you can call ./runModel as follows

```bash
$ ./runModel <SAVE_FLD> <SAVE_IDX> <P_TYPE_STATE_TRANSITION> <DEATH_PROBABILITY_FACTOR> <KILL_PROBABILITY_FACTOR>
```
As an example:
```bash
$ ./runModel modelPredictions 0 3 2 2
```
