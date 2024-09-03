## Boolean-TME-ABM: restructured

***NOTE:***
This model is under re-construction to more strictly follow OOP principles. Because the folks working on this use Apple machines, the restructuring will be first completed with Apple Silicon compatability and then for UNIX compatability. 



## Prerequisites

This repository contains a C++ program that builds using CMake and Make.

Before running the program, ensure you have the following installed on your system:

- CMake (https://cmake.org/download/)
- Make (Usually available by default on Unix-like systems)
- C++ compiler (such as g++ for Linux or macOS, or Visual C++ for Windows)
- A Python Interpreter 

If you are using **Apple silicon**, please ensure you have `libomp` installed. You can install it via Homebrew:

```bash
$ brew install libomp
```
## Build and Run Instructions

Follow these steps to build and run the program:

1. Clone this repository to your local machine:
2. If you are using Apple silicon please work with the model contained in _APPLE or _UNIX if you are using an alternative system. 


```bash
$ git clone https://github.com/FinleyLabUSC/Boolean-TME-ABM.git
$ cd <_APPLE or _UNIX_>
$ cmake .
$ make 
$ pip install -r requirements.txt 
```
 
3. This will build an executable runModel and install the necessary python packages to generate the parameter file. To run, you can call ./runModel as follows

```bash
$ ./runModel <SAVE_FLD> <SAVE_DIR> <P_TYPE_STATE_TRANSITION> <DEATH_PROBABILITY_FACTOR> <KILL_PROBABILITY_FACTOR>
```
As an example:
```bash
$ ./runModel modelPredictions 0 3 2 2
```

## Notes
As of now, the entirety of the build is not seeded and, as such, running the model with the same parameter set will produce different results. 