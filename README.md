## Boolean-TME-ABM: restructured

***NOTE:***
AThis model is under re-construction to more strictly follow OOP principles. Because the folks working on this use Apple machines, the restructuring will be first completed with Apple Silicon compatability and then for UNIX compatability. 



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


## Model Updates 
Cell Classes
- Made seperate header and cpp files for Cancer Cells, CD4 Cells, CD8 Cells and Macrophage Cells
    - Made new constructors for each type of cell
    - Moved all cell specific fuctions and parameters to the cpp files and header files 
- All the specific cell classes inherited all properties of the Cell class 
    - Included general functions and parameters that applied to all cells
        - ex: force functions, overlap functions, physical properties like radius, compressed, neighbors, etc
- All functions that were cell type specific were renamed when put in the new file
    - ex: prolifState() became cancer_prolifState()
    - others like cd4_differentiation only needed to be moved
- functions that were moved to specific cells were:
    - proliferate --> cancer, cd4, cd8, macrophage
    - age --> cancer, cd4, cd8, macrophage
    - prolifeState --> cancer, cd8
    - indirectInteractions --> cancer, cd8
    - directInteractions --> cancer, cd8
    - directInteractionProperties --> cd4, cd8, macrophage
    - differentiate --> cd4, macrophage
    - addChemotaxos --> cd8
    - dieFromCD8 --> cancer
    - gainPDL1 --> cancer
    - pdl1Inhibition --> cd8
- there was also specific variables that were added to cell header files
    - cancer --> pdl1Shift, pdl1WhenExpressed
    - cd4 --> pdl1WhenExpressed, probTh, kTr;
    - cd8 --> killProb, baseKillProb, infScale, t_cell_phenotype_Trajectory, pTypeStateTransition
    - macrophage --> kTr, kM1, kM2, plasticity, pdl1WhenExpressed

    ** all had mt added to their header file
- added getter and setters for the functions that needed them in the cpp/header files
    - getters and setters
        - CD4
            - pdl1WhenExpressed
            - probTh
            - kTr
        - CD8
            - setKillProb
            - baseKillProb
            - infScale
            - t_cell_phenotype_Trajectory
        - Macrophage
            - kTr
            - kM1
            - kM2
            - plasticity
            - pdl1WhenExpressed
Environment 
- In environment header file:
    - Split internalCellFunctions into:
        - internalCancerCellFunctions
        - internalCD4CellFunctions
        - internalCD8CellFunctions
        - internalMacrophageCellFunctions
    - Split cell_list into
        - cancer_list
        - cd4_list
        - cd8_list
        - macrophage_list
- In environmentInfo.cpp:
    - printStep
        - changed from iterating through the cell list and type checking to just to iterating through cancer_list, cd4_list, cd8_list and macrophage_list
    - updateTimeSeries
        - changed from iterating through the cell list and type checking to just iterating through cancer_list, cd4_list, cd8_list and macrophage_list
- In environmentLoadSave.cpp
    - save
        - made for loops so function could iterate through cancer_list, cd4_list, cd8_list and macrophage_list
            - used to just iterate through cell_list
- In environmentMain.cpp
    - changed to cancer_list instead of cell_list() with a type check
- In environmentMisc.cpp
    - initializeCells 
        - changed to pushback cancer cells onto cancer_list
        - the type was prevuiolsy 0 on the cell push backed on the cell list
    - recruitImmuneCells
        - added else if statements to the while loop to push back macrophages and cd4 cells onto their respective lists  
    - tumorSize
        - changed for loop to go through cancer_list not cell list
    - necrosis
        - changed for loop to go through cancer_list not cell list
- In environmentRunCells.cpp
    - neighborInfluenceInteractions
        - indirect
            - changed iterating through cell_lists to cancer_list, cd4_list, cd8_list and macrophage_list
            - then within each iteration looked at all the other cell lists again to look at all cells
            - only cancer cells and cd8 cells have indirectInteraction fucntions 
                - so the chunks of code that went from 0 to size of macrophage list or cd4 list were commented out 
        - direct
            - similar to above, changed iterating through cell_lists to cancer_list, cd4_list, cd8_list and macrophage_list and looked at cell types neighbors 
            - only cancer cells and cd8 cells have directInteraction fucntions 
                - so the chunks of code that went from 0 to size of macrophage list or cd4 list were commented out 
        - differentiate
            - similar to above, changed iterating through cell_lists to cancer_list, cd4_list, cd8_list and macrophage_list
            - only cd4 and macrophage have differentiation functions
                - commented out cancer and cd8 
    - calculateForces
        - migrate
            - split cell_list into cancer_list, cd4_list, cd8_list and macrophage_list
        - calculateForces
            - redid all for loops so there was a for loop for cancer_list, cd4_list, cd8_list and macrophage_list
            - then within each iteration looked at all the other cell lists again to look at all cells
        - resolveForces
            - split cell_list into ancer_list, cd4_list, cd8_list and macrophage_list
        - calculateOverlap
            - split cell_list into ancer_list, cd4_list, cd8_list and macrophage_list
    - internalCellFunctions
        - split up internalCellFunction to each of its specific cell types so that they could call their specific fucntion and lists 
    - made a dead vector for each type of cell

       