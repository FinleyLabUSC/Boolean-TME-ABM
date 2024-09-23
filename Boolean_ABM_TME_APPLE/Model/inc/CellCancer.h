#ifndef IMMUNE_MODEL_CELL_CANCER_H
#define IMMUNE_MODEL_CELL_CANCER_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include "Cell.h"



class CellCancer: public Cell{
public:

    // constructor 
    CellCancer(std::vector<std::vector<double>> &cellParams, size_t init_tstamp=0);

    //functions
    std::array<double, 3> cancer_proliferate(double dt);
    void cancer_prolifState();
    void cancer_dieFromCD8(std::array<double, 2> otherX, double otherRadius, double kp, double dt);
    void cancer_gainPDL1(double dt);
    void cancer_age(double dt, size_t step_count);
    void cancer_indirectInteractions(double tstep);
    void cancer_directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep);
    std::vector<double> cancer_directInteractionProperties(int interactingState, size_t step_count);

  





     // variables
    double pdl1Shift;
    
private:
    std::mt19937 mt;
};
    

#endif 