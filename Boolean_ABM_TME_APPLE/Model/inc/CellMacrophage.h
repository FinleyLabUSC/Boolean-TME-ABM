#ifndef IMMUNE_MODEL_CELL_MACROPHAGE_H
#define IMMUNE_MODEL_CELL_MACROPHAGE_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include "Cell.h"


class CellMacrophage: public Cell{
public:

    //constructor 
    CellMacrophage(std::vector<std::vector<double>> &cellParams, size_t init_tstamp=0);

    //functions
    void macrophage_differentiation(double dt);
    void macrophage_age(double dt, size_t step_count);
    std::vector<double> macrophage_directInteractionProperties(int interactingState, size_t step_count);



private:
    std::mt19937 mt;
};


#endif 