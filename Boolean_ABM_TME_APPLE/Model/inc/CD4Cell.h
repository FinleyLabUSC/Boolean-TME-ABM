#ifndef IMMUNE_MODEL_CELL_CD4_H
#define IMMUNE_MODEL_CELL_CD4_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include "Cell.h"


class CD4Cell: public Cell{
public:

    //constructor
    CD4Cell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp, std::array<double, 2> loc, int idx, 
            int cellType);
    
    // functions
    void cd4_differentiation(double dt);
    void cd4_age(double dt, size_t step_count);
    std::vector<double> cd4_directInteractionProperties(int interactingState, size_t step_count);
    void set_pdl1WhenExpressed(double value);
    double get_pdl1WhenExpressed();
    void set_probTh(double value);
    double get_probTh();
    void set_kTr(double value);
    double get_kTr();


private:
    std::mt19937 mt;
    double pdl1WhenExpressed; 
    double probTh; 
    double kTr;
};


#endif 