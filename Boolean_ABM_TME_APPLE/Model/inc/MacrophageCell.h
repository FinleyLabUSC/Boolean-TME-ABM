#ifndef IMMUNE_MODEL_CELL_MACROPHAGE_H
#define IMMUNE_MODEL_CELL_MACROPHAGE_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include "Cell.h"


class MacrophageCell: public Cell{
public:

    //constructor 
    MacrophageCell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp, std::array<double, 2> loc, 
          int cellType);

    //functions
    std::array<double, 3> macrophage_proliferate(double dt);
    void macrophage_differentiation(double dt);
    void macrophage_age(double dt, size_t step_count);
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count) override; 
    
    void set_kTr(double value);
    double get_kTr();
    void set_kM1(double value);
    double get_kM1();
    void set_kM2(double value);
    double get_kM2();
    void set_plasticity(double value);
    double get_plasticity();
    void set_pdl1WhenExpressed(double value);
    double get_pdl1WhenExpressed();

     
    
    



private:
    std::mt19937 mt;
    
    // macrophage differentiation
    double kTr;
    double kM1;
    double kM2;
    double plasticity;
    double pdl1WhenExpressed;
};


#endif 