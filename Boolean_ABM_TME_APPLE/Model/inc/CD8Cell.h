#ifndef IMMUNE_MODEL_CELL_CD8_H
#define IMMUNE_MODEL_CELL_CD8_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include "Cell.h"

class CD8Cell: public Cell{
public:

    //constructor
    CD8Cell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp, std::array<double, 2> loc, int idx, 
            int cellType,std::vector<std::string> phenotypeTrajectory);

    //variables
    int pTypeStateTransition; // for CD8+ T cells only

    // functions
    // CD8 specific
    void cd8_addChemotaxis(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void cd8_setKillProb();
    void cd8_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt);
    void cd8_age(double dt, size_t step_count);
    void cd8_indirectInteractions(double tstep);
    void cd8_directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep);
    std::vector<double> cd8_directInteractionProperties(int interactingState, size_t step_count);

    void set_killProb(double value);
    double get_killProb();

    void set_baseKillProb(double value);
    double get_baseKillProb();

    void set_infScale(double value);
    double get_infScale();

    void set_t_cell_phenotype_Trajectory(std::vector<std::string> value);
    std::vector<std::string> get_t_cell_phenotype_Trajectory();

    



private:
    std::mt19937 mt;

    // T cell killing (cd8)
    double killProb;
    double baseKillProb;
    double infScale;
    std::vector<std::string> t_cell_phenotype_Trajectory; 


};

#endif 
