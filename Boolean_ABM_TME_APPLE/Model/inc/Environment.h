#ifndef IMMUNE_MODEL_ENVIRONMENT_H
#define IMMUNE_MODEL_ENVIRONMENT_H

#include <vector>
#include <algorithm>
#include <random>
#include "Cell.h"
#include <iostream>
#include <iomanip>      
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem> 
#include <array>
#include "/opt/homebrew/opt/libomp/include/omp.h"




class Environment{

public:
    Environment(std::string folder, std::string set, std::string tCellTrajectoryPath);
    //destructor needed
    void simulate(double tstep);

private:
    void runCells(double tstep, size_t step_count);
    void neighborInfluenceInteractions(double tstep, size_t step_count);
    void internalCellFunctions(double tstep, size_t step_count);
    void recruitImmuneCells(double tstep, size_t step_count);
    std::array<double, 2> recruitmentLocation();
    void tumorSize();
    void necrosis(double tstep);
    double calculateDiffusibles(std::array<double, 2> x);

    void save(double tstep, double tstamp);
    void loadParams();

    void initializeCells();
    void calculateForces(double tstep);
    void updateTimeSeries();

    void printStep(double time);
    
    double dt;

    // cell lists
    std::vector<Cell> cell_list;

    // time courses
    std::vector<int> cancerTS;
    std::vector<int> cd8TS;
    std::vector<int> cd4TS;
    std::vector<int> m0TS;
    std::vector<int> m1TS;
    std::vector<int> m2TS;
    std::vector<int> radiusTS;

    /*
    t cell trajectory matrix: we can either represent as a vector of chars 
    where a char maps to a phenotypic state or an int where the int maps to 
    a phenotypic state
    */
    std::vector<std::string> tCellPhenotypeTrajectory_1;

    std::vector<std::vector<std::string>> tCellPhenotypeTrajectory; 

    

    // parameter lists
    std::vector<std::vector<double>> cellParams;
    std::vector<double> recParams;
    std::vector<double> envParams;

    std::string saveDir;
    int steps;
    std::vector<double> immuneCellRecRates;
    std::vector<int> immuneCellRecTypes;
    std::vector<double> immuneCells2rec;
    double recDist;
    double maxRecCytoConc;
    double tumorRadius;
    double necroticGrowth;
    double necroticRadius;
    double necroticForce;
    double necroticLimit;
    std::array<double, 2> tumorCenter;
    double recruitmentDelay;

    // environment params
    double simulationDuration;
    int day;

    std::mt19937 mt;
};

#endif //IMMUNE_MODEL_ENVIRONMENT_H
