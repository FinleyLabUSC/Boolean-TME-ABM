#ifndef IMMUNE_MODEL_CELL_H
#define IMMUNE_MODEL_CELL_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include "ModelUtil.h"

class Cell{
public:
    /*
     * FUNCTIONS
     */

    // initialization
    Cell(std::array<double, 2> loc, int idx, std::vector<std::vector<double>> &cellParams, int cellType,
    std::vector<std::string> tCellPhenotypeTrajectory, size_t init_tstamp=0);
   
    // force functions
    std::array<double, 2> attractiveForce(std::array<double, 2> dx, double otherRadius);
    std::array<double, 2> repulsiveForce(std::array<double, 2> dx, double otherRadius);
    void calculateForces(std::array<double, 2> otherX, double otherRadius, int &otherType);
    void resolveForces(double dt, std::array<double, 2> &tumorCenter, double &necroticRadius, double &necroticForce);
    void resetForces();
    void neighboringCells(std::array<double, 2> otherX, int otherID);

    // overlap functions
    void calculateOverlap(std::array<double, 2> otherX, double otherRadius);
    void resetOverlap();
    void isCompressed();

    // cell behavior functions
    
    void inherit(std::vector<double> properties);
    std::vector<double> inheritanceProperties();
    void age(double dt, size_t step_count);
    void migrate(double dt, std::array<double, 2> tumorCenter);
    void indirectInteractions(double tstep);
    void directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep);
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count);

    // differentiation
    void differentiate(double dt);

    // cell influences
    void addInfluence(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void clearInfluence(); 

    // other functions
    double calcDistance(std::array<double, 2> otherX);
    double calcInfDistance(double dist, double xth);
    static double calcNorm(std::array<double, 2> dx);
    static std::array<double, 2> unitVector(std::array<double, 2> v);
    void updateID(int idx);

    /*
     * PARAMETERS
     */

    // location
    std::array<double, 2> x;

    // physical properties
    double radius;
    bool compressed;
    double currentOverlap;
    std::vector<int> neighbors;

    // age, division, and lifespan
    double divProb;
    double divProb_base;
    double deathProb;
    bool canProlif;

    // force properties
    double mu;
    double kc;
    double damping;
    double maxOverlap;
    double rmax;
    std::array<double, 2> currentForces;

    // migration
    double migrationSpeed;
    double migrationBias;


    // interactions with other cells
    double influenceRadius; 
    double pdl1;  
    std::array<double, 8> influences; 
    std::array<double, 8> chemotaxVals;
    

    // identification
    int id;
    int type;
    int state;

    //lifespan
    size_t init_time; 

private:
    std::mt19937 mt;
};

#endif //IMMUNE_MODEL_CELL_H
