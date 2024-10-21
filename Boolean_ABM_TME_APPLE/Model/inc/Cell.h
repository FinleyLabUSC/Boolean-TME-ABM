#ifndef IMMUNE_MODEL_CELL_H
#define IMMUNE_MODEL_CELL_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>

class Cell{
public:
    /*
     * FUNCTIONS
     */

    // initialization
    Cell(std::array<double, 2> loc, int idx, std::vector<std::vector<double>> &cellParams, int cellType,
    std::vector<std::string> tCellPhenotypeTrajectory, size_t init_tstamp=0);
    void initializeCancerCell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp=0);
    void initializeCD8Cell(std::vector<std::vector<double>> &cellParams, std::vector<std::string> phenotypeTrajectory, size_t init_tstamp);
    void initializeCD4Cell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp=0);
    void initializeMacrophageCell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp=0);

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
    std::array<double, 3> proliferate(double dt);
    void prolifState();
    void inherit(std::vector<double> properties);
    std::vector<double> inheritanceProperties();
    void age(double dt, size_t step_count);
    void migrate(double dt, std::array<double, 2> tumorCenter, double tumorRadius);
    void indirectInteractions(double tstep, size_t step_count);
    void directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep);
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count);

    // differentiation
    void differentiate(double dt);

    // cell influences
    void addInfluence(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void addChemotaxis(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void clearInfluence();

    // macrophage
    void macrophage_differentiation(double dt);

    // CD4 specific
    void cd4_differentiation(double dt);

    // CD8 specific
    void cd8_setKillProb(size_t step_count);
    void cd8_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt);

    // cancer specific
    void cancer_dieFromCD8(std::array<double, 2> otherX, double otherRadius, double kp, double dt);
    void cancer_gainPDL1(double dt);

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
    int pTypeStateTransition; // for CD8+ T cells only

    // force properties
    double mu;
    double kc;
    double damping;
    double maxOverlap;
    double rmax;
    std::array<double, 2> currentForces;

    // migration
    double migrationSpeed;
    double migrationSpeed_inTumor; 
    double migrationBias;
    double migrationBias_inTumor; 

    // cancer properties
    double pdl1Shift;

    // interactions with other cells
    double influenceRadius;
    double pdl1;
    double pdl1WhenExpressed;
    std::array<double, 8> influences;
    std::array<double, 8> chemotaxVals;
    double probTh;

    // differentiation
    double kTr;
    double kM1;
    double kM2;
    double plasticity;

    // T cell killing
    double killProb;
    double baseKillProb;
    double infScale;
    std::vector<std::string> t_cell_phenotype_Trajectory; 

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
