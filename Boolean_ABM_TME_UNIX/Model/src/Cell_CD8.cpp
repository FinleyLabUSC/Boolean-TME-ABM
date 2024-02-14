#include "Cell.h"
#include "ModelUtil.h"

void Cell::initializeCD8Cell(std::vector<std::vector<double> > &cellParams, std::vector<std::string> phenotypeTrajectory, size_t init_tstamp) {
    state = 6;

    mu = cellParams[0][2];
    kc = cellParams[1][2];
    damping = cellParams[2][2];
    maxOverlap = cellParams[3][2]*cellParams[4][2];
    radius = cellParams[4][2]/2.0;
    deathProb = cellParams[5][2];
    migrationSpeed = cellParams[6][2];
    killProb = cellParams[7][2];
    infScale = cellParams[8][2];
    influenceRadius = cellParams[9][2];
    migrationBias = cellParams[10][2];
    divProb_base = cellParams[11][2];
    pTypeStateTransition = cellParams[12][2]; 
    rmax = 1.5*radius*2;

    if(phenotypeTrajectory.size() == 0 || phenotypeTrajectory.empty()){
        std::cerr << "WARNING CONSTRUCTOR: t_cell_phenotype_Trajectory is empty!" << std::endl; 
    }
    
    
    t_cell_phenotype_Trajectory = phenotypeTrajectory; 
    init_time = init_tstamp;
}

void Cell::cd8_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt) {
    // inhibition via direct contact

    if(state != 6){return;}

    double distance = calcDistance(otherX);
    if(distance <= radius+otherRadius){
        std::uniform_real_distribution<double> dis(0.0,1.0);
        if(dis(mt) < otherpdl1){
            state = 7;
            killProb = 0;
            migrationSpeed = 0.0;
        }
    }
}

void Cell::cd8_setKillProb(){
    // set kill prob based on influence
    // essentially effects of cytokines
    // Petty and Yang, Tumor-associated macrophages: implications in cancer immunotherapy, 2017
    //  - cytokines suppress T cell function

    if(state != 6){return;}

    // posInfluence is M1 + Th
    // negInfluence is M2 + Treg
    double posInfluence = 1 - (1 - influences[1])*(1 - influences[4]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[5]);

    double scale = posInfluence - negInfluence;
    // killProb = baseKillProb*pow(infScale, scale);
}