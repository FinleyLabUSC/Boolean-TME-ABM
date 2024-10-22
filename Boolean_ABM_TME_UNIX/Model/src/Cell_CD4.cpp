#include "Cell.h"

void Cell::initializeCD4Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp) {
    state = 4;

    mu = cellParams[0][1];
    kc = cellParams[1][1];
    damping = cellParams[2][1];
    maxOverlap = cellParams[3][1]*cellParams[4][1];
    radius = cellParams[4][1]/2.0;
    deathProb = cellParams[5][1];
    migrationSpeed = cellParams[6][1];
    kTr = cellParams[7][1];
    influenceRadius = cellParams[8][1];
    migrationBias = cellParams[9][1];
    pdl1WhenExpressed = cellParams[10][1];
    migrationBias_inTumor = cellParams[11][1]; 
    migrationSpeed_inTumor = cellParams[12][1]; 

    rmax = 1.5*radius*2;
}

void Cell::cd4_differentiation(double dt) {
    if(state != 4){return;}

    std::uniform_real_distribution<double> dis(0.0,1.0);
    // negInfuence is M2 + alive cancer
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[3]);
    if(dis(mt) < kTr*negInfluence){
        state = 5;
        pdl1 = pdl1WhenExpressed;
    }
}
