#include "Cell.h"

void Cell::initializeCancerCell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp) {
    state = 3;
    canProlif = true;

    mu = cellParams[0][0];
    kc = cellParams[1][0];
    damping = cellParams[2][0];
    maxOverlap = cellParams[3][0]*cellParams[4][0];
    radius = cellParams[4][0]/2.0;
    divProb = cellParams[5][0];
    deathProb = cellParams[6][0];
    influenceRadius = cellParams[7][0];
    pdl1WhenExpressed = cellParams[8][0];
    pdl1Shift = cellParams[9][0];

    rmax = 1.5*radius*2;

    float p_ferroptosis_sensitive = cellParams[10][0]; //probability a cell is ferroptosis sensitive
    std::uniform_real_distribution<double> dis(0.0, 1.0); //instantiate uniform distribution

    if(dis(mt) <= p_ferroptosis_sensitive){
        ferroptosis_sensitive = 1; 
    }
    else{
        ferroptosis_sensitive = 0; 
    }



    
}

void Cell::cancer_dieFromCD8(std::array<double, 2> otherX, double otherRadius, double kp, double dt) {
    /*
     * die from CTL based on a probability
     * contact required
     */
    if(type != 0){return;}

    if(calcDistance(otherX) <= radius+otherRadius){
        std::uniform_real_distribution<double> dis(0.0,1.0); 
        if(dis(mt) < kp){
            state = -1;
        }
    }
}

void Cell::cancer_gainPDL1(double dt) {
    /*
     * shift pdl1 value based on influence from CTL and Th
     *  ifn-g is shown to increase PD-L1 expression
     */
    if(type != 0 || pdl1 > 0.0){return;}

    // induced by ifn-g secreting cells
    // posInfluence is Th + CD8
    double posInfluence = 1 - (1 - influences[4])*(1 - influences[6]);
    std::uniform_real_distribution<double> dis(0.0,1.0);
    if(dis(mt) < posInfluence*pdl1Shift){
        pdl1 = pdl1WhenExpressed;
    }
}
