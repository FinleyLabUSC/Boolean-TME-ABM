#include "CD4Cell.h"
#include "Cell.h"

CD4Cell::CD4Cell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp, std::array<double, 2> loc, int idx, 
            int cellType)
            :Cell(loc, idx, cellParams, init_tstamp)
            {
    
    type = cellType; 
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

    rmax = 1.5*radius*2;
}

/*
M2 cells promote differentiation of cd4 into regulatory state (and also promotes M0 into M2)
CD4 in this model comes in as helper state and can be converted to regulatory state based on enviromental factors ^^

*/
void CD4Cell::cd4_differentiation(double dt) {
    if(state != 4){return;}

    std::uniform_real_distribution<double> dis(0.0,1.0);
    // negInfuence is M2 + alive cancer
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[3]);
    if(dis(mt) < kTr*negInfluence){
        state = 5;
        pdl1 = pdl1WhenExpressed;
    }
}

void CD4Cell::cd4_age(double dt, size_t step_count) {
    /*
     * cells die based on a probability equal to 1/lifespan
     */

    std::uniform_real_distribution<double> dis(0.0, 1.0);

    if(dis(mt) < deathProb){
        state = -1;
    }
}


/*
In the helper state, CD4+ cells promote M0 differentiation into the M1 state 
As regulatory cells, they express CTLA-4 (has same function as PD-L1) and promote M0 differentiation into the M2 state 
*/
std::vector<double> CD4Cell::cd4_directInteractionProperties(int interactingState, size_t step_count) {

    /*
     * returns the properties that go into Cell::directInteractions
     */
    if (state == 4){
        // CD4 helper
        return {};
    }
    else if (state == 5){
        // CD4 regulatory
        if(interactingState == 6){
            return {radius, pdl1};
        }
        return {};
    }
    return {};
}
void CD4Cell::set_pdl1WhenExpressed(double value) {pdl1WhenExpressed = value;}
double CD4Cell::get_pdl1WhenExpressed() {return pdl1WhenExpressed;}
void CD4Cell::set_probTh(double value) {probTh = value;}
double CD4Cell::get_probTh() {return probTh;}
void CD4Cell::set_kTr(double value) {kTr = value;}
double CD4Cell::get_kTr() {return kTr;}