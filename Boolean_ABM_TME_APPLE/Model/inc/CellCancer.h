// CELL BEHAVIOR FUNCTIONS
#include "CellCancer.h"



// CellCancer::CellCancer(std::vector<std::vector<double>> &cellParams, size_t init_tstamp) : Cell::Cell(args){ 

CellCancer::CellCancer(std::vector<std::vector<double>> &cellParams, size_t init_tstamp): Cell::Cell(std::array<double, 2> loc, int idx, std::vector<std::vector<double>> &cellParams, int cellType, std::vector<std::string> tCellPhenotypeTrajectory, size_t init_tstamp): mt((std::random_device())()){
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
    
}

void CellCancer::cancer_dieFromCD8(std::array<double, 2> otherX, double otherRadius, double kp, double dt) {
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

void CellCancer::cancer_gainPDL1(double dt) {
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

std::array<double, 3> CellCancer::cancer_proliferate(double dt) {
    // positions 0 and 1 are cell location
    // position 2 is boolean didProliferate?
    if(!canProlif){return {0,0,0};}

    std::uniform_real_distribution<double> dis(0.0, 1.0);
    if(dis(mt) < divProb){
        // place daughter cell a random angle away from the mother cell at a distance of radius
        std::normal_distribution<double> rd(0.0, 1.0);
        std::array<double, 2> dx = {rd(mt),
                                      rd(mt)};
        double norm = calcNorm(dx);
        return{radius*(dx[0]/norm)+x[0],
               radius*(dx[1]/norm)+x[1],
               1};
    } else{
        return {0,0,0};
    }
}

void CellCancer::cancer_prolifState() {
    /*
     * cancer cells and CD8 can proliferate
     * right now, CD8 proliferation prob is set to 0, however leaving it in for future changes
     */
    if(type == 0){
        canProlif = !(state == -1 || compressed);
    } else if(type == 3){
        // CTLs -> presence of Th promotes their proliferation, M2 and Treg decrease it
        // assume CTLs need IL-2 from Th to proliferate
        canProlif = !(state == 7 || compressed);
        double posInfluence = influences[4];
        double negInfluence = 1 - (1 - influences[2])*(1 - influences[5]);

        double scale = posInfluence - negInfluence;
        //divProb = divProb_base*pow(infScale, scale);
        divProb = scale*divProb_base;
    } else {
        canProlif = false;
    }
}

void CellCancer::cancer_age(double dt, size_t step_count) {
    /*
     * cells die based on a probability equal to 1/lifespan
     */

    std::uniform_real_distribution<double> dis(0.0, 1.0);

    if(dis(mt) < deathProb){
        state = -1;
    }
}

void CellCancer::cancer_indirectInteractions(double tstep) {
    cancer_gainPDL1(tstep);
    return;
}

void CellCancer::cancer_directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep) {
    if(interactingState == 6){
        // interactionProperties = {radius, killProb}
        cancer_dieFromCD8(interactingX, interactionProperties[0], interactionProperties[1], tstep);
        }
        return;
}

std::vector<double> CellCancer::cancer_directInteractionProperties(int interactingState, size_t step_count) {

    /*
     * returns the properties that go into Cell::directInteractions
     */
    if (state == 3){
        // cancer
        if(interactingState == 6){
            return {radius, pdl1};
        }
        return {};
    } 
    return {};
}


