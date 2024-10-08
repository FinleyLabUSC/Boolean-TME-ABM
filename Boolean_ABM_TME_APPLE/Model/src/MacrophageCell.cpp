#include "MacrophageCell.h"

MacrophageCell::MacrophageCell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp, std::array<double, 2> loc, int idx, 
          int cellType)
    : Cell(loc, idx, cellParams, init_tstamp){
    state = 0;

    mu = cellParams[0][3];
    kc = cellParams[1][3];
    damping = cellParams[2][3];
    maxOverlap = cellParams[3][3]*cellParams[4][3];
    radius = cellParams[4][3]/2.0;
    deathProb = cellParams[5][3];
    migrationSpeed = cellParams[6][3];
    kM1 = cellParams[7][3];
    kM2 = cellParams[8][3];
    influenceRadius = cellParams[9][3];
    migrationBias = cellParams[10][3];
    pdl1WhenExpressed = cellParams[11][3];

    rmax = 1.5*radius*2;
}

void MacrophageCell::macrophage_differentiation(double dt) {
    /*
     * differentiate
     * probability of not differentiating is constant (before scaling)
     *  the closer a macrophage is to other cells, the more likely to differentiate
     * inf-g secreting cells (CTL and Th) promote M1
     * M2, cancer, and Treg promote M2
     */
    if(type != 1){return;}

    // posInfluence is active CD8 + Th
    // negInfluence is M2 + alive cancer + Treg
    double posInfluence = 1 - (1 - influences[6])*(1 - influences[4]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[3])*(1 - influences[5]);
    double p1 = kM1*posInfluence;
    double p2 = kM2*negInfluence;
    auto p0 = static_cast<double>(state == 0);
    double sum = p0 + p1 + p2;

    std::array<double, 3> probs = {p0/sum,
                                   (p0 + p1)/sum,
                                   (p0 + p1 + p2)/sum};
    int choice = 0;
    std::uniform_real_distribution<double> dis(0.0,1.0);
    double p = dis(mt);
    for(int i=0; i<3; ++i){
        if(p > probs[i]){choice++;}
    }
    state = choice;
    if(state == 1){
        pdl1 = 0;
    } else if(state == 2){
        pdl1 = pdl1WhenExpressed;
    }
}

void MacrophageCell::macrophage_age(double dt, size_t step_count) {
    /*
     * cells die based on a probability equal to 1/lifespan
     */

    std::uniform_real_distribution<double> dis(0.0, 1.0);

    if(dis(mt) < deathProb){
        state = -1;
    }
}

std::vector<double> MacrophageCell::macrophage_directInteractionProperties(int interactingState, size_t step_count) {
    /*
     * returns the properties that go into Cell::directInteractions
     */
    if(state == 0){
        // M0 macrophages
        return {};
    } else if (state == 1){
        // M1 macrophages
        return {};
    } else if (state == 2){
        // M2 macrophages
        if(interactingState == 6){
            return {radius, pdl1};
        }
        return {};
    }
    return {};
}

void MacrophageCell::set_kTr(double value) {kTr = value;}
double MacrophageCell::get_kTr() {return kTr;}
void MacrophageCell::set_kM1(double value) {kTr = value;}
double MacrophageCell::get_kM1() {return kTr;}
void MacrophageCell::set_kM2(double value) {kTr = value;}
double MacrophageCell::get_kM2() {return kTr;}
void MacrophageCell::set_plasticity(double value) {kTr = value;}
double MacrophageCell::get_plasticity() {return kTr;}
void MacrophageCell::set_pdl1WhenExpressed(double value) {pdl1WhenExpressed = value;}
double MacrophageCell::get_pdl1WhenExpressed() {return pdl1WhenExpressed;}




