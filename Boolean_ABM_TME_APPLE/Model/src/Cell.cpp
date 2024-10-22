#include "Cell.h"

/*
 * CELL TYPES
 * 0 - cancer
 * 1 - macrophage
 * 2 - CD4
 * 3 - CD8
 *
 * CELL STATES
 * -1 - dead
 * 0 - M0
 * 1 - M1
 * 2 - M2
 * 3 - alive (cancer)
 * 4 - Th
 * 5 - Treg
 * 6 - active (CD8)
 * 7 - suppressed (CD8)
 *
 */

// INITIALIZE CELL TYPE
Cell::Cell(std::array<double, 2> loc, int idx, std::vector<std::vector<double>> &cellParams, size_t init_tstamp): mt((std::random_device())()) {
    /*
     * initialize all parameters to 0
     * set parameters based on cellType
     */

    x = loc;
    id = idx;

    radius = 0;
    compressed = false;
    currentOverlap = 0;
    divProb = 0;
    divProb_base = 0;
    deathProb = 0;
    canProlif = false;
    mu = 0;
    kc = 0;
    damping = 0;
    maxOverlap = 0;
    rmax = 0;
    currentForces = {0,0};
    migrationSpeed = 0;
    migrationBias = 0;
    influenceRadius = 0;
    pdl1 = 0;
    // pdl1WhenExpressed = 0;
    for(int i=0; i<influences.size(); ++i){
        influences[i] = 0;
    }
    for(int i=0; i<chemotaxVals.size(); ++i){
        chemotaxVals[i] = 0;
    }
    // kTr = 0;
    // kM1 = 0;
    // kM2 = 0;
    // plasticity = 0;
    // killProb = 0;
    // baseKillProb = 0;
    // infScale = 0;

    state = 0;

    // for influence distance, assume a soft-cutoff where p(distance) = probTh
    // probTh = 0.01;
}
/*
FORCE FUNCTIONS
from Osborne 2017

cells are modeled using a center-based approach which considers each cell as a point and radius 
(used to calculate physical force between cells)
*/
std::array<double, 2> Cell::attractiveForce(std::array<double, 2> dx, double otherRadius) {
    double dxNorm = calcNorm(dx);
    std::array<double, 2> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm};
    double sij = radius + otherRadius;

    double scaleFactor = mu*(dxNorm - sij)*exp(-kc*(dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;

    return {F0, F1};
}

std::array<double, 2> Cell::repulsiveForce(std::array<double, 2> dx, double otherRadius) {
    double dxNorm = calcNorm(dx);
    std::array<double, 2> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm};
    double sij = radius + otherRadius;

    double scaleFactor = mu*sij*log10(1 + (dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;

    return {F0, F1};
}

void Cell::calculateForces(std::array<double, 2> otherX, double otherRadius, int &otherType) {
    /*
     * assume attractive force only between cancer cells
     */
    double distance = calcDistance(otherX);
    if(distance < rmax){
        std::array<double, 2> dx = {(otherX[0]-x[0]),
                                    (otherX[1]-x[1])};
        if(distance < (radius + otherRadius)){
            std::array<double, 2> force = repulsiveForce(dx, otherRadius);
            currentForces[0] += force[0];
            currentForces[1] += force[1];
        } else if(type == 0 && otherType == 0){ // attraction if both are cancer cells
            std::array<double, 2> force = attractiveForce(dx, otherRadius);
            currentForces[0] += force[0];
            currentForces[1] += force[1];
        }
    }
}

void Cell::resolveForces(double dt, std::array<double, 2> &tumorCenter, double &necroticRadius, double &necroticForce) {
    /*
     * if the cell is touching the necrotic core, push it outward
     * numerically solve forces
     */
    if(calcDistance(tumorCenter) < necroticRadius+radius){
        std::array<double, 2> dx = {x[0] - tumorCenter[0],
                                    x[1] - tumorCenter[1]};
        dx = unitVector(dx);
        if(!std::isnan(dx[0])) {
            currentForces[0] += necroticForce * dx[0];
            currentForces[1] += necroticForce * dx[1];
        }
    }

    x[0] += (dt/damping)*currentForces[0];
    x[1] += (dt/damping)*currentForces[1];

    resetForces();
}

void Cell::resetForces() {
    /*
     * resets forces with a slight randomizing factor
     */
    double D = 1;
    std::uniform_real_distribution<double> dis(-D, D);
    currentForces = {dis(mt),dis(mt)};
}

void Cell::neighboringCells(std::array<double, 2> otherX, int otherID){
    /*
     * determine which cells are within 2*maximum interaction distance
     * stores the index in cell_list (in Environment) of the neighboring cells
     */
    double dis = calcDistance(otherX);
    if(dis <= 10*rmax){
        neighbors.push_back(otherID);
    }
}

// OVERLAP FUNCTIONS
void Cell::calculateOverlap(std::array<double, 2> otherX, double otherRadius) {
    /*
     * sum up the total overlap distance between the cell and all other cells
     */
    double distance = calcDistance(otherX);
    if(distance < radius + otherRadius){
        currentOverlap += radius + otherRadius - distance;
    }
}

void Cell::resetOverlap() {
    currentOverlap = 0;
}

void Cell::isCompressed() {
    compressed = currentOverlap > maxOverlap;
    resetOverlap();
}




void Cell::migrate(double dt, std::array<double,2> tumorCenter) {
    /*
     * biased random-walk towards tumor center
     *
     * commented out code for migrating up a pseudo-chemotaxix gradient
     * it produces weird spatial behaviors and makes it difficult to recruit immune cells around the tumor
     */
    if(type == 0 || state == -1 || state == 7){return;} // cancer cells, dead cells, suppressed CD8

    /*std::uniform_int_distribution<int> choose_cv(0, chemotaxVals.size()-1);
    int cv = choose_cv(mt);
    for(int i=0; i<chemotaxVals.size(); ++i){
        if(chemotaxVals[i] > chemotaxVals[cv]){cv = i;}
    }
    int z = 0;
    std::array<double, 2> dx_direction = {0,0};
    for(int i=-1; i<2; ++i){
        for(int j=-1; j<2; ++j){
            if(i==0 && j==0){continue;}
            if(z == cv){
                dx_direction = {static_cast<double>(i),
                                static_cast<double>(j)};
            }
            z++;
        }
    }*/
    std::array<double, 2> dx_direction = {tumorCenter[0] - x[0],
                                          tumorCenter[1] - x[1]};
    std::normal_distribution<double> vect(0.0, 1.0);
    std::array<double, 2> dx_random = {vect(mt), vect(mt)};

    dx_direction = unitVector(dx_direction);
    dx_random = unitVector(dx_random);
    std::array<double, 2> dx_movement = {0,0};
    for(int i=0; i<2; ++i){
        dx_movement[i] = migrationBias*dx_direction[i] + (1- migrationBias)*dx_random[i];
    }
    dx_movement = unitVector(dx_movement);

    for(int i=0; i<x.size(); ++i){
        x[i] += dt*migrationSpeed*dx_movement[i];
        if(std::isnan(x[i])){
            throw std::runtime_error("migration NaN");
        }
    }
}



void Cell::inherit(std::vector<double> properties) {
    /*
     * daughter cells inherit properties from mother
     * even though this is applicable only for cancer cells, I left it this way for ease of running cell functions
     */
    if(state == 0){
        // M0 macrophages
        return;
    } else if (state == 1){
        // M1 macrophages
        return;
    } else if (state == 2){
        // M2 macrophages
        return;
    } else if (state == 3){
        // cancer
        pdl1 = properties[0];
        return;
    } else if (state == 4){
        // CD4 helper
        return;
    } else if (state == 5){
        // CD4 regulatory
        return;
    } else if (state == 6){
        // CD8 active
        return;
    } else if (state == 7){
        // CD8 suppressed
        return;
    }
}

std::vector<double> Cell::inheritanceProperties() {
    /*
     * returns the properties that go into Cell::inherit
     */
    if(state == 0){
        // M0 macrophages
        return {};
    } else if (state == 1){
        // M1 macrophages
        return {};
    } else if (state == 2){
        // M2 macrophages
        return {};
    } else if (state == 3){
        // cancer
        return {pdl1};
    } else if (state == 4){
        // CD4 helper
        return {};
    } else if (state == 5){
        // CD4 regulatory
        return {};
    } else if (state == 6){
        // CD8 active
        return {};
    } else if (state == 7){
        // CD8 suppressed
        return {};
    }

    return {};
}

// CELL INFLUENCE
void Cell::addInfluence(std::array<double, 2> otherX, double otherInfluence, int otherState) {
    /*
     * determine influence based on distance for each cell state
     * total influence is 1
     */
    if(otherState == -1){return;}

    influences[otherState] = 1 - (1 - influences[otherState])*(1 - calcInfDistance(calcDistance(otherX), otherInfluence));
}

void Cell::clearInfluence() {
    for(int i=0; i<influences.size(); ++i){
        influences[i] = 0;
    }
    for(int i=0; i<chemotaxVals.size(); ++i){
        chemotaxVals[i] = 0;
    }
}



// OTHER FUNCTIONS
double Cell::calcDistance(std::array<double, 2> otherX) {
    double d0 = (otherX[0] - x[0]);
    double d1 = (otherX[1] - x[1]);

    return sqrt(d0*d0 + d1*d1);
}

double Cell::calcInfDistance(double dist, double xth) {
    /*
     * calculate influence using an exponential decay based on distance from cell center
     */
    double alpha = -log2(1);//probTh
    double lambda = alpha*0.693/xth;

    return exp(-lambda*dist);
}

double Cell::calcNorm(std::array<double, 2> dx){
    return sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
}

std::array<double, 2> Cell::unitVector(std::array<double, 2> v) {
    double norm = calcNorm(v);
    return {v[0]/norm, v[1]/norm};
}

void Cell::updateID(int idx) {
    id = idx;
}