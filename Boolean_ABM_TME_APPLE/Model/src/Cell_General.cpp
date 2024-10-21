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
Cell::Cell(std::array<double, 2> loc, int idx, std::vector<std::vector<double>> &cellParams, int cellType, std::vector<std::string> tCellPhenotypeTrajectory, size_t init_tstamp): mt((std::random_device())()) {
    /*
     * initialize all parameters to 0
     * set parameters based on cellType
     */

    type = cellType;

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
    pdl1Shift = 0;
    influenceRadius = 0;
    pdl1 = 0;
    pdl1WhenExpressed = 0;
    for(int i=0; i<influences.size(); ++i){
        influences[i] = 0;
    }
    for(int i=0; i<chemotaxVals.size(); ++i){
        chemotaxVals[i] = 0;
    }
    kTr = 0;
    kM1 = 0;
    kM2 = 0;
    plasticity = 0;
    killProb = 0;
    baseKillProb = 0;
    infScale = 0;

    state = 0;

    // for influence distance, assume a soft-cutoff where p(distance) = probTh
    probTh = 0.01;

    if(cellType == 0){
        initializeCancerCell(cellParams);
    } else if(cellType == 1){
        initializeMacrophageCell(cellParams);
    } else if(cellType == 2){
        initializeCD4Cell(cellParams);
    } else if(cellType == 3){
        initializeCD8Cell(cellParams, tCellPhenotypeTrajectory, init_tstamp);
    } else{
        throw std::runtime_error("Cell::Cell -> unavailable cell type");
    }
}

// FORCE FUNCTIONS
// from Osborne 2017
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

// CELL BEHAVIOR FUNCTIONS
std::array<double, 3> Cell::proliferate(double dt) {
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

void Cell::age(double dt, size_t step_count) {
    /*
     * cells die based on a probability equal to 1/lifespan
     */

    
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    //CD8 cells only 
    if(type == 3){
        size_t step_alive = step_count -  init_time; 
                    
        if(step_alive < 0){
            std::cout << "ERROR in OPENMP REGION: NEGATIVE LIFESPAN" << std::endl; 
        }

        //either it is in one of 7 end states or still progressing through state of gene evolution 
        if((step_alive*pTypeStateTransition-1) < t_cell_phenotype_Trajectory.size()){

            std::string phenotype = t_cell_phenotype_Trajectory[step_alive*pTypeStateTransition - 1];
            char phenotype_char = phenotype[0]; 
            
            switch(phenotype_char){
            case 'N': 
                if(dis(mt) < deathProb){
                    state = -1;
                }
                break; 
            case 'M': 
                if(dis(mt) < (deathProb/2)){
                    state = -1;
                }
                break; 
            default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate     
                if(dis(mt) < deathProb){
                    state = -1;
                }
                
        }
        }
        else {

            char phenotype_char; 
            if (t_cell_phenotype_Trajectory.empty() || (t_cell_phenotype_Trajectory.size() == 0)){
                std::cerr << "WARNING age: t_cell_phenotype_Trajectory is empty!" << std::endl;
                // suppress bad access with Exhausted phenotype
                phenotype_char = 'E'; 
            }
            else{
                phenotype_char = t_cell_phenotype_Trajectory.back()[0];         
            }
            //char phenotype_char = 'E'; //assume all cells out of their timecourse are exhausted 
            switch(phenotype_char){
            case 'N': 
                if(dis(mt) < deathProb){
                    state = -1;
                }
                break; 
            case 'M': 
                if(dis(mt) < deathProb/2){
                    state = -1;
                }
                break; 
            default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate
                if(dis(mt) < deathProb){
                    state = -1;
                }
            }
        }
    }
    if(dis(mt) < deathProb){
        state = -1;
    }
}

void Cell::migrate(double dt, std::array<double,2> tumorCenter, double tumorRadius) {
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

    if(calcDistance(tumorCenter) <= tumorRadius){
        //in tumor
        for(int i=0; i<2; ++i){
            dx_movement[i] = migrationBias_inTumor*dx_direction[i] + (1- migrationBias_inTumor)*dx_random[i];
        }
        dx_movement = unitVector(dx_movement);
            for(int i=0; i<x.size(); ++i){
                x[i] += dt*migrationSpeed_inTumor*dx_movement[i];
                if(std::isnan(x[i])){
                    throw std::runtime_error("migration NaN");
            }
        }
    }
    else{
        //out of tumor
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
    // dx_movement = unitVector(dx_movement);

    // for(int i=0; i<x.size(); ++i){
    //     x[i] += dt*migrationSpeed*dx_movement[i];
    //     if(std::isnan(x[i])){
    //         throw std::runtime_error("migration NaN");
    //     }
    // }
}

void Cell::prolifState() {
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

void Cell::indirectInteractions(double tstep, size_t step_count) {
    /*
     * after determining total influences on the cell, run the indirect interaction functions
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
        cancer_gainPDL1(tstep);
        return;
    } else if (state == 4){
        // CD4 helper
        return;
    } else if (state == 5){
        // CD4 regulatory
        return;
    } else if (state == 6){
        // CD8 active
        cd8_setKillProb(step_count);
        return;
    } else if (state == 7){
        // CD8 suppressed
        return;
    }
}

void Cell::directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep) {
    /*
     * when the cell touches another cell, run direct interactions
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
        if(interactingState == 6){
            // interactionProperties = {radius, killProb}
            cancer_dieFromCD8(interactingX, interactionProperties[0], interactionProperties[1], tstep);
        }
        return;
    } else if (state == 4){
        // CD4 helper
        return;
    } else if (state == 5){
        // CD4 regulatory
        return;
    } else if (state == 6){
        // CD8 active
        if(interactingState == 2 || interactingState == 3 || interactingState == 5){
            // interactionProperties = {radius, pdl1}
            cd8_pdl1Inhibition(interactingX, interactionProperties[0], interactionProperties[1], tstep);
        }
        return;
    } else if (state == 7){
        // CD8 suppressed
        return;
    }
}

std::vector<double> Cell::directInteractionProperties(int interactingState, size_t step_count) {
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
    } else if (state == 3){
        // cancer
        if(interactingState == 6){
            return {radius, pdl1};
        }
        return {};
    } else if (state == 4){
        // CD4 helper
        return {};
    } else if (state == 5){
        // CD4 regulatory
        if(interactingState == 6){
            return {radius, pdl1};
        }
        return {};
    } else if (state == 6){
        // CD8 active
        if(interactingState == 3){
            return {radius, killProb};
        }
        
    } else if (state == 7){
        // CD8 suppressed
        return {};
    }

    return {};
}

// DIFFERENTIATION
void Cell::differentiate(double dt) {
    /*
     * runs either macrophage differentiation or cd4 differentiation
     */
    if(type == 1){
        macrophage_differentiation(dt);
    } else if(type == 2){
        cd4_differentiation(dt);
    }
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

void Cell::addChemotaxis(std::array<double, 2> otherX, double otherInfluence, int otherType) {
    /*
     * assume immune cells chemotax only towards cancer cells
     * cancer cells do not move
     *
     * determines influence at 8 points around the cell, so that it will chemotax towards the highest
     */
    if(type == 0 || otherType != 0){return;}

    std::array<double, 2> chemotax_x = {0.0,0.0};
    std::array<double, 2> dx = {0.0,0.0};

    int z = 0;
    for(int i=-1; i<2; ++i){
        for(int j=-1; j<2; ++j){
            if(i == 0 && j == 0){continue;}

            dx = {static_cast<double>(i), static_cast<double>(j)};
            dx = unitVector(dx);
            chemotax_x = {x[0] + radius*dx[0], x[1] + radius*dx[1]};
            double ct = calcInfDistance(calcNorm({otherX[0] - chemotax_x[0], otherX[1] - chemotax_x[1]}),
                                        otherInfluence);
            chemotaxVals[z] = 1 - (1 - chemotaxVals[z])*(1 - ct);
            z += 1;
        }
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
    double alpha = -log2(probTh);
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