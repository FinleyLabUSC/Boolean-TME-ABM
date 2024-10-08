#include "CD8Cell.h"

CD8Cell::CD8Cell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp, std::array<double, 2> loc, int idx, 
          int cellType,std::vector<std::string> phenotypeTrajectory) 
    : Cell(loc, idx, cellParams, init_tstamp)
    {
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

    t_cell_phenotype_Trajectory = phenotypeTrajectory; 
    init_time = init_tstamp;
}

void CD8Cell::cd8_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt) {
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

void CD8Cell::cd8_setKillProb(){
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
    killProb = baseKillProb*pow(infScale, scale);
}

void CD8Cell::cd8_addChemotaxis(std::array<double, 2> otherX, double otherInfluence, int otherType) {
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


void CD8Cell::cd8_age(double dt, size_t step_count) {
    /*
     * cells die based on a probability equal to 1/lifespan
     */

    std::uniform_real_distribution<double> dis(0.0, 1.0);
    
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
            
        case 'M': 
            if(dis(mt) < (deathProb/2)){
                state = -1;
            }
            
        default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate     
            if(dis(mt) < deathProb){
                state = -1;
            }
    }
    }
    else {
        std::string phenotype = t_cell_phenotype_Trajectory[t_cell_phenotype_Trajectory.size() - 1];
        char phenotype_char = phenotype[0]; 
        switch(phenotype_char){
        case 'N': 
            if(dis(mt) < deathProb){
                state = -1;
            }
        case 'M': 
            if(dis(mt) < deathProb){
                state = -1;
            }
            
        default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate
            if(dis(mt) < deathProb){
                state = -1;
            }
        }
    }  
}

void CD8Cell::cd8_indirectInteractions(double tstep) {
    cd8_setKillProb();
    return;
}

void CD8Cell::cd8_directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep) {
    if(interactingState == 2 || interactingState == 3 || interactingState == 5){
        // interactionProperties = {radius, pdl1}
        cd8_pdl1Inhibition(interactingX, interactionProperties[0], interactionProperties[1], tstep);
        }
        return;
}

std::vector<double> CD8Cell::cd8_directInteractionProperties(int interactingState, size_t step_count) {
    /*
     * returns the properties that go into Cell::directInteractions
     */
    
    if (state == 6) {
        // CD8 active
        if(interactingState == 3) {    
            size_t step_alive = step_count -  init_time; 
        
            if((step_alive*pTypeStateTransition-1) < t_cell_phenotype_Trajectory.size()){
                
                std::string phenotype = t_cell_phenotype_Trajectory[step_alive*6 - 1];
                char phenotype_char = phenotype[0]; 
                
                switch(phenotype_char){
                case 'N': 
                    return {radius, killProb};
                case 'M': 
                    return {radius, killProb*2};
                default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate
                    return {radius, killProb/10};
                }
            }
            else { //we assume the t cell stays at the end of its trajectory until it dies
                std::string phenotype = t_cell_phenotype_Trajectory[
                    t_cell_phenotype_Trajectory.size() - 1
                ];
                
                char phenotype_char = phenotype[0]; 
                switch(phenotype_char){
                case 'N': 
                    return {radius, killProb};
                case 'M': 
                    return {radius, killProb*2};
                default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate
                    return {radius, killProb/10};
                }
            }
            return {radius, killProb};
        }
    } 
    else if (state == 7){
        // CD8 suppressed
        return {};
    }
    
    return {};
}
void CD8Cell::set_killProb(double value) {killProb = value;}
double CD8Cell::get_killProb() {return killProb;}
void CD8Cell::set_baseKillProb(double value) {baseKillProb = value;}
double CD8Cell::get_baseKillProb() {return baseKillProb;}
void CD8Cell::set_infScale(double value) {infScale = value;}
double CD8Cell::get_infScale() {return infScale;}
void CD8Cell::set_t_cell_phenotype_Trajectory(std::vector<std::string> value) {t_cell_phenotype_Trajectory = value;}
std::vector<std::string>  CD8Cell::get_t_cell_phenotype_Trajectory() {return t_cell_phenotype_Trajectory;}
