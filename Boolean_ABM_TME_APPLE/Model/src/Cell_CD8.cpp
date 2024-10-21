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
    baseKillProb = cellParams[7][2];
    infScale = cellParams[8][2];
    influenceRadius = cellParams[9][2];
    migrationBias = cellParams[10][2];
    divProb_base = cellParams[11][2];
    pTypeStateTransition = cellParams[12][2]; 
    migrationBias_inTumor = cellParams[13][2]; 
    migrationSpeed_inTumor = cellParams[14][2]; 

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

void Cell::cd8_setKillProb(size_t step_count){
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

    //index from the boolean model output to find phenotype and modify 
    size_t step_alive = step_count -  init_time;

    bool flag = false; 
    if((step_alive*pTypeStateTransition-1) < t_cell_phenotype_Trajectory.size()){
        std::string phenotype = t_cell_phenotype_Trajectory[step_alive*pTypeStateTransition - 1];
        char phenotype_char = phenotype[0]; 
        switch(phenotype_char){
        case 'N':  
            killProb = baseKillProb;
            break; 
        case 'M': 
            killProb = baseKillProb * 2;
            flag = true; 
            break; 
        default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate
            killProb = baseKillProb / 10;
        }
    }
    else{ //we assume the t cell stays at the end of its trajectory until it dies, this simulates the behavior of terminally differentiated cells
        char phenotype_char; 
        if (t_cell_phenotype_Trajectory.empty() || (t_cell_phenotype_Trajectory.size() == 0)){
            std::cerr << "WARNING directInteractionProperties: t_cell_phenotype_Trajectory is empty!" << std::endl;
            //handle any bada alloc error by assuming exhausted state...will debug this, very rare and not fatal 
            phenotype_char = 'E'; 
        }
        else{
            phenotype_char = t_cell_phenotype_Trajectory.back()[0];         
        }
        // can set to E
        switch(phenotype_char){
        case 'N': 
            killProb = baseKillProb;
            break; 
        case 'M': 
            killProb = baseKillProb * 2;
            break; 
        default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate
            killProb = baseKillProb / 10;
        }
    }
    killProb = killProb*pow(infScale, scale); //realize the impact of pos and negative influence in coordination with T cell boolean network 
}