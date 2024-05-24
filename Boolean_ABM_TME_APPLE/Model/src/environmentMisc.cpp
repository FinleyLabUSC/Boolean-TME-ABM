#include "Environment.h"
#include "ModelUtil.h"

void Environment::initializeCells() {
    /*
     * places the initial tumor
     */

    double radiiCells = envParams[0];
    cell_list.push_back(Cell({0.0,0.0}, 0, cellParams, 0, tCellPhenotypeTrajectory_1));
    int q = 1;
    for(int i=1; i<radiiCells; ++i){
        double circumfrence = 2*i*cellParams[4][0]*3.1415;
        double nCells = circumfrence/cellParams[4][0];
        for(int j=0; j<nCells; ++j){
            double x = i * cellParams[4][0] * cos(2 * 3.1415 * j / nCells);
            double y = i * cellParams[4][0] * sin(2 * 3.1415 * j / nCells);
            cell_list.push_back(Cell({x, y}, q, cellParams, 0, tCellPhenotypeTrajectory_1));
            q++;
        }
    }
}

void Environment::initializeVessels(){

    // we take the initial tumor radius, pad it with the vesselRadius 
    //then we place the vessels somewhere outside this initial radius at
    //a density dictated by vesDens, we sample from a uniform distribution 
    //to uniformly place the number of vesels.
    double some_radius = envParams[0]+vesRad; //this needs to be confirmed
    
    int Nves = M_PI*pow(some_radius, 2)/(vesDens); //type coerce to int

    std::uniform_real_distribution<double> spacing(vesDens-100,vesDens+100);
    std::uniform_real_distribution<double> vesDis(100,some_radius-100);
    
    

    for(int q=0; q<Nves; ++q){
        double i = vesDis(mt);
        double j = vesDis(mt);
        double space = spacing(mt);
        bool badSpace = true;
        int tries = 0;
        while(badSpace){
            tries++;
            if(tries > 10000){
                // if can't find a suitable location, just choose a random one
                badSpace = false;
            }
            double minDis = 250;
            for(auto & vessel : vessel_list){
                double di = i - vessel.x[0];
                double dj = j - vessel.x[1];
                double distance = sqrt(di*di + dj*dj);
                minDis = std::min(distance, minDis);
            }
            if(minDis < space){
                i = vesDis(mt);
                j = vesDis(mt);
            } else{
                badSpace = false;
            }
        }
        
        vessel_list.push_back(Vessel({i, j}, vesRad, envParams[2], vesMu, vesDec, 0, vesInfluenceDistance));
        vessel_list[vessel_list.size()-1].generalLocation(0);
    }
    
}

void Environment::recruitImmuneCells(double tstep,  size_t step_count) {
    // recruitment is scaled by number of cancer cells
    auto day = static_cast<double>(tstep*steps/24.0);
    if(day < recruitmentDelay){return;}

    int numC = cancerTS[steps - recruitmentDelay*24/tstep];

    for(int i=0; i<immuneCellRecRates.size(); ++i){
        double recRate = immuneCellRecRates[i]*static_cast<double>(numC);
        immuneCells2rec[i] += tstep*recRate;
        while(immuneCells2rec[i] >= 1){
            std::array<double, 2> recLoc = recruitmentLocation();
            
            //i==0 represents the idx associated with initializing a cd8 t cell
            if(i == 0){
                
                size_t phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 

                std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                    std::cerr << "WARNING recruitImmuneCells: t_cell_phenotype_Trajectory is empty!" << std::endl;
                    std::cerr << phenotypeIdx << std::endl; 
                    std::cerr << "ATTEMPTING FORCE RESOLVE... " << std::endl; 
                    phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 
                    trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                    if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                        std::cerr << "FORCE RESOLVE FAILED... " << std::endl; 
                    }
                    else{
                        std::cerr << "FORCE RESOLVE SUCCESSFUL... " << std::endl; 
                    }
                }
                

                cell_list.push_back(Cell(recLoc, 
                                    static_cast<int>(cell_list.size()), 
                                    cellParams, 
                                    immuneCellRecTypes[i],
                                    trajec_phenotype, 
                                    step_count)); 
            }
            else{
                cell_list.push_back(Cell(recLoc,
                                     static_cast<int>(cell_list.size()),
                                     cellParams,
                                     immuneCellRecTypes[i], tCellPhenotypeTrajectory_1));

            }
            
            immuneCells2rec[i] -= 1;
        }
    }
}

void Environment::recruitVessels(double tstep, size_t step_count){

    auto day = static_cast<double>(tstep*steps/24.0);
    if(day < recruitmentDelay){return;}

    float vesselRecRate = 0.01; 
    float pThreshold = 0.05; 
    int numRecruit = static_cast<int>(vesselRecRate * cell_list.size()); 

    std::uniform_real_distribution<float> p_recruit(0.0, 1.0); 
    float a = p_recruit(mt); 
    if(a <= pThreshold){
        for(int i =0; i < 1; ++i){

        std::normal_distribution<double> angle(0.0,1.0);
        std::array<double, 2> dx = {angle(mt),
                                    angle(mt)};
        double norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
    
        std::uniform_real_distribution<double> loc(100, recDist+100); //test make sure vessel
        double distance = loc(mt)+ tumorRadius; 
    
        float u_i = dx[0]/norm; 
        float u_j = dx[1]/norm; 
        float loc_x = float(distance)*u_i;
        float loc_y = float(distance)*u_j; 

        vessel_list.push_back(Vessel({loc_x, loc_y}, vesRad, envParams[2], vesMu, vesDec, 0, vesInfluenceDistance));
        vessel_list[vessel_list.size()-1].generalLocation(0);
        
    }
    }
    return; 
    

}
std::array<double, 2> Environment::recruitmentLocation() {
    /*
     * cells enter a random distance away from the tumor radius
     * cells enter at a random angle from the tumor center
     */
    //std::uniform_real_distribution<double> angle(-1.0, 1.0);
    std::normal_distribution<double> angle(0.0,1.0);
    std::array<double, 2> dx = {angle(mt),
                                angle(mt)};
    double norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

    std::uniform_real_distribution<double> loc(0.0, recDist);
    double distance = loc(mt) + tumorRadius;

    return {distance*(dx[0]/norm), distance*(dx[1]/norm)};

    /*
     * recruit based on cytokine "concentrations". Recruit if a location is in a certain range
     */
    /*std::array<double, 2> x = {0.0, 0.0};
    double d = 0.0;
    std::normal_distribution<double> angle(0.0,1.0);
    std::uniform_real_distribution<double> distance_from_center(0.5*tumorRadius, 1.5*tumorRadius);
    std::uniform_real_distribution<double> recruitmentProbability(0.0, 1.0);
    std::array<double, 2> dx = {0.0, 0.0};
    double norm = 0;
    double distance = 0.0;
    bool recruited = false;
    while(!recruited){
        dx = {angle(mt),
              angle(mt)};
        norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
        distance = distance_from_center(mt);
        x = {distance*dx[0]/norm, distance*dx[1]/norm};
        d = calculateDiffusibles(x);
        if(d <= maxRecCytoConc && recruitmentProbability(mt) < d){
            recruited = true;
        }
    }*/

    /*std::array<double, 2> x = {0.0, 0.0};
    double minD = 1e6;
    double maxD = 200.0;
    std::normal_distribution<double> angle(0.0,1.0);
    std::uniform_real_distribution<double> distance_from_center(0.5*tumorRadius, tumorRadius + 200);
    std::array<double, 2> dx = {0.0, 0.0};
    double norm = 0;
    double distance = 0.0;
    bool recruited = false;
    while(!recruited){
        dx = {angle(mt),
              angle(mt)};
        norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
        distance = distance_from_center(mt);
        x = {distance*dx[0]/norm, distance*dx[1]/norm};
        for(auto & c : cell_list){
            if(c.state == 3){
                double dist = c.calcDistance(x);
                minD = std::min(dist, minD);
                if(minD < 50){
                    break;
                }

            }
        }
        if(minD >= 50.0 && minD < 200.0){
            recruited = true;
        } else{
            minD = 1e6;
        }
    }

    return x;*/
}

void Environment::tumorSize(){
    tumorCenter = {0,0};
    double avgX = 0;
    double avgY = 0;
    double numC = 0;
    for(auto &c : cell_list){
        if(c.type == 0) {
            avgX += c.x[0];
            avgY += c.x[1];
            numC += 1;
        }
    }
    avgX /= numC;
    avgY /= numC;

    tumorCenter = {avgX, avgY};

    tumorRadius = 0;
    for(auto& c : cell_list){
        if(c.type == 0){
            tumorRadius = std::max(tumorRadius, c.calcDistance(tumorCenter));
        }
    }

    //originally commented out

    // tumorCenter = {0.0, 0.0};
    // double avgX = 0;
    // double avgY = 0;
    // double numC = 0;
    // std::vector<int> mainTumorMass;
    // for(auto & c : cell_list){
    //     if(c.type == 0){
    //         int numNeighbors = 0;
    //         for(auto &c2 : cell_list){
    //             if(c2.type == 0 && c.calcDistance(c2.x) < 50.0){
    //                 ++numNeighbors;
    //             }
    //         }
    //         if(numNeighbors > 10){
    //             numC++;
    //             avgX += c.x[0];
    //             avgY += c.x[1];
    //             mainTumorMass.push_back(1);
    //         } else{
    //             mainTumorMass.push_back(0);
    //         }
    //     } else{
    //         mainTumorMass.push_back(0);
    //     }
    // }

    // avgX /= numC;
    // avgY /= numC;

    // tumorCenter = {avgX, avgY};

    //----END originally commented out

    // tumorRadius = 0.0;
    // for(int i=0; i<cell_list.size(); ++i){
    //     if(mainTumorMass[i] == 1){
    //         if(cell_list[i].type == 0){
    //             tumorRadius = std::max(tumorRadius, cell_list[i].calcDistance(tumorCenter));
    //         }
    //     }
    // }

}

void Environment::necrosis(double tstep) {
    int nCancer = 0;
    for(auto &c : cell_list){
        if(c.type == 0){
            ++nCancer;
        }
    }

    necroticRadius += tstep*nCancer*necroticGrowth;
    necroticRadius = std::max(std::min(necroticRadius, tumorRadius-necroticLimit), 0.0);
}

double Environment::calculateDiffusibles(std::array<double, 2> x) {
    double d = 0;
    for(auto & c: cell_list){
        if(c.state == 3){
            // is a cancer cell
            double alpha = -log2(c.probTh);
            double lambda = alpha*0.693/c.influenceRadius;
            double dist = c.calcDistance(x);
            d = 1 - (1 - d)*(1 - exp(-lambda*dist));
        }
    }

    return d;
}