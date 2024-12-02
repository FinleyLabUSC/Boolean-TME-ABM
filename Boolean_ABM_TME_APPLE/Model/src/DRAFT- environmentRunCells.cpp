#include "Environment.h"
#include "ModelUtil.h"

void Environment::neighborInfluenceInteractions(double tstep, size_t step_count) {

    /*
     * FIRST LOOP
     * - determine neighbors
     * - determine influences on a cell
     * - perform indirect interactions
     *
     * SECOND LOOP
     * - direct interactions
     *
     * THIRD LOOP
     * - differentiate
     */

#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        // reset neighborhood and influence
        cell_list[i].neighbors.clear();
        cell_list[i].clearInfluence();
        for(auto &c : cell_list){
            // assume that a cell cannot influence itself
            if(cell_list[i].id != c.id){
                cell_list[i].neighboringCells(c.x, c.id);
                cell_list[i].addInfluence(c.x, c.influenceRadius, c.state);
                //cell_list[i].addChemotaxis(c.x, c.influenceRadius, c.type);
            }
        }
        cell_list[i].indirectInteractions(tstep);
    }

#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        for(auto &c : cell_list[i].neighbors){
            cell_list[i].directInteractions(cell_list[c].state,
                                            cell_list[c].x,
                                            cell_list[c].directInteractionProperties(cell_list[i].state, step_count),
                                            tstep);
        }
    }

#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        cell_list[i].differentiate(tstep);
    }
}

void Environment::calculateForces(double tstep) {
    /*
     * 1. Calculate total force vector for each cell
     * 2. Resolve forces on each cell
     * 3. Determine current overlap for each cell
     * 4. Determine if each cell is compressed
     */

    // divide tstep into smaller steps for solving
    // only solve forces between neighboring cells to improve computation time
    int Nsteps = static_cast<int>(tstep/dt);

    // iterate through Nsteps, calculating and resolving forces between neighbors
    // also includes migration
    for(int q=0; q<Nsteps; ++q){
        // migrate first
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            cell_list[i].migrate(dt, tumorCenter);
        }

        // calc forces
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            for(auto &c : cell_list[i].neighbors){
                cell_list[i].calculateForces(cell_list[c].x, cell_list[c].radius, cell_list[c].type);
            }
        }

        // resolve forces
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            cell_list[i].resolveForces(dt, tumorCenter, necroticRadius, necroticForce);
        }
    }

    // calculate overlaps and proliferation states
#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        for(auto &c : cell_list[i].neighbors){
            cell_list[i].calculateOverlap(cell_list[c].x, cell_list[c].radius);
        }
        cell_list[i].isCompressed();
    }
}


void Environment::internalCellFunctions(double tstep, size_t step_count) { // from git hub

    /*
     * cell death via aging
     */

    for(int i=0; i< cancer_list.size(); ++i){
        cancer_list[i].cancer_age(tstep, step_count);
        // if in necrotic core, die
        if(cancer_list[i].calcDistance(tumorCenter) < necroticRadius){
            cancer_list[i].state = -1;
        }
    }
    
    for(int i=0; i< cd4_list.size(); ++i){
        cd4_list[i].cd4_age(tstep, step_count);
        // if in necrotic core, die
        if(cd4_list[i].calcDistance(tumorCenter) < necroticRadius){
            cd4_list[i].state = -1;
        }
    }
    
    for(int i=0; i< cd8_list.size(); ++i){
        cd8_list[i].cd8_age(tstep, step_count);
        // if in necrotic core, die
        if(cd8_list[i].calcDistance(tumorCenter) < necroticRadius){
            cd8_list[i].state = -1;
        }
    }

    for(int i=0; i< macrophage_list.size(); ++i){
        macrophage_list[i].macrophage_age(tstep, step_count);
        // if in necrotic core, die
        if(macrophage_list[i].calcDistance(tumorCenter) < necroticRadius){
            macrophage_list[i].state = -1;
        }
    }
    
    
    /*
     * cell proliferation
     * remove cell if out of bounds
     */

    cell_list[i].prolifState();
    
    
    
    std::array<double, 3> newLoc = cell_list[i].proliferate(tstep);
    
    if(newLoc[2] == 1){
        if(cell_list[i].type == 3){
            int phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 
            std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
            if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                std::cerr << "WARNING INTERNAL CELL FUNCTIONS: t_cell_phenotype_Trajectory is empty!" << std::endl; 
            }
            cell_list.push_back(Cell({newLoc[0], newLoc[1]},
                                    cell_list.size(),
                                    cellParams,
                                    cell_list[i].type, trajec_phenotype, step_count));
        }
        else{
            cell_list.push_back(Cell({newLoc[0], newLoc[1]},
                                    cell_list.size(),
                                    cellParams,
                                    cell_list[i].type, tCellPhenotypeTrajectory_1));

        }
        cell_list[cell_list.size() - 1].inherit(cell_list[i].inheritanceProperties());
    }






















    int numCells = (cancer_list.size() + cd4_list.size() + cd8_list.size() + macrophage_list.size());

    for(int i=0; i<numCells; ++i){
        //cell_list[i].prolifState();
        //std::array<double, 3> newLoc = cell_list[i].proliferate(tstep);

        cancer_list[i].cancer_prolifState();
        std::array<double, 3> newLoc = cancer_list[i].cancer_proliferate(tstep);



        if(newLoc[2] == 1){
            for (int i = 0; i < cd8_list.size(); i++) {
                size_t phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size());
                //int phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 
                std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                    std::cerr << "WARNING INTERNAL CELL FUNCTIONS: t_cell_phenotype_Trajectory is empty!" << std::endl; 
                }
                cd8_list.push_back(CD8Cell(cellParams, size_t init_tstamp, {newLoc[0], newLoc[1]}, int idx, 3, trajec_phenotype));
                // need to figure out init_tstamp and idx
            }

            for (int i = 0; i < cancer_list.size(); i++) {
                cancer_list.push_back(CancerCell(cellParams, size_t init_tstamp, {newLoc[0], newLoc[1]}, int idx, 0));
                // need to figure out init_tstamp and idx
            }


            for (int i = 0; i < cd4_list.size(); i++) {
                cd4_list.push_back(CD4Cell(cellParams, size_t init_tstamp, {newLoc[0], newLoc[1]}, int idx, 2));
                // need to figure out init_tstamp and idx
            }

            for (int i = 0; i < macrophage_list.size(); i++) {
                macrophage_list.push_back(MacrophageCell(cellParams, size_t init_tstamp, {newLoc[0], newLoc[1]}, int idx, 1));
                // need to figure out init_tstamp and idx
            }

            cell_list[cell_list.size() - 1].inherit(cell_list[i].inheritanceProperties());
        }
    




            

          

    



                





            if(cell_list[i].type == 3){
                int phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 
                std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                    std::cerr << "WARNING INTERNAL CELL FUNCTIONS: t_cell_phenotype_Trajectory is empty!" << std::endl; 
                }
                cell_list.push_back(Cell({newLoc[0], newLoc[1]},
                                     cell_list.size(),
                                     cellParams,
                                     cell_list[i].type, trajec_phenotype, step_count));
            }
            else{
                

            }
            cell_list[cell_list.size() - 1].inherit(cell_list[i].inheritanceProperties());
        }
    

    // remove dead cells
    std::vector<int> dead;
    for(int i=0; i<cell_list.size(); ++i){
        if(cell_list[i].state == -1){
            dead.push_back(i);
        }
    }
    std::reverse(dead.begin(), dead.end());
    for(auto &i : dead){
        cell_list.erase(cell_list.begin()+i);
    }

    for(int i=0; i<cell_list.size(); ++i){
        cell_list[i].updateID(i);
        if(cell_list[i].state == -1){
            throw std::runtime_error("Environment::internalCellFunctions -> dead cell not removed");
        }
    }
}



void Environment::internalCellFunctions(double tstep, size_t step_count) {
    /*
     * cell death via aging
     * cell proliferation
     * remove cell if out of bounds
     */
    int numCells = (cancer_list.size() + cd4_list.size() + cd8_list.size() + macrophage_list.size());
    for(int i=0; i<numCells; ++i){
        cell_list[i].age(tstep, step_count);
        // if in necrotic core, die
        
        
        if(cell_list[i].calcDistance(tumorCenter) < necroticRadius){
            cell_list[i].state = -1;
        }

        cancer_list[i].prolifState();
        std::array<double, 3> newLoc = cancer_list[i].proliferate(tstep);
        if(newLoc[2] == 1){
            
            if(cell_list[i].type == 3){

                int phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 
                std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                cd8_list.push_back(CD8Cell({newLoc[0], newLoc[1]},
                                     cell_list.size(),
                                     cellParams,
                                     cell_list[i].type, trajec_phenotype, step_count));
            }
            else{ 
                cd8_list.push_back(CD8Cell({newLoc[0], newLoc[1]},
                                     cell_list.size(),
                                     cellParams,
                                     cell_list[i].type, tCellPhenotypeTrajectory_1));

            }
            cell_list[cell_list.size() - 1].inherit(cell_list[i].inheritanceProperties());
        }
    }
}


void Environment::internalCancerCellFunctions(double tstep, size_t step_count) {
    /*
     * cell death via aging
     * cell proliferation
     * remove cell if out of bounds
     */
    int numCells = (cancer_list.size() + cd4_list.size() + cd8_list.size() + macrophage_list.size());
    for(int i=0; i<numCells; ++i){
        cancer_list[i].cancer_age(tstep, step_count);
        // if in necrotic core, die
        if(cancer_list[i].calcDistance(tumorCenter) < necroticRadius){
            cancer_list[i].state = -1;
        }
        
        cancer_list[i].cancer_prolifState();
        std::array<double, 3> newLoc = cancer_list[i].cancer_proliferate(tstep);
        
        if(newLoc[2] == 1){
            cd8_list.push_back(CD8Cell({newLoc[0], newLoc[1]}, cell_list.size(), cellParams, cell_list[i].type, tCellPhenotypeTrajectory_1));
            
            
            
            
            cancer_list.push_back(CancerCell(cellParams, size_t init_tstamp, {newLoc[0], newLoc[1]}, int idx, 0));
        }

        cancer_list[cancer_list.size() - 1].inherit(cancer_list[i].inheritanceProperties());

    }

    // remove dead cells
    std::vector<int> dead;
    for(int i=0; i<cancer_list.size(); ++i){
        if(cancer_list[i].state == -1){
            dead.push_back(i);
        }
    }

    std::reverse(dead.begin(), dead.end());
    for(auto &i : dead){
        cancer_list.erase(cancer_list.begin()+i);
    }

    for(int i=0; i<cancer_list.size(); ++i){
        cancer_list[i].updateID(i);
        if(cancer_list[i].state == -1){
            throw std::runtime_error("Environment::internalCancerCellFunctions -> dead cell not removed");
        }
    }
}

void Environment::internalCD4CellFunctions(double tstep, size_t step_count) {

    // aging, out of bounds removal 

    int numCells = (cancer_list.size() + cd4_list.size() + cd8_list.size() + macrophage_list.size());
    for(int i=0; i<numCells; ++i){
        cd4_list[i].cd4_age(tstep, step_count);
        // if in necrotic core, die
        if(cd4_list[i].calcDistance(tumorCenter) < necroticRadius){
            cd4_list[i].state = -1;
        }
        /*cd4_list[i].prolifState();
        std::array<double, 3> newLoc = cd4_list[i].proliferate(tstep)
        if(newLoc[2] == 1){
            cd4_list.push_back(CancerCell(cellParams, size_t init_tstamp, {newLoc[0], newLoc[1]}, int idx, 0));
        }
        */
        cd4_list[cd4_list.size() - 1].inherit(cd4_list[i].inheritanceProperties());

    }
    



    // remove dead cells
    std::vector<int> dead;
    for(int i=0; i<cd4_list.size(); ++i){
        if(cd4_list[i].state == -1){
            dead.push_back(i);
        }
    }

    std::reverse(dead.begin(), dead.end());
    for(auto &i : dead){
        cd4_list.erase(cd4_list.begin()+i);
    }

    for(int i=0; i<cd4_list.size(); ++i){
        cd4_list[i].updateID(i);
        if(cd4_list[i].state == -1){
            throw std::runtime_error("Environment::internalCD4CellFunctions -> dead cell not removed");
        }
    }
}

void Environment::internalCD8CellFunctions(double tstep, size_t step_count) {

    // add in aging, out of bounds removal 

    int numCells = (cancer_list.size() + cd4_list.size() + cd8_list.size() + macrophage_list.size());
    for(int i=0; i<numCells; ++i){
        cd8_list[i].age(tstep, step_count);
        // if in necrotic core, die
        if(cell_list[i].calcDistance(tumorCenter) < necroticRadius){
            cell_list[i].state = -1;
        }
        cancer_list[i].prolifState();
        std::array<double, 3> newLoc = cancer_list[i].proliferate(tstep);
        if(newLoc[2] == 1){
            
            if(cell_list[i].type == 3){

                int phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size()); 
                std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                cd8_list.push_back(CD8Cell({newLoc[0], newLoc[1]},
                                     cell_list.size(),
                                     cellParams,
                                     cell_list[i].type, trajec_phenotype, step_count));
            }
            else{ 
                cd8_list.push_back(CD8Cell({newLoc[0], newLoc[1]},
                                     cell_list.size(),
                                     cellParams,
                                     cell_list[i].type, tCellPhenotypeTrajectory_1));

            }
            cell_list[cell_list.size() - 1].inherit(cell_list[i].inheritanceProperties());
        }
    }


    // remove dead cells
    std::vector<int> dead;
    for(int i=0; i<cd8_list.size(); ++i){
        if(cd8_list[i].state == -1){
            dead.push_back(i);
        }
    }

    std::reverse(dead.begin(), dead.end());
    for(auto &i : dead){
        cd8_list.erase(cd8_list.begin()+i);
    }

    for(int i=0; i<cd8_list.size(); ++i){
        cd8_list[i].updateID(i);
        if(cd8_list[i].state == -1){
            throw std::runtime_error("Environment::internalCD8CellFunctions -> dead cell not removed");
        }
    }
}

void Environment::internalMacrophageCellFunctions(double tstep, size_t step_count) {

    // add in aging, proliferation, out of bounds removal 


    // remove dead cells
    std::vector<int> dead;
    for(int i=0; i<macrophage_list.size(); ++i){
        if(macrophage_list[i].state == -1){
            dead.push_back(i);
        }
    }

    std::reverse(dead.begin(), dead.end());
    for(auto &i : dead){
        macrophage_list.erase(macrophage_list.begin()+i);
    }

    for(int i=0; i<macrophage_list.size(); ++i){
        macrophage_list[i].updateID(i);
        if(macrophage_list[i].state == -1){
            throw std::runtime_error("Environment::internalMacrophageCellFunctions -> dead cell not removed");
        }
    }
}



void Environment::runCells(double tstep, size_t step_count) {
    neighborInfluenceInteractions(tstep, step_count);
    calculateForces(tstep);
    //internalCellFunctions(tstep, step_count);
    internalCancerCellFunctions(tstep, step_count);
    internalCD4CellFunctions(tstep, step_count);
    internalCD8CellFunctions(tstep, step_count);
    internalMacrophageCellFunctions(tstep, step_count);
}
