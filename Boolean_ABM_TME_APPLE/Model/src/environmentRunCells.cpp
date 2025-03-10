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

        // reset neighborhood and influence
        for(auto &cancer_cell: cancer_list){
            //clear neighbor list
            cancer_cell.cell_neighbors.clear();
            cancer_cell.clearInfluence(); 
            //neighbors and influence for cancer cells
            for(auto &other_cc: cancer_list){
                if(cancer_cell.getId() != other_cc.getId()){
                    cancer_cell.neighbor_updated(&other_cc, other_cc.type, other_cc.x);
                    cancer_cell.addInfluence(other_cc.x, other_cc.influenceRadius, other_cc.state);
                }
            }
            for(auto &other_mac: macrophage_list){
                if(cancer_cell.getId() != other_mac.getId()){
                    cancer_cell.neighbor_updated(&other_mac, other_mac.type, other_mac.x);
                    cancer_cell.addInfluence(other_mac.x, other_mac.influenceRadius, other_mac.state);
                }
            }
            for(auto &other_c4: cd4_list){
                if(cancer_cell.getId() != other_c4.getId()){
                    cancer_cell.neighbor_updated(&other_c4, other_c4.type, other_c4.x);
                    cancer_cell.addInfluence(other_c4.x, other_c4.influenceRadius, other_c4.state);
                }
            }
            for(auto &other_c8: cd8_list){
                if(cancer_cell.getId() != other_c8.getId()){
                    cancer_cell.neighbor_updated(&other_c8, other_c8.type, other_c8.x);
                    cancer_cell.addInfluence(other_c8.x, other_c8.influenceRadius, other_c8.state);
                }
            }
            cancer_cell.cancer_indirectInteractions(tstep); 
        }
        /*
        for(auto &mac_cell: macrophage_list){
            //clear neighbor list
            mac_cell.cell_neighbors.clear();
            mac_cell.clearInfluence(); 
            //neighbors and influence for cancer cells
            for(auto &other_cc: cancer_list){
                if(mac_cell.getId() != other_cc.getId()){
                    mac_cell.neighbor_updated(&other_cc, other_cc.type, other_cc.x);
                    mac_cell.addInfluence(other_cc.x, other_cc.influenceRadius, other_cc.state);
                }
            }
            for(auto &other_mac: macrophage_list){
                if(mac_cell.getId() != other_mac.getId()){
                    mac_cell.neighbor_updated(&other_mac, other_mac.type, other_mac.x);
                    mac_cell.addInfluence(other_mac.x, other_mac.influenceRadius, other_mac.state);
                }
            }
            for(auto &other_c4: cd4_list){
                if(mac_cell.getId() != other_c4.getId()){
                    mac_cell.neighbor_updated(&other_c4, other_c4.type, other_c4.x);
                    mac_cell.addInfluence(other_c4.x, other_c4.influenceRadius, other_c4.state);
                }
            }
            for(auto &other_c8: cd8_list){
                if(mac_cell.getId() != other_c8.getId()){
                    mac_cell.neighbor_updated(&other_c8, other_c8.type, other_c8.x);
                    mac_cell.addInfluence(other_c8.x, other_c8.influenceRadius, other_c8.state);
                }
            }
            //can add macrophage indirect interactions
        }

        for(auto &cd4_cell: cd4_list){
            //clear neighbor list
            cd4_cell.cell_neighbors.clear();
            cd4_cell.clearInfluence(); 
            //neighbors and influence for cancer cells
            for(auto &other_cc: cancer_list){
                if(cd4_cell.getId() != other_cc.getId()){
                    cd4_cell.neighbor_updated(&other_cc, other_cc.type, other_cc.x);
                    cd4_cell.addInfluence(other_cc.x, other_cc.influenceRadius, other_cc.state);
                }
            }
            for(auto &other_mac: macrophage_list){
                if(cd4_cell.getId() != other_mac.getId()){
                    cd4_cell.neighbor_updated(&other_mac, other_mac.type, other_mac.x);
                    cd4_cell.addInfluence(other_mac.x, other_mac.influenceRadius, other_mac.state);
                }
            }
            for(auto &other_c4: cd4_list){
                if(cd4_cell.getId() != other_c4.getId()){
                    cd4_cell.neighbor_updated(&other_c4, other_c4.type, other_c4.x);
                    cd4_cell.addInfluence(other_c4.x, other_c4.influenceRadius, other_c4.state);
                }
            }
            for(auto &other_c8: cd8_list){
                if(cd4_cell.getId() != other_c8.getId()){
                    cd4_cell.neighbor_updated(&other_c8, other_c8.type, other_c8.x);
                    cd4_cell.addInfluence(other_c8.x, other_c8.influenceRadius, other_c8.state);
                }
            }
            //can add cd4 indirect interactions
        }
        */

        for(auto &cd8_cell: cd8_list){
            //clear neighbor list
            cd8_cell.cell_neighbors.clear();
            cd8_cell.clearInfluence(); 
            //neighbors and influence for cancer cells
            for(auto &other_cc: cancer_list){
                if(cd8_cell.getId() != other_cc.getId()){
                    cd8_cell.neighbor_updated(&other_cc, other_cc.type, other_cc.x);
                    cd8_cell.addInfluence(other_cc.x, other_cc.influenceRadius, other_cc.state);
                }
            }
            for(auto &other_mac: macrophage_list){
                if(cd8_cell.getId() != other_mac.getId()){
                    cd8_cell.neighbor_updated(&other_mac, other_mac.type, other_mac.x);
                    cd8_cell.addInfluence(other_mac.x, other_mac.influenceRadius, other_mac.state);
                }
            }
            for(auto &other_c4: cd4_list){
                if(cd8_cell.getId() != other_c4.getId()){
                    cd8_cell.neighbor_updated(&other_c4, other_c4.type, other_c4.x);
                    cd8_cell.addInfluence(other_c4.x, other_c4.influenceRadius, other_c4.state);
                }
            }
            for(auto &other_c8: cd8_list){
                if(cd8_cell.getId() != other_c8.getId()){
                    cd8_cell.neighbor_updated(&other_c8, other_c8.type, other_c8.x);
                    cd8_cell.addInfluence(other_c8.x, other_c8.influenceRadius, other_c8.state);
                }
            }
            cd8_cell.cd8_indirectInteractions(tstep); 
        }

/*
TODO
*/
#pragma omp parallel for

    for(auto& cc: cancer_list){
        for(auto cell_ptr: cc.cell_neighbors){
            cc.cancer_directInteractions(cell_ptr->state, cell_ptr->x, cc.directInteractionProperties(cell_ptr->state, step_count), tstep); 
        }
    }
    for(auto& cd8: cd8_list){
        for(auto cell_ptr: cd8.cell_neighbors){
            cd8.cd8_directInteractions(cell_ptr->state, cell_ptr->x, cd8.directInteractionProperties(cell_ptr->state, step_count), tstep);
        }
    }

#pragma omp parallel for

    for(int i=0; i<macrophage_list.size(); ++i){
        macrophage_list[i].macrophage_differentiation(tstep);
    }

    for(int i=0; i<cd4_list.size(); ++i){
        cd4_list[i].cd4_differentiation(tstep);
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
        for(int i=0; i<cancer_list.size(); ++i){
            cancer_list[i].migrate(dt, tumorCenter);
        }

        for(int i=0; i<macrophage_list.size(); ++i){
            macrophage_list[i].migrate(dt, tumorCenter);
        }

        for(int i=0; i<cd4_list.size(); ++i){

            cd4_list[i].migrate(dt, tumorCenter);
        }

        for(int i=0; i<cd8_list.size(); ++i){
            cd8_list[i].migrate(dt, tumorCenter);
        }

        /*
        for(int i=0; i<cell_list.size(); ++i){
            cell_list[i].migrate(dt, tumorCenter);
        }
        */
        

        // calc forces
#pragma omp parallel for
        
        for(auto& cc: cancer_list){
            for(auto cell_ptr: cc.cell_neighbors){
                cc.calculateForces(cell_ptr->x, cell_ptr->radius, cell_ptr->type); 
            }
        }
        for(auto& mac: macrophage_list){
            for(auto cell_ptr: mac.cell_neighbors){
                mac.calculateForces(cell_ptr->x, cell_ptr->radius, cell_ptr->type); 
            }
        }
        for(auto& cd4: cd4_list){
            for(auto cell_ptr: cd4.cell_neighbors){
                cd4.calculateForces(cell_ptr->x, cell_ptr->radius, cell_ptr->type); 
            }
        }
        for(auto& cd8: cd8_list){
            for(auto cell_ptr: cd8.cell_neighbors){
                cd8.calculateForces(cell_ptr->x, cell_ptr->radius, cell_ptr->type); 
            }
        }
        

        // resolve forces
#pragma omp parallel for
        
        for(int i=0; i<cancer_list.size(); ++i){
            cancer_list[i].resolveForces(dt, tumorCenter, necroticRadius, necroticForce);
        }

        for(int i=0; i<macrophage_list.size(); ++i){
            macrophage_list[i].resolveForces(dt, tumorCenter, necroticRadius, necroticForce);
        }

        for(int i=0; i<cd4_list.size(); ++i){
            cd4_list[i].resolveForces(dt, tumorCenter, necroticRadius, necroticForce);
        }

        for(int i=0; i<cd8_list.size(); ++i){
            cd8_list[i].resolveForces(dt, tumorCenter, necroticRadius, necroticForce);
        }
       
    }

    // calculate overlaps and proliferation states
#pragma omp parallel for
    
    for(auto& cc: cancer_list){
        for(auto cell_ptr: cc.cell_neighbors){
            cc.calculateOverlap(cell_ptr->x, cell_ptr->radius); 
        }
        cc.isCompressed();
    }
    for(auto& mac: macrophage_list){
        for(auto cell_ptr: mac.cell_neighbors){
            mac.calculateOverlap(cell_ptr->x, cell_ptr->radius); 
        }
        mac.isCompressed(); 
    }
    for(auto& cd4: cd4_list){
        for(auto cell_ptr: cd4.cell_neighbors){
            cd4.calculateOverlap(cell_ptr->x, cell_ptr->radius); 
        }
        cd4.isCompressed(); 
    }

    for(auto& cd8: cd8_list){
        for(auto cell_ptr: cd8.cell_neighbors){
            cd8.calculateOverlap(cell_ptr->x, cell_ptr->radius); 
        }
        cd8.isCompressed();
    }
        
    }



void Environment::internalCancerCellFunctions(double tstep, size_t step_count) {
    /*
     * cancer cell death via aging
     * cancer cell proliferation
     * remove cell if out of bounds
     */

    for(int i=0; i< cancer_list.size(); ++i){
        cancer_list[i].cancer_age(tstep, step_count);
        // if in necrotic core, die
        if(cancer_list[i].calcDistance(tumorCenter) < necroticRadius){
            cancer_list[i].state = -1;
        }

        cancer_list[i].cancer_prolifState();
        std::array<double, 3> newLoc = cancer_list[i].cancer_proliferate(tstep);

        if(newLoc[2] == 1){
            cancer_list.push_back(CancerCell(cellParams, step_count, {newLoc[0], newLoc[1]}, 0));
            cancer_list[cancer_list.size() - 1].inherit(cancer_list[i].inheritanceProperties());
        }
    }

    // remove dead cancer cells
    std::vector<int> cancer_dead;
    
    for(int i=0; i<cancer_list.size(); ++i){
        if(cancer_list[i].state == -1){
            cancer_dead.push_back(i);
        }
    }
    std::reverse(cancer_dead.begin(), cancer_dead.end());
    for(auto &i : cancer_dead){
        cancer_list.erase(cancer_list.begin()+i);
    }

    for(int i=0; i<cancer_list.size(); ++i){
        cancer_list[i].updateID(i);
        if(cancer_list[i].state == -1){
            throw std::runtime_error("Environment::internalCellFunctions -> dead cancer cell not removed");
        }
    }
}
        
void Environment::internalCD4CellFunctions(double tstep, size_t step_count) {
    /*
     * cd4 cell death via aging
     * cd4 cell proliferation
     * remove cell if out of bounds
     */

    for(int i=0; i< cd4_list.size(); ++i){
        cd4_list[i].cd4_age(tstep, step_count);
        // if in necrotic core, die
        if(cd4_list[i].calcDistance(tumorCenter) < necroticRadius){
            cd4_list[i].state = -1;
        }

        std::array<double, 3> newLoc = cd4_list[i].cd4_proliferate(tstep);
        if(newLoc[2] == 1){
            cd4_list.push_back(CD4Cell(cellParams, step_count , {newLoc[0], newLoc[1]}, 2));
            cd4_list[cd4_list.size() - 1].inherit(cd4_list[i].inheritanceProperties());
        }
    }

    // remove dead cd4 cells
    std::vector<int> cd4_dead;

    for(int i=0; i<cd4_list.size(); ++i){
        if(cd4_list[i].state == -1){
            cd4_dead.push_back(i);
        }
    }
    std::reverse(cd4_dead.begin(), cd4_dead.end());
    for(auto &i : cd4_dead){
        cd4_list.erase(cd4_list.begin()+i);
    }

    for(int i=0; i<cd4_list.size(); ++i){
        cd4_list[i].updateID(i);
        if(cd4_list[i].state == -1){
            throw std::runtime_error("Environment::internalCellFunctions -> dead cd4 cell not removed");
        }
    }
}

void Environment::internalCD8CellFunctions(double tstep, size_t step_count) {
    /*
     * cd8 cell death via aging
     * cd8 cell proliferation
     * remove cell if out of bounds
     */

    for(int i=0; i< cd8_list.size(); ++i){ 
        cd8_list[i].cd8_age(tstep, step_count);
        // if in necrotic core, die
        if(cd8_list[i].calcDistance(tumorCenter) < necroticRadius){
            cd8_list[i].state = -1;
        }

        cd8_list[i].cd8_prolifState();
        std::array<double, 3> newLoc = cd8_list[i].cd8_proliferate(tstep);

        if(newLoc[2] == 1){
            int phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size() - 1); 
            std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
            
            if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                std::cout << "ATTEMPTING TO ACCESS TRAJECTORY" << phenotypeIdx << std::endl;  
                std::cerr << "WARNING INTERNAL CELL FUNCTIONS: t_cell_phenotype_Trajectory is empty!" << std::endl; 
            }

            cd8_list.push_back(CD8Cell(cellParams, step_count, {newLoc[0], newLoc[1]}, 3, trajec_phenotype));
            cd8_list[cd8_list.size() - 1].inherit(cd8_list[i].inheritanceProperties());
        }
    }


    // remove dead cd8 cells
    std::vector<int> cd8_dead;

    for(int i=0; i<cd8_list.size(); ++i){
        if(cd8_list[i].state == -1){
            cd8_dead.push_back(i);
        }
    }
    std::reverse(cd8_dead.begin(), cd8_dead.end());
    for(auto &i : cd8_dead){
        cd8_list.erase(cd8_list.begin()+i);
    }

    for(int i=0; i<cd8_list.size(); ++i){
        cd8_list[i].updateID(i);
        if(cd8_list[i].state == -1){
            throw std::runtime_error("Environment::internalCellFunctions -> dead cd8 cell not removed");
        }
    }
}

void Environment::internalMacrophageCellFunctions(double tstep, size_t step_count) {

    /*
     * macrophage death via aging
     * macrophage proliferation
     * remove cell if out of bounds
     */

    for(int i=0; i< macrophage_list.size(); ++i){
        macrophage_list[i].macrophage_age(tstep, step_count);
        // if in necrotic core, die
        if(macrophage_list[i].calcDistance(tumorCenter) < necroticRadius){
            macrophage_list[i].state = -1;
        }

        std::array<double, 3> newLoc = macrophage_list[i].macrophage_proliferate(tstep);
        if(newLoc[2] == 1){
            macrophage_list.push_back(MacrophageCell(cellParams, step_count, {newLoc[0], newLoc[1]}, 1));
            macrophage_list[macrophage_list.size() - 1].inherit(macrophage_list[i].inheritanceProperties());
        }
    }


    // remove dead macrophage cells
    std::vector<int> macrophage_dead;

    for(int i=0; i<macrophage_list.size(); ++i){
        if(macrophage_list[i].state == -1){
            macrophage_dead.push_back(i);
        }
    }
    std::reverse(macrophage_dead.begin(), macrophage_dead.end());
    for(auto &i : macrophage_dead){
        macrophage_list.erase(macrophage_list.begin()+i);
    }

    for(int i=0; i<macrophage_list.size(); ++i){
        macrophage_list[i].updateID(i);
        if(macrophage_list[i].state == -1){
            throw std::runtime_error("Environment::internalCellFunctions -> dead macrophage cell not removed");
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


