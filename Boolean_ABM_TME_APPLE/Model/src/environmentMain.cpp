#include "Environment.h"
#include "ModelUtil.h"


namespace fs = std::__fs::filesystem; 


void read_trajectory_csv(std::vector<std::vector<std::string>>& csv_data, std::string fpath){
    std::ifstream myfile(fpath);

    //try to read in the csv file
    if(!myfile.is_open()){
        std::cerr << "Error: Unable to open file at " << fpath << std::endl; 
        return; 
    }
    char ch;
    std::string cell;
    std::vector<std::string> current_row;
    bool inside_quotes = false;
    while (myfile.get(ch)) {
        if (ch == '"') {
            inside_quotes = !inside_quotes;
        } else if (ch == ',' && !inside_quotes) {
            current_row.push_back(cell);
            cell.clear();
        } else if (ch == '\n' && !inside_quotes) {
            current_row.push_back(cell);
            cell.clear();
            csv_data.push_back(current_row);
            current_row.clear();
        } else {
            cell += ch;
        }
    }
    // Close the file
    myfile.close(); 
}

std::vector<std::string> get_phenotype(std::vector<std::vector<std::string>>& trajectory_2d_vec){ 
    std::vector<std::string> res; 

    size_t max_col = trajectory_2d_vec[0].size(); 

    for(size_t i = 1; i < trajectory_2d_vec.size(); ++i){ 
        //start at i=1 because the first value is the title, change this if need be
        res.push_back(trajectory_2d_vec[i][max_col - 1]); 
    }

    return res; 

}


Environment::Environment(std::string folder, std::string set, std::string tCellTrajectoryPath): mt((std::random_device())()) {
    /*
     * initialize a simulation environment
     * -----------------------------------
     * loads the three parameter files
     *  cell parameters
     *  environment parameters
     *  recruitment parameters
     *
     * set environment variables to their respective values
     */

    saveDir = "./"+folder+"/set_"+set;
    loadParams();

    // cd8RecRate = recParams[0];
    // mRecRate = recParams[1];
    // cd4RecRate = recParams[2];
    for(int i=0; i<3; ++i){
        immuneCellRecRates.push_back(recParams[i]);
        immuneCells2rec.push_back(0.0);
    }
    immuneCellRecTypes = {3, 1, 2}; // {CD8, macrophage, CD4} -> same order as RecRates

    recDist = recParams[3];
    maxRecCytoConc = recParams[3];
    recruitmentDelay = recParams[4];

    simulationDuration = envParams[1];
    necroticGrowth = envParams[2];
    necroticForce = envParams[3];
    necroticLimit = envParams[4];

    tumorCenter = {0,0};
    tumorRadius = 0;
    necroticRadius = 0;


    
    vesRec = vesselParams[0];
    vesRad = vesselParams[1] / 2; 
    vesDens = vesselParams[2];
    vesAge = vesselParams[4]; 
    vesMu = vesselParams[5]; 
    vesDec = vesselParams[6]; 
    vesInfluenceDistance = vesselParams[7]; 

    steps = 0;
    day = 0;

    dt = 0.005;


    // attempt to load in t cell trajectory files by iterating over files in tCellTrajectoryPath
    printf("Attempting to read trajectories located in %s...\n", tCellTrajectoryPath.c_str()); 

    
    for (const auto& entry : std::__fs::filesystem::directory_iterator(tCellTrajectoryPath)){
        //make sure we are only reading csv files 
        if(entry.path().extension().string() == ".csv"){
            
            std::vector<std::vector<std::string>> trajec_csv;

            read_trajectory_csv(trajec_csv, entry.path().string()); 
            std::vector<std::string> ptype = get_phenotype(trajec_csv); 
            tCellPhenotypeTrajectory.push_back(ptype); 
        } 
        else{
            std::cout << entry.path().string() << " cannot be read as a trajectory, if this is unexpected please confirm the integrity of the .csv file \n"; 
        }
    }

    printf("Finished reading trajectories located in %s...\n", tCellTrajectoryPath.c_str()); 
    
    //attempt to load in t cell trajectory file 
    std::string trajecPath =  "t_cell_trajectory/1Tcell_Sim_ABM.csv"; 
    // std::string trajecPath2 = tCellTrajectoryPath + "/2Tcell_Sim_ABM.csv"; 
    // std::string trajecPath3 = tCellTrajectoryPath + "/3Tcell_Sim_ABM.csv"; 

    std::vector<std::vector<std::string>> trajec_csv_1;
    
    // read_trajectory_csv(trajec_csv_1, trajecPath1); 
    // // if we only care about the phenotype we can simply create a std::vector<std::string or char> of phenotypes 
    std::vector<std::string> phenotype_trajec_1 = get_phenotype(trajec_csv_1); 
    tCellPhenotypeTrajectory_1 = phenotype_trajec_1; 
    

}

void Environment::simulate(double tstep) {
    /*
     * initializes and runs a simulation
     * ---------------------------------
     * place initial tumor
     * run simulation loop
     *  recruit immune cells
     *  run cell functions
     * ends once time limit is reached or there are no more cancer cells
     */
    initializeVessels(); 
    initializeCells();
    tumorSize();

    save(tstep, steps*tstep);
    updateTimeSeries();
    std::cout << "starting simulations...\n";
    while(tstep*steps/24 < simulationDuration) {
        recruitImmuneCells(tstep, tstep*steps); 
        recruitVessels(tstep, tstep*steps);
        runCells(tstep, tstep*steps);
        tumorSize();
        necrosis(tstep);

        steps += 1;
        printStep(steps * tstep);
        updateTimeSeries();
        if (fmod(steps * tstep, 24) == 0) {
            // save every simulation day
            save(tstep, steps*tstep);
        }

        int numC = 0;
        for (auto &c: cell_list) {
            if (c.type == 0) {
                numC++;
            }
        }
        if (numC == 0) {
            save(tstep, steps*tstep);
            break;
        }
    }
}
