#include "Environment.h"

void Environment::loadParams() {
    std::ifstream dataCP(saveDir+"/params/cellParams.csv");
    std::string line;
    while(std::getline(dataCP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        cellParams.push_back(parsedRow);
    }
    dataCP.close();

    std::ifstream dataRP(saveDir+"/params/recParams.csv");
    while(std::getline(dataRP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        recParams.push_back(parsedRow[0]);
    }
    dataRP.close();

    std::ifstream dataEP(saveDir+"/params/envParams.csv");
    while(std::getline(dataEP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        envParams.push_back(parsedRow[0]);
    }
    dataEP.close();
}

void Environment::save(double tstep, double tstamp) {

    std::ofstream myfile;
    std::string day_dir = saveDir + "/cellLists/day_" + std::to_string(day);
    std::string str = "mkdir -p " + day_dir;
    const char *command = str.c_str();
    std::system(command);

    myfile.open(saveDir+"/necroticRadius.csv");
    myfile << necroticRadius << std::endl;
    myfile.close();

    myfile.open(day_dir+"/cells.csv");
    for(auto &cell : cell_list){
        //logging cell location, type and state 
        if(cell.type == 3){ //cd8 t cell
            size_t idx = (tstamp - cell.init_time)*cell.pTypeStateTransition; 
            if(idx > cell.t_cell_phenotype_Trajectory.size() - 1){
                //we are outside of the array and want to get the last 
                std::string phenotype_char; 
                if (cell.t_cell_phenotype_Trajectory.empty() || (cell.t_cell_phenotype_Trajectory.size() == 0)){
                    std::cerr << "WARNING LOGGING: t_cell_phenotype_Trajectory is empty!" << std::endl;
                    phenotype_char = 'E'; 
                }
                else{
                    phenotype_char = cell.t_cell_phenotype_Trajectory.back();         
                }

                myfile << cell.type << ","
                    << cell.x[0] << ","
                    << cell.x[1] << ","
                    << cell.radius << ","
                    << phenotype_char << ","
                    << cell.pdl1 << std::endl;
            }
            else{
                //we are in the array and can index
                std::string pType = cell.t_cell_phenotype_Trajectory[idx]; 
                myfile << cell.type << ","
                    << cell.x[0] << ","
                    << cell.x[1] << ","
                    << cell.radius << ","
                    << pType << ","
                    << cell.pdl1 << std::endl;
            }
        }
        else{
            myfile << cell.type << ","
               << cell.x[0] << ","
               << cell.x[1] << ","
               << cell.radius << ","
               << cell.state << ","
               << cell.pdl1 << std::endl;
        }
    }
    myfile.close();

    /*myfile.open(day_dir+"/cancerCells.csv");
    for(auto &cell : cell_list){
        if(cell.type == 0) {
            myfile << cell.x[0] << "," << cell.x[1] << "," << cell.radius << "," << cell.migrationSpeed << "," << cell.pdl1 << "," << static_cast<double>(cell.compressed) << std::endl;
        }
    }
    myfile.close();

    myfile.open(day_dir+"/cd8Cells.csv");
    for(auto &cell : cell_list){
        if(cell.type == 3) {
            int state = 0;
            if (cell.state == 6) {
                state = 0;
            } else if (cell.state == 7) {
                state = 1;
            }
            myfile << cell.x[0] << "," << cell.x[1] << "," << cell.radius << "," << state << std::endl;
        }
    }
    myfile.close();

    myfile.open(day_dir+"/cd4Cells.csv");
    for(auto &cell : cell_list){
        if(cell.type == 2) {
            int state = 0;
            if (cell.state == 4) {
                state = 0;
            } else if (cell.state == 5) {
                state = 1;
            }
            myfile << cell.x[0] << "," << cell.x[1] << "," << cell.radius << "," << state << std::endl;
        }
    }
    myfile.close();

    myfile.open(day_dir+"/mCells.csv");
    for(auto &cell : cell_list){
        if(cell.type == 1) {
            int state = 0;
            if(cell.state == 0){state = 0;}
            if(cell.state == 1){state = 1;}
            if(cell.state == 2){state = 2;}
            myfile << cell.x[0] << "," << cell.x[1] << "," << cell.radius << "," << state << std::endl;
        }
    }
    myfile.close();*/

    if(day == 0){
        day++;
        return;}

    myfile.open(saveDir+"/cancerTS.csv");
    myfile << cancerTS[0];
    for(int i=1; i<cancerTS.size(); ++i){
        myfile << "," << cancerTS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/cd8TS.csv");
    myfile << cd8TS[0];
    for(int i=1; i<cd8TS.size(); ++i){
        myfile << "," << cd8TS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/cd4TS.csv");
    myfile << cd4TS[0];
    for(int i=1; i<cd4TS.size(); ++i){
        myfile << "," << cd4TS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/m0TS.csv");
    myfile << m0TS[0];
    for(int i=1; i<m0TS.size(); ++i){
        myfile << "," << m0TS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/m1TS.csv");
    myfile << m1TS[0];
    for(int i=1; i<m1TS.size(); ++i){
        myfile << "," << m1TS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/m2TS.csv");
    myfile << m2TS[0];
    for(int i=1; i<m2TS.size(); ++i){
        myfile << "," << m2TS[i];
    }
    myfile << std::endl;
    myfile.close();

    myfile.open(saveDir+"/radiusTS.csv");
    myfile << radiusTS[0];
    for(int i=1; i<radiusTS.size(); ++i){
        myfile << "," << radiusTS[i];
    }
    myfile << std::endl;
    myfile.close();

    ++day;

    /*int cd4idx = 0;
    int cd8idx = 0;
    int midx = 0;
    for(auto& c : cell_list){
        if(c.type == 1){
            myfile.open(saveDir+"/tracking/m/cell_"+std::to_string(midx)+".csv");
            for(auto & x : c.locHist){
                myfile << x[0] << "," << x[1] << std::endl;
            }
            myfile.close();
            ++midx;
        } else if(c.type == 2){
            myfile.open(saveDir+"/tracking/cd4/cell_"+std::to_string(cd4idx)+".csv");
            for(auto & x : c.locHist){
                myfile << x[0] << "," << x[1] << std::endl;
            }
            myfile.close();
            ++cd4idx;
        } else if(c.type == 3){
            myfile.open(saveDir+"/tracking/cd8/cell_"+std::to_string(cd8idx)+".csv");
            for(auto & x : c.locHist){
                myfile << x[0] << "," << x[1] << std::endl;
            }
            myfile.close();
            ++cd8idx;
        }
    }*/
}
