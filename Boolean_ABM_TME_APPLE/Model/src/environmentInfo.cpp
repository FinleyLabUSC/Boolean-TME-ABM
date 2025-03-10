#include "Environment.h"

void Environment::printStep(double time) {
    int numM = macrophage_list.size();
    int numT8 = cd8_list.size();
    int numT4 = cd4_list.size();
    int numC = cancer_list.size();


   std::cout << std::fixed << std::setprecision(5);
    std::cout << "Time: " << std::setw(10) << (time / 24) << " | cancer: " << std::setw(10) << numC 
          << " | cd8: " << std::setw(10) << numT8 << " | cd4: " << std::setw(10) << numT4 
          << " | macrophage: " << std::setw(10) << numM << std::endl;
}

void Environment::updateTimeSeries() {
    int numT8 = 0;
    int numT4 = 0;
    int numC = 0;

    for(auto &cell : cd8_list){
        numT8++;
    }

    for(auto &cell : cd4_list){
        numT4++;
    }

    for(auto &cell : cancer_list){
        numC++;
    }

    /*
    for(auto &cell : cell_list){
        if(cell.type == 3){
            numT8++;
        } else if(cell.type == 2){
            numT4++;
        } else if(cell.type == 0){
            numC++;
        }
    }
    */

    cancerTS.push_back(numC);
    cd8TS.push_back(numT8);
    cd4TS.push_back(numT4);

    int m0 = 0;
    int m1 = 0;
    int m2 = 0;

    for(auto &m : macrophage_list){
        if (m.state == 0) { m0++; }
        if (m.state == 1) { m1++; }
        if (m.state == 2) { m2++; }
    }

    /*
    for(auto &c : cell_list){
        if(c.type == 1) {
            if (c.state == 0) { m0++; }
            if (c.state == 1) { m1++; }
            if (c.state == 2) { m2++; }
        }
    }
    */

    m0TS.push_back(m0);
    m1TS.push_back(m1);
    m2TS.push_back(m2);

    radiusTS.push_back(tumorRadius);
}
