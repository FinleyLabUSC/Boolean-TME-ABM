#ifndef IMMUNE_MODEL_VESSEL_H
#define IMMUNE_MODEL_VESSEL_H

#include "Cell.h"
#include "ModelUtil.h"
#include <vector>
#include <cmath> 
#include <random> 
#include <string> 
#include <array> 




class Cell; 

class Vessel{

    public: 
        Vessel(std::array<double, 2> loc, double CD, double mol, double hm, double hd, int recruited, double influenceDistance);
        //constructor args 
        //(loc, diameter, vessel max overlap vessel mu, vessel decrease factor, has it been recruited)


        void treatment(bool treatmentOn);

        // overlap functions
        void calculateOverlap(std::array<double, 2> otherX, double otherRadius);
        void resetOverlap();
        void isCompressed();
        double calcDistance(std::array<double, 2> otherX);
        void neighboringCells(Cell &cell); 
        void generalLocation(double dx);
        void neighborhoodInteractions(); 


        double radius;
        std::array<double, 2> x;
        std::array<int, 2> generalX;
        double maxOverlap;
        double currentOverlap;
        int state; //0 dead, 1 alive
        double age;
        double mu;

        std::vector<std::array<double, 2>> neighborhoodLoc; 
        std::vector<Cell> neighborhood; 

        double healthyMu;
        double health;

        double influenceDistance; 

}; 





#endif

