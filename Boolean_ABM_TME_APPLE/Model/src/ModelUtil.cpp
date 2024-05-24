#include "ModelUtil.h"

namespace fs = std::__fs::filesystem; 

int getRandomNumber(int maxNumber) {
    //takes in an integer maxNumber and returns a random integer from 0 to maxNumber-1
    std::random_device rd;  // Obtain a random seed from the OS entropy device
    std::mt19937 gen(rd()); // Mersenne Twister 32-bit PRNG using seed from random device
    std::uniform_int_distribution<> distrib(0, maxNumber-1); // Define the range

    return distrib(gen); // Generate and return a random number within the specified range
}

int countNumFiles(std::string dirPath){
    //counts the number of csv files in a directory
    //input: std::string directory path
    //output: int counting number of files with a .csv extension

    int csv_file_count = 0;

    for (const auto& entry : std::__fs::filesystem::directory_iterator(dirPath)) {
        if (entry.path().extension().string() == ".csv") {
            ++csv_file_count;
        }
    }
    return csv_file_count; 
}

std::vector<std::string> get2dvecrow(std::vector<std::vector<std::string>>& vec, size_t row_idx){
    std::vector<std::string> res; 
    for(size_t i = 0; i < vec[row_idx].size(); ++i){
        
        res.push_back(vec[row_idx][i]); 
    }

    return res; 

}

std::array<double, 2> getBrownianUpdate(std::array<double, 2>& old_loc){
    //simulate brownian motion, gets old location, adds some dx and dy and 
    //returns the new location

    double x = old_loc[0]; 
    double y = old_loc[1]; 

    float mu = 0; 
    float sigma = 5; 

    static std::mt19937 gen(std::random_device{}());
    std::normal_distribution<double> dist(mu, sigma);

    double dx = dist(gen); 
    double dy = dist(gen);

    std::array<double, 2> res = {x+dx, y+dy}; 
    return res; 
}

std::vector<double> getDiff(std::vector<double>& a, std::vector<double>& b){
    std::vector<double> res; 
    for(int i =0; i < a.size(); ++i){
        res.emplace_back(b[i] - a[i]);
    }

    return res; 

}

double getScaleFactor(double dp, double distance, double influenceLimit){
    return exp((log2(dp) * log(2) * distance) / influenceLimit);
} 