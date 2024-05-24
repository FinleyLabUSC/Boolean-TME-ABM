#ifndef MODELUTIL_H
#define MODELUTIL_H

#include <random>
#include <filesystem> 
#include <string> 
#include <vector>
#include <iostream>
#include <cmath> 

int getRandomNumber(int maxNumber);

int countNumFiles(std::string dirPath);

std::vector<std::string> get2dvecrow(std::vector<std::vector<std::string>>& vec, size_t row_idx); 
std::array<double, 2> getBrownianUpdate(std::array<double, 2>& old_loc); 
std::vector<double> getDiff(std::vector<double>& a, std::vector<double>& b);

double getScaleFactor(double dp, double distance, double influenceLimit); 



#endif /* RANDOMUTIL_H */