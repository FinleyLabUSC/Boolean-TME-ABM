#ifndef MODELUTIL_H
#define MODELUTIL_H

#include <random>
#include <filesystem> 
#include <string> 
#include <vector>
#include <iostream>

int getRandomNumber(int maxNumber);

int countNumFiles(std::string dirPath);

std::vector<std::string> get2dvecrow(std::vector<std::vector<std::string>>& vec, size_t row_idx); 
std::array<double, 2> getBrownianUpdate(std::array<double, 2>& old_loc, float mu=0, float sigma=1); 


#endif /* RANDOMUTIL_H */