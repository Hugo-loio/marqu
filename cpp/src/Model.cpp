#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "Model.h"

marqu::Model::Model(int N) : N(N){
}

int marqu::Model::findNSites(const std::string & name){
  std::string filename = basePath + name + "/props.csv";
  std::ifstream file(filename);
  if (!file.is_open()){
    throw std::runtime_error("Error: Could not open file: " + filename);
  }

  std::string line;
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string key, val;

    if (!std::getline(ss, key, ',')) continue;

    if (key == "nsites") {
      if (!std::getline(ss, val, ',')) {
	throw std::runtime_error("Error: Missing value for key: " + key);
      }

      try {
	return std::stoi(val);
      } catch (...) {
	throw std::runtime_error("Error: Invalid int for value: " + val);
      }
    }
  }

  throw std::runtime_error("Error: site number not found");

}

void marqu::Model::addRateMatrix(
    const std::string & name, 
    const std::vector<std::vector<int>> & sites)
{
  this->sites.push_back(sites);
  int nSites = findNSites(name);

  std::string filename = basePath + name + "/M.csv";
  std::ifstream file(filename);
  if (!file.is_open()){
    throw std::runtime_error("Error: Could not open file: " + filename);
  }

  int size = 1;
  for(int i = 0; i < nSites; i++) size *= 6;

  localMMinus.emplace_back();
  localMPlus.emplace_back();

  std::string line;
  for (int i = 0; i < size; i++) {
    if (!getline(file, line)) {
      throw std::runtime_error("Error: Not enough rows in CSV file");
    }

    std::stringstream ss(line);
    std::string value;
    double cumulativePlusRate = 0;
    double cumulativeMinusRate = 0;
    localMMinus.back().emplace_back();
    localMPlus.back().emplace_back();

    for (int j = 0; j < size; j++) {
      if (!getline(ss, value, ',')) {
	throw std::runtime_error(
	    "Error: Not enough columns in row " + std::to_string(i));
      }

      try {
	double val = stod(value);
	if(val == 0) continue;
	double & rate = (val < 0) ? cumulativeMinusRate : cumulativePlusRate;
	auto & localM = (val < 0) ? localMMinus : localMPlus;
	rate += val;
	localM.back()[i].emplace_back(marqu::Configuration(j, nSites), rate);
      } catch (...) {
	throw std::runtime_error("Error: Invalid float value at row " +
	    std::to_string(i) + ", col " + std::to_string(j));
      }
    }
  }
}
