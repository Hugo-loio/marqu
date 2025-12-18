#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include "Model.h"

marqu::Model::Model(std::size_t N) : N(N){
}

std::size_t marqu::Model::findNSites(const std::string & name) const{
  std::string filename = basePath + name + "/props.csv";
  std::ifstream file(filename);
  if (!file.is_open()){
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::string line;
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string key, val;

    if (!std::getline(ss, key, ',')) continue;

    if (key == "nsites") {
      if (!std::getline(ss, val, ',')) {
	throw std::runtime_error("Missing value for key: " + key);
      }

      try {
	return std::stoi(val);
      } catch (...) {
	throw std::runtime_error("Invalid int for value: " + val);
      }
    }
  }

  throw std::runtime_error("site number not found");

}

void marqu::Model::checkSites(const std::vector<std::vector<std::size_t>> & sites) const{
  for(const std::vector<std::size_t> & group : sites){
    for(std::size_t site : group){
      if(site < 0 || site >= N){
	throw std::invalid_argument("site " + std::to_string(site) + " not valid");
      }
    }
  }
}

void marqu::Model::addRateMatrix(
    const std::string & name, 
    const std::vector<std::vector<std::size_t>> & sites)
{
  checkSites(sites);
  siteCollection.push_back(sites);
  std::size_t nSites = findNSites(name);

  std::string filename = basePath + name + "/M.csv";
  std::ifstream file(filename);
  if (!file.is_open()){
    throw std::runtime_error("Could not open file: " + filename);
  }

  std::size_t size = 1;
  for(std::size_t i = 0; i < nSites; i++) size *= 6;

  localMMinus.emplace_back();
  localMPlus.emplace_back();

  std::string line;
  for (std::size_t i = 0; i < size; i++) {
    if (!getline(file, line)) {
      throw std::runtime_error("Not enough rows in CSV file");
    }

    std::stringstream ss(line);
    std::string value;
    double cumulativePlusRate = 0;
    double cumulativeMinusRate = 0;
    localMMinus.back().emplace_back();
    localMPlus.back().emplace_back();

    for (std::size_t j = 0; j < size; j++) {
      if (!getline(ss, value, ',')) {
	throw std::runtime_error(
	    "Not enough columns in row " + std::to_string(i));
      }

      try {
	double val = stod(value);
	if(val == 0) continue;
	double & rate = (val < 0) ? cumulativeMinusRate : cumulativePlusRate;
	auto & localM = (val < 0) ? localMMinus : localMPlus;
	rate += std::abs(val);
	localM.back()[i].emplace_back(marqu::Configuration(j, nSites), rate);
      } catch (...) {
	throw std::runtime_error("Invalid float value at row " +
	    std::to_string(i) + ", col " + std::to_string(j));
      }
    }
  }
}

void marqu::Model::clear(){
  siteCollection.clear();
  localMPlus.clear();
  localMMinus.clear();
}
