#include <stdexcept>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

#include "BellPairStateSampler.h"

marqu::BellPairStateSampler::BellPairStateSampler(
    const std::string & pairsPath,
    const std::vector<double> & phases, std::size_t N, 
    const std::pair<Sign, Axis> & unpaired): 
  phases(phases), N(N), unpaired(unpaired)
{
  std::ifstream file(pairsPath);

  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + pairsPath);
  }

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) continue; 

    std::stringstream ss(line);
    std::string val1, val2;

    if (std::getline(ss, val1, ',') && std::getline(ss, val2, ',')) {
      size_t first = std::stoull(val1);
      size_t second = std::stoull(val2);
      pairs.emplace_back(first, second);
    }
  }

  file.close();

  commonInit();
}


marqu::BellPairStateSampler::BellPairStateSampler(
    const std::vector<std::pair<std::size_t, std::size_t>> & pairs,
    const std::vector<double> & phases, std::size_t N, 
    const std::pair<Sign, Axis> & unpaired): 
  pairs(pairs), phases(phases), N(N), unpaired(unpaired)
{
  commonInit();
}

marqu::BellPairStateSampler::BellPairStateSampler(
    const BellPairStateSampler & other) : BaseStateSampler(other), 
  pairs(other.pairs), unpairedSites(other.unpairedSites), phases(other.phases), 
  unpaired(other.unpaired), N(other.N), prob(other.prob), 
  uni_dist(other.uni_dist), bool_dist(other.bool_dist), int_dist(other.int_dist) 
{
  this->orientations = new std::pair<Sign, Axis>[N];
  for (std::size_t i = 0; i < N; ++i) {
    this->orientations[i] = other.orientations[i];
  }
}

marqu::BellPairStateSampler & marqu::BellPairStateSampler::operator=(
    const BellPairStateSampler & other) {
  if (this != &other) {
    delete[] orientations;

    pairs = other.pairs;
    unpairedSites = other.unpairedSites;
    phases = other.phases;
    unpaired = other.unpaired;
    N = other.N;
    prob = other.prob;
    uni_dist = other.uni_dist;
    bool_dist = other.bool_dist;
    int_dist = other.int_dist;

    orientations = new std::pair<Sign, Axis>[N];
    for (std::size_t i = 0; i < N; ++i) {
      orientations[i] = other.orientations[i];
    }
  }
  return *this;
}

marqu::BellPairStateSampler::~BellPairStateSampler(){
  delete[] orientations;
}

void marqu::BellPairStateSampler::commonInit(){
  this->orientations = new std::pair<Sign, Axis>[N];

  //std::cout << N << std::endl;
  std::vector<bool> isPaired(N, false);
  for(const auto & pair : pairs) {
    //std::cout << pair.first << " " << pair.second << std::endl;
    isPaired.at(pair.first) = true;
    isPaired.at(pair.second) = true;
  }

  for(std::size_t i = 0; i < N; ++i) {
    if(!isPaired[i]) unpairedSites.push_back(i);
  }
}

void marqu::BellPairStateSampler::setPair(
    const std::pair<std::size_t, std::size_t> & pair, 
    Sign s1, Axis a1, Sign s2, Axis a2){
  orientations[pair.first] = std::make_pair(s1, a1);
  orientations[pair.second] = std::make_pair(s2, a2);
}

marqu::Configuration marqu::BellPairStateSampler::sample(std::mt19937 & gen){
  for (std::size_t j : unpairedSites) {
    if(uni_dist(gen) < prob){
      orientations[j] = std::make_pair(unpaired.first, unpaired.second);
    }
    else{
      Sign sign = static_cast<Sign>(bool_dist(gen));
      Axis axis = static_cast<Axis>(
	  (static_cast<int>(unpaired.second) + bool_dist(gen) + 1) % 3
	  );
      orientations[j] = std::make_pair(sign, axis);
    }
  }

  Sign p = Sign::plus, m = Sign::minus;
  Axis x = Axis::x, y = Axis::y, z = Axis::z;
  double a = 1.0 / 9.0;

  for(std::size_t i = 0; i < pairs.size(); i++) {
    double r = uni_dist(gen);

    if (r < a + std::cos(phases[i])/9.0) {
      switch (int_dist(gen)) {
	case 0: setPair(pairs[i], p, x, p, x); break;
	case 1: setPair(pairs[i], m, x, m, x); break;
	case 2: setPair(pairs[i], p, y, m, y); break;
	case 3: setPair(pairs[i], m, y, p, y); break;
      }
    }
    else if (r < 2*a) {
      switch (int_dist(gen)) {
	case 0: setPair(pairs[i], p, x, m, x); break;
	case 1: setPair(pairs[i], m, x, p, x); break;
	case 2: setPair(pairs[i], p, y, p, y); break;
	case 3: setPair(pairs[i], m, y, m, y); break;
      }
    }
    else if (r < 3*a + std::sin(phases[i])/9.0) {
      switch (int_dist(gen)) {
	case 0: setPair(pairs[i], p, x, p, y); break;
	case 1: setPair(pairs[i], p, y, p, x); break;
	case 2: setPair(pairs[i], m, x, m, y); break;
	case 3: setPair(pairs[i], m, y, m, x); break;
      }
    }
    else if (r < 4*a) {
      switch (int_dist(gen)) {
	case 0: setPair(pairs[i], p, x, m, y); break;
	case 1: setPair(pairs[i], m, y, p, x); break;
	case 2: setPair(pairs[i], m, x, p, y); break;
	case 3: setPair(pairs[i], p, y, m, x); break;
      }
    }
    else if (r < 5*a) {
      if (bool_dist(gen)) setPair(pairs[i], p, z, p, z);
      else setPair(pairs[i], m, z, m, z);
    }
    else {
      Sign s1 = (bool_dist(gen)) ? m : p, s2 = (bool_dist(gen)) ? m : p;
      switch (int_dist(gen)) {
	case 0: setPair(pairs[i], s1, x, s2, z); break;
	case 1: setPair(pairs[i], s1, z, s2, x); break;
	case 2: setPair(pairs[i], s1, y, s2, z); break;
	case 3: setPair(pairs[i], s1, z, s2, y); break;
      }
    }
  }

  return Configuration(orientations, N);
}

std::unique_ptr<marqu::BaseStateSampler> 
marqu::BellPairStateSampler::clone() const{
  return std::make_unique<BellPairStateSampler>(*this);
}
