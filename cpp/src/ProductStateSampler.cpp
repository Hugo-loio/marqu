#include "ProductStateSampler.h"

marqu::ProductStateSampler::ProductStateSampler(const std::string & orientations)
  : config(orientations), N(config.getN()){
    this->orientations = new std::pair<Sign, Axis>[N];
  }

marqu::ProductStateSampler::ProductStateSampler(const ProductStateSampler & other)
  : config(other.config), N(other.N), prob(other.prob), 
  uni_dist(other.uni_dist), bool_dist(other.bool_dist) 
{
  this->orientations = new std::pair<Sign, Axis>[N];
  for (int i = 0; i < N; ++i) {
    this->orientations[i] = other.orientations[i];
  }
}

marqu::ProductStateSampler & marqu::ProductStateSampler::operator=(
    const ProductStateSampler & other) {
  if (this != &other) {
    delete[] orientations; 

    config = other.config;
    N = other.N;
    prob = other.prob;
    uni_dist = other.uni_dist;
    bool_dist = other.bool_dist;

    orientations = new std::pair<Sign, Axis>[N];
    for (int i = 0; i < N; ++i) {
      orientations[i] = other.orientations[i];
    }
  }
  return *this;
}

marqu::ProductStateSampler::~ProductStateSampler(){
  delete[] orientations;
}

marqu::Configuration marqu::ProductStateSampler::sample(std::mt19937 & gen){
  for (int j = 0; j < N; ++j) {
    if(uni_dist(gen) < prob){
      orientations[j] = std::make_pair(config.sign(j), config.axis(j));
    }
    else{
      Sign sign = static_cast<Sign>(bool_dist(gen));
      Axis axis = static_cast<Axis>(
	  (static_cast<int>(config.axis(j)) + bool_dist(gen) + 1) % 3
	  );
      orientations[j] = std::make_pair(sign, axis);
    }
  }

  return Configuration(orientations, N);
}

std::unique_ptr<marqu::BaseStateSampler> 
marqu::ProductStateSampler::clone() const{
  return std::make_unique<ProductStateSampler>(*this);
}
