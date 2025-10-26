#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "BaseParticleSimulator.h"

void marqu::BaseParticleSimulator::setInitialState(const std::string & orientations){
  initConfig = new Configuration(orientations);
}

marqu::BaseParticleSimulator::~BaseParticleSimulator(){
  if(initConfig != nullptr){
    delete initConfig;
  }
}

int marqu::BaseParticleSimulator::initialize(int particleNumber, bool removeStatic){
  clearParticles();
  if(initConfig == nullptr){
    throw std::runtime_error("Cannot initialize the simulator without setting the initial state!");
  }
  int N = initConfig->getN();
  int nStatic = 0;
  std::bernoulli_distribution bool_dist(0.5);

  for(int i = 0; i < particleNumber; ++i){
    std::pair<Sign, Axis> * orientations = new std::pair<Sign, Axis>[N];
    for (int j = 0; j < N; ++j) {
      if(uni_dist(gen) < (1.0/3)){
	orientations[j] = std::make_pair(initConfig->sign(j), initConfig->axis(j));
      }
      else{
	Sign sign = static_cast<Sign>(bool_dist(gen));
	Axis axis = static_cast<Axis>(
	    (static_cast<int>(initConfig->axis(j)) + bool_dist(gen) + 1) % 3
	    );
	orientations[j] = std::make_pair(sign, axis);
      }
    }

    Configuration config(orientations, N);
    Particle particle(config, true, eventRate(config));

    // Do not accept particles with no dynamics
    if(removeStatic and particle.eventRate < 1e-9){
      updateObservable(particle, true);
      nStatic++;
      --i;
    }
    else{
      addParticle(particle);
    }

    delete[] orientations;
  }

  initParticleNumber = particleNumber + nStatic;
  //for(double & obs : observableTracker){
  //  obs /= initParticleNumber;;
  //}

  return nStatic;
}

void marqu::BaseParticleSimulator::updateObservable(const Particle & particle, bool add){
  std::vector<double> res = observables(particle.configuration);
  int sign = (add ^ particle.type) ? -1 : 1;
  //std::cout << std::setprecision(8) << "add " << add << " sign " << sign << " res " << res[0] << " " << res[0]/initParticleNumber << " " << initParticleNumber << std::endl;
  if(observableTracker.size() == 0){
    for(const double & val : res){
      if(! add) throw std::invalid_argument("Cannot subtract to empty tracker");
      //observableTracker.push_back(sign*val/initParticleNumber);
      observableTracker.push_back(sign*val);
    }
  }
  else{
    if(observableTracker.size() != res.size()){
	throw std::runtime_error("Number of observables changed mid simulation!");
    }
    for(int i = 0; i < observableTracker.size(); i++){
      //observableTracker[i] += sign*res[i]/initParticleNumber;
      observableTracker[i] += sign*res[i];
    }
  }
}

void marqu::BaseParticleSimulator::trackClear(){
  observableTracker = std::vector<double>();
  particleNumber = 0;
  initParticleNumber = 1;
  totalRate = 0;
}

void marqu::BaseParticleSimulator::trackAdd(const Particle & particle){
  particleNumber++;
  totalRate += particle.eventRate;
  updateObservable(particle, true);
}

void marqu::BaseParticleSimulator::trackRemove(const Particle & particle){
  particleNumber--;
  totalRate -= particle.eventRate;
  updateObservable(particle, false);
}

void marqu::BaseParticleSimulator::trackMove(const Particle & pIn, const Particle & pOut){
  totalRate += pOut.eventRate - pIn.eventRate;
  updateObservable(pIn, false);
  updateObservable(pOut, true);
}

void marqu::BaseParticleSimulator::trackAnnihilate(int nAnnihilations, double configRate){
  //std::cout << 1 << " " << particleNumber << " " << totalRate << std::endl;
  particleNumber -= 2*nAnnihilations;
  totalRate -= 2*nAnnihilations*configRate;
  //std::cout << 2 << " " << particleNumber << " " << totalRate << std::endl;
}


//int marqu::BaseParticleSimulator::occupation(const std::string & orientations){
//  Configuration config(orientations);
//  return occupation(config.flattened());
//}
//
//std::vector<int> marqu::BaseParticleSimulator::occupations(const std::vector<std::string> & orientationsVec){
//  std::vector<int> configurations;
//  for(const std::string & orientations : orientationsVec){
//    Configuration config(orientations);
//    configurations.push_back(config.flattened());
//  }
//  return occupations(configurations);
//}
