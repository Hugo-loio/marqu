#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "BaseParticleSimulator.h"

void marqu::BaseParticleSimulator::setInitialState(const std::string & orientations){
  if(initConfig != nullptr){
    delete initConfig;
  }
  initConfig = new Configuration(orientations);
}

marqu::BaseParticleSimulator::~BaseParticleSimulator(){
  if(initConfig != nullptr){
    delete initConfig;
  }
  if(model != nullptr){
    delete model;
  }
}

double marqu::BaseParticleSimulator::eventRate(const Configuration & configuration) const {
  if(model == nullptr){
    throw std::runtime_error("Rates not implemented, model not found");
  }

  double res = 0;
  for(std::size_t i = 0; i < model->siteCollection.size(); i++){
    for(const auto & sites : model->siteCollection[i]){
      int col = configuration.subFlattened(sites);
      auto & colMinus = model->localMMinus[i][col];
      auto & colPlus = model->localMPlus[i][col];
      if(! colMinus.empty()) res += colMinus.back().second;
      if(! colPlus.empty()) res += colPlus.back().second;
    }
  }
  return res;
}

std::pair<marqu::Configuration, marqu::Sign> marqu::BaseParticleSimulator::randomEvent(const Particle & particle){
  if(model == nullptr){
    throw std::runtime_error("Transitions not implemented, model not found");
  }

  Configuration config(particle.configuration);
  double rate = uni_dist(gen)*particle.eventRate;
  for(std::size_t i = 0; i < model->siteCollection.size(); i++){
    for(const auto & sites : model->siteCollection[i]){
      int col = config.subFlattened(sites);
      auto & colMinus = model->localMMinus[i][col];
      auto & colPlus = model->localMPlus[i][col];

      if(!colMinus.empty()){
	if(rate <= colMinus.back().second){
	  for(const auto & transition : colMinus){
	    if(rate <= transition.second){
	      config.subSet(transition.first, sites);
	      return std::make_pair(config, Sign::minus);
	    }
	  }
	}
	rate -= colMinus.back().second;
      }

      if(!colPlus.empty()){
	if(rate <= colPlus.back().second){
	  for(const auto & transition : colPlus){
	    if(rate <= transition.second){
	      config.subSet(transition.first, sites);
	      return std::make_pair(config, Sign::plus);
	    }
	  }
	}
	rate -= colPlus.back().second;
      }
    }
  }

  throw std::runtime_error("Rate accumulator exceeded total event rate");
}

int marqu::BaseParticleSimulator::initialize(int particleNumber, bool removeStatic){
  clearParticles();
  observableTracker = std::vector<double>(getObservableCount(), 0.0);
  observableBuffer = std::vector<double>(getObservableCount(), 0.0);

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

  return nStatic;
}

void marqu::BaseParticleSimulator::setModel(const Model & model){
  checkModel(model);
  if(this->model != nullptr){
    delete this->model;
  }
  this->model = new Model(model);
}

void marqu::BaseParticleSimulator::setModel(Model && model){
  checkModel(model);
  if(this->model != nullptr){
    delete this->model;
  }
  this->model = new Model(std::move(model));
}

void marqu::BaseParticleSimulator::checkModel(const Model & model) const{
  if(model.isClassical()){
    std::cout << "The model is classical, consider switching for the ClassicalParticleSimulator to increase efficiency" << std::endl; 
  }
}

void marqu::BaseParticleSimulator::updateObservable(const Particle & particle, bool add){
  observables(particle.configuration, observableBuffer);
  int sign = (add ^ particle.type) ? -1 : 1;
  for(std::size_t i = 0; i < observableTracker.size(); i++){
    observableTracker[i] += sign*observableBuffer[i];
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
  particleNumber -= 2*nAnnihilations;
  totalRate -= 2*nAnnihilations*configRate;
}
