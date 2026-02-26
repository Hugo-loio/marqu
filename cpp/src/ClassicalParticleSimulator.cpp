#include <stdexcept>
#include <iostream>

#include "ClassicalParticleSimulator.h"

void marqu::ClassicalParticleSimulator::discreteTimeStep(double dt){
  if(uni_dist(gen) < rate*dt) markovStep();
}

double marqu::ClassicalParticleSimulator::gillespieTimeStep(){
  double dt = std::exponential_distribution<>{rate}(gen);
  markovStep();
  return dt;
}

void marqu::ClassicalParticleSimulator::displayParticles() const{
  std::cout << "Active configuration: " << config->toString() << std::endl;
}

void marqu::ClassicalParticleSimulator::checkModel(const Model & model) const{
  if(! model.isClassical()){
    throw std::invalid_argument("The model is not Classical");
  } 
}

void marqu::ClassicalParticleSimulator::markovStep(){
  if(model == nullptr){
    throw std::runtime_error("Transitions not implemented, model not found");
  }

  double randRate = uni_dist(gen)*rate;
  for(std::size_t i = 0; i < model->siteCollection.size(); i++){
    for(const auto & sites : model->siteCollection[i]){
      int col = config->subFlattened(sites);
      const auto & colPlus = model->localMPlus[i][col];
      if(colPlus.empty()) continue;
      if(randRate <= colPlus.back().second){
	for(const auto & transition : colPlus){
	  if(randRate > transition.second) continue;
	  config->subSet(transition.first, sites);
	  observables(*config, observableTracker);
	  for(std::size_t i = 0; i < observableTracker.size(); i++){
	    observableTracker[i] += observableBuffer[i];
	  }
	  //This could possibly be improved due to only
	  //local changes to config, but it doesn't seem that simple
	  rate = eventRate(*config); 
	  return;
	}
      }
      randRate -= colPlus.back().second;
    }
  }
  throw std::runtime_error("Rate accumulator exceeded total event rate");
}

void marqu::ClassicalParticleSimulator::addParticle(const Particle & particle){
  config = std::make_unique<Configuration>(particle.configuration);
  observables(*config, observableBuffer);
  for(std::size_t i = 0; i < observableTracker.size(); i++){
    observableTracker[i] += observableBuffer[i];
  }
  rate = particle.eventRate;
}

int marqu::ClassicalParticleSimulator::initialize(
    int particleNumber, bool removeStatic){
  if(particleNumber != 1) throw std::invalid_argument("Only single particle simulations accepted for this simulator");
  if(! removeStatic) throw std::invalid_argument("Static particles not accepted for this simulator");
  int nStatic = BaseParticleSimulator::initialize();
  observableBuffer = observableTracker;
  observables(*config, observableTracker);
  for(std::size_t i = 0; i < observableTracker.size(); i++){
    observableBuffer[i] -= observableTracker[i];
  }
  return nStatic;
}
