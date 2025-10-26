#include "ListParticleSimulator.h"
#include <iostream>

int marqu::ListParticleSimulator::occupation(int configuration) const {
  int result = 0;
  for (const Particle & particle : particles) {
    if (particle.configuration == configuration) {
      if(particle.type){
	result++;
      }
      else{
	result--;
      }
    }
  }
  return result;
}

std::vector<int> marqu::ListParticleSimulator::occupations(const std::vector<int>& configurations) const {
  std::vector<int> results(configurations.size(), 0);
  for (const Particle & particle : particles) {
    for (size_t i = 0; i < configurations.size(); i++) {
      if(particle.configuration == configurations[i]){
	if(particle.type){
	  results[i]++;
	}
	else{
	  results[i]--;
	}
      }
    }
  }
  return results;
}

void marqu::ListParticleSimulator::discreteTimeStep(double dt) {
  for (auto it = particles.begin(); it != particles.end(); it++) {
    if (uni_dist(gen) < it->eventRate * dt) {
      markovStep(it);
    }
  }
}

double marqu::ListParticleSimulator::gillespieTimeStep() {
  double trate = totalRate();
  std::exponential_distribution<> exp_dist(trate);
  double dt = exp_dist(gen);
  double randRate = uni_dist(gen)*trate;
  double integratedRate = 0;
  for (auto it = particles.begin(); it != particles.end(); it++) {
    integratedRate += it->eventRate;
    if(integratedRate >= randRate){
      markovStep(it);
      break;
    }
  }
  return dt;
}

// More efficient version for uniform event rates
//double marqu::ListParticleSimulator::gillespieTimeStep() {
//  std::exponential_distribution<> exp_dist(totalRate());
//  double dt = exp_dist(gen);
//
//  int index = int(uni_dist(gen) * particleNumber);
//  auto it = particles.begin();
//  std::advance(it, index);
//
//  markovStep(it);
//  return dt;
//}

void marqu::ListParticleSimulator::displayParticles() const{
  int index = 0;
  for (auto it = particles.begin(); it != particles.end(); ++it, ++index) {
    std::cout << "Particle " << index << ": configuration " << it->configuration << ", type " << it->type << std::endl;
  }
}

void marqu::ListParticleSimulator::clearParticles() {
  particles.clear();
  particleNumber = 0;
}

void marqu::ListParticleSimulator::removeParticle(std::list<marqu::Particle>::iterator it) {
  particles.erase(it);
  particleNumber--;
}

void marqu::ListParticleSimulator::addParticle(int configuration, bool type) {
  particles.push_back({configuration, type, eventRate(configuration)});
  particleNumber++;
}

void marqu::ListParticleSimulator::addParticle(Particle particle) {
  particles.push_back(particle);
  particleNumber++;
}

void marqu::ListParticleSimulator::moveParticle(std::list<marqu::Particle>::iterator it, int configuration, bool switchType) {
  Particle & particle = *it;
  particle.configuration = configuration;
  particle.eventRate = eventRate(configuration);
  if(switchType){
    particle.type = ! particle.type;
  }
}

std::list<marqu::Particle>::iterator marqu::ListParticleSimulator::findFirst(const Particle & particle){
  for (auto it = particles.begin(); it != particles.end(); ++it) {
    if (particle == *it) {
      return it;
    }
  }
  return particles.end();
}

std::list<marqu::Particle>::iterator marqu::ListParticleSimulator::findLast(const Particle & particle){
  for (auto rit = particles.rbegin(); rit != particles.rend(); ++rit) {
    if (particle == *rit) {
      auto it = rit.base(); 
      return --it;         
    }
  }
  return particles.end(); 
}
