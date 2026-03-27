#include "VectorParticleSimulator.h"
#include <iostream>

void marqu::VectorParticleSimulator::discreteTimeStep(double dt){
  for(int i = 0; i < particleNumber; i++){
    if(uni_dist(gen) < particles[i].eventRate*dt){
      markovStep(i);
    }
  }
}

double marqu::VectorParticleSimulator::gillespieTimeStep(){
  if(totalRate <= 0) {
    throw std::runtime_error("All particles are stationary");
  }
  std::exponential_distribution<> exp_dist(totalRate);
  double dt = exp_dist(gen);
  double randRate = uni_dist(gen)*totalRate;
  double integratedRate = 0;
  for(int i = 0; i < particleNumber; i++){
    integratedRate += particles[i].eventRate;
    if(integratedRate >= randRate){
      markovStep(i);
      break;
    }
  }
  return dt;
}

// This one is more efficient for uniform event rates, consider implementing in derived classes
//double marqu::VectorParticleSimulator::gillespieTimeStep(){
//  std::exponential_distribution<> exp_dist(totalRate());
//  double dt = exp_dist(gen);
//  markovStep(int(uni_dist(gen)*particleNumber)); //Only works for uniform eventRates
//  return dt;
//}

void marqu::VectorParticleSimulator::displayParticles() const{
  for (int i = 0; i < particleNumber; i++){
    const Particle & particle = particles[i];
    std::cout << "Particle " << i << ": configuration " << particle.configuration.toString() << ", type " << particle.type << std::endl;
  }
}

void marqu::VectorParticleSimulator::clearParticles(){
  particles = std::vector<Particle>();
  trackClear();
}

void marqu::VectorParticleSimulator::removeParticle(int index){
  trackRemove(particles[index]);
  particles.erase(particles.begin() + index);
}

void marqu::VectorParticleSimulator::addParticle(Particle && particle){
  trackAdd(particle);
  particles.push_back(std::move(particle));
}

void marqu::VectorParticleSimulator::addParticle(const Particle & particle){
  trackAdd(particle);
  particles.push_back(particle);
}

void marqu::VectorParticleSimulator::moveParticle(int index, Particle && particle){
  trackMove(particles[index], particle);
  particles[index] = std::move(particle);
}

void marqu::VectorParticleSimulator::moveParticle(int index, const Particle & particle){
  trackMove(particles[index], particle);
  particles[index] = particle;
}

int marqu::VectorParticleSimulator::findFirst(const Particle & particle) const{
  for (int i = 0; i < particleNumber; i++){
    if(particle == particles[i]){
      return i;
    }
  }
  return -1;
}

int marqu::VectorParticleSimulator::findLast(const Particle & particle) const{
  for (int i = particleNumber - 1; i >= 0; i--){
    if(particle == particles[i]){
      return i;
    }
  }
  return -1;
}

/* Markov Step from the initial interpretation of the master equation
   void marqu::VectorParticleSimulator::markovStep(int particleIndex){
   Particle & particle = particles[particleIndex];
   std::pair<int, Sign> res = randomEvent(particle);
   int newConfig = res.first;
   Sign sign = res.second;
   Particle newParticle = {
   newConfig, 
   (sign == Sign::plus) ? particle.type : ! particle.type, 
   eventRate(newConfig)
   };
   Particle newAntiParticle = opposite(newParticle);
   int index = findLast(newAntiParticle);

   if(sign == Sign::plus){ 
   if(index >= 0){
   if(index > particleIndex){
   removeParticle(index);
   removeParticle(particleIndex);
   }
   else{
   removeParticle(particleIndex);
   removeParticle(index);
   }
   }
   else{
   moveParticle(particleIndex, newConfig);
   }
   }
   else{
   addParticle(particle);
   if(index >= 0){
   removeParticle(index);
   }
   else{
   addParticle(newParticle);
   }
   }
   }
   */

// After proof sum M- = sum M+ = Vc/2
void marqu::VectorParticleSimulator::markovStep(int particleIndex){
  Particle & particle = particles[particleIndex];
  std::pair<Configuration, Sign> res = randomEvent(particle);
  bool type = (res.second == Sign::plus) ? particle.type : ! particle.type;
  double rate = eventRate(res.first);
  Particle newParticle(std::move(res.first), type, rate);
  Particle newAntiParticle = opposite(newParticle);
  int index = findLast(newAntiParticle);

  if(index >= 0){
    removeParticle(index);
  }
  else{
    addParticle(newParticle);
  }
}
