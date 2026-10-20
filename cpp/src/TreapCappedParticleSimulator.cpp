#include "TreapCappedParticleSimulator.h"

void marqu::TreapCappedParticleSimulator::markovStep(const Particle &particle) {
  std::pair<Configuration, Sign> res = randomEvent(particle);
  bool type = (res.second == Sign::plus) ? particle.type : !particle.type;
  Node *config = findConfig(root, res.first);
  double rate = config ? config->rate : eventRate(res.first);

  Particle newParticle(std::move(res.first), type, rate);
  if (res.second == Sign::plus) {
    moveParticle(particle, newParticle);
  } else {
    moveParticle(particle, newParticle);
    if (particleNumber + 1 < maxParticleNumber) {
      addParticle(particle);
      addParticle(particle);
    }
    // if (particleNumber + 1 < maxParticleNumber) {
    //   addParticle(particle);
    //   addParticle(newParticle);
    // } else {
    //   moveParticle(particle, newParticle);
    // }
  }
}
