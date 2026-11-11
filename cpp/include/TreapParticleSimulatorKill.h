#ifndef TREAPPARTICLESIMULATORKILL_H
#define TREAPPARTICLESIMULATORKILL_H

#include "Particle.h"
#include "TreapParticleSimulator.h"

namespace marqu {

class TreapParticleSimulatorKill : public TreapParticleSimulator {
public:
  TreapParticleSimulatorKill(double killRate = 0.0);
  virtual ~TreapParticleSimulatorKill() = default;

  void setKillRate(double rate);
  double getKillRate() const { return killRate; }

  double gillespieTimeStep() override;
  void deleteRandomParticle();

protected:
  double killRate = 0.0;

  // Selects the k-th particle uniformly across the treap in O(log K) time
  Particle getParticleByIndex(Node *tree, int count);
};

} // namespace marqu

#endif
