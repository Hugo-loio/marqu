#ifndef TREAPCAPPEDPARTICLESIMULATOR2_H
#define TREAPCAPPEDPARTICLESIMULATOR2_H

#include "Particle.h"
#include "TreapParticleSimulator.h"

namespace marqu {
class TreapCappedParticleSimulator2 : public TreapParticleSimulator {
public:
  TreapCappedParticleSimulator2(int n) : maxParticleNumber(n) {};
  ~TreapCappedParticleSimulator2() = default;

  void setCap(int n) { maxParticleNumber = n; };
  int getCap() const { return maxParticleNumber; };

  void deleteRandomParticle();

protected:
  int maxParticleNumber;
  virtual void markovStep(const Particle &particle) override;
  Particle getParticleByIndex(Node *tree, int count);
};
} // namespace marqu

#endif
