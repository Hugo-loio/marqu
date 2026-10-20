#ifndef TREAPCAPPEDPARTICLESIMULATOR_H
#define TREAPCAPPEDPARTICLESIMULATOR_H

#include "Particle.h"
#include "TreapParticleSimulator.h"

namespace marqu {
class TreapCappedParticleSimulator : public TreapParticleSimulator {
public:
  TreapCappedParticleSimulator(int n) : maxParticleNumber(n) {};
  ~TreapCappedParticleSimulator() = default;

  void setCap(int n) { maxParticleNumber = n; };

protected:
  int maxParticleNumber;
  virtual void markovStep(const Particle &particle) override;
};
} // namespace marqu

#endif
