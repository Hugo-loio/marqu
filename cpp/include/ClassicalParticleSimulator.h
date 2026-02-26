#ifndef CLASSICALPARTICLESIMULATOR_H
#define CLASSICALPARTICLESIMULATOR_H

#include <memory>

#include "BaseParticleSimulator.h"
#include "Particle.h"

namespace marqu{
  class ClassicalParticleSimulator : public BaseParticleSimulator{
    public:
      //ClassicalParticleSimulator() {};
      //~ClassicalParticleSimulator() {};

      void discreteTimeStep(double dt) override;
      double gillespieTimeStep() override;

      void displayParticles() const override;

      int initialize(int particleNumber = 1, bool removeStatic = true) override; 

    protected:
      std::unique_ptr<Configuration> config;
      double rate;

      void checkModel(const Model & model) const override;
      virtual void markovStep(); 
      void updateObservable(const Configuration & configuration);

      // Only used for initialization
      void clearParticles() override {};
      void addParticle(const Particle & particle) override;
      void addParticle(Particle && particle) override {addParticle(particle);};

      //bool warnDiscrete = true;
  };
}

#endif
