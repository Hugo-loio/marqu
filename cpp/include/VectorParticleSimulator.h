#ifndef VECTORPARTICLESIMULATOR_H
#define VECTORPARTICLESIMULATOR_H

#include "BaseParticleSimulator.h"
#include "Particle.h"

namespace marqu{

  class VectorParticleSimulator : public BaseParticleSimulator{
    public:
      void discreteTimeStep(double dt);
      double gillespieTimeStep();

      void displayParticles() const;

    protected:
      void clearParticles();
      void addParticle(Particle && particle);
      void addParticle(const Particle & particle);

      //Not in the Base class
      virtual void removeParticle(int index);
      virtual void moveParticle(int index, Particle &&);
      virtual void moveParticle(int index, const Particle &);
      int findFirst(const Particle &) const;
      int findLast(const Particle &) const;

      virtual void markovStep(int particleIndex); 

      std::vector<Particle> particles; 
  };
}

#endif
