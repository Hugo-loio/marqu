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

/*
Notes:
 
1. It could be faster to save integrated event rates which take O(particleNumber)
time to update on average, whenever a particle is removed or moved and O(1)
when a particle is added.
The advantage is that choosing the particle in the gillespie time step would be
O(log(particleNumber)) with a binary search.
For the discrete time step this is not efficient.
Actually in principle segment trees would allow you to integrate in log(particleNumber) time, so with a binary search you could do the gillespie in log^2(particleNumber) time

*/
#endif
