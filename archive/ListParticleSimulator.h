#ifndef LISTPARTICLECONTAINER_H
#define LISTPARTICLECONTAINER_H

#include <list>
#include "BaseParticleSimulator.h"
#include "Particle.h"

namespace marqu{

  class ListParticleSimulator : public BaseParticleSimulator{
    public:
      int occupation(int configuration) const;
      std::vector<int> occupations(const std::vector<int> & configurations) const;

      void discreteTimeStep(double dt);
      double gillespieTimeStep();

      void displayParticles() const;

    protected:
      void clearParticles();
      void addParticle(int configuration, bool type);
      void addParticle(Particle particle);

      //Not in the Base class
      void removeParticle(std::list<Particle>::iterator it);
      void moveParticle(std::list<Particle>::iterator it, int configuration, bool switchType = false);
      std::list<Particle>::iterator findFirst(const Particle &);
      std::list<Particle>::iterator findLast(const Particle &);

      virtual void markovStep(std::list<Particle>::iterator it) = 0; //User implementation
      std::list<Particle> particles;
  };
}

#endif
