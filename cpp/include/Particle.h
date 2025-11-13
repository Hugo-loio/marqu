#ifndef PARTICLE_H
#define PARTICLE_H

#include "Configuration.h"

namespace marqu{
  class Particle{
    public:
      Particle(Configuration && configuration, bool type, double eventRate) noexcept;
      Particle(const Configuration & configuration, bool type, double eventRate);

      Configuration configuration;
      bool type; //true particle, false antiparticle
      double eventRate;
  }; 

  bool operator==(const Particle &, const Particle &);
  Particle opposite(const Particle & particle);
}

#endif
