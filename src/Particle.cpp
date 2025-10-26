#include "Particle.h"

marqu::Particle::Particle(marqu::Configuration && configuration, bool type, double eventRate) noexcept
: configuration(std::move(configuration)),
  type(type),
  eventRate(eventRate)
{
}

marqu::Particle::Particle(const marqu::Configuration & configuration, bool type, double eventRate):
  configuration(configuration),
  type(type),
  eventRate(eventRate)
{
}

bool marqu::operator==(const marqu::Particle & p1, const marqu::Particle & p2){
  return (p1.configuration.flattened() == p2.configuration.flattened()) 
    && (p1.type == p2.type);
}

marqu::Particle marqu::opposite(const Particle & particle){
  Particle antiParticle(particle);
  antiParticle.type = ! antiParticle.type;
  return antiParticle;
}
