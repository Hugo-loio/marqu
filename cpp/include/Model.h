#ifndef MODEL_H
#define MODEL_H

#include <string>

#include "Particle.h"

namespace marqu{
  class Model{
    public:
      Model(int N);
      //~BaseModel();

      //void addRateMatrix(std::string path &, const std::vector<std::vector<int>> sites); 
      //std::pair<Configuration, Sign> randomEvent(const Particle & particle);

    protected:
      int N;
      //std::vector<std::vector<std::pair<marqu::Configuration, double>>> localMplus;
      //std::vector<std::vector<std::pair<marqu::Configuration, double>>> localMminus;

  };
}

#endif
