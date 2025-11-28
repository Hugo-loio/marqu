#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <vector>

#include "Particle.h"

namespace marqu{
  // Output configuration, accumulated row escape rate
  using Transition = std::pair<marqu::Configuration, double>;
  using TransitionRow = std::vector<Transition>; //Only non-zero elements
  using RateMatrix = std::vector<TransitionRow>;

  using InteractionSites = std::vector<int>;
  using InteractionSitesGroup = std::vector<std::vector<int>>;

  class Model{
    public:
      Model(int N);
      //~Model();

      // name : Rate matrix name, created in the python lib 
      // sites : Set of sites to apply the matrix
      void addRateMatrix(const std::string & name,
	  const std::vector<std::vector<int>> & sites); 

      //std::pair<Configuration, Sign> randomEvent(const Particle & particle);

    protected:
      const int N;
      //Groups of groups of sites with the same corresponding rate matrix
      std::vector<InteractionSitesGroup> sites;

      // [site set group][input config][output index]
      std::vector<RateMatrix> localMPlus;
      std::vector<RateMatrix> localMMinus;

    private:
      std::string basePath = "marqu_data/rate_matrix/";
      int findNSites(const std::string & name);
  };
}

#endif
