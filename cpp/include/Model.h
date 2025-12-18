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

  using InteractionSites = std::vector<std::size_t>;
  using InteractionSitesGroup = std::vector<std::vector<std::size_t>>;

  class Model{
    public:
      Model(std::size_t N);

      ~Model() = default;
      Model(const Model& other) = default;
      Model& operator=(const Model& other) = default;
      Model(Model&& other) noexcept = default;
      Model& operator=(Model&& other) noexcept = default;

      // name : Rate matrix name, created in the python lib 
      // sites : Set of sites to apply the matrix
      void addRateMatrix(const std::string & name,
	  const std::vector<std::vector<std::size_t>> & sites); 

      void clear();

      //Groups of groups of sites with the same corresponding rate matrix
      std::vector<InteractionSitesGroup> siteCollection;

      // [site set group][input config][output index]
      std::vector<RateMatrix> localMPlus;
      std::vector<RateMatrix> localMMinus;

    protected:
      const std::size_t N;
      std::string basePath = "marqu_data/rate_matrix/";

    private:
      std::size_t findNSites(const std::string & name) const;
      void checkSites(const std::vector<std::vector<std::size_t>> & sites) const;
  };
}

#endif
