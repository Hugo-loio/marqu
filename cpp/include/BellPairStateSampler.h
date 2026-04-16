#ifndef BELLPAIRSTATESAMPLER_H
#define BELLPAIRSTATESAMPLER_H

#include <random>
#include <vector>

#include "BaseStateSampler.h"

namespace marqu{
  class BellPairStateSampler : public BaseStateSampler{
    public:
      BellPairStateSampler(
	  const std::vector<std::pair<std::size_t, std::size_t>> & pairs,
	  const std::vector<double> & phases, std::size_t N, 
	  const std::pair<Sign, Axis> & unpaired = {Sign::plus, Axis::z});
      BellPairStateSampler(
	  const std::string & pairsPath,
	  const std::vector<double> & phases, std::size_t N, 
	  const std::pair<Sign, Axis> & unpaired = {Sign::plus, Axis::z});
      BellPairStateSampler(const BellPairStateSampler & other);
      BellPairStateSampler & operator=(const BellPairStateSampler & other);
      ~BellPairStateSampler();

      void commonInit();

      Configuration sample(std::mt19937 & gen) override;
      std::unique_ptr<BaseStateSampler> clone() const override;

      std::vector<std::pair<std::size_t, std::size_t>> getPairs() const{
	return pairs;
      }


    private:
      void setPair(const std::pair<std::size_t, std::size_t> & pair,
	  Sign, Axis, Sign, Axis);

    protected:
      std::vector<std::pair<std::size_t, std::size_t>> pairs;
      std::vector<std::size_t> unpairedSites;
      std::vector<double> phases;
      std::pair<Sign, Axis> unpaired;
      std::pair<Sign, Axis> * orientations;
      std::size_t N;
      double prob = 1.0/3.0;
      std::uniform_real_distribution<double> uni_dist{0,1}; 
      std::bernoulli_distribution bool_dist{0.5};
      std::uniform_int_distribution<int> int_dist{0, 3}; 
  };
}

#endif
