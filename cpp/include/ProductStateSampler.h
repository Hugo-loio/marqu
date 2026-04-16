#ifndef PRODUCTSTATESAMPLER_H
#define PRODUCTSTATESAMPLER_H

#include <random>
#include <string>

#include "BaseStateSampler.h"

namespace marqu{
  class ProductStateSampler : public BaseStateSampler{
    public:
      ProductStateSampler(const std::string & orientations);
      ProductStateSampler(const ProductStateSampler & other);
      ProductStateSampler & operator=(const ProductStateSampler & other);
      ~ProductStateSampler();

      Configuration sample(std::mt19937 & gen) override;
      std::unique_ptr<BaseStateSampler> clone() const override;

    protected:
      std::pair<Sign, Axis> * orientations;
      Configuration config;
      int N;
      double prob = 1.0/3.0;
      std::uniform_real_distribution<double> uni_dist{0,1}; 
      std::bernoulli_distribution bool_dist{0.5};
  };
}

#endif
