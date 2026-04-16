#ifndef BASESTATESAMPLER_H
#define BASESTATESAMPLER_H

#include <random>
#include <memory>

#include "Configuration.h"

namespace marqu{
  class BaseStateSampler{
    public:
      virtual ~BaseStateSampler() = default;
      BaseStateSampler() = default;
      BaseStateSampler(const BaseStateSampler&) = default;
      BaseStateSampler& operator=(const BaseStateSampler&) = default;

      virtual Configuration sample(std::mt19937 & gen) = 0;
      virtual std::unique_ptr<BaseStateSampler> clone() const = 0;
  };
}

#endif
