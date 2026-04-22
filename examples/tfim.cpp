#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

#include "MarQu/BellPairStateSampler.h"

#include "tfim.h"

models::TFIM::TFIM(std::size_t N, std::size_t dim, double J, double h, 
    double gamma) : gamma(gamma), N(N), dim(dim), J(J), h(h){
  nsites = N;
  for(size_t i = 1; i < dim; i++){
    nsites *= N;
  }

  std::string pairsPath = "pairs/dim" + std::to_string(dim); 
  pairsPath += "_N" + std::to_string(N) + ".csv";
  std::vector<double> phases(pairsPath.length(), M_PI/4);
  marqu::BellPairStateSampler sampler(pairsPath, phases, nsites);
  pairs = sampler.getPairs();
  setInitialStateSampler(sampler);

  std::vector<std::vector<std::size_t>> sites;
  for(std::size_t i = 0; i < nsites; i++){
    std::size_t multiplier = 1;
    for(std::size_t dir = 0; dir < dim; dir++){
      if(dir > 0) multiplier *= N;
      std::size_t coord = (i / multiplier) % N;
      std::size_t new_coord = (coord + 1) % N;
      std::size_t neighbor = i - (coord * multiplier) + (new_coord * multiplier);
      sites.push_back({i, neighbor});
    }
  }

  marqu::Model model(nsites);

  std::ostringstream oss;
  oss << std::setprecision(8) << std::noshowpoint << "tfim_dim" << dim << "_J" << J << "_h" << h << "_gamma" << gamma;
  model.addRateMatrix(oss.str(), sites);
  setModel(std::move(model));

  std::vector<size_t> short_pair{pairs[0].first, pairs[0].second};
  std::vector<size_t> long_pair{pairs[1].first, pairs[1].second};

  yx_short.reset(new marqu::PauliString({y, x}, short_pair));
  yy_short.reset(new marqu::PauliString({y, y}, short_pair));
  yx_long.reset(new marqu::PauliString({y, x}, long_pair));
  yy_long.reset(new marqu::PauliString({y, y}, long_pair));
}

void models::TFIM::observables(
    const marqu::Configuration & configuration, std::vector<double> & out){
  out[0] = yx_short->estimate(configuration);
  out[0] += yy_short->estimate(configuration);
  out[1] = marqu::pauliOperator(y, pairs[0].first, configuration);
  out[2] = marqu::pauliOperator(x, pairs[0].second, configuration);
  out[2] += marqu::pauliOperator(y, pairs[0].second, configuration);

  out[3] = yx_long->estimate(configuration);
  out[3] += yy_long->estimate(configuration);
  out[4] = marqu::pauliOperator(y, pairs[1].first, configuration);
  out[5] = marqu::pauliOperator(x, pairs[1].second, configuration);
  out[5] += marqu::pauliOperator(y, pairs[1].second, configuration);
}
