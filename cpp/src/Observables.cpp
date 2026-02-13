#include "Observables.h"
#include <cmath>

double marqu::pauliOperator(marqu::Axis axis, int site, const marqu::Configuration & configuration){
  if(configuration.axis(site) == axis){
    return configuration.sign(site)*3;
  }
  return 0;
}

double marqu::pauliString(const std::vector<Axis> & axes, 
    const std::vector<std::size_t> & sites, const Configuration & configuration){
  std::size_t length = axes.size();
  double res = std::pow(3, length);
  for(std::size_t i = 0; i < length; i++){
    if(configuration.axis(sites[i]) == axes[i]){
      res = res * configuration.sign(sites[i]); 
    }
    else return 0;
  }
  return res;
}

double marqu::magnetization(marqu::Axis axis, const marqu::Configuration & configuration){
  double res = 0;
  for(int site = 0; site < configuration.getN(); site++){
    res += pauliOperator(axis, site, configuration);
  }
  return res;
}
