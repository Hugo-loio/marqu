#include "Observables.h"

double marqu::pauliOperator(marqu::Axis axis, int site, const marqu::Configuration & configuration){
  if(configuration.axis(site) == axis){
    return configuration.sign(site)*3;
  }
  return 0;
}

double marqu::magnetization(marqu::Axis axis, const marqu::Configuration & configuration){
  double res = 0;
  for(int site = 0; site < configuration.getN(); site++){
    res += pauliOperator(axis, site, configuration);
  }
  return res;
}
