#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

#include "tfim.h"
#include "MarQu/Runner.h"

int main(int argc, char * argv[]){
  if(argc != 8){
    throw std::runtime_error("Wrong number of main arguments!");
  }

  int nsamples = std::stoi(argv[1]);
  int N = std::stoi(argv[2]);
  int dim = std::stoi(argv[3]);
  double J = std::stod(argv[4]);
  double h = std::stod(argv[5]);
  double gamma = std::stod(argv[6]);
  double T = std::stod(argv[7]);

  int nsavedts = 100;

  models::TFIM simulator(N, dim, J, h, gamma);
  marqu::Runner runner(simulator);
  runner.run(T, nsamples);

  const std::vector<double> & a = runner.avgHistObservables[0];
  const std::vector<double> & a1 = runner.avgHistObservables[1];
  const std::vector<double> & a2 = runner.avgHistObservables[2];
  const std::vector<double> & b = runner.avgHistObservables[3];
  const std::vector<double> & b1 = runner.avgHistObservables[4];
  const std::vector<double> & b2 = runner.avgHistObservables[5];

  for(std::size_t i = 0; i < a.size(); i++){
    std::cout << "\nt = " << runner.times[i]  << std::setprecision(8)
      << "\n\tShort correlation = " << a[i] - a1[i] * a2[i] 
      << "\n\tLong correlation = " << b[i] - b1[i] * b2[i]  << std::endl;
  }

  return 0;
}
