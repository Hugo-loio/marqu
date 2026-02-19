#include "Configuration.h"

namespace marqu{
  double pauliOperator(Axis, int site, const Configuration &);

  double pauliString(const std::vector<Axis> &, 
      const std::vector<std::size_t> & sites, const Configuration &);

  class PauliString{
    public :
      PauliString(std::vector<Axis> axes, std::vector<std::size_t> sites);
      double estimate(const Configuration & configuration);

    protected:
      std::vector<Axis> axes;
      std::vector<std::size_t> sites;
      std::size_t length;
      double norm;
  };

  double magnetization(Axis, const Configuration &);
}
