#include "Configuration.h"

namespace marqu{
  double pauliOperator(Axis, int site, const Configuration &);

  double pauliString(const std::vector<Axis> &, 
      const std::vector<std::size_t> & sites, const Configuration &);

  double magnetization(Axis, const Configuration &);
}
