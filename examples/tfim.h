#ifndef TFIM_H
#define TFIM_H

#include <memory>

#include "MarQu/TreapParticleSimulator.h"
#include "MarQu/Configuration.h"
#include "MarQu/Observables.h"

namespace models {
  class TFIM: public marqu::TreapParticleSimulator{
    public:
      TFIM(std::size_t N, std::size_t dim, double J, 
	  double h, double gamma);

      std::size_t getObservableCount() const {return 6;};
      void observables(const marqu::Configuration &, std::vector<double> & out);

    private:
      std::size_t N;  // Size length
      std::size_t dim; // Dimensions
      std::size_t nsites; // Total spins
      double J;
      double h;
      double gamma; //Damping rate
      std::vector<std::pair<std::size_t, std::size_t>> pairs;
      std::unique_ptr<marqu::PauliString> yx_short;
      std::unique_ptr<marqu::PauliString> yy_short;
      std::unique_ptr<marqu::PauliString> yx_long;
      std::unique_ptr<marqu::PauliString> yy_long;

      marqu::Axis x = marqu::Axis::x;
      marqu::Axis y = marqu::Axis::y;
  };
}

#endif
