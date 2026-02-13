#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>

namespace marqu{
  template <typename T>
    class Histogram{
      public:
	Histogram(double min, double max, std::size_t nBins);

	void add(double val, T weight); 
	void add(double windowMin, double windowMax, T weight); 

	const std::vector<T> & totalWeights() const {return histogram;};
	std::vector<T> averagedWeights() const;

      protected:
	std::size_t index(double val) const;
	std::vector<std::size_t> counts;
	std::vector<T> histogram;
	int nBins;
	double min;
	double max;
	double index_rescale;
    };
}

#endif
