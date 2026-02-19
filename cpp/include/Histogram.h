#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>

//TODO: possible future optimization, same bin counts for multiple observables

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
	std::size_t nBins;
	double min;
	double max;
	double binWidth;
	double indexRescale;

	std::pair<double, double> lastBin = {0,0}; 
	std::size_t lastIndex = 0;
	void updateLast(std::size_t index);
    };
}

#endif
