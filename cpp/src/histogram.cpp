#include "Histogram.h"
#include <complex>

template <typename T>
marqu::Histogram<T>::Histogram(double min, double max, std::size_t nBins) : 
  //histogram(nBins, T{}), counts(nBins, 0), nBins(nBins), min(min), max(max) {
  histogram(nBins, T{}), nBins(nBins), min(min), max(max) {
    binWidth = (max - min)/nBins;
    indexRescale = 1/binWidth;
    for(std::size_t i = 0; i < nBins; i++){
      binBounds.emplace_back(min + binWidth*i, min + binWidth*(i + 1));
    }
  }

template <typename T>
std::size_t marqu::Histogram<T>::index(double val) const{
  if(val == max) return nBins-1;
  return std::size_t(indexRescale*(val - min));
}

template <typename T>
void marqu::Histogram<T>::add(double windowMin, double windowMax, T weight){
  if(windowMin >= binBounds[lastIndex].first && 
      windowMax <= binBounds[lastIndex].second){
    histogram[lastIndex] += weight * (windowMax - windowMin);
    return;
  }

  std::size_t indexMin = index(windowMin);
  std::size_t indexMax = index(windowMax);

  if(indexMin == indexMax){
    histogram[indexMin] += weight * (windowMax - windowMin);
    lastIndex = indexMin;
    return;
  }

  for(std::size_t i = indexMin; i <= indexMax; i++){
    double overlap = std::min(windowMax, binBounds[i].second) -
      std::max(windowMin, binBounds[i].first);
    histogram[i] += weight * overlap;
  }
}

template <typename T>
std::vector<T> marqu::Histogram<T>::averagedWeights() const{
  std::vector<T> res(histogram);
  for(std::size_t i = 0; i < nBins; i++){
    res[i] *= indexRescale;
  }
  return res;
}

template <typename T>
void marqu::Histogram<T>::clear(){
  std::fill(histogram.begin(), histogram.end(), T{});
  lastIndex = 0;
}

template class marqu::Histogram<std::size_t>;
template class marqu::Histogram<int>;
template class marqu::Histogram<double>;
template class marqu::Histogram<std::complex<double>>;
