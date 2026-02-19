#include "Histogram.h"
#include <complex>


template <typename T>
marqu::Histogram<T>::Histogram(double min, double max, std::size_t nBins) : 
  histogram(nBins, T{}), counts(nBins, 0), nBins(nBins), min(min), max(max) {
    binWidth = (max - min)/nBins;
    indexRescale = 1/binWidth;
  }

template <typename T>
std::size_t marqu::Histogram<T>::index(double val) const{
  if(val == max) return nBins-1;
  return std::size_t(indexRescale*(val - min));
}

template <typename T>
void marqu::Histogram<T>::updateLast(std::size_t index){
  lastIndex = index;
  lastBin.first = min + binWidth*index;
  lastBin.second = min + binWidth*(index + 1);
}

template <typename T>
void marqu::Histogram<T>::add(double val, T weight){
  if(val >= lastBin.first && val <= lastBin.second){
    histogram[lastIndex] += weight;
    counts[lastIndex]++;
    return;
  }
  std::size_t i = index(val);
  histogram[i] += weight;
  counts[i]++;
  updateLast(i);
}

template <typename T>
void marqu::Histogram<T>::add(double windowMin, double windowMax, T weight){
  if(windowMin >= lastBin.first && windowMax <= lastBin.second){
    histogram[lastIndex] += weight;
    counts[lastIndex]++;
    return;
  }

  std::size_t indexMin = index(windowMin);
  std::size_t indexMax = index(windowMax);

  if(indexMin == indexMax){
    histogram[indexMin] += weight;
    counts[indexMin]++;
    updateLast(indexMin);
    return;
  }

  for(std::size_t i = indexMin; i <= indexMax; i++){
    histogram[i] += weight;
    counts[i]++;
  }
}

template <typename T>
std::vector<T> marqu::Histogram<T>::averagedWeights() const{
  std::vector<T> res(histogram);
  for(std::size_t i = 0; i < nBins; i++){
    res[i] /= counts[i];
  }
  return res;
}

template class marqu::Histogram<std::size_t>;
template class marqu::Histogram<int>;
template class marqu::Histogram<double>;
template class marqu::Histogram<std::complex<double>>;
