#include "utils.h"

#include <stdexcept>

unsigned long long marqu::binomialCoefficient(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;

    // (n choose k) = n * (n-1 choose k-1) / k
    return (n * marqu::binomialCoefficient(n-1, k-1))/k;
} 

int marqu::binarySearch(const std::vector<double> & vec, double val){
  int left = 0;
  int right = vec.size();
  while(left < right){
    int mid = left + (right - left)/2;
    if(vec[mid] < val) left = mid + 1;
    else right = mid;
  }
  return left;
}

std::string marqu::subString(
    const std::string & str, const std::vector<int> & inds){
  std::string res;
  res.reserve(inds.size());
  for(int index : inds){
    res += str[index];
  }
  return res;
}
