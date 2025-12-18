#include <vector>
#include <string>

namespace marqu{
  // (n choose k)
  unsigned long long binomialCoefficient(int n, int k);

  int binarySearch(const std::vector<double> & vec, double val);

  std::string subString(const std::string & str, const std::vector<int>  & inds);
}
