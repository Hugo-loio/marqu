#include "utils.h"

unsigned long long marqu::binomialCoefficient(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;

    // (n choose k) = n * (n-1 choose k-1) / k
    return (n * marqu::binomialCoefficient(n-1, k-1))/k;
} 
