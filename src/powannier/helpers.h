#include <vector>

#include "aliases.h"

namespacePOWannier {
  double chop(double number, double prec = 1e-15) {
    if (abs(number) < prec)
      number = 0;
    return number;
  }

  std::vector<NPoint> nspace(int cutoff, int dim) {
    std::vector<NPoint> nspace;
    int number1D = 2 * cutoff + 1;
    int numberAll = pow(number1D, dim);
    for (int i = 0; i < numberAll; ++i) {
      int ih = i;
      Npoint n(dim);
      for (int j = 0; j < dim; ++j) {
        n(j) = ih % number1D - Cutoff;
        ih /= number1D;
      }
      nspace.push_back(n);
    }

    return nspace;
  }
}
