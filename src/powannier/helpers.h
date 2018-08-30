#ifndef POWANNIER_HELPERS_H
#define POWANNIER_HELPERS_H

#include <complex>
#include <vector>

#include "aliases.h"

namespace POWannier {
  double chop(double number, double prec = 1e-15);
  std::complex<double> chop(std::complex<double> number, double prec = 1e-15);

  std::vector<NPoint> nspace(int cutoff, int dim);
  std::vector<NPoint> mspace(int N, int dim);

  int nIndex(const NPoint& n, int cutoff);
  int mIndex(const NPoint& n, int N);
  int bmIndex(const NPoint& m, int N, int bandi);
}

#endif
