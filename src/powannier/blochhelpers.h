#ifndef POWANNIER_BLOCHHELPERS_H
#define POWANNIER_BLOCHHELPERS_H

#include <complex>
#include <vector>

#include "aliases.h"
#include "blochsystem.h"
#include "bloch.h"

namespace POWannier {
  class BlochSystem;
  enum class DrawRange {
    WholeLattice,
    ElementaryCell,
    BrillouinZone
  };

  void draw(const BlochSystem& bs, auto&& func, std::string filename, int density, Position beg, Position end) {
    arma::vec range = arma::linspace<arma::vec>(0, 1, density);
    std::ofstream outfile(filename);
    outfile.precision(14);

    std::function<void(Position,int)> nDloop;
    nDloop = [&] (Position alpha, int n) {
      if (n == bs.dim()) {
        Position r = alpha * bs.latticeBasis();
        auto result = func(r);
        arma::rowvec toPrint(bs.dim() + result.n_elem);
        toPrint.head(bs.dim()) = r;
        toPrint.tail(result.n_elem) = result;
        toPrint.raw_print(outfile);
      } else {
        for (const auto& x : range) {
          alpha(n) = beg(n) + x * (end(n) - beg(n));
          nDloop(alpha, n+1);
        }
      }
    };

    Position alpha(bs.dim());
    nDloop(alpha, 0);
  }

  void draw(const BlochSystem& bs, auto&& func, std::string filename, int density, ReciprocalVector beg, ReciprocalVector end) {
    arma::vec range = arma::linspace<arma::vec>(0, 1, density);
    std::ofstream outfile(filename);
    outfile.precision(14);

    std::function<void(ReciprocalVector,int)> nDloop;
    nDloop = [&] (ReciprocalVector k, int n) {
      if (n == bs.dim()) {
        //auto result = arma::conv_to<arma::vec>::from(func(k));
        arma::vec result = func(k);
        arma::vec toPrint(bs.dim() + result.n_elem);
        toPrint.head(bs.dim()) = k;
        toPrint.tail(result.n_elem) = result;
        toPrint.t().raw_print(outfile);
      } else {
        for (const auto& x : range) {
          k(n) = beg(n) + x * (end(n) - beg(n));
          nDloop(k, n+1);
        }
      }
    };

    ReciprocalVector k(bs.dim());
    nDloop(k, 0);
  }

  void draw(const BlochSystem& bs, auto&& func, std::string filename, int density, DrawRange drawRange) {
    Position beg, end;
    if (drawRange == DrawRange::WholeLattice) {
      double ibeg = -bs.N/2.0;
      double iend =  bs.N/2.0;
      if (bs.N %2 == 0) {
        ibeg = -(bs.N - 1) / 2.0;
        iend =  (bs.N + 1) / 2.0;
      }
      beg = ibeg * Position(bs.dim(), arma::fill::ones);
      end = iend * Position(bs.dim(), arma::fill::ones);
      draw(bs, func, filename, density, beg, end);
    } else if (drawRange == DrawRange::ElementaryCell) {
      beg = Position(bs.dim(), arma::fill::zeros);
      end = Position(bs.dim(), arma::fill::ones);
      draw(bs, func, filename, density, beg, end);
    } else {
      throw std::runtime_error("invalid draw range!");
    }
  }

  void drawBands(int numberOfBands, const Potential& V, double kl, int cutoff, double s, std::string fileName, int density, DrawRange drawRange);
}

#endif