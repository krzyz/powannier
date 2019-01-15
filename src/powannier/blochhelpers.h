#ifndef POWANNIER_BLOCHHELPERS_H
#define POWANNIER_BLOCHHELPERS_H

#include <complex>
#include <vector>

#include "aliases.h"
#include "blochsystem.h"
#include "bloch.h"

namespace POWannier {
  class BlochSystem;

  /**
   * @enum DrawRange
   * Enum that describes the range in which the function values 
   * are to be calculated.
   */
  enum class DrawRange {
    /// The whole lattice for which the problem is defined (real space)
    WholeLattice,
    /// Elementary cell (spanned by lattice basis vectors \f$\vec{a_i}\f$, real space)
    ElementaryCell,
    /**
     * not an actual Brillouin zone - just the space spanned by reciprocal 
     * lattice basis vectors \f$ \vec{b_i} \f$ (reciprocal space)
     */
    BrillouinZone
  };

  /**
   * @brief Calculate and write to a file the values of a given function 
   * in the provided region in real space.
   * 
   * @details
   *   The function calculates values in the given range, which is given
   *   by two vectors (@p beg and @p end) in the lattice basis vectors 
   *   (\f$ \vec{a_i} \f$) space. The range is simply \f$ n\f$-dimensional cuboid
   *   (where \f$ n \f$ is the dimension of the problem) in that space. In real space,
   *   the outermost edges are: \f$ r_{\textrm{beg}} = \sum_i \textrm{beg}_i 
   *   \vec{a_i} \f$ and \f$ r_{\textrm{end}} = \sum_i \textrm{end}_i \vec{a_i} \f$
   *   and the range is \f$ n\f$-dimensional parallelepiped. \n
   *   The values of the function are calculated for points that are uniformly
   *   distributed in the range and there is \f$\textrm{density}^n\f$) of them in total
   *   (that is, \f$ \textrm{density} \f$ in each direction).
   * 
   * @tparam Function
   *   Callable type.
   * 
   * @param latticeBasis
   *   The lattice basis vectors of the system.
   * @param func
   *   Function that is calculated at each considered point.
   * @param filename
   *   Name of the output file.
   * @param density
   *   The number of points the value of the function is calculated for 
   *   in each direction.
   * @param beg
   *   The outermost edge of the range (given in the lattice basis vectors space).
   * @param end
   *   The outermost edge of the range (opposite to @p beg, also given in the 
   *   lattice basis vectors space).
   * 
   * @pre
   *   @p func must take the argument of type Position and return real number.
   */
  template <class Function>
  void draw(const LatticeBasis& latticeBasis, Function&& func, std::string filename, int density, Position beg, Position end) {
    arma::vec range = arma::linspace<arma::vec>(0, 1, density);
    std::ofstream outfile(filename);
    outfile.precision(14);
    int dim = latticeBasis.n_rows;

    std::function<void(Position,int)> nDloop;
    nDloop = [&] (Position alpha, int n) {
      if (n == dim) {
        Position r = alpha * latticeBasis;
        auto result = func(r);
        arma::rowvec toPrint(dim + result.n_elem);
        toPrint.head(dim) = r;
        toPrint.tail(result.n_elem) = result;
        toPrint.raw_print(outfile);
      } else {
        for (const auto& x : range) {
          alpha(n) = beg(n) + x * (end(n) - beg(n));
          nDloop(alpha, n+1);
        }
      }
    };

    Position alpha(dim);
    nDloop(alpha, 0);
  }

  /**
   * @brief Calculate and write to a file the values of a given function 
   * in the provided region in reciprocal space.
   * 
   * @details
   *   The function calculates values in the given range, which is given
   *   by two vectors (@p beg and @p end) in the reciprocal space. The range
   *   is the \f$ n\f$-dimensional cuboid in that space.
   *   The values of the function are calculated for points that are uniformly
   *   distributed in the range and there is \f$\textrm{density}^n\f$) of them 
   *   in total (that is, \f$ \textrm{density} \f$ in each direction).
   * 
   * @tparam Function
   *   Callable type.
   * @param latticeBasis
   *   The lattice basis vectors of the system.
   * @param func
   *   Function that is calculated at each considered point.
   * @param filename
   *   Name of the output file.
   * @param density
   *   The number of points the value of the function is calculated for 
   *   in each direction.
   * @param beg
   *   The outermost edge of the range.
   * @param end
   *   The outermost edge of the range (opposite to @p beg).
   * @pre
   *   @p func must take the argument of type ReciprocalVector and return real number.
   */

  template <class Function>
  void draw(const LatticeBasis& latticeBasis, Function&& func, std::string filename, int density, ReciprocalVector beg, ReciprocalVector end) {
    arma::vec range = arma::linspace<arma::vec>(0, 1, density);
    std::ofstream outfile(filename);
    outfile.precision(14);
    int dim = latticeBasis.n_rows;

    std::function<void(ReciprocalVector,int)> nDloop;
    nDloop = [&] (ReciprocalVector k, int n) {
      if (n == dim) {
        arma::vec result = func(k);
        arma::vec toPrint(dim + result.n_elem);
        toPrint.head(dim) = k;
        toPrint.tail(result.n_elem) = result;
        toPrint.t().raw_print(outfile);
      } else {
        for (const auto& x : range) {
          k(n) = beg(n) + x * (end(n) - beg(n));
          nDloop(k, n+1);
        }
      }
    };

    ReciprocalVector k(dim);
    nDloop(k, 0);
  }

  /**
   * @brief Calculate and write to a file the values of a given function 
   * over one of the predefined (by enum DrawRange) regions.
   * 
   * @param bs
   *   The BlochSystem type object for the problem.
   * @param func
   *   Function that is calculated at each considered point.
   * @param filename
   *   Name of the output file.
   * @param density
   *   The number of points the value of the function is calculated for 
   *   in each direction.
   * @param drawRange
   *   One of the draw ranges defined in the @ref DrawRange enum class.
   * 
   * @pre
   *   @p func must take the argument of type appropriate for the given 
   *   @p drawRange and return real number.
   */
  template <class Function>
  void draw(const BlochSystem& bs, Function&& func, std::string filename, int density, DrawRange drawRange) {
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
      draw(bs.latticeBasis(), func, filename, density, beg, end);
    } else if (drawRange == DrawRange::ElementaryCell) {
      beg = Position(bs.dim(), arma::fill::zeros);
      end = Position(bs.dim(), arma::fill::ones);
      draw(bs.latticeBasis(), func, filename, density, beg, end);
    } else {
      throw std::runtime_error("invalid draw range!");
    }
  }

  /**
   * @brief Calculate and write to a file the energy spectrum of the 
   * band structure over one of the predefined (by enum DrawRange) regions.
   * 
   * @param numberOfBands
   *   Number of bands for which to calculate the band structure (starting with
   *   the lowest band).
   * @param V
   *   Given Potential. 
   * @param kl
   *   Laser wavelength (see Potential class).
   * @param cutoff
   *   Cutoff in Fourier series approximation.
   * @param depth
   *   Potential depth.
   * @param fileName
   *   Name of the output file.
   * @param density
   *   The number of points the value of the function is calculated for 
   *   in each direction.
   * @param drawRange
   *   One of the draw ranges defined in the @ref DrawRange enum class.
   */

  void drawBands(int numberOfBands, const Potential& V, double kl, int cutoff, double depth, std::string fileName, int density, DrawRange drawRange);
}

#endif