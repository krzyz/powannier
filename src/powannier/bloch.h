#ifndef POWANNIER_BLOCH_H
#define POWANNIER_BLOCH_H
#include "blochspec.h"

namespace POWannier {
  /**
   * @brief %Bloch function calculated for a given reciprocal lattice vector.
   * 
   * @details
   * %Bloch function calculated for a given potential \f$ V(\vec{r}) \f$ and 
   * reciprocal lattice vector \f$ \vec{k} \f$. It has the form: 
   * \f$ \phi_{i,\vec{k}} (\vec{r}) = \frac{1}{\sqrt{\Omega}} e^{i \vec{k} 
   * \vec{r}} u(\vec{r}) \f$ such that:
   *   - \f$ \Omega \f$ is the volume of the system (in the case considered here,
   *     that is finite system with periodic boundaries),
   *   - \f$ n \f$ is the band index
   *   - \f$ u(\vec{r}) \f$ is a periodic function with the same period as 
   *     \f$ V(\vec{r}) \f$ (accordingly, \f$ u(\vec{r}) \f$ can be expanded
   *     in a Fourier series - called also plane wave expansion: \f$ u(\vec{r})
   *     = \sum_n u_{i,n} e^{i (\sum_j n_j \vec{b_j}) \vec{r}} \f$. In practice, the range
   *     of @p n indices has to be finite with some @p cutoff value)
   *   - it is a solution to the time-independent Schrödinger equation:
   *     \f$ \hat{H} \phi_{i,\vec{k}} (\vec{r}) = E_k \phi_{i,\vec{k}}(\vec{r}) \f$,
   *       where \f$ \hat{H} = -\frac{\hbar^2}{2m} \nabla^2 + V(\vec{r}) \f$.
   * 
   * The energies (accessible by function energies() ) are given in units of recoil energy,
   * \f$ E_R = \frac{\hbar^2 k_l^2}{2m} \f$, assuming that the lattice potential 
   * is created with the laser of wavelength \f$ k_l \f$.
   * 
   * @see
   *   BlochSystem class.
   */
  class Bloch {
    public:
      /**
       * Constructor.
       * 
       * @param k A reciprocal lattice vector.
       * @param potential A potential of the same dimension as @p k vector.
       * @param laserWavelength Wavewector of the laser \f$ k_l \f$.
       *  used only to set energy scale of the problem (see the detailed description).
       * @param fourierCutoff Cutoff of plane wave expansion.
       * @param potentialDepth The scale of the potential
       *   (that is the parameter by which the function provided in
       *    construction of @p potential is multiplied).
       */
      Bloch(ReciprocalVector k, const Potential& potential,
        double laserWavelength, int fourierCutoff, double potentialDepth = 1);

      /**
       * @returns
       *   the reciprocal vector \f$ \vec{k} \f$ for which the %Bloch function is calculated.
       */
      ReciprocalVector k() const {
        return _k;
      }

      /**
       * Energies of the bloch functions, given in units of \f$ E_R = \frac{\hbar^2 k_l^2}{2m} \f$.
       * @returns
       *   Vector containing the energies, where \f$ i \f$th position of the vector is the energy for the \f$ i \f$th band.
       */
      arma::vec energies() const {
        return _energies;
      }

      /**
       * List of sets of indices for which the Fourier series coefficients are calculated.
       */
      std::vector<POWannier::NPoint> ns() const {
        return _blochSpec.ns;
      };

      /**
       * @param r
       *   Position in the lattice \f$ \vec{r} \f$.
       * @param band
       *   The band index \f$ i \f$.
       * 
       * @returns
       *   The value of the %Bloch function \f$ \phi_{n,\vec{k}}(\vec{r}) \f$.
       */
      std::complex<double> operator()(Position r, int band);

      /**
       * @returns
       *   The lattice basis for a given problem.
       */
      const LatticeBasis& latticeBasis() const {
        return _blochSpec.latticeBasis;
      }

    private:
      const BlochSpec _blochSpec;
      const ReciprocalVector _k;
      arma::vec _energies;
      arma::cx_mat _eigenvectors;
  };
}

#endif