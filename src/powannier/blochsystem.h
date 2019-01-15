#ifndef POWANNIER_BLOCHSYSTEM_H
#define POWANNIER_BLOCHSYSTEM_H

#include "blochspec.h"

namespace POWannier { 
  /**
   * @brief %Bloch functions calculated for a given lattice system.
   * 
   * BlochSystem class stores all the %Bloch functions (or rather, 
   * the coefficients of their plane wave expansions) that are to be used
   * in Wannier function calculations for a given system. \n 
   * 
   * The number of the stored functions depends on the system size @p N, 
   * as the relevant quasimomentum wavewectors take the form 
   * \f$ \vec{k} = \frac{\sum_i m_i \vec{b}_i}{N} \f$, where \f$ \vec{b}_i \f$
   * are the reciprocal lattice basis vectors and \f$ \{ m_1, m_2,
   * \dots, m_{dim} \}\f$, are sets of indices such that \f$ m_i \in 
   * \{-(N-1)/2, \dots, (N+1)/2\} \f$ (in case of odd @p N) and \f$ m_i \in
   * \{-N/2+1, \dots, N/2\} \f$ (in case of even @p N). The indices actually
   * used in this library (which are returned by mspace() helper function 
   * and stored in @p ms member) are actually \f$ m_i \in \{0, \dots, N \}\f$,
   * as only their differences (see RSystem) are used in the calculations.
   * The functions are approximated as per @p cutoff parameter. \n
   * 
   * Construction of the BlochSystem object is done through the use of a builder:
   * @code
   * auto bs = POWannier::BlochSystem::createWithPotential(V)
   *  .setLaserWavelength(pi)
   *  .setFourierCutoff(5)
   *  .setSystemSize(20)
   *  .setPotentialDepth(8)
   *  .complete();
   * @endcode
   * 
   * The call to createWithPotential() initializes the building process,
   * while the call to complete() returns the final object. All the @p set
   * calls are optional and if not used, the default values are used for 
   * corresponding members.
   * 
   * @see
   *   Bloch 
   *   BlochSystemBuilder
   */
  class BlochSystem {
    public:
      /**
       * @brief Builder for BlochSystem class.
       * 
       * Usage is described in BlochSystem class.
       */
      class BlochSystemBuilder {
        public:
          /**
           * @brief Constructor.
           * 
           * @param potential
           *   Potential for a given system.
           */
          BlochSystemBuilder(const Potential& potential) :
            _potential(potential) {}

          /**
           * @brief Basic setter.
           * The default value is \f$ \pi \f$.
           * @param kl
           *   Laser wavelength \f$ k_l \f$ (see Bloch)
           * @returns
           *   This instantce of BlochSystemBuilder with modified laser wavelength.
           */
          BlochSystemBuilder& setLaserWavelength(double kl) {
            _laserWavelength = kl;

            return *this;
          }

          /**
           * @brief Basic setter.
           * The default value is 2.
           * @param cutoff
           *   Cutoff in plane wave expansion of the %Bloch function (see Bloch).
           * @returns
           *   This instance of BlochSystemBuilder with modified cutoff.
           */
          BlochSystemBuilder& setFourierCutoff(int cutoff) {
            _fourierCutoff = cutoff;

            return *this;
          }

          /**
           * @brief Basic setter.
           * The default value is 5.
           * @param N
           *   System size (the number of lattice cells in each direction).
           * @returns
           *   This instance of BlochSystemBuilder with modified system size.
           */
          BlochSystemBuilder& setSystemSize(int N) {
            _systemSize = N;

            return *this;
          }

          /**
           * @brief Basic setter.
           * The default value is 8.
           * 
           * @param potentialDepth The scale of the potential
           *   (that is the parameter by which the function provided in
           *   construction of @p potential is multiplied).
           * @returns
           *   This instance of BlochSystemBuilder with modified potential depth.
           */
          BlochSystemBuilder& setPotentialDepth(int potentialDepth) {
            _potentialDepth = potentialDepth;

            return *this;
          }

          /**
           * @brief Finish building the BlochSystem and return the result.
           * @returns
           *   Finished object of BlochSystem class.
           */
          BlochSystem complete() {
            return BlochSystem(_potential, _laserWavelength, _fourierCutoff, _systemSize, _potentialDepth);
          }


        private:
          const Potential _potential;
          double _laserWavelength = pi;
          int _fourierCutoff = 2;
          int _systemSize = 5;
          int _potentialDepth = 8;
      };

      /// Number of the cells in the lattice in each direction.
      const int N;
      /// Indices corresponding to quasimomentum wavevectors.
      const std::vector<POWannier::NPoint> ms;

      /**
       * @brief Initialization of the builder
       * @param potential
       *   Potential for which the bloch functions are to be calculated.
       * @returns
       *   Builder, using which further characterization can be made.
       */
      static BlochSystemBuilder createWithPotential(const Potential& potential);

      /**
       * @brief Energies of the %Bloch functions (given in \f$E_R\f$ units, see Bloch)
       * @param m
       *   Index corresponding to some quasimomentum \f$\vec{k}\f$.
       * @returns
       *   List of energies of %Bloch functions \f$\phi_{i,\vec{k}}(\vec{r})\f$.
       */
      const arma::vec& energies(NPoint m) const;

      /**
       * @brief Coefficients of the bloch functions plane wave expansion.
       * @param m
       *   Index corresponding to some quasimomentum \f$\vec{k}\f$.
       * @returns
       *   A matrix containing the coefficients \f$ u_{i,n}\f$ (see Bloch) of
       *   the %Bloch functions \f$\phi_{i,\vec{k}} (\vec{r})\f$.
       *   Different columns correspond to different band indices @p i,
       *   while the rows correspond correspond to indices @p ms in order.
       */
      const arma::cx_mat& eigenvectors(NPoint m) const;

      /**
       * @brief Calculate %Bloch function for a given @p m indices.
       * 
       * @param m
       *   Index corresponding to some quasimomentum \f$\vec{k}\f$.
       * @param r
       *   %Position.
       * @param band
       *   Band index @p i.
       * @returns
       *   value of the Bloch function \f$ \phi_{i,\vec{k}} (\vec{r}) \f$.
       */
      std::complex<double> bloch(NPoint m, Position r, int band = 0) const;

      /**
       * @brief Single coefficient of the bloch function plane wave expansion.
       * 
       * @param m
       *   Index corresponding to some quasimomentum \f$\vec{k}\f$.
       * @param n 
       *   Index of the plane wave (Fourier series coefficient), see Bloch.
       * @param band
       *   Band index @p i.
       * @returns
       *   Value of the \f$ u_{i,n} \f$ for quasimomentum \f$ \vec{k} \f$.
       */
      std::complex<double> blochC(NPoint m, NPoint n, int band = 0) const;

      /**
       * @brief Quasimomentum corresponding to a given index.
       * 
       * @param m
       *   Index corresponding to some quasimomentum \f$\vec{k}\f$.
       * @returns
       *   The quasimomentum \f$ \vec{k} \f$.
       */
      ReciprocalVector kFromM(NPoint m) const;

      /**
       * @brief Simple getter.
       * @return Laser wavelength \f$ k_l \f$.
       */
      double kl() const {
        return _blochSpec.kl;
      }

       /**
       * @brief Simple getter.
       * @return Cutoff of the plane wave expansion.
       */
      int cutoff() const {
        return _blochSpec.cutoff;
      }

      /**
       * @brief Simple getter.
       * @return The dimension of the system.
       */
      int dim() const {
        return _blochSpec.dim;
      }

      /**
       * @brief Simple getter.
       * @return Volume of the elementary cell.
       */
      double elementaryCellVolume() const {
        return _elementaryCellVolume;
      }

      /**
       * @brief Simple getter.
       * @return The lattice basis vectors.
       */
      const LatticeBasis& latticeBasis() const {
        return _blochSpec.latticeBasis;
      }

    private:
      const double _elementaryCellVolume;
      const BlochSpec _blochSpec;
      std::vector<arma::vec> _energies;
      std::vector<arma::cx_mat> _eigenvectors;

      void generateAll(const Potential& potential);
      BlochSystem(const Potential& potential, 
        double kl, int cutoff, int N, double s = 1);

    friend class BlochSystemBuilder;
  };
}

#endif
