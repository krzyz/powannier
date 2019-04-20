#ifndef POWANNIER_RSYSTEM_H
#define POWANNIER_RSYSTEM_H

#include <armadillo>
#include <memory>
#include "blochsystem.h"
#include "wannier.h"
#include "wannierpositions.h"
#include "../util/cache.h"

namespace POWannier {
  /**
   * @brief Object responsible for position operators diagonalization 
   * and subsequent Wannier function creation.
   * 
   * @par Introduction
   * Most lattice models have an effective (small) range of interactions
   * between the particles present therein. This assumption naturally
   * leads to the <a href="https://en.wikipedia.org/wiki/Tight_binding">
   * tight binding model </a> in which the particles are represented 
   * in the basis of local (exponentially decaying) functions. An example
   * of such basis are the %Wannier functions, which for simple systems
   * can be calculated directly from the %Bloch functions (see Bloch): \n
   * 
   * \f$ W_{\alpha,\vec{R}} (\vec{r}) = \frac{1}{\sqrt{N^{dim}}} 
   * \sum_{\vec{k}} e^{-i\vec{k}\vec{R}}\phi_{i,\vec{k}}(\vec{r}). \f$ \n
   * 
   * In the case of more complicated lattice potentials the above formula
   * does not work (that is, the resulting basis functions are not local).
   * Here a simple approach to generating Wannier functions is taken,
   * namely the diagonalization of the position operators in the %Bloch
   * function bases @cite kivelson82, which has already been succesfully 
   * applied to optical lattice systems @cite bissbort2012. 
   * 
   * @par Position operators diagonalization method details
   * The position operators here are defined as: \f$ \hat{\lambda}_i = 
   * \frac{\vec{b}_i}{2\pi} \f$. In order to diagonalize these operators,
   * a Hamiltonian matrix may be easily constructed: (here, \f$ p \f$ and
   * \f$ q \f$ are the indices combining \f$ \vec{k} \f$ (in the form of @p m)
   * and the band index): \n
   * \f$ (\Lambda_i)_{p,q} = \langle \phi_{\alpha, \vec{k}} | \hat{\lambda_i} | 
   * \phi_{\beta, \vec{k}'} \rangle = \frac{1}{\Omega} \int_{\Omega}
   * d\vec{r} \lambda_i \sum_{m, n, m', n'} u^*_{\alpha,n} u_{\alpha,n'}
   * e^{i \sum_j (n_j' - n_j + (m_j' - m_j)/N) \vec{b}_j \vec{r} } \f$ \n
   * 
   * In order to evaluate this integral, the following change of variables
   * is used: \f$ \vec{r}_i = N \vec{a_i} \lambda_i + c \f$, where the constant
   * \f$ c \f$ depends on whether @p N is odd or even: \f$ c_{odd} = 0 \f$,
   * \f$ c_{even} = \vec{a_i}/2 \f$. This change has an associated Jacobian
   * determinant equal to \f$ \Omega \f$ and, assuming that \f$\vec{r} = \vec{0}\f$
   * is in the middle of the \f$ (N+1)/2 \f$-th cell in the odd @p N case (and 
   * \f$ (N/2) \f$-th cell in th even @p N case), changes the boundaries of 
   * integration (now in dimensionless \f$ d\lambda_i \f$) to {-1/2, 1/2}.
   * After additional change of order between integration and summations and
   * noticing that the integrals in each (\f$ \lambda_i \f$) direction can 
   * now be sepereated, we finally get (in the even @p N case): \n
   * \f$ (\Lambda_i)_{p,q} = \sum_{m, n, m', n'} I_{i, m_j, n_j, m_j', n_j'} 
   * \prod_{j\neq i} I_{m_j, n_j, m_j', n_j'} \f$,
   * where: \n
   * \f$ I_{j, m_j, n_j, m_j', n_j'} = \int_{-1/2}^{1/2} d\lambda_j 
   * \lambda_j e^{2\pi i s_j \lambda_j} \f$, \n
   * \f$ I_{m_j, n_j, m_j', n_j'} = \int_{-1/2}^{1/2} d\lambda_j 
   * e^{2\pi i s_j \lambda_j} \f$, \n
   * and \f$ s_i = (N(n_i' - n_i) + m_i' - m_i) \f$. \n \n
   * These integrals can be done analytically, resulting in: \n
   * \f$ I_{j, m_j, n_j, m_j', n_j'} = \sinc(\pi s_j) = 
   * \delta_{s_j, 0} \f$, \n
   * \f$ I_{m_j, n_j, m_j', n_j'} = \frac{i (-s_j \pi \cos(\pi s_j)) + \sin(\pi 
   * s_j)}{2 \pi^2 s_j^2} = \frac{i}{2\pi s_j}(-1)^{s_j+1} (1-\delta_{s_j,0}) \f$, \n
   * where the second equations in both lines follow from the fact that
   * \f$ s_j \in \mathbb{Z} \f$. We also have \f$s_j = 0 \iff m_j = m_j' \wedge
   * n_j = n_j' \f$, because \f$ 0 \leq \max(m_j' - m_j) - \min(m_i' - m_i) < N \f$.
   * Because of this fact and the delta functions in the result, each \f$ \Lambda_i \f$
   * can be represented in block diagonal form, where for each block the following
   * is true: \f$ \forall_{j \neq i} m_j = m_j' \wedge n_j = n_j' \f$, which can
   * significantly speed up the diagonalization process. \n \n 
   * The final procedure of finding the %Wannier function centered at position
   * \f$ \vec{R} = \{R_1, R_2, \dots, R_{dim}\} \f$ consists of the following loop:
   * <ul>
   * <li>diagonalize \f$\Lambda_i\f$ in the subspace of eigenfunctions, 
   *   whose eigenvalues correspond to positions \f$ R_1, R_2, \dots, R_{i-1} 
   *   \f$ </li>
   * <li>select the eigenfunctions of \f$ \Lambda_i \f$ that form the 
   *   degenerated subspace with the eigenvalue corresponding to the position 
   *   \f$ R_i \f$, as that subspace will be used in the next iteration of 
   *   the loop.</li></ul>
   * Before the first (@p i = 1) step, the current subspace is the whole space 
   * spanned by all  \f$ \phi_{\alpha,\vec{k}} \f$ %Bloch functions. After the 
   * final (@p i = @p dim) step, the result will be a localized function 
   * centered at position  \f$\vec{R}\f$, that is \f$W_{\vec{R}} (\vec{r})\f$.
   * Note that here the %Wannier function does not have a band index -- in 
   * general, this function may be spanning multiple bands (for example, this
   * happens when there are multiple potential wells of similar depth in a
   * single elementary cell).
   * }
   * 
   * @par %Wannier function selection process
   * The selection of \f$\vec{R}\f$ in this class is done in two steps: first,
   * the cell is selected (using NPoint @p cellLocation parameter -- 
   * \f$ \vec{0} \f$ correspond to the middle of the lattice) and after that,
   * using NPoint @p wanLocation an appropiate %Wannier function can be chosen.
   * NPoint @p wanLocation is set up to cover all cases and the process can
   * be simplified in some of them:
   * - If there is only one %Wannier function in the cell for the energy 
   *   given by the band energy, @p wanLocation may be omitted
   * - If the eigenvalues of \f$\Lambda_1\f$ is nondegenerate, it is enough
   *   to provide only single index as @p wanLocation, @p {n}, and by doing so
   *   %Wannier corresponding the <tt>n</tt>-th \f$ \Lambda_1 \f$ eigenvalue
   *   will be given.
   * - If some of eigenvalues of \f$\Lambda_1\f$ are degenerate, it is
   *   necessary to provide information on this degeneracy while constructing
   *   the RSystem (using WannierLocation object) and then each index 
   *   \f$ \textrm{wanLocation}_i \f$ specifies which of the eigenvalues is to 
   *   be chosen at <tt>i</tt> step of the diagonalization.
   *   loop.
   * 
   * @par Example
   * @code
   *   // Trivial case: one-dimensional sin^2(k_l x) potential
   *   double kl = pi;
   * 
   *   auto basis = std::vector<POWannier::Vector> { {pi/kl} };
   * 
   *   auto sineSquare = [] (Position r) {
   *     // Only r[0] is needed since the potential is in 1d
   *     return std::pow(std::sin(kl * r[0]), 2);
   *   };
   * 
   *   auto V = POWannier::Potential::create()
   *     .setBasis(basis)
   *     .evaluateFromFunction(func)
   *     .setCutoff(1)  // enough to capture this potential exactly
   *     .complete();
   *
   *   auto bs = std::make_shared<POWannier::BlochSystem>(
   *     POWannier::BlochSystem::createWithPotential(V)
   *       .setLaserWavelength(kl)
   *       .setFourierCutoff(20)
   *       .setSystemSize(8)
   *       .setPotentialDepth(8)
   *       .complete()
   *   );
   *
   *   POWannier::RSystem rs(bs, {0});
   * 
   * 
   *   POWannier::NPoint n({0}); // Wannier in the middle of the lattice
   *   auto wan = rs.getWannier(n);
   *   wan.draw("wan1d.dat", 100, POWannier::DrawRange::WholeLattice);
   *
   * @endcode
   * 
   * @see Bloch
   * @see BlochSystem
   * @see Wannier
   * @see WannierPositions
   */
  class RSystem {
    public:
      /// System size (number of cells in each direction).
      const int N;
      /// Dimension of the system
      const int dim;
      /// Dimension of the quasimomentum basis.
      const int mdim;
      /// Dimension of the quasimomentum+band indices basis.
      const int bmdim;
      /// Bands used in the calculations.
      const std::vector<int> bands;
      /// Cutoff in plane wave expansion of %Bloch and %Wannier functions.
      const int cutoff;

      /**
       * @brief Basic constructor.
       * 
       * @param bs
       *   BlochSystem containing %Bloch functions of the system.
       * @param bands
       *   Specify from which bands %Bloch functions are to be used.
       */
      RSystem(std::shared_ptr<BlochSystem> bs,
        std::vector<int> bands);

      /**
       * @brief Cnstructor for the case of degenerated \f$ \hat{\lambda_1}
       *  \f$ eigenvalues
       * 
       * @param bs
       *   std::shared_ptr to BlochSystem containing %Bloch functions of the system.
       * @param bands
       *   Specify from which bands %Bloch functions are to be used.
       * @param positions
       *   Characterization of the eigenvalues degeneracies for a single 
       *   elementary cell (see WannierPositions).
       * @pre
       *   @p positions must contain total positions number (leafs number)
       *   equal to the number of bands supplied (as @p bands).
       * @throws
       *   std::runtime_error if the condition on positions number is not met.
       * @see WannierPositions
       */
 
      RSystem(std::shared_ptr<BlochSystem> bs,
        std::vector<int> bands,
        WannierPositions positions);

      /**
       * @brief Construct a %Wannier function at a given call (single band only)
       * 
       * @param cellLocation
       *   Location of the cell containing selected %Wannier function. 
       *   \f$ \vec{0} \f$ is the middle of the lattice.
       * @returns
       *   %Wannier function \f$ W_{\vec{R}}(\vec{r}) \f$, where \f$ \vec{R} \f$
       *   corresponds to a given $@ cellLocation.
       */
      std::unique_ptr<POWannier> getWannier(NPoint cellLocation);

      /**
       * @brief Construct a %Wannier function at a given call (general case)
       * 
       * @param cellLocation
       *   Location of the cell containing selected %Wannier function. 
       *   \f$ \vec{0} \f$ is the middle of the lattice.
       * @param wanLocation
       *   Location of the %Wannier function inside of the cell chosen by
       *   @p cellLocation. The \f$ \textrm{wanLocation}_i \f$ index corresponds
       *   to which eigenvalue of \f$\hat{\lambda}_i\f$ should be chosen.
       * @returns
       *   %Wannier function \f$ W_{\vec{R}}(\vec{r}) \f$, where \f$ \vec{R} \f$
       *   corresponds to a given @p cellLocation.
       */
      std::unique_ptr<POWannier> getWannier(NPoint cellLocation, NPoint wanLocation);

      /**
       * @brief Get the eigensystem of \f$\hat{\lambda_n}\f$.
       * @param n
       *   The chosen dimension.
       * @returns
       *   The Eigensystem of \f$\hat{\lambda_n}\f$ operator in the whole space.
       */
      const Eigensystem& getEigensystem(int n);

      /**
       * @brief Get the eigensystem \f$\hat{\lambda_n}\f$ in <tt>n</tt>-th subspace.
       * @param n
       *   The chosen dimension.
       * @param cellLocation
       *   See getWannier(NPoint).
       * @param wanLocation
       *   See getWannier(NPoint,NPoint).
       * @returns
       *   The Eigensystem of \f$\hat{\lambda_n}\f$ operator in the subspace 
       *   used in the <tt>i</tt>-th step of the diagonalization loop, while 
       *   using @p cellLocation and @p wanLocation.
       * @note
       *   Only the first @p n of indices in both @p cellLocation and
       *   @p wanLocation are used.
       */
      const Eigensystem& getSubEigensystem(int n, NPoint cellLocation, NPoint wanLocation);

      /**
       * @brief Get the eigenvalues of \f$\hat{\lambda_n}\f$.
       * @see getEigensystem()
       * @returns
       *   The eigenvalues of the \f$ \hat{\lambda}_n \f$ operator.
       */
      const arma::vec& getEigenvalues(int n) {
        return std::get<0>(getEigensystem(n));
      }

      /**
       * @brief Get the eigenvectors of \f$\hat{\lambda_n}\f$.
       * @see getEigensystem()
       * @returns
       *   The eigenvectors of the \f$ \hat{\lambda}_n \f$ operator.
       */
      const arma::cx_mat& getEigenvectors(int n) {
        return std::get<1>(getEigensystem(n));
      }

      /**
       * @brief Get the eigenvalues \f$\hat{\lambda_n}\f$ in <tt>n</tt>-th subspace.
       * @see getSubEigensystem()
       * @returns
       *   The eigenvalues of the <tt>n</tt>-th subspace.
       */
      const arma::vec& getSubEigenvalues(int n, NPoint cellLocation, NPoint wanLocation) {
        return std::get<0>(getSubEigensystem(n, cellLocation, wanLocation));
      }

      /**
       * @brief Get the eigenvectors \f$\hat{\lambda_n}\f$ in <tt>n</tt>-th subspace.
       * @see getSubEigensystem()
       * @returns
       *   The eigenvectors of the <tt>n</tt>-th subspace.
       */
      const arma::cx_mat& getSubEigenvectors(int n, NPoint cellLocation, NPoint wanLocation) {
        return std::get<1>(getSubEigensystem(n, cellLocation, wanLocation));
      }

    private:
      // here std::vector<int> is used instead of NPoint (i.e. arma::ivec),
      // because armadillo vectors do not play well with being map keys as 
      // a part of tuple (operator< has to return bool but when comparing
      // arma vectors it rather returns a vector of bools)
      using Ipp = std::tuple<int,std::vector<int>,std::vector<int>>;

      std::shared_ptr<BlochSystem> _bs;
      WannierPositions _wannierPositions;

      // caches for storing already computed eigensystems
      Cache<int,Eigensystem> _eigensystemCache;
      Cache<Ipp,Eigensystem> _subEigensystemCache;
      Cache<Ipp,arma::cx_mat> _subspaceCache;

      arma::uvec transformToBm() const;
      arma::uvec transformFromInner(int inner) const;
      arma::uvec transformFromInnerToBm(int inner) const;

      void sortEigensystem(arma::vec& eigenvalues, arma::cx_mat& eigenvectors) const;

      Eigensystem calculateEigensystem(int inner) const;
      Eigensystem calculateSubEigensystem(int inner, const arma::cx_mat& subspace);
      Eigensystem calculateSubEigensystem(Ipp ipp);
      arma::cx_mat calculateSubspace(Ipp ipp);
      std::tuple<int,int> getSubspaceRange(int n, NPoint cellLocation, NPoint wanLocation, bool last = false);
      const arma::cx_mat& getSubspace(int n, NPoint cellLocation, NPoint wanLocation);
      int wannierPosition(int n, NPoint cellLocation, NPoint wanLocation);

     arma::cx_mat rInnerMatrix(int inner, NPoint m, int bi, int bj) const;

      arma::cx_double sumOverNs(int inner, NPoint m, NPoint mp, int bi, int bj) const;
  };
}

#endif
