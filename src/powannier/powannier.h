#ifndef POWANNIER_POWANNIER_H
#define POWANNIER_POWANNIER_H

#include <memory>
#include "blochsystem.h"
#include "blochhelpers.h"

namespace POWannier {
  /**
   * @brief %POWannier function.
   * 
   * POWannier functions are localized functions defined in the %Bloch
   * functions space (or rather, its subspace for a given set of bands
   * \f$ \{ \alpha \} \f$) in the system with periodic potential:
   * \f$ W_{\vec{R}} (\vec{r}) = \sum_{\alpha, \vec{k}} 
   * w_{\vec{R},\alpha,\vec{k}} \phi_{\vec{k}} (\vec{r}) \f$.
   * In this library, the \f$ w_{\vec{R},\alpha,\vec{k}} \f$ coefficients are 
   * obtained by position operators diagonalization.
   * 
   * The creation of the POWannier is done through RSystem's getWannier()
   * method.
   * @see Bloch
   * @see RSystem
   */
  class POWannier : public Wannier {
    public:
      friend double tunneling(const POWannier& wan1, const POWannier& wan2);

      /// Default constructor
      POWannier() = default;

      /**
       * @brief Basic constructor.
       * @param bs
       *   BlochSystem containing %Bloch functions of the system.
       * @param coefficients
       *   Coefficients \f$ w_{\vec{R},\alpha,\vec{k}} \f$ that define
       *   the %POWannier function.
       * @param bands
       *   Specify from which bands %Bloch functions are to be used.
       *   Default value is {0}.
       * @see coefficients()
       */
      POWannier(std::shared_ptr<BlochSystem> bs, 
              const arma::cx_vec& coefficients,
              std::vector<int> bands = {0});

      /**
       * @brief Value of the %POWannier function at a given position.
       * @param r
       *   Given position \f$ \vec{r} \f$.
       * @returns
       *   Value of \f$ W_{\vec{R}}(\vec{r}) \f$.
       */
 
      double operator()(Position r);

      /**
       * @brief Simple getter.
       * @returns The coefficients \f$ w_{\vec{R},\alpha,\vec{k}} \f$.
       * The coefficients are given in the order as given by bmIndex().
       */
      const arma::cx_vec& coefficients() {
        return _coefficients;
      }
      
      const arma::cx_double coefficient(int mi, int band) {
        return _coefficients[bandsi * _bs->ms.size() + mi];
      }

      /**
       * @brief Write values of the %POWannier functions over a given range to the specified file.
       * @see POPOWannier::draw(const LatticeBasis&,Function&&,std::string,int,Position,Position)
       */
      void draw(std::string filename, int density, Position beg, Position end);

      /**
       * @brief Write values of the %POWannier functions over a given range to the specified file.
       * @see POPOWannier::draw(const BlochSystem&,Function&&,std::string,int,DrawRange)
       */
      void draw(std::string filename, int density, DrawRange drawRange);

      /**
       * @brief Simple getter. 
       * @returns  Band indices \f$ \{ \alpha \} \f$.
       */
      const std::vector<int>& bands();

      int N() {
        return _bs->N;
      }

    private:
      Complex complexValue(Position r);
      std::vector<int> _bands;
      std::shared_ptr<BlochSystem> _bs;
      arma::cx_vec _coefficients;
  };
}
#endif