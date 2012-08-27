//Classes for handling 2dFGRS constraints on the cosmological
// parameters.

#ifndef __w2df__
#define __w2df__

#include <cosgrids.h>


/*!
  \brief Code for handling 2dFGRS growth of structure constraints

  \ingroup other_const
*/

namespace w2dFGRS {

  /*!
    \brief Mini class containing relevant cosmological parameters
  */
  class w2df_params {
  public:
    double om; //!< \f$\Omega_m\f$
    double w0; //!< w at z=0
    double wa; //!< Derivative of w with a.  Specifically, \f$w\left(a\right) = w0 + wa \left(1 - a\right)\f$
    
    w2df_params(double _om, double _w0, double _wa) {
      om = _om; w0 = _w0; wa = _wa;
    } //!< Constructor
  };
  
  /*!
    \brief Class for actually calculating probability surface
  */
  double fcalc(double,const w2df_params&,int); //!< Integrator of growth factor
  void rk4(double[2], double[2], double, double, double[2], const w2df_params&); //!< Rk4 integrator from Numerical Recipes
  void derivs(double, double[2], double[2], const w2df_params&); //!< Equation derivatives

  double get_growth_factor(double,double,double,double); //!< Return the growth factor at a given set of cosmological parameters
  
  cosgrid2D make_prob_grid(double,double,double,int,double,double,int,double,double); //!< Returns a grid of probability values from the 2dFGRS data
  cosgrid2D make_prob_grid(double,double,double,const cosgrid2D&); //!< Returns a grid of probability values from the 2dFGRS data
  
}

#endif
