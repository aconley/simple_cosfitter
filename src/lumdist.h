//lumdist.h

#ifndef __lumdist__
#define __lumdist__

#include <vector>

#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>

#include <snedata.h>

/*!
  \defgroup lumdist Luminosity distance integrators

  Set of functions to calcluate luminosity distance integrals.
*/

/*!
  \brief Luminosity distance calculator

  \ingroup lumdist

  This is a class rather than a namespace so that the auxillary memory
  structures used by the gsl can be pre-allocated, rather than being
  allocated on each call.

  There are a number of routines avaiable.  Which is the fastest depends
  somewhat on the number of SN and how they are arranged.  For a single
  luminosity distance, the singleLumDist functions should always be
  the fastest.  For a large number of SN, one is better off sorting
  them by zcmb and using one of the Sum integrators.  In most cases,
  the rk4 method should be the best, but for a large enough number of
  SN (~600+) the trapezoidal rule is actually faster because it has
  smaller overheads.

  getLumDist attempts to choose between these routines for you.
  If you have a small number of SN in a \f$\Omega_{\Lambda}\f$
  universe, then it calls the single luminosity distance function
  repeatedly (this is fast because it can be written in terms
  of elliptic integrals).  If this isn't the case, or you have
  more SN, it switches to a RK4 scheme.  Then, for even more SN,
  it switches to a simple trapezoidal rule.  

  Experiments using splines to hold the distances proved rather
  unsatisfactory, as they seemed to cause considerable loss in
  precision.

*/
class lumdist {
 private:

  double param_tol; //!< How close param must be to special values to be treated as such (i.e., if abs(ocurv) < param_tol, we treat it as flat)

  //Form of w_a
  bool useKomatsuForm; //!< Use the Komatsu form for w(a)
  double atrans; //!< If using Komatsu form, transition scale factor \f$a_t\f$

  static const double prec_lum; //!< Controls desired precision of results

  double dz; //< Step size for trapezoidal rule

  //ODE stuff
  gsl_odeiv_step_type *T45; //!< Controls type of integration scheme
  gsl_odeiv_step *s_ode; //!< GSL ODE Stepper
  gsl_odeiv_control *c_ode; //!< GSL ODE control
  gsl_odeiv_evolve *e_ode; //!< GSL ODE evolution
  
  //Adaptive integrator
  gsl_integration_workspace *intspace; //!< GSL adaptive integration workspace

  mutable std::vector<double> zcmb; //!< Stores CMB frame redshifts

  //Working space
  mutable std::vector<double> integral_value; //!< Stores value of integral

  void init_gsl(); //<! Initializes GSL structs

  //Luminosity distance integrators

  int singleLumDistInt(double, double, double, double&) const; //!< Single luminosity distance, w=-1
  int singleLumDist(const SNeDataEntry&, double, double, double&) const; //!< Single luminosity distance, w=-1

  int lumDistSumInt(std::vector<double>& intval, double om, double ol) const; //!< Does integral bit, w=-1, ordered distances, trapezoidal rule
  int lumDistSum(const SNeData&, std::vector<double>& dl, double om, 
		 double ol) const; //!< luminosity distance, w=-1 multiple ordered distances, trapezoidal rule
  
  int lumDistSumInt_rk4(std::vector<double>& intval, double om, 
			double ol) const; //!< Does integral bit, w=-1, ordered distances, RK45 integration
  int lumDistSum_rk4(const SNeData&, std::vector<double>& dl, double, double) const; //!< luminosity distance, w=-1 multiple ordered distances, RK45 integration

  int singleLumDistW(const SNeDataEntry&, double, double, double,
		     double&) const; //!< Single luminosity distance, w=const

  int lumDistSumWInt(std::vector<double>& intval, double, 
		     double, double) const; //!< luminosity distance, w=const multiple ordered distances, trapezoidal rule, integral only
  int lumDistSumW(const SNeData&, std::vector<double>& dl, double, 
		  double, double) const; //!< luminosity distance, w=const multiple ordered distances, trapezoidal rule

  int lumDistSumWInt_rk4(std::vector<double>& intval, double, 
			 double, double) const; //!< luminosity distance, w=const multiple ordered distances, RK45 integration, integral only
  int lumDistSumW_rk4(const SNeData&, std::vector<double>& dl, double, 
		      double, double) const; //!< luminosity distance, w=const multiple ordered distances, RK45 integration

  int singleLumDistW0WA(const SNeDataEntry&, double, double, double, double,
			double&) const; //!< Single luminosity distance, \f$w_0, w_a\f$

  int lumDistSumW0WAInt(std::vector<double>& intval, double, double,
			double, double) const; //!< luminosity distance, $w_0, w_a$ multiple ordered distances, trapezoidal rule, integral only
  int lumDistSumW0WA(const SNeData&, std::vector<double>& dl, double, double,
		     double, double) const; //!< luminosity distance, $w_0, w_a$ multiple ordered distances, trapezoidal rule

  int lumDistSumW0WAInt_rk4(std::vector<double>& intval, double, double,
			 double, double) const; //!< luminosity distance, $w_0, w_a$ multiple ordered distances, RK45 integration, integral only
  int lumDistSumW0WA_rk4(const SNeData&, std::vector<double>& dl, double, 
			 double, double, double) const; //!< luminosity distance, $w_0, w_a$ multiple ordered distances, RK45 integration


  /*! \brief Value inside integral, w=-1 case */
  inline double calcVal(double z,double om,double ol) const {
    double opz=1.0+z;
    return 1.0/sqrt( opz*opz*(1+om*z)-z*(2+z)*ol );
  }
  /*! \brief Value inside integral, general w case.  
   wfac is -3(1+w)*/
  double calcValW(double z,double om, double ow, double ocurv, 
		  double wfac) const {
    double opz = 1.0 + z;
    return 1.0/sqrt( (opz*om + ocurv)*opz*opz + ow*pow(opz,wfac));
  }
  /*! \brief Value inside integral, \f$w_0, w_a\f$ case. 
   wfac is -3(1+w_0+w_a)*/
  double calcValW0WA(double z,double om, double ow, double ocurv, 
		     double wfac, double wa) const {
    double opz = 1.0 + z;
    return 1.0/sqrt( (opz*om + ocurv)*opz*opz +
		     ow * pow(opz, wfac) * exp( - 3.0*wa*z/opz) );
  }
  
  /*! \brief Value inside integral, \f$w_0, w_a\f$, Komatsu form */
  double calcValW0WAKom(double z, double om, double ow, double ocurv,
			double w0, double wa) const {
    double opz = 1.0 + z;
    double opat = 1.0 + atrans;
    // 3*(1+weff) = 3( (1-a)wa/ln a + (1+w0+(1+at)wa)*ln( (a+at)/(1+at) )/ln a)
    // a^3(1+weff) = exp( 3 [ (1-a) wa + 
    //  ( 1 + w0 + (1+at) wa) ln (a+at)/(1+at)] )
    double defac = exp( -3.0* ( z*wa/opz + (1+w0+opat*wa)*
				log( (1/opz+atrans)/opat ) ) );
    return 1.0/sqrt( (opz*om + ocurv)*opz*opz + ow*defac );
  }

 public:

  lumdist(); //!< Default constructor
  lumdist(const std::vector<double>&); //!< Assign zvals
  ~lumdist(); //!< Destructor

  enum int_method { AUTO, SINGLE, TRAPEZOIDAL, RK4 }; //!< Method to use

  double getParamTol() const { return param_tol; } //!< Get param_tol for user
  void setParamTol(double val) { param_tol = val; } //!< Set param_tol

  void setZcmb(const std::vector<double>&) const; //!< Set zcmb
  void setZcmb(const SNeData&) const; //!< Set zcmb

  double getDz() const { return dz; } //!< Return trapezoidal dz
  void setDz(const double val) { dz = val; } //!< Set dz for trapezoidal rule

  bool usingKomatsuForm() const { return useKomatsuForm; } //!< Are we using the Komatsu form
  void setUseKomatsuForm() { useKomatsuForm=true; }
  void unsetUseKomatsuForm() { useKomatsuForm = false; }
  double getAtrans() const { return atrans; }
  void setAtrans(double val) { atrans = val; }

  //Nice user interface that branch cuts between all the others
  // as needed
  /*! \brief User interface to luminosity distances */
  int getLumDist( const SNeData&, std::vector<double>& dl, double om,
		  double ode, double w0=-1.0, double wa=0.0, 
		  int_method method = AUTO) const; 
  int getLumDist( const SNeData&, std::vector<double>& dl, 
		  const std::vector<double>& params,
		  int_method method = AUTO) const; //!< Vector argument version
  int getLumDist( const SNeDataEntry&, double& dl, double om,
		  double ode, double w0=-1.0, double wa=0.0 ) const; //!< Single distance interface

};

#endif
