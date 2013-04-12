#ifndef __auxconstraint__
#define __auxconstraint__

#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "cosgrids.h"
/*!
  \defgroup other_constraint  Constraints from other measurements
*/

/*!
  \brief Auxillar constraints
  
  \ingroup other_constraint

  A namespace for other types of constraints, like BAO and WMAP.
*/

namespace auxconstraint {
  //PDG (2004) value
  //const ogammah2 = 2.471e-5; //!< \f$\Omega_{\gamma} h^2\f$

  //Komatsu et al. value
  const double ogammah2 = 2.469e-5; //!< \f$\Omega_{\gamma} h^2\f$

  const double neff = 3.04; //!< Number of Neutrino species
  const double oradh2 = ogammah2*(1.0+0.2271*neff); //!< \f$\Omega_{r} h^2\f$

  double one_over_E( double z, void * params );  //!< \f$1/E, w=-1\f$
  double one_over_E_w( double z, void * params ); //!< \f$1/E\f$ w=const
  double one_over_E_w0wa( double z, void * params ); //!< \f$1/E , w_0, w_a\f$
  double one_over_E_w0waKom( double z, void * params ); //!< \f$1/E , w_0, w_a, Komatsu form\f$

  double one_over_E_rad( double z, void * params );  //!< \f$1/E, w=-1\f$ with radiation
  double one_over_E_w_rad( double z, void * params ); //!< \f$1/E\f$ w=const, with radiation
  double one_over_E_w0wa_rad( double z, void * params ); //!< \f$1/E , w_0, w_a\f$, with radiation
  double one_over_E_w0waKom_rad( double z, void * params ); //!< \f$1/E , w_0, w_a\f$, Komatsu form with radiation
  
  double rsint( double a, void * params ); //!< \f$r_s \left( z_{cmb} \right)\f$ integral
  double rsint_w( double a, void * params ); //!< \f$r_s \left( z_{cmb} \right)\f$ integral with \f$w\f$
  double rsint_w0wa( double a, void * params ); //!< \f$r_s \left( z_{cmb} \right)\f$ integral with \f$w_0, w_a\f$.
  double rsint_w0waKom( double a, void * params ); //!< \f$r_s \left( z_{cmb} \right)\f$ integral with \f$w_0, w_a\f$, Komatsu form

  double rsHovercint( double a, void * params); //<! \f$r_s H_0 / c\f$ integral
  double rsHovercint_w( double a, void * params); //<! \f$r_s H_0 / c\f$ integral with \f$w\f$
  double rsHovercint_w0wa( double a, void * params); //<! \f$r_s H_0 / c\f$ integral with \f$w_0, w_a\f$
  double rsHovercint_w0waKom( double a, void * params); //<! \f$r_s H_0 / c\f$ integral with \f$w_0, w_a\f$, Komatsu form

  bool isRebounding(double, double); //!< Check for rebounding Universe if \f$w=-1\f$

  //Can't be part of member class because we are interfacing to the GSL,
  // which doesn't handle member functions
  double percival_chisq(const gsl_vector *v, void *params); //<! \f$\chi^2\f$ of Percival BAO fits
  double wmap7_shift_chisq(const gsl_vector *v, void *params); //<! \f$\chi^2\f$ of WMAP7 shift parameters
  double percival_wmap7_chisq(const gsl_vector *v, void *params); //<! \f$\chi^2\f$ of WMAP7 shift parameters and Percival BAO fits

  //We do the constraints as classes because we need to initialize
  // internal work space for the integrator, and it's nice to 
  // re-use that allocation

  /*!
    \brief Helper class for some functions common to WMAP shift and BAO
  */
  class distance_helper {
  private:
    gsl_integration_workspace *work; //!< Workspace for numeric integratio
    bool useKomatsuForm; //!< Use Komatsu et al. (2008) form for \f$w\left(a\right)\f$
    double atrans; //!< Transition scale factor if Komatsu form used
  public:
    distance_helper(bool fixcurve=false, double ocurve=0.0); //!< Constructor
    ~distance_helper(); //!< Destructor

    double GetAtrans() const { return atrans; } //!< Returns transition a
    void SetAtrans(double val) { atrans = val; } //!< Sets transition a
    bool UsingKomatsuForm() const { return useKomatsuForm; } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { useKomatsuForm=true; } //!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { useKomatsuForm = false; } //!< Turns off the Komatsu form

    double GetK1( double, double, double, double, double, double, double, bool& ) const; //<! Get \f$K_{1} = H_0 r_{s} \left( z  \right) / c\f$
    double GetK2( double, double, double, double, double, double, bool& ) const; //<! Get \f$K_{w} = \left(1 + z_{\star}\right) H_0 D_{A} \left( z \right) / c\f$
    double GetOneOverE(double, double, double, double, double, double, bool& ) const; //!< Get \f$1/E = H_0 / H\f$

    double GetZstar( double, double ) const; //<! \f$z_{\star}\f$
    double GetZdrag( double, double ) const; //<! \f$z_d\f$

  };

  /*!
    \brief Abstract base class for auxillary constraints

    This is not pure abstract, since some of the actual code is common
  */
  //There are a number of functions here that I would like to make const,
  // but can't because of interface issues with the GSL for WMAP7 and 
  // baoP09
  class base_auxconstraint {
  protected:
    void makeProbSurface1D(const std::map< param_tags::paramcodes,
			   param_struct >&, const std::string&,
			   bool binaryout=false); //!< 1D version
    void makeProbSurface2D(const std::map< param_tags::paramcodes,
			   param_struct >&, const std::string&,
			   bool binaryout=false); //!< 2D version
    void makeProbSurface3D(const std::map< param_tags::paramcodes,
			   param_struct >&, const std::string&,
			   bool binaryout=false); //!< 3D version

    void updateStandardCosVec( param_tags::paramcodes, std::vector<double>&,
			       double ) const;
    bool fixcurv; //!< Fixed curvature set
    double ocurv; //!< \f$\Omega_k\f$, if fixcurv set

  public:

    base_auxconstraint(bool fixcurve=false, double ocurv=0.0); //!< Constructor
    virtual ~base_auxconstraint(); //!< Destructor

    void SetFixCurv(double val) { fixcurv=true; ocurv=val; } //!< Set fixed curvature
    bool IsFixedCurv() const { return fixcurv; } //!< Do we have fixed curvature
    void UnsetFixCurv() { fixcurv = false; } //!< Unset fixed curvature

    virtual double GetChiSq(double,double,double,double) = 0; //!< Returns \f$\chi^2\f$.
    double GetChiSq(const std::vector<double>&); //!< Returns \f$\chi^2\f$.

    void makeProbSurface(const std::map< param_tags::paramcodes,
			 param_struct >&, const std::string&); //!< Makes probability surface
  };


  /*!
    \brief Class for Baryon Acoustic Peak from the SDSS LRGs using the
    approach of Eisenstein et. al (2005)
  */
  class baoE05 : public base_auxconstraint {
  private:
    mutable gsl_integration_workspace *work; //!< Workspace for numeric integration
    mutable bool qag_alloc; //!< Is qag workspace allocated

    double z;  //!< Redshift of galaxy sample
    double a_meas; //! Measured value of A (def: 0.469)
    double a_err;  //!< Error in A (def: 0.017)
    double ns; //!< Matter spectral index
    double fnu; //!< Neutrino mass fraction \f$\Omega_{\nu}/\Omega_m\f$

    bool useKomatsuForm; //!< Use Komatsu et al. (2008) form for \f$w\left(a\right)\f$
    double atrans; //!< Transition scale factor if Komatsu form used

    void allocateQAGWork() const; //Allocate the qag workspace for adaptive integration

  public:
    baoE05(double z=0.35, double a_meas=0.469, double a_err=0.017,
	   double ns=0.98, double fnu=0.0, bool fixcurv=false,
	   double ocurv=0.0 ); //!< Constructor
    ~baoE05(); //!< Destructor
    
    void SetZ(double val) { z = val; } //!< Set $z$ of measurement
    void SetA(double val) { a_meas = val; } //!< Set measured A to compare with
    void SetAErr(double val) { a_err = val; } //!< Set \f$\sigma_A\f$
    void SetNs(double val) { ns = val; } //!< Set Matter spectral index
    void SetFnu(double val) { fnu = val; } //!< Set Neutrino mass fraction

    double GetAtrans() const { return atrans; } //!< Returns transition a
    void SetAtrans(double val) { atrans = val; } //!< Sets transition a
    bool UsingKomatsuForm() const { return useKomatsuForm; } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { useKomatsuForm=true; } //!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { useKomatsuForm = false; } //!< Turns off the Komatsu form

    double aval(double,double,double,double) const; //!< Return calculated of A parameter
    double aval(const std::vector<double>&) const; //!< Return calculated A

    double GetChiSq(double,double,double,double); //!< Returns \f$\chi^2\f$.
  };


  /*!
    \brief Class for Percival et al. (2009) SDSS+2dF BAO constraints
  */
  class baoP09 : public base_auxconstraint {
  private:

    //GSL stuff
    static const gsl_multimin_fminimizer_type *T; //!< Type of fitter
    gsl_vector * v; //!< Internal use vector for simplex minimization
    gsl_vector *ss; //!< Step sizes for minimization
    gsl_multimin_fminimizer *s;
    gsl_multimin_function minex_func;

    distance_helper dhelp; //!< Class for doing integrals
  public:
    baoP09(bool fixcurve=false, double ocurv=0.0); //!< Constructor
    ~baoP09(); //!< Destructor

    double GetAtrans() const { return dhelp.GetAtrans(); } //!< Returns transition a
    void SetAtrans(double val) { dhelp.SetAtrans(val); } //!< Sets transition a
    bool UsingKomatsuForm() const { return dhelp.UsingKomatsuForm(); } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { dhelp.SetUseKomatsuForm(); }//!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { dhelp.UnsetUseKomatsuForm(); } //!< Turns off the Komatsu form

    double GetChiSq(double,double,double,double); //!< Returns \f$\chi^2\f$.

    double GetZdrag( double a, double b ) const { return dhelp.GetZdrag(a,b); } //<! \f$z_{d}\f$

  };


  /*!
    \brief Class for WMAP 3rd year distance to last scattering

    \f$H_0\f$ free form from Wang & Mukherjee (2007)
  */
  class wmap3yr_dls : public base_auxconstraint {
  private:
    gsl_integration_workspace *work; //!< Workspace for numeric integration

    bool usela; //!< If set, then the second shift parameter is also used.

    double h; //!< Prior on the Hubble constant
    double h_error; //!< Uncertainty in Hubble constant prior

    unsigned int nhsteps; //!< Number of steps in h for numerical marginalization, only used if \f$l_a\f$ is
    unsigned int nombh2steps; //!< Number of steps in \f$\Omega_b h^2\f$ for numerical marginalization, only used if \f$l_a\f$ is

    bool useKomatsuForm; //!< Use Komatsu et al. (2008) form for \f$w\left(a\right)\f$
    double atrans; //!< Transition scale factor if Komatsu form used

    double GetChiSq_Ronly( double, double, double, double ) const; //!< \f$\chi^2\f$ using only the first shift parameter
    double GetChiSq_Rla( double, double, double, double ) const; //!< \f$\chi^2\f$ using both shift parameters and the \f$h\f$ prior.

  public:
    wmap3yr_dls(bool fixcurve=false, double ocurv=0.0); //!< Constructor
    ~wmap3yr_dls(); //!< Destructor

    void SetUseLa() { usela = true; } //!< Set to use \f$l_a\f$ 
    void UnsetUseLa() { usela = false; } //!< Set to not use \f$l_a\f$

    void SetHubblePrior( double H, double H_ERROR ) {
      h = H; H_ERROR = h_error;
    }  //<! Set the Hubble prior.  Only used if the second shift param is used

    void SetHubblePriorMean( double H ) { h = H; } //!< Set mean value of Hubble prior
    void SetHubblePriorSigma( double HERROR ) { h_error = HERROR; } //!< Set sigma of hubble prior

    std::pair<double,double> GetHubblePrior() const; //!< Return the mean and sigma for the Hubble prior

    double GetAtrans() const { return atrans; } //!< Returns transition a
    void SetAtrans(double val) { atrans = val; } //!< Sets transition a
    bool UsingKomatsuForm() const { return useKomatsuForm; } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { useKomatsuForm=true; } //!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { useKomatsuForm = false; } //!< Turns off the Komatsu form

    double rval(double,double,double,double) const; //!< Value of R parameter
    double rval(const std::vector<double>&) const; //!< Value of R parameter

    double laval(double,double,double,double,double,double) const; //!< Value of \f$l_a\f$ parameter
    double laval(const std::vector<double>&, double, double) const; //!< Value of \f$l_a\f$ parameter
    double laval( double, double, double, double, double, double, double ) const; //!< Value of \f$l_a\f$ if you've already calculated R
    double laval( double, const std::vector<double>&, double, double ) const; //!< Value of \f$l_a\f$ if you've already calculated R

    double GetChiSq(double,double,double,double); //!< Returns \f$\chi^2\f$.
  };


  /*!
    \brief Class for WMAP 5th year distance to last scattering

    Prescription due to Komatsu et al. (2008) 
  */
  class wmap7yr_dls : public base_auxconstraint {
  private:

    //GSL stuff
    static const gsl_multimin_fminimizer_type *T; //!< Type of fitter
    gsl_vector * v; //!< Internal use vector for simplex minimization
    gsl_vector *ss; //!< Step sizes for minimization
    gsl_multimin_fminimizer *s;
    gsl_multimin_function minex_func;

    distance_helper dhelp; //!< Class for doing integrals
  public:
    wmap7yr_dls(bool fixcurve=false, double ocurv=0.0); //!< Constructor
    ~wmap7yr_dls(); //!< Destructor

    double GetAtrans() const { return dhelp.GetAtrans(); } //!< Returns transition a
    void SetAtrans(double val) { dhelp.SetAtrans(val); } //!< Sets transition a
    bool UsingKomatsuForm() const { return dhelp.UsingKomatsuForm(); } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { dhelp.SetUseKomatsuForm(); }//!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { dhelp.UnsetUseKomatsuForm(); } //!< Turns off the Komatsu form

    double GetChiSq(double,double,double,double); //!< Returns \f$\chi^2\f$.
  };


  /*!
    \brief Class for joint WMAP 5th year distance to last scattering
    and Percival et al. (2009) SDSS+2dF BAO constraints

    This is useful because both fits have to minimize over the
    same nuisance variables, so it should be faster to do both at once.
  */
  class baoP09_wmap7yr_dls : public base_auxconstraint {
  private:
    //GSL stuff
    static const gsl_multimin_fminimizer_type *T; //!< Type of fitter
    gsl_vector * v; //!< Internal use vector for simplex minimization
    gsl_vector *ss; //!< Step sizes for minimization
    gsl_multimin_fminimizer *s;
    gsl_multimin_function minex_func;

    distance_helper dhelp; //!< Class for doing integrals
  public:
    baoP09_wmap7yr_dls(bool fixcurve=false, double ocurv=0.0); //!< Constructor
    ~baoP09_wmap7yr_dls(); //!< Destructor

    double GetAtrans() const { return dhelp.GetAtrans(); } //!< Returns transition a
    void SetAtrans(double val) { dhelp.SetAtrans(val); } //!< Sets transition a
    bool UsingKomatsuForm() const { return dhelp.UsingKomatsuForm(); } //!< Is the Komatsu form for \f$w \left(a\right)\f$ in use?
    void SetUseKomatsuForm() { dhelp.SetUseKomatsuForm(); }//!< Turns on Komatsu form
    void UnsetUseKomatsuForm() { dhelp.UnsetUseKomatsuForm(); } //!< Turns off the Komatsu form

    double GetChiSq(double,double,double,double); //!< Returns \f$\chi^2\f$.
  };


  //Gaussian prior on Om
  class omprior : public base_auxconstraint {
  private:
    double ommean; //!<Mean om of prior
    double omerror; //!< Error on om for prior
  public:
    omprior(double om=0.28, double omerr=0.05, bool fixcurv=false,
	    double ok=0.0); //!< Constructor
    omprior(std::pair<double,double> vals, bool fixcurv=false,
	    double ok=0.0); //!< Constructor

    std::pair<double,double> GetPriorVals() const {
      return std::pair<double,double>(ommean,omerror); } //!< Return prior vals

    void SetPriorVals(std::pair<double,double> vals) { ommean = vals.first;
      omerror = vals.second; } //!< Set values as pair
    
    double GetChiSq(double,double,double,double); //!< Return \f$\chi^2\f$
  };

}

#endif
