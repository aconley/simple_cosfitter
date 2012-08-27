//fitter.h

#ifndef __fitter__
#define __fitter__

#include <string>
#include <ostream>
#include <vector>
#include <map>

#include <ctime>

#include <paramfile.h>
#include <cosgrids.h>
#include <snedata.h>
#include <lumdist.h>
#include <lumdist_array.h>
#include <covmatrix.h>
#include <param_results.h>

/*!
  \defgroup fitter Cosmological Fitter
  \brief Handles the actual luminosity distance fitting

  This grouping consists of routines which control the actual
  luminosity distance fits to supernovae data.
*/


/*! 
  \brief Class which performs cosmological fits. 

  \ingroup fitter

  This should always be prepped (by calling cosfitter::prep()) before
  running a fit with dofit.  
  \code
  cosfitter ft;
  ft.prep( paramfilename );
  ft.dofit();
  \endcode
  The prep code isn't part of the constructor because it has the
  potential to throw exceptions, and you really don't want that inside
  a constructor.  Maybe someday I'll write a safer version.
*/
class cosfitter {
 private:

#ifdef TIMING
  mutable clock_t tot_inv_time, tot_mult_time, tot_nondia_time;
#endif

  enum cosmo_fittype { flat_omegam, flat_w0, omegam_omegade, flat_omegam_w0,
		       omegam_omegade_w0, flat_omegam_w0_wa }; //!< Types of cosmology fit we support

  static const double inv_twopi;  //!< \f$1 / 2 \pi\f$

  bool isprepped;  //!< Are we ready to run a fit?

  std::string paramfilename; //!< Name of parameter file

  SNeData sne; //!< Supernova data
  unsigned int nsn; //!< Number of supernovae

  bool have_fixed_alpha; //!< Is \f$\alpha\f$ present and fixed
  bool have_fixed_beta; //!< Is \f$\beta\f$ present and fixed

  std::vector<double> pre_vars; //!< Variance array for non-parameter dependent portion.  Will hold full errors if \f$\alpha\f$ and \f$\beta\f$ are fixed.

  //Convenience variables for inverse covariance matrix
  mutable covMatrix invcovmatrix; //!< Holds current inverse covariance matrix \f$V^{-1}\f$.  Keeps us from having to continually re-allocate

  mutable std::vector<double> A1; //!< Has 1 where in scriptmset 1
  mutable bool has_A2; //!< We have some elements in A2
  mutable std::vector<double> A2; //!< Has 2 where in scriptmset 2
  mutable std::vector<double> invcovmatrix_A1; //!< Holds the sums of the rows of the inverse covariance matrix \f$V^{-1} \cdot \vec{A1}\f$
  mutable std::vector<double> invcovmatrix_A2; //!< Holds the sums of the rows of the inverse covariance matrix \f$V^{-1} \cdot \vec{A2}\f$

  //Analytic marginalization convenience variables
  mutable double amarg_D;  //!< \f$ \vec{A_1}^T \cdot  V^{-1}  \cdot \vec{A_2} \f$, useful for analytic marginalization
  mutable double amarg_E; //!< \f$ \vec{A_1}^T \cdot  V^{-1}  \cdot \vec{A_1} \f$, useful for analytic marginalization
  mutable double amarg_F; //!< \f$ \vec{A_2}^T \cdot  V^{-1}  \cdot \vec{A_2} \f$, useful for analytic marginalization

  //Convenience arrays for non-diagonal case
  mutable double* diffarr;  //Difference array between SN mags and cosmo prediction, non-diagonal case
  mutable double* prod; //!< Internal convenience array for products in \f$\chi^2\f$ evaluation, non-diagonal case

  lumdist lm;
  mutable std::vector<double> dl; //!< Luminosity distances working variable
  mutable std::vector<double> working_vec; //!< General purpose working variable

  void calcPreErrs(); //!< Calculates non-parameter dependent part of errors

  void readParamFile(); //!< Initializes fit parameters from paramfilename

  void print_results( std::ostream&, const std::string&, 
		      const std::map< param_tags::paramcodes,
		      param_results >&) const; //!< Print extended results to output stream
  void write_paramsummary(const std::map< param_tags::paramcodes,
			  param_results >&) const; //!< Writes short parameter summary to file

  //Actual fitters, organized by the ordering of the cosmological
  // paramter loops and number of cosmologicalp parameters
  void cosmofit0D_diagonal( std::map< param_tags::paramcodes,
			 param_results>&, double& ) const; //!< 0D fit, cosmological parameters fixed, diagonal errors
  void cosmofit0D_nondiagonal( std::map< param_tags::paramcodes,
			 param_results>&, double& ) const; //!< 0D fit, cosmological parameters fixed, non diagonal errors
  void cosmofit1D_outer( cosfitter::cosmo_fittype cosmo_fit, 
			 std::map< param_tags::paramcodes,
			 param_results>&, double& ) const; //!< 1D fit, cosmological params as outer loop
  void cosmofit1D_inner( cosfitter::cosmo_fittype cosmo_fit, 
			 std::map< param_tags::paramcodes,
			 param_results>&, double& ) const; //!< 1D fit, cosmological params as inner loop
  void cosmofit2D_outer( cosfitter::cosmo_fittype cosmo_fit,
			 std::map< param_tags::paramcodes,
			 param_results>&, double& ) const; //!< 2D fit, cosmological params as outer loop
  void cosmofit2D_inner( cosfitter::cosmo_fittype cosmo_fit,
			 std::map< param_tags::paramcodes,
			 param_results>&, double& ) const; //!< 2D fit, cosmological params as inner loop
  void cosmofit3D_outer( cosfitter::cosmo_fittype cosmo_fit,
			 std::map< param_tags::paramcodes,
			 param_results>&, double& ) const; //!< 3D fit, cosmological params as outer loop
  void cosmofit3D_inner( cosfitter::cosmo_fittype cosmo_fit,
			 std::map< param_tags::paramcodes,
			 param_results>&, double& ) const; //!< 3D fit, cosmological params as inner loop

  //Likelihood calcluators.  We have one for each combination of
  // diagonal/non-diagonal errors and nuisance parameter types
  double calcLikelihood() const; //!< No \f$\alpha, \beta\f$, diagonal errors
  void calcLikelihood_la(cosgrid1D&) const; //!< No \f$\beta\f$, diagonal errors
  void calcLikelihood_lb(cosgrid1D&) const; //!< No \f$\alpha\f$, diagonal errors
  void calcLikelihood_lab(cosgrid2D&) const; //!< \f$\alpha, \beta\f$, diagonal errors

  double calcLikelihood_nondia( double, double, const lumdist_array0D& ) const; //!< 0D, non-diagonal case
  void calcLikelihood_nondia( cosgrid1D&, double, double,
			      const lumdist_array1D& ) const; //!< 1D, non-diagonal case
  void calcLikelihood_nondia( cosgrid2D&, double, double, 
			      const lumdist_array2D& ) const; //!< 2D, non-diagonal case
  void calcLikelihood_nondia( cosgrid3D&, double, double,
			      const lumdist_array3D& ) const; //!< 3D, non-diagonal case

  //Fit drivers
  void cosmofit0D( std::map< param_tags::paramcodes, param_results>&, 
		   double&, std::string& ) const; //!< cosmological params fixed
  void cosmofit1D( std::map< param_tags::paramcodes, param_results>&,
		   double&, std::string& ) const; //!< Driver for fit to one cosmological parameter
  void cosmofit2D( std::map< param_tags::paramcodes, param_results>&,
		   double&, std::string& ) const; //!< Driver for fit to two cosmological parameters
  void cosmofit3D( std::map< param_tags::paramcodes, param_results>&,
		   double&, std::string& ) const; //!< Driver for fit to three cosmological parameters

  //Covariance matrix 
  void computeInvCovMatrix( double alpha=0.0, double beta=0.0 ) const; //!< Calculates inverse covariance matrix

 public:

  cosfitter(); //!< Default constructor
  ~cosfitter(); //!< Destructor

  fitparam fparam; //!< Parameters of fit

  unsigned int getNsn() const { return sne.size(); } //!< Return number of SNe
  std::set<unsigned int> getDataSetList() const { 
    return sne.getDataSetList(); 
  } //!< Return the set of dataset numbers

  void prep(const std::string&); //!< Reads in parameter file and supernova data
  std::string getParamFileName() const { return paramfilename; } //!< Returns name of parameter file

  void dofit() const; //!< Global fit driver.

  void indivchi(const std::map< param_tags::paramcodes, param_results >&,
		const std::map< unsigned int, double >&) const; //!< Prints out individual \f$ \chi^2 \f$ contributions

  std::pair<double,double> estimate_scriptm(const std::map< param_tags::paramcodes, param_results >& ) const; //!< Given the other parameters, gives the best value of \f${\mathcal M}\f$


};

#endif
