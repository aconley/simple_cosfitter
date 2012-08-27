//paramfile.h

#ifndef __paramfile__
#define __paramfile__

#include <vector>
#include <map>

#include "param_tags.h"

/*!
  \brief Information about parameters we might be interested in

  \ingroup fitter
*/
struct param_struct {
  enum fittypes { not_set, dependent, fixed, analytic, loop }; //!< Ways we will treat this parameter

  std::string name; //!< Name of parameter, from param_tags
  param_tags::paramspec param_spec; //!< Code and type
  fittypes fit; //!< How we will treat this parameter
  double fixval; //!< Fixed value, if fittype is fixed
  double min; //!< Minimum value, if fittype is loop
  double max; //!< Maximum value, if fittype is loop
  unsigned int n; //!< Number of bins to use, if fittype is loop
  double dval; //!< Step size, if fittype is loop

  param_struct(); //!< Constructor

  

};

/*!
  \brief Sorting class for param_struct based on param_spec.param_codes
  \ingroup fitter
*/
class param_struct_SortByCode {
 public :
  /*! \brief Comparison operator on paramcodes */
  int operator()(const param_struct& p1, const param_struct& p2) {
    return p1.param_spec.first < p2.param_spec.first;
  }
};

/*!
  \brief Data structure holding parameters of fit.

  \ingroup fitter
*/
struct fitparam {

  std::string datafilename; //!<Data file to read from
  std::string outputfilename; //!< File to output marginalized prob
  std::string albetaoutputfilename; //!< File to output alpha/beta marg prob to
  std::string extendedoutfilename; //!< File to output extended fit residual info to
  std::string paramsummaryfilename; //!< File to write short parameter summary to
  std::string mag_covfilename; //!< File containing covariance matrix
  std::string width_covfilename; //!< File containing width covariance matrix
  std::string colour_covfilename; //!< File containing colour covariance matrix
  std::string magwidth_covfilename; //!< File containing mag-width covariance matrix
  std::string magcolour_covfilename; //!< File containing mag-colour covariance matrix
  std::string widthcolour_covfilename; //!< File containing width-colour covariance matrix

  std::string woodbury0_covfilename; //!< File containing woodbury style covariance matrix for \f$m_B\f$
  std::string woodburya_covfilename; //!< File containing woodbury style covariance matrix for \f$s\f$
  std::string woodburyb_covfilename; //!< File containing woodbury style covariance matrix for \f$\mathcal C\f$

  //Information about variables
  std::map< param_tags::paramcodes, param_struct > params; //!< Information about params

  double ocurv; //!< \f$\Omega_{curv}\f$, if being done with fixed value

  //How the nuisance params will be handled
  enum nuisance_loop { unknown, cosmo_outer, cosmo_inner };
  enum nuisance_strategy { not_set, m_analytic, m_analytic_a_loop,
			   m_analytic_b_loop, m_analytic_ab_loop };
  nuisance_loop loop_order; //!< Which way to loop on nuisance params
  nuisance_strategy strategy; //!< How are we going to deal with the nuisance params
  std::string strategy_string; //!< Description of nuisance strategy

  //Error propagation
  bool errsfixed; //!< Both \f$\alpha\f$ and \f$\beta\f$ are effectively fixed, so the errors can be totally pre-computed

  unsigned int nintrinsicdisp; //!< Number of independent intrinsic dispersions, in magnitudes (def:0, which means there is only one for all data)
  std::map< unsigned int, double > intrinsicdisp; //<! Map from dataset numbers to intrinsic dispersions, in magnitudes

  double  pecz;  //!< Peculiar velocity errors in redshift units (def: 0.001)
  double dzint;  //!< Size of trapezoidal integration step in lumdist integral (def: 0.0005). Not used currently.

  bool usekomatsuform; //!< Komatsu et al. form for w(a)
  double atrans; //!< Transition scale factor for Komatsu form

  bool albetaout; //!< Ouput marginalized \f$\alpha, \beta\f$.  (def: false)
  bool extendedout; //!< Print extended fit information to file (def: false)
  bool paramsummary; //!< Print short parameter summary to file (def: false)
  bool binaryout; //!< Output probability files in binary format (def: false)
  bool any_covfile; //!< Is there any covaraince matrix at all? (def: false)
  bool mag_covfileset; //!< Is a mag covariance matrix provided (def: false)
  bool width_covfileset; //!< Is a width covariance matrix provided (def: false)  
  bool colour_covfileset; //!< Is a colour covariance matrix provided (def: false)
  bool magwidth_covfileset; //!< Is a mag-width covariance matrix provided (def: false)
  bool magcolour_covfileset; //!< Is a mag-colour covariance matrix provided (def: false)
  bool widthcolour_covfileset; //!< Is a width-colour covariance matrix provided (def: false)

  bool woodbury_covfileset; //!< Is a Woodbury covariance matrix provided (def: false)

  bool flatonlyfit;  //!< Assume flat (i.e, \f$\Omega_{curv}=0\f$) (def: false)
  bool fixcurv; //!< Do fit with fixed \f$\Omega_{curv}\f$ (def: false)
  bool verbose;  //!< If set, more verbose messages are printed as it runs.
  bool showprogbar; //!< Shows the progress meter

  double scriptmcut; //!< Cut on third parameter to decide which scriptmset

  fitparam(); //!< Default constructor.  Sets default values above

  void readFile(const std::string& filename, bool silent=false); //!< Reads in from parameter file.

};


#endif
