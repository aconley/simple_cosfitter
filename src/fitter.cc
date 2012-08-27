//fitter.cc

//Cosmology fitting program, assuming gaussian errors

#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cstdio>

#include <param_tags.h>
#include <param_results.h>
#include <fitter.h>
#include <cosfitterexcept.h>
#include <cosgrids.h>
#include <snedata.h>
#include <utility.h>

#if USEMKL
#include <mkl_cblas.h>
#elif USEATLAS
extern "C" {
  #include <cblas.h>
}
#elif USEACCELERATE
#include <Accelerate/Accelerate.h>
#endif

using namespace std;

const double cosfitter::inv_twopi = 
  0.5/3.14159265358979323846264338327950288419716939937510582;

/////////////////////////////////////////////////
//              cosfitter                      //
/////////////////////////////////////////////////

cosfitter::cosfitter() {
  nsn = 0;
  amarg_E = 0.0;
  pre_vars.clear();
  have_fixed_alpha = false;
  have_fixed_beta = false;
  isprepped = false;
  prod = diffarr = NULL;
}

cosfitter::~cosfitter() {
  if (prod != NULL) delete[] prod;
  if (diffarr != NULL) delete[] diffarr;
}

/*!
  Performs the part of the error propagation that does not depend on
 the parameters being fit by incorporating the intrinsic dispersion
 and redshift errors with the magnitude errors.  Thus, we don't need
 to keep recalculating these.  If we are doing a fit where \f$\alpha\f$
 and \f$\beta\f$ have been completely pre-calculated, then we
 go even further and calculate the full errors for each SN.

 Also, if the nuisance parameters are fixed, and we have a covariance
 matrix we can pre-invert it for later use.

 Note that we are ignoring the fact that the redshift errors converted
 into magnitude space actually depend on the cosmology.  This is a
 very good approximation as long as you have spectroscopic redshifts
 for your SN.  Someday, when people are using photo-zs, this will have
 to be changed.

 If we have fixed values of \f$\alpha\f$ and \f$\beta\f$ at this point
 we will incorporate those parts of the errors.  Note that if you
 have only one fixed, there is a cross term that has to be evaluated
 later.
*/
void cosfitter::calcPreErrs() {

  //dz multiplicative factor
  const double zfacsq = 25.0 / ( log(10.0)*log(10.0) );

  bool fixedintrinsic; //!< True if there is only one for all sne
  double currintrinsicsq;
  unsigned int i;

  if ( fparam.nintrinsicdisp == 0 ) {
    fixedintrinsic = 1;
    currintrinsicsq = fparam.intrinsicdisp[0] * fparam.intrinsicdisp[0];
  } else {
    fixedintrinsic = 0;
    currintrinsicsq = 0.0;
  }

  double dzerrsq, emptyfac;

  double curr_alpha, curr_beta;
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  it = fparam.params.find( param_tags::alpha );
  if (it != fparam.params.end() && it->second.fit == param_struct::fixed) {
    have_fixed_alpha = true; curr_alpha = it->second.fixval;
  } else { have_fixed_alpha = false; curr_alpha = 0.0; }
  it = fparam.params.find( param_tags::beta );
  if (it != fparam.params.end() && it->second.fit == param_struct::fixed) {
    have_fixed_beta = true; curr_beta = it->second.fixval;
  } else { have_fixed_beta = false; curr_beta = 0.0; }

  if ( pre_vars.size() != nsn )
    pre_vars.resize(nsn);

  if (fparam.errsfixed) {
    std::map< std::string, param_struct >::const_iterator it;
    //All the way, meaning we are using a single fixed values for alpha/beta
    // or aren't using them at all.

    //Get diagonal terms in covariance matrix
    //We can always pre-calculate this part
    //Now uses empty universe formula
    if (fixedintrinsic) {
      for (i = 0; i < nsn; ++i) {
	emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
	dzerrsq = fparam.pecz*fparam.pecz + sne[i].var_z;
	dzerrsq *= zfacsq*emptyfac*emptyfac;
	pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
      }
    } else {
      std::map<unsigned int, double>::iterator pos;
      for (i = 0; i < nsn; ++i) {
	pos = fparam.intrinsicdisp.find( sne[i].dataset );
	if (pos == fparam.intrinsicdisp.end()) {
	  std::stringstream errstr;
	  errstr << "Unknown dataset number " << sne[i].dataset;
	  errstr << " for sn: " << sne[i].name;
	  throw CosFitterExcept("cosfitter","calcPreErrs",
				errstr.str(),4);
	}
	currintrinsicsq = pos->second * pos->second;

	emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
	dzerrsq = fparam.pecz*fparam.pecz + sne[i].var_z;
	dzerrsq *= zfacsq*emptyfac*emptyfac;
	pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
      }
    }

    //We may or may not be able to do these, depending on whether
    // we have single values of them
    if (have_fixed_alpha) {
      double curr_alpha_squared = curr_alpha * curr_alpha;
      for (i = 0; i < nsn; ++i) 
	pre_vars[i] += curr_alpha_squared * sne[i].var_widthpar
	  + (2.0 * curr_alpha * sne[i].cov_mag_widthpar);
    }
    
    if (have_fixed_beta) {
      double curr_beta_squared = curr_beta * curr_beta;
      if (have_fixed_alpha) {
	double curr_alpha_beta = curr_alpha * curr_beta;
	for (i = 0; i < nsn; ++i)
	  pre_vars[i] += curr_beta_squared * sne[i].var_colourpar 
	    - (2.0 * curr_beta * sne[i].cov_mag_colourpar)
	    - (2.0 * curr_alpha_beta * sne[i].cov_widthpar_colourpar);
      } else {
	for (i = 0; i < nsn; ++i)
	  pre_vars[i] += curr_beta_squared * sne[i].var_colourpar 
	    - (2.0 * curr_beta * sne[i].cov_mag_colourpar);
      }
    }

    //Pre-invert covaraince matrix
    if ( !sne.areErrorsDiagonal() ) {

      if (fparam.verbose) 
	cout << "Pre-inverting covariance matrix as errs are fixed...";

      //This also fills in the row sums and amarg_E
      computeInvCovMatrix( curr_alpha, curr_beta );

      if ( amarg_E == 0.0 )
	throw CosFitterExcept("cosfitter","calcPreErrs",
			      "Sum of inverse cov matrix entries == 0.0",8);

      if (fparam.verbose) cout << "done." << std::endl;

    } else {
      //Cov matrix is diagonal. 
      //Calculate the auto-marginalization parameter, which we can
      // do because the errors are fixed
      double current_invvar;
      current_invvar = 1.0/pre_vars[0];
      amarg_E = current_invvar;

      for (unsigned int i = 1; i < nsn; ++i) {
	current_invvar = 1.0 / pre_vars[i];
	amarg_E += current_invvar;
      }
      if ( amarg_E == 0.0 )
	throw CosFitterExcept("cosfitter","calcPreErrs",
			      "Sum of inverse cov matrix entries == 0.0",8);

    }    

  } else {
    //Errors aren't fixed, so only do intrinsic and redshift errors,
    //plus alpha beta terms if they have single values
    //Note that they can't both be single values or we would have
    // already tripped into fparam.errsfixed.  So we don't need to
    // evaluate the alpha/beta cross term

    if (fixedintrinsic) {
      for (unsigned int i = 0; i < nsn; ++i) {
	emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
	dzerrsq = fparam.pecz*fparam.pecz + sne[i].var_z;
	dzerrsq *= zfacsq*emptyfac*emptyfac;
	pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
      }
    } else {
      std::map<unsigned int, double>::iterator pos;
      for (i = 0; i < nsn; ++i) {
	pos = fparam.intrinsicdisp.find( sne[i].dataset );
	if (pos == fparam.intrinsicdisp.end()) {
	  std::stringstream errstr;
	  errstr << "Unknown dataset number " << sne[i].dataset;
	  errstr << " for sn: " << sne[i].name;
	  throw CosFitterExcept("cosfitter","calcPreErrs",
				errstr.str(),4);
	}
	currintrinsicsq = pos->second * pos->second;
	emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
	dzerrsq = fparam.pecz*fparam.pecz + sne[i].var_z;
	dzerrsq *= zfacsq*emptyfac*emptyfac;
	pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
      }
    }

    for (i = 0; i < nsn; ++i) {
      //z error term -- only true at low z, but should be largely
      // negligible for high z Sne anyways
      emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
      dzerrsq = fparam.pecz*fparam.pecz + sne[i].var_z;
      dzerrsq *= zfacsq*emptyfac*emptyfac;
      
      if (! fixedintrinsic ) {
	std::map<unsigned int, double>::iterator pos;
	pos = fparam.intrinsicdisp.find( sne[i].dataset );
	if (pos == fparam.intrinsicdisp.end()) {
	  std::stringstream errstr;
	  errstr << "Unknown dataset number " << sne[i].dataset;
	  errstr << " for sn: " << sne[i].name;
	  throw CosFitterExcept("cosfitter","calcPreErrs",
				errstr.str(),4);
	}
	currintrinsicsq = pos->second * pos->second;
      }

      pre_vars[i] = sne[i].var_mag + dzerrsq + currintrinsicsq;
      
      if ( have_fixed_alpha ) 
	pre_vars[i] += curr_alpha * curr_alpha *
	  sne[i].var_widthpar + 2.0 * curr_alpha *
	  sne[i].cov_mag_widthpar;
      
      if ( have_fixed_beta )
	pre_vars[i] += curr_beta * curr_beta *
	  sne[i].var_colourpar - 2.0 * curr_beta *
	  sne[i].cov_mag_colourpar;
    }
    
  }

}


/*!
  Updates the inverse covariance matrix and fills in the row
  sums and analytic marginalization E parameter.

  \param[in] alpha \f$\alpha\f$
  \param[in] beta \f$\beta\f$

  If the requisite \f$\alpha, \beta\f$ covariance matrix
  isn't available, alpha and/or beta are ignored.  
*/
void cosfitter::computeInvCovMatrix( double alpha, 
				     double beta ) const {

  //Now update the diagonal errors
  //We could check to see if they are non-zero, but in reality the
  // cost of this function is so dominated by the inversion it
  // just doesn't matter
  working_vec.resize( nsn );

  double asq = alpha * alpha;
  double bsq = beta * beta;
  double ab = alpha * beta;

  if ( have_fixed_alpha && have_fixed_beta ) {
    //Then the alpha and beta terms were already included
    // in pre_vars by calcPreErrs
    for (unsigned int i = 0; i < nsn; ++i)
      working_vec[i] = pre_vars[i]; //Needed because they are diff types
  } else if ( have_fixed_alpha ) {
    //The alpha term was included, but the beta and cross term wasn't
    for (unsigned int i = 0; i < nsn; ++i)
      working_vec[i] = pre_vars[i] + bsq * sne[i].var_colourpar 
	- 2.0 * beta  * sne[i].cov_mag_colourpar
	- 2.0 * ab * sne[i].cov_widthpar_colourpar;
  } else if ( have_fixed_beta ) {
    //Likewise, but for beta
    for (unsigned int i = 0; i < nsn; ++i)
      working_vec[i] = pre_vars[i] + asq * sne[i].var_widthpar 
	+ 2.0 * alpha  * sne[i].cov_mag_widthpar
	- 2.0 * ab * sne[i].cov_widthpar_colourpar;
  } else {
    //Neither
    for (unsigned int i = 0; i < nsn; ++i)
      working_vec[i] = pre_vars[i] + asq * sne[i].var_widthpar 
	+ bsq * sne[i].var_colourpar 
	+ 2.0 * alpha * sne[i].cov_mag_widthpar
	- 2.0 * beta  * sne[i].cov_mag_colourpar
	- 2.0 * ab * sne[i].cov_widthpar_colourpar;
  }

  //And invert
  int status;
  status = 0;
  if ( sne.haveFullCovMatrix() ) {
    //std::cerr << "Doing inversion with " << alpha << " " << beta << std::endl;

    //We need to do a full, N^3 inversion
    sne.getCovMatrixSNeRef().getCombinedCovMatrix( invcovmatrix, alpha, beta );
    for (unsigned int i = 0; i < nsn; ++i)
      invcovmatrix[i][i] += working_vec[i];

    //ofstream of("covmatrix.txt");
    //of << invcovmatrix << std::endl;
    //of.close();

    invcovmatrix.invert(status);
    //std::cout << status << std::endl;
  } else {
    //Woodbury form -- should be much faster
    sne.getWoodburyCovMatrixSNeRef().getInvCovMatrix( invcovmatrix, 
						      working_vec,
						      alpha, beta, status);
  }
  if (status) {
#ifdef DEBUG
    //Write out the one we crashed on
    if (sne.haveFullCovMatrix() ) {
      invcovmatrix = sne.getCombinedCovMatrix( alpha, beta );
      for (unsigned int i = 0; i < nsn; ++i)
	invcovmatrix[i][i] += working_vec[i];
      std::cerr << std::endl << 
	"Cosfitter inversion error status: " << status << 
	" for alpha: " << alpha << 
	" beta: " << beta << std::endl;
      std::cerr << "Writing cov matrix as cov_matrix_failed.txt" <<
	std::endl;
      std::ofstream os("cov_matrix_failed.txt");
      os << invcovmatrix;
      os.close();
    } else {
      sne.getWoodburyCovMatrixRefSNe().getCovMatrix( invcovmatrix, 
						     working_vec,
						     alpha, beta, status);
      std::cerr << std::endl << 
	"Cosfitter Woodbury inversion error status: " << status << 
	" for alpha: " << alpha << 
	" beta: " << beta << std::endl;
      std::cerr << "Writing cov matrix as cov_matrix_failed.txt" <<
	std::endl;
      std::ofstream os("cov_matrix_failed.txt");
      os << invcovmatrix;
      os.close();
    }
#endif
    throw CosFitterExcept("cosfitter","computeInvCovMatrix",
			  "Error inverting",1);
  }

  //Get the total and row-sums for later
  amarg_E = invcovmatrix.getRowTotals( invcovmatrix_rowsums );

}

/*!
  Driver for 0D cosmological fits -- i.e., only the nuisance
  parameters are fit.  It is assumed that the user specifies
  at least \f$\Omega_m\f$ and a flat fit.

  \param[out] results_map  Results of fit are added to this structure
  \param[out] maxlikelihood Maximum likelihood of fit
  \param[out] fittag String description of type of fit performed
*/
void cosfitter::cosmofit0D( std::map< param_tags::paramcodes,
			    param_results>& results_map,
			    double& maxlikelihood,
			    std::string& fittag) const {
  fittag = "Cosmological parameters fixed";

  if (fparam.loop_order == fitparam::cosmo_outer) {
    cosmofit0D_diagonal( results_map, maxlikelihood );
  } else if (fparam.loop_order == fitparam::cosmo_inner) {
    cosmofit0D_nondiagonal( results_map, maxlikelihood );
  } else {
    throw CosFitterExcept("cosfitter","cosmofit0D",
			  "Unknown loop ordering",2);
  }
}

/*!
  Driver for 1D cosmological fits.  Currently this only supports
  fixed curvature \f$\Omega_m\f$ fits.

  \param[out] results_map  Results of fit are added to this structure
  \param[out] maxlikelihood Maximum likelihood of fit
  \param[out] fittag String description of type of fit performed
*/
void cosfitter::cosmofit1D( std::map< param_tags::paramcodes,
			    param_results>& results_map,
			    double& maxlikelihood,
			    std::string& fittag) const {
  //This is the easiest case because we only support 
  // fixed curvature, omega_m fits in 1D

  fittag = "Flat Universe, Omega_m with fixed curvature";

  if (fparam.outputfilename == "") {
    std::string errstrng = "An output file must be specified";
    throw CosFitterExcept("cosfitter","cosmofit1D",errstrng,16);
  }

  if (! fparam.fixcurv ) {
    std::string errstr = "All supported 1D fits have fixed curvature";
    throw CosFitterExcept("cosfitter","cosmofit1D",errstr,1);
  }

  //Figure out what type of fit we are doing
  cosfitter::cosmo_fittype cosmo_fit;
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  it = fparam.params.find( param_tags::omegam );
  if (it == fparam.params.end()) {
    std::string errstr = "All supported 1D fits use Omega_m";
    throw CosFitterExcept("cosfitter","cosmofit1D",errstr,1);
  }  
  if (it->second.fit == param_struct::loop) {
    fittag = "Omega_m with fixed curvature";
    cosmo_fit = cosfitter::flat_omegam;
  } else {
    it = fparam.params.find( param_tags::w0 );
    if (it == fparam.params.end() || it->second.fit != param_struct::loop) { 
      std::string errstr = "Can't determine 1D fit type";
      throw CosFitterExcept("cosfitter","cosmofit1D",errstr,1);
    }
    fittag = "w_0 with fixed curvature and fixed Omega_m";
    cosmo_fit = cosfitter::flat_w0;
  }

  if (fparam.loop_order == fitparam::cosmo_outer) {
    cosmofit1D_outer( cosmo_fit, results_map, maxlikelihood );
  } else if (fparam.loop_order == fitparam::cosmo_inner) {
    cosmofit1D_inner( cosmo_fit, results_map, maxlikelihood );
  } else {
    throw CosFitterExcept("cosfitter","cosmofit1D",
			  "Unknown loop ordering",2);
  }
}

/*!
  Driver for 2D cosmological fits.  Currently this only supports
  \f$\Omega_m, \Omega_{DE}\f$ and fixed curvature \f$\Omega_m, w_0\f$ fits.

  \param[out] results_map  Results of fit are added to this structure
  \param[out] maxlikelihood Maximum likelihood of fit
  \param[out] fittag String description of type of fit performed
*/
void cosfitter::cosmofit2D( std::map< param_tags::paramcodes,
			    param_results>& results_map,
			    double& maxlikelihood,
			    std::string& fittag) const {

  if (fparam.outputfilename == "") {
    std::string errstrng = "An output file must be specified";
    throw CosFitterExcept("cosfitter","cosmofit2D",errstrng,16);
  }

  //This one is a bit more complicated.  Determine what type of
  // fit is being performed
  cosfitter::cosmo_fittype cosmo_fit;
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  it = fparam.params.find( param_tags::omegam );
  //Make sure fit type makes sense
  if (it == fparam.params.end() || it->second.fit != param_struct::loop ) {
    std::string errstr = "All supported 2D fit types use Omega_m ";
    throw CosFitterExcept("cosfitter","cosmofit2D",errstr,1);
  }

  it = fparam.params.find( param_tags::omegade );
  if (it != fparam.params.end() ) {
    //We know there were 2 cosmo params with loops, and we already
    // made sure omega_m was present
    if (fparam.fixcurv) 
      throw CosFitterExcept("cosfitter","cosmofit2D",
			    "Have Omega_m and Omega_de in 2D, but fixed curvature",2);
    cosmo_fit = cosfitter::omegam_omegade;
    fittag = "Omega_m, Omega_Lambda fit";
  } else {
    it = fparam.params.find( param_tags::w0 );
    if (it == fparam.params.end() ) {
      it = fparam.params.find( param_tags::wa );
      if (it == fparam.params.end() ) {
	std::string errstr = "Can't determine fit type: have omega_m, ";
	errstr += "no omega_de, no w0";
	throw CosFitterExcept("cosfitter","cosmofit2D",
			      errstr,4);
      } else {
	throw CosFitterExcept("cosfitter","cosmofit2D",
			      "2D fit to Omega_m, wa not supported",8);
      }
    }
    if (!fparam.fixcurv)
      throw CosFitterExcept("cosfitter","cosmofit2D",
			    "Need fixed curvature to 2D fit omega_m, w0",16);
    fittag = "Flat Omega_m, Omega_DE fit with constant equation of state";
    cosmo_fit = cosfitter::flat_omegam_w0;
  }

  if (fparam.loop_order == fitparam::cosmo_outer) {
    cosmofit2D_outer( cosmo_fit, results_map, maxlikelihood );
  } else if (fparam.loop_order == fitparam::cosmo_inner) {
    cosmofit2D_inner( cosmo_fit, results_map, maxlikelihood );
  } else {
    throw CosFitterExcept("cosfitter","cosmofit2D",
			  "Unknown loop ordering",2);
  }

}

/*!
  Driver for 3D cosmological fits.  Currently this only supports
  \f$\Omega_m, \Omega_{DE}, w_0\f$ and 
  fixed curvature \f$\Omega_m, w_0, w_a\f$ fits.

  \param[out] results_map  Results of fit are added to this structure
  \param[out] maxlikelihood Maximum likelihood of fit
  \param[out] fittag String description of type of fit performed
*/
void cosfitter::cosmofit3D( std::map< param_tags::paramcodes,
			    param_results>& results_map,
			    double& maxlikelihood,
			    std::string& fittag) const {

  if (fparam.outputfilename == "") {
    std::string errstrng = "An output file must be specified";
    throw CosFitterExcept("cosfitter","cosmofit3D",errstrng,16);
  }

  //Determine what type of fit is being performed
  cosfitter::cosmo_fittype cosmo_fit;
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  it = fparam.params.find( param_tags::omegam );
  //Make sure fit type makes sense
  if (it == fparam.params.end() || it->second.fit != param_struct::loop ) {
    std::string errstr = "All supported 3D fit types use Omega_m ";
    throw CosFitterExcept("cosfitter","cosmofit3D",errstr,1);
  }

  it = fparam.params.find( param_tags::omegade );
  if (it != fparam.params.end() ) {
    //We know there were 2 cosmo params with loops, and we already
    // made sure omega_m was present
    if (fparam.fixcurv) 
      throw CosFitterExcept("cosfitter","cosmofit3D",
			    "Have Omega_m and Omega_de in 2D, but fixed curvature",2);
    fittag = "Omega_m, Omega_de, w0 fit";
    cosmo_fit = cosfitter::omegam_omegade_w0;
  } else {
    it = fparam.params.find( param_tags::w0 );
    if (it == fparam.params.end() ) {
      std::string errstr = "Can't determine fit type: have omega_m, ";
      errstr += "no omega_de, no w0";
      throw CosFitterExcept("cosfitter","cosmofit3D",
			    errstr,4);
    }
    it = fparam.params.find( param_tags::wa );
    if (it == fparam.params.end() ) {
      std::string errstr = "Can't determine fit type: have omega_m, ";
      errstr += "no omega_de, no wa";
      throw CosFitterExcept("cosfitter","cosmofit3D",
			    errstr,4);
    }
    if (!fparam.fixcurv)
      throw CosFitterExcept("cosfitter","cosmofit3D",
			    "Need fixed curvature to 3D fit omega_m, w0, wa",
			    16);
    fittag = "Flat Omega_m, Omega_de, w0, wa fit";
    if (fparam.usekomatsuform) fittag += ", Komatsu w(a) form";
    cosmo_fit = cosfitter::flat_omegam_w0_wa;
  }


  if (fparam.loop_order == fitparam::cosmo_outer) {
    cosmofit3D_outer( cosmo_fit, results_map, maxlikelihood );
  } else if (fparam.loop_order == fitparam::cosmo_inner) {
    cosmofit3D_inner( cosmo_fit, results_map, maxlikelihood );
  } else {
    throw CosFitterExcept("cosfitter","cosmofit3D",
			  "Unknown loop ordering",2);
  }
}


/*!
  Handles 0D cosmological fit in case where the errors are diagonal

  \param[out] results_map  A map that will be filled with the
                   results of the fit
  \param[out] maxlikelihood The maximum likelihood

  As a side effect, this may write out a cosgrid as
  a file containing the probability surface of the nuisance parameters.
*/
void cosfitter::cosmofit0D_diagonal( std::map< param_tags::paramcodes,
				     param_results >& results_map, 
				     double& maxlikelihood ) const {

  if ( fparam.loop_order != fitparam::cosmo_outer )
    throw CosFitterExcept("cosfitter","cosmofit0D_diagonal",
			  "Can't call this with loop_order set to outer",1);

  std::map< param_tags::paramcodes, param_struct >::const_iterator it;

  //Figuring out om, ode is a bit tricky -- the user may have specified
  // one and fixed curvature, or both and not fixed curvature
  bool haveom, haveode;
  double om, ode;
  om = 0.3; ode = 0.7;  //To stop compiler warning
  it = fparam.params.find( param_tags::omegam );
  haveom = (it != fparam.params.end() );
  if (haveom) om = it->second.fixval;
  it = fparam.params.find( param_tags::omegade );
  haveode = (it != fparam.params.end() );
  if (haveode) ode = it->second.fixval;
  if ( ! ( haveom || haveode ) )
	throw CosFitterExcept("cosfitter","cosmofit0D_diagonal",
			      "Need one of om or ode",2);
  if ( fparam.fixcurv ) {
    const double errtol = 1e-6;
    if ( haveom && haveode ) {
      //Make sure they add up correctly
      if ( fabs( 1.0 - om - ode - fparam.ocurv ) >= errtol ) 
	throw CosFitterExcept("cosfitter","cosmofit0D_diagonal",
			      "om+ode+ocurv != 1.0 for fixed curvature",4);
    }
    if (haveom) ode = 1.0 - fparam.ocurv - om; else
      om = 1.0 - fparam.ocurv - ode;
  } else if ( ! ( haveom && haveode ) ) {
	throw CosFitterExcept("cosfitter","cosmofit0D_diagonal",
			      "Must have om and ode if not fixed curv set",8);
  }

  //Check on w0, wa: if present and fixed use, otherwise use defaults
  double w0, wa;
  it = fparam.params.find( param_tags::w0 );
  if ( it != fparam.params.end() ) w0 = it->second.fixval; else w0 = -1.0;
  it = fparam.params.find( param_tags::wa );
  if ( it != fparam.params.end() ) wa = it->second.fixval; else wa = 0.0;

  param_struct::fittypes alpha_fittype, beta_fittype;
  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  italpha = fparam.params.find( param_tags::alpha );
  if ( italpha == fparam.params.end() ) 
    alpha_fittype = param_struct::not_set; else
    alpha_fittype = italpha->second.fit;
  itbeta = fparam.params.find( param_tags::beta );
  if ( itbeta == fparam.params.end() ) 
    beta_fittype = param_struct::not_set; else
    beta_fittype = itbeta->second.fit;

  //Get lumdist
  int st;
  dl.resize(nsn);
  st = lm.getLumDist( sne, dl, om, ode, w0, wa );
  if (st != 0)
    throw CosFitterExcept("cosfitter","cosmofit0D_diagonal",
			  "Can't get distance for fixed cosmological params",2);
  if ( alpha_fittype == param_struct::loop && 
       beta_fittype == param_struct::loop ) {
    //Loop on alpha and beta
    cosgrid2D nuisance( italpha->second, itbeta->second );
    nuisance = 0.0;

    calcLikelihood_lab( nuisance );
    maxlikelihood = nuisance.getTotal();
    nuisance.normalize();

    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);

  } else if ( alpha_fittype == param_struct::loop ) {
    //Loop on alpha, analytic
    // scriptm and either no beta or a fixed value
    cosgrid1D nuisance( italpha->second );
    nuisance = 0.0;
    calcLikelihood_la( nuisance );
    maxlikelihood = nuisance.getTotal();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout) {
      //Well, write alpha, although we won't have beta
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
    }
  } else if ( beta_fittype == param_struct::loop ) {
    //Loop on beta, analytic
    // scriptm and either no alpha or a fixed value
    cosgrid1D nuisance( itbeta->second );
    nuisance = 0.0;
    calcLikelihood_lb( nuisance );
    maxlikelihood = nuisance.getTotal();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout) {
      //Well, write alpha, although we won't have beta
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
    }
  } else {
    //Analytic scriptm, and
    // some combination of no alpha or beta or fixed values
    //Not sure what we are doing here, but whatever.  Nothing
    // to fit, after all
    maxlikelihood = calcLikelihood();
  }

}

/*!
  Handles 0D cosmological fit (i.e., cosmo parameters fixed)
  with non-diagonal errors

  \param[out] results_map  A map that will be filled with the
                   results of the fit
  \param[out] maxlikelihood The maximum likelihood

  As a side effect, this will write out a cosgrid as
  a file containing the probability surface.		   
*/
void cosfitter::cosmofit0D_nondiagonal( std::map< param_tags::paramcodes,
				     param_results >& results_map, 
				     double& maxlikelihood ) const {
  if ( sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","cosmofit0D_nondiagonal",
			  "Can't call this with diagonal errors",1);
  
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;

  //Figuring out om, ode is a bit tricky -- the user may have specified
  // one and fixed curvature, or both and not fixed curvature
  bool haveom, haveode;
  double om, ode;
  om = 0.25; ode = 0.75; //To shut up the compiler warnging
  it = fparam.params.find( param_tags::omegam );
  haveom = (it != fparam.params.end() );
  if (haveom) om = it->second.fixval;
  it = fparam.params.find( param_tags::omegade );
  haveode = (it != fparam.params.end() );
  if (haveode) ode = it->second.fixval;
  if ( ! ( haveom || haveode ) )
	throw CosFitterExcept("cosfitter","cosmofit0D_diagonal",
			      "Need one of om or ode",2);
  if ( fparam.fixcurv ) {
    const double errtol = 1e-6;
    if ( haveom && haveode ) {
      //Make sure they add up correctly
      if ( 1.0 - fabs( om + ode + fparam.ocurv ) >= errtol ) 
	throw CosFitterExcept("cosfitter","cosmofit0D_diagonal",
			      "om+ode+ocurv != 1.0 for fixed curvature",4);
    }
    if (haveom) ode = 1.0 - fparam.ocurv - om; else
      om = 1.0 - fparam.ocurv - ode;
  } else if ( ! ( haveom && haveode ) ) {
	throw CosFitterExcept("cosfitter","cosmofit0D_diagonal",
			      "Must have om and ode if not fixed curv set",8);
  }
  
  //Check on w0, wa: if present and fixed use, otherwise use defaults
  double w0, wa;
  it = fparam.params.find( param_tags::w0 );
  if ( it != fparam.params.end() ) w0 = it->second.fixval; else w0 = -1.0;
  it = fparam.params.find( param_tags::wa );
  if ( it != fparam.params.end() ) wa = it->second.fixval; else wa = 0.0;

  //Calculate the luminosity distances
  lumdist_array0D lar;
  lar.loadData( om, ode, w0, wa, sne, fparam.verbose );

  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  param_struct::fittypes alpha_fittype, beta_fittype;
  italpha = fparam.params.find( param_tags::alpha );
  if ( italpha == fparam.params.end() ) 
    alpha_fittype = param_struct::not_set; else
    alpha_fittype = italpha->second.fit;
  itbeta = fparam.params.find( param_tags::beta );
  if ( itbeta == fparam.params.end() ) 
    beta_fittype = param_struct::not_set; else
    beta_fittype = itbeta->second.fit;

  //Allocate working arrays
  diffarr = new double[ nsn ];
  prod = new double[ nsn ];

  double curr_alpha, curr_beta;
  if ( alpha_fittype == param_struct::loop && 
       beta_fittype == param_struct::loop ) {
    //Loop on alpha/beta
    cosgrid2D nuisance( italpha->second, itbeta->second );
    unsigned int nsteps = nuisance.getAxisN(0);
    if (fparam.verbose) std::cout << "Calculating Likelihoods" << std::endl;
    for (unsigned int i = 0; i < nsteps; ++i) {
      curr_alpha = nuisance.getAxisVal(0,i);
      for (unsigned int j = 0; j < nuisance.getAxisN(1); ++j) {
	curr_beta = nuisance.getAxisVal(1,j);
	nuisance[i][j] = calcLikelihood_nondia( curr_alpha, curr_beta, lar );
      }
    }
    unsigned int i1,i2;
    maxlikelihood = nuisance.getMaximum(i1,i2);

    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);

  } else if ( alpha_fittype == param_struct::loop ) {
    //Loop on alpha
    cosgrid1D nuisance( italpha->second );
    if ( beta_fittype == param_struct::not_set ) {
      if ( fparam.colour_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit0D_outer",
			      "Have Colour cov matrix, but not treating",8);
      curr_beta = 0.0;
    } else if ( beta_fittype == param_struct::fixed ) {
      curr_beta = itbeta->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit0D_outer",
				 "Unknown Beta fit type",16);
    unsigned int nsteps = nuisance.getAxisN();
    if (fparam.verbose) std::cout << "Calculating Likelihoods" << std::endl;
    for (unsigned int i = 0; i < nsteps; ++i) {
      curr_alpha = nuisance.getAxisVal(i);
      nuisance[i] = calcLikelihood_nondia( curr_alpha, curr_beta, lar );
    }
    unsigned int i1;
    maxlikelihood = nuisance.getMaximum(i1);
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
  } else if ( beta_fittype == param_struct::loop ) {
    //Loop on beta
    cosgrid1D nuisance( itbeta->second );
    if ( alpha_fittype == param_struct::not_set ) {
      if ( fparam.width_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit0D_outer",
			      "Have Alpha cov matrix, but not treating",8);
      curr_alpha = 0.0;
    } else if ( alpha_fittype == param_struct::fixed ) {
      curr_alpha = italpha->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit0D_outer",
				 "Unknown Alpha fit type",16);
    
    unsigned int nsteps = nuisance.getAxisN();
    if (fparam.verbose) std::cout << "Calculating Likelihoods" << std::endl;
    for (unsigned int i = 0; i < nsteps; ++i) {
      curr_beta = nuisance.getAxisVal(i);
      nuisance[i] = calcLikelihood_nondia( curr_alpha, curr_beta, lar );
    }
    unsigned int i1;
    maxlikelihood = nuisance.getMaximum(i1);
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
  } else {
    //Alpha and beta are both fixed or missing
    //This makes things pretty simple
    if ( alpha_fittype == param_struct::not_set ) {
      if ( fparam.width_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit0D_outer",
			      "Have Alpha cov matrix, but not treating",8);
      curr_alpha = 0.0;
    } else if ( alpha_fittype == param_struct::fixed ) {
      curr_alpha = italpha->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit0D_outer",
				 "Unknown Alpha fit type",16);
    if ( beta_fittype == param_struct::not_set ) {
      if ( fparam.colour_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit0D_outer",
			      "Have Beta cov matrix, but not treating",8);
      curr_beta = 0.0;
    } else if ( beta_fittype == param_struct::fixed ) {
      curr_beta = itbeta->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit0D_outer",
				 "Unknown Beta fit type",16);
    maxlikelihood = calcLikelihood_nondia( curr_alpha, curr_beta, lar );
  }
  lar.free();

  delete[] prod;
  delete[] diffarr;
  prod = diffarr = NULL;

}

/*!
  Handles 1D cosmological fit in case where the outer loop is
  the cosmological parameters (i.e., when the luminosity distance
  is the most expensive thing to calcluate).

  \param[in] cosmo_fit Type of fit to perform
  \param[out] results_map  A map that will be filled with the
                   results of the fit
  \param[out] maxlikelihood The maximum likelihood

  As a side effect, this will write out a cosgrid as
  a file containing the probability surface.		   
*/
void cosfitter::cosmofit1D_outer(cosfitter::cosmo_fittype cosmo_fit,  
				 std::map< param_tags::paramcodes,
				 param_results >& results_map, 
				 double& maxlikelihood ) const {

  if ( fparam.loop_order != fitparam::cosmo_outer )
    throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			  "Can't call this with loop_order set to outer",1);

  //Holds current params in order om, ode, w0, wa, with defaults
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

  //Set fixed values
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for (it = fparam.params.begin(); it != fparam.params.end(); ++it) 
    if ( it->second.fit == param_struct::fixed )
      utility::updateStandardCosVec( it->first, fparam, current_cosparams, 
				     it->second.fixval );

  //Next point it at the loop variable
  if ( cosmo_fit == cosfitter::flat_omegam )
    it = fparam.params.find( param_tags::omegam ); 
  else if (cosmo_fit == cosfitter::flat_w0 )
    it = fparam.params.find( param_tags::w0 ); 
  else 
    throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			  "Unsupported 1D fit type",1);
  if ( it == fparam.params.end() )
    throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			  "Can't find loop variable",2);

  cosgrid1D margprob1D( it->second );
  margprob1D = 0.0;

  param_struct::fittypes alpha_fittype, beta_fittype;
  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  italpha = fparam.params.find( param_tags::alpha );
  if ( italpha == fparam.params.end() ) 
    alpha_fittype = param_struct::not_set; else
    alpha_fittype = italpha->second.fit;
  itbeta = fparam.params.find( param_tags::beta );
  if ( itbeta == fparam.params.end() ) 
    beta_fittype = param_struct::not_set; else
    beta_fittype = itbeta->second.fit;

  int st;
  double cos0;
  dl.resize(nsn);
  param_tags::paramcodes code0 = margprob1D.getAxisSpec().first;
  if ( alpha_fittype == param_struct::loop && 
       beta_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on alpha/beta and analytic scriptm
    cosgrid2D nuisance( italpha->second, itbeta->second );
    nuisance = 0.0;
    cosgrid2D working( italpha->second, itbeta->second );
    unsigned int nsteps = margprob1D.getAxisN();
    unsigned int progbar_stepsize = nsteps / 10;
    std::stringstream progbar;
    progbar << "Progress 0%";
    fprintf(stderr,"%s",progbar.str().c_str());
    for (unsigned int i=0; i < margprob1D.getAxisN(); ++i) {
      if ( i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      cos0 = margprob1D.getAxisVal( i );
      utility::updateStandardCosVec(code0,fparam,current_cosparams,cos0);
      st = lm.getLumDist( sne, dl, current_cosparams);
      if (st == 0) {
	calcLikelihood_lab( working );
	nuisance += working; //Update alpha/beta
	margprob1D[i] = working.getTotal();
      } else margprob1D[i] = 0.0;  //Bad lumdist
    }
    fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);

  } else if ( alpha_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on alpha, analytic
    // scriptm and either no beta or a fixed value
    cosgrid1D nuisance( italpha->second );
    nuisance = 0.0;
    cosgrid1D working( italpha->second );
    unsigned int nsteps = margprob1D.getAxisN();
    unsigned int progbar_stepsize = nsteps / 10;
    std::stringstream progbar;
    progbar << "Progress 0%";
    fprintf(stderr,"%s",progbar.str().c_str());
    for (unsigned int i=0; i < nsteps; ++i) {
      if ( i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      cos0 = margprob1D.getAxisVal( i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      st = lm.getLumDist( sne, dl, current_cosparams );
      if (st == 0) {
	calcLikelihood_la( working );
	margprob1D[i] = working.getTotal();
	nuisance += working; //Update alpha
      } else margprob1D[i] = 0.0;  //Bad lumdist
    }
    fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout) {
      //Well, write alpha, although we won't have beta
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
    }
  } else if ( beta_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on beta, analytic
    // scriptm and either no alpha or a fixed value
    cosgrid1D nuisance( itbeta->second );
    nuisance = 0.0;
    cosgrid1D working( itbeta->second );
    unsigned int nsteps = margprob1D.getAxisN();
    unsigned int progbar_stepsize = nsteps / 10;
    std::stringstream progbar;
    progbar << "Progress 0%";
    fprintf(stderr,"%s",progbar.str().c_str());
    for (unsigned int i=0; i < margprob1D.getAxisN(); ++i) {
      if ( i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      cos0 = margprob1D.getAxisVal( i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      st = lm.getLumDist( sne, dl, current_cosparams );
      if (st == 0) {
	calcLikelihood_lb( working );
	margprob1D[i] = working.getTotal();
	nuisance += working; //Update beta
      } else margprob1D[i] = 0.0;  //Bad lumdist
    }
    fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout) {
      //Well, write alpha, although we won't have beta
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
    }
  } else {
    //Omega_m with fixed curvature, analytic scriptm, and
    // some combination of no alpha or beta or fixed values
    for (unsigned int i=0; i < margprob1D.getAxisN(); ++i) {
      cos0 = margprob1D.getAxisVal( i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      st = lm.getLumDist( sne, dl, current_cosparams );
      if (st == 0) {
	margprob1D[i] = calcLikelihood();
      } else margprob1D[i] = 0.0;  //Bad lumdist
    }
  }
  unsigned int iom;
  maxlikelihood = margprob1D.getMaximum(iom);

  //Get info about loop variable
  margprob1D.normalize();
  GetParamMap( margprob1D, results_map );

  //Write probability distribution
  margprob1D.writeFile( fparam.outputfilename, fparam.binaryout );

}

/*!
  Handles 1D cosmological fit in case where the inner loop is
  the cosmological parameters (i.e., when the covariance matrix
  inversion is the most expensive thing to calcluate).

  \param[in] cosmo_fit Type of fit to perform
  \param[out] results_map  A map that will be filled with the
                   results of the fit
  \param[out] maxlikelihood The maximum likelihood

  As a side effect, this will write out a cosgrid as
  a file containing the probability surface.		   
*/
//This one is relatively simple because the only fit we support
// for one-D is fixed curvature, omega_m
void cosfitter::cosmofit1D_inner( cosfitter::cosmo_fittype cosmo_fit, 
				  std::map< param_tags::paramcodes,
				  param_results >& results_map, 
				  double& maxlikelihood ) const {
  if ( sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","cosmofit1D_inner",
			  "Can't call this with diagonal errors",1);
  
  if ( fparam.loop_order != fitparam::cosmo_inner )
    throw CosFitterExcept("cosfitter","cosmofit1D_inner",
			  "Can't call this with loop_order set to inner",2);

  

  //Point it at our cosmological parameter
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  if ( cosmo_fit == cosfitter::flat_omegam )
    it = fparam.params.find( param_tags::omegam ); 
  else if (cosmo_fit == cosfitter::flat_w0 )
    it = fparam.params.find( param_tags::w0 ); 
  else 
    throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			  "Unsupported 1D fit type",1);
  if ( it == fparam.params.end() )
    throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			  "Can't find loop variable",2);
  cosgrid1D margprob1D( it->second );
  margprob1D = 0.0;
  
  //Calculate the luminosity distances
  lumdist_array1D lar;
  lar.loadData( margprob1D, sne, fparam, fparam.showprogbar );

  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  param_struct::fittypes alpha_fittype, beta_fittype;
  italpha = fparam.params.find( param_tags::alpha );
  if ( italpha == fparam.params.end() ) 
    alpha_fittype = param_struct::not_set; else
    alpha_fittype = italpha->second.fit;
  itbeta = fparam.params.find( param_tags::beta );
  if ( itbeta == fparam.params.end() ) 
    beta_fittype = param_struct::not_set; else
    beta_fittype = itbeta->second.fit;

  //Allocate working arrays
  unsigned int n0 = margprob1D.getAxisN();
  diffarr = new double[ nsn * n0 ];
  prod = new double[ nsn * n0 ];

#ifdef TIMING
  clock_t starttime;
#endif
  double curr_alpha, curr_beta;
  if ( alpha_fittype == param_struct::loop && 
       beta_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on alpha/beta and analytic scriptm
    cosgrid1D working( it->second );
    cosgrid2D nuisance( italpha->second, itbeta->second );
    unsigned int nsteps = nuisance.getAxisN(0);
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.verbose) std::cout << "Calculating Likelihoods" << std::endl;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
#ifdef TIMING
    tot_mult_time = tot_inv_time = tot_nondia_time = 0;
#endif
    for (unsigned int i = 0; i < nsteps; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_alpha = nuisance.getAxisVal(0,i);
      for (unsigned int j = 0; j < nuisance.getAxisN(1); ++j) {
	curr_beta = nuisance.getAxisVal(1,j);
#ifdef TIMING
	starttime = clock();
#endif
	calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
#ifdef TIMING
	tot_nondia_time += clock() - starttime;
#endif
	margprob1D += working;
	nuisance[i][j] = working.getTotal();
      }
    }

    if (fparam.showprogbar) fprintf(stderr,"\n");

#ifdef TIMING
    std::cout << "Total time inverting: " << tot_inv_time << " ticks= "
	      << tot_inv_time / static_cast<double>(CLOCKS_PER_SEC) 
	      << " sec" << std::endl;
    std::cout << "Total time doing post-inversion multiply: " 
	      << tot_mult_time << " ticks= "
	      << tot_mult_time / static_cast<double>(CLOCKS_PER_SEC) 
	      << " sec" << std::endl;
    std::cout << "Total time likelihooding: " << tot_nondia_time << " ticks= "
	      << tot_nondia_time / static_cast<double>(CLOCKS_PER_SEC) 
	      << " sec" << std::endl;
#endif
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);

  } else if ( alpha_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on alpha, analytic
    // scriptm and either no beta or a fixed value
    cosgrid1D working( it->second );
    cosgrid1D nuisance( italpha->second );
    if ( beta_fittype == param_struct::not_set ) {
      if ( fparam.colour_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			      "Have Colour cov matrix, but not treating",8);
      curr_beta = 0.0;
    } else if ( beta_fittype == param_struct::fixed ) {
      curr_beta = itbeta->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit1D_outer",
				 "Unknown Beta fit type",16);
    unsigned int nsteps = nuisance.getAxisN();
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.verbose) std::cout << "Calculating Likelihoods" << std::endl;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i = 0; i < nuisance.getAxisN(); ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_alpha = nuisance.getAxisVal(i);
      calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
      margprob1D += working;
      nuisance[i] = working.getTotal();
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
  } else if ( beta_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on beta, analytic
    // scriptm and either no alpha or a fixed value
    cosgrid1D working( it->second );
    cosgrid1D nuisance( itbeta->second );
    if ( alpha_fittype == param_struct::not_set ) {
      if ( fparam.width_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			      "Have Alpha cov matrix, but not treating",8);
      curr_alpha = 0.0;
    } else if ( alpha_fittype == param_struct::fixed ) {
      curr_alpha = italpha->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit1D_outer",
				 "Unknown Alpha fit type",16);
    
    unsigned int nsteps = nuisance.getAxisN();
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.verbose) std::cout << "Calculating Likelihoods" << std::endl;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i = 0; i < nsteps; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_beta = nuisance.getAxisVal(i);
      calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
      margprob1D += working;
      nuisance[i] = working.getTotal();
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
  } else {
    //Alpha and beta are both fixed or missing
    //This makes things pretty simple
    if ( alpha_fittype == param_struct::not_set ) {
      if ( fparam.width_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			      "Have Alpha cov matrix, but not treating",8);
      curr_alpha = 0.0;
    } else if ( alpha_fittype == param_struct::fixed ) {
      curr_alpha = italpha->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit1D_outer",
				 "Unknown Alpha fit type",16);
    if ( beta_fittype == param_struct::not_set ) {
      if ( fparam.colour_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit1D_outer",
			      "Have Beta cov matrix, but not treating",8);
      curr_beta = 0.0;
    } else if ( beta_fittype == param_struct::fixed ) {
      curr_beta = itbeta->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit1D_outer",
				 "Unknown Beta fit type",16);
    calcLikelihood_nondia( margprob1D, curr_alpha, curr_beta, lar );
  }
  lar.free();

  delete[] prod;
  delete[] diffarr;
  prod = diffarr = NULL;

  unsigned int iom;
  maxlikelihood = margprob1D.getMaximum(iom);

  //Get info about Omega_m
  margprob1D.normalize();
  GetParamMap( margprob1D, results_map );

  //Write probability distribution
  margprob1D.writeFile( fparam.outputfilename, fparam.binaryout );

}


/*!
  Handles 2D cosmological fit in case where the outer loop is
  the cosmological parameters (i.e., when the luminosity distance
  is the most expensive thing to calcluate).  Only two types of
  fits are currently supported: omega_m/omega_de and flat omega_m/ w0

  \param[in] cosmo_fit  The type of fit to do
  \param[out] results_map  A map that will be filled with the
                   results of the fit
  \param[out] maxlikelihood The maximum likelihood

  As a side effect, this will write out a cosgrid as
  a file containing the probability surface.		   
*/
void cosfitter::cosmofit2D_outer( cosfitter::cosmo_fittype cosmo_fit,
				  std::map< param_tags::paramcodes,
				  param_results >& results_map, 
				  double& maxlikelihood ) const {
  if ( ! sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			  "Errors must be diagonal",1);
  if ( fparam.loop_order != fitparam::cosmo_outer )
    throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			  "Can't call this with loop_order set to outer",2);

  std::map< param_tags::paramcodes, param_struct >::const_iterator it1, it2;

  //Holds current params in order om, ode, w0, wa, with defaults
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

  //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for (it = fparam.params.begin(); it != fparam.params.end(); ++it) 
    if ( it->second.fit == param_struct::fixed )
      utility::updateStandardCosVec( it->first, fparam, current_cosparams, 
				     it->second.fixval );

  //Now point it1, it2 at the vars we will loop over
  it1 = fparam.params.find( param_tags::omegam );
  if ( it1 == fparam.params.end() )
    throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			  "Can't find Omega_m--all supported 2D fits have it",
			  16);
  if ( it1->second.fit != param_struct::loop )
    throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			  "Omega_m must be loop fit type",32);

  if ( cosmo_fit == cosfitter::omegam_omegade ) {
    it2 = fparam.params.find( param_tags::omegade );
    if (it2 == fparam.params.end() )
      throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			    "Can't find Omega_de despite specified fit type",
			    64);
    if (it2->second.fit != param_struct::loop )
      throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			    "Omega_DE must have loop fit type for om/ode fit",
			    128);
  } else if ( cosmo_fit == cosfitter::flat_omegam_w0 ) {
    it2 = fparam.params.find( param_tags::w0 );
    if (it2 == fparam.params.end() )
      throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			    "Can't find w0 despite specified fit type",
			    32);
    if (it2->second.fit != param_struct::loop )
      throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			    "w0 must have loop fit type for flat om/w fit",
			    64);
  } else {
      throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			    "Unknown fit type",8);
  }

  cosgrid2D margprob2D( it1->second, it2->second );
  margprob2D = 0.0;

  param_struct::fittypes alpha_fittype, beta_fittype;
  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  italpha = fparam.params.find( param_tags::alpha );
  if ( italpha == fparam.params.end() ) 
    alpha_fittype = param_struct::not_set; else
    alpha_fittype = italpha->second.fit;
  itbeta = fparam.params.find( param_tags::beta );
  if ( itbeta == fparam.params.end() ) 
    beta_fittype = param_struct::not_set; else
    beta_fittype = itbeta->second.fit;

  int st;
  double cos0, cos1;
  param_tags::paramcodes code0 = margprob2D.getAxisSpec(0).first;
  param_tags::paramcodes code1 = margprob2D.getAxisSpec(1).first;
  unsigned int n0 = margprob2D.getAxisN(0);
  unsigned int n1 = margprob2D.getAxisN(1);
  dl.resize(nsn);
  if (fparam.verbose) std::cout << "Doing likelihood loop" << std::endl;
  if ( alpha_fittype == param_struct::loop && 
       beta_fittype == param_struct::loop ) {
    cosgrid2D nuisance( italpha->second, itbeta->second );
    nuisance = 0.0;
    cosgrid2D working( italpha->second, itbeta->second );
    unsigned int progbar_stepsize = n0 / 10;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i=0; i < n0; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / n0)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      cos0 = margprob2D.getAxisVal( 0, i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      for (unsigned int j=0; j < n1; ++j) {
	cos1 = margprob2D.getAxisVal( 1, j );
	utility::updateStandardCosVec( code1, fparam, current_cosparams, 
				       cos1 );
	st = lm.getLumDist(sne,dl,current_cosparams);
	if (st == 0) {
	  calcLikelihood_lab( working );
	  nuisance += working; //Update alpha/beta
	  margprob2D[i][j] = working.getTotal();
	} else margprob2D[i][j] = 0.0;  //Bad lumdist
      }
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);

  } else if ( alpha_fittype == param_struct::loop ) {
    //loop on alpha, analytic scriptm and either no beta or a fixed value
    cosgrid1D nuisance( italpha->second );
    nuisance = 0.0;
    cosgrid1D working( italpha->second );
    for (unsigned int i=0; i < n0; ++i) {
      cos0 = margprob2D.getAxisVal( 0, i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      for (unsigned int j=0; j < n1; ++j) {
	cos1 = margprob2D.getAxisVal( 1, j );
	utility::updateStandardCosVec( code1, fparam, 
				       current_cosparams, cos1 );
	st = lm.getLumDist(sne,dl,current_cosparams);
	if (st == 0) {
	  calcLikelihood_la( working );
	  margprob2D[i][j] = working.getTotal();
	  nuisance += working; //Update alpha
	} else margprob2D[i][j] = 0.0;  //Bad lumdist
      }
    }
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout) {
      //Well, write alpha, although we won't have beta
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
    }
  } else if ( beta_fittype == param_struct::loop ) {
    //loop on beta, analytic scriptm and either no alpha or a fixed value
    cosgrid1D nuisance( itbeta->second );
    nuisance = 0.0;
    cosgrid1D working( itbeta->second );
    for (unsigned int i=0; i < n0; ++i) {
      cos0 = margprob2D.getAxisVal( 0, i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      for (unsigned int j=0; j < n1; ++j) {
	cos1 = margprob2D.getAxisVal( 1, j );
	utility::updateStandardCosVec( code1, fparam, 
				       current_cosparams, cos1 );
	st = lm.getLumDist(sne,dl,current_cosparams);
	if (st == 0) {
	  calcLikelihood_la( working );
	  margprob2D[i][j] = working.getTotal();
	  nuisance += working; //Update alpha
	} else margprob2D[i][j] = 0.0;  //Bad lumdist
      }
    }
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout) {
      //Well, write alpha, although we won't have beta
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
    }
  } else {
    //analytic scriptm, and some combination of no alpha or beta or 
    // fixed values
    for (unsigned int i=0; i < n0; ++i) {
      cos0 = margprob2D.getAxisVal( 0, i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      for (unsigned int j=0; j < n1; ++j) {
	cos1 = margprob2D.getAxisVal( 1, j );
	utility::updateStandardCosVec( code1, fparam, current_cosparams, 
				       cos1 );
	st = lm.getLumDist(sne,dl,current_cosparams);
	if (st == 0) {
	  margprob2D[i][j] = calcLikelihood();
	} else margprob2D[i][j] = 0.0;  //Bad lumdist
      }
    }
  }

  unsigned int i1,i2;
  maxlikelihood = margprob2D.getMaximum(i1,i2);

  //Get info about cosmological params
  margprob2D.normalize();
  GetParamMap( margprob2D, results_map );

  //Write probability distribution
  margprob2D.writeFile( fparam.outputfilename, fparam.binaryout );

}

/*!
  Handles 2D cosmological fit in case where the inner loop is
  the cosmoogical parameters (i.e., when the covariance matrix
  inversion is the most expensive thing to calcluate).  Only two types of
  fits are currently supported: omega_m/omega_de and flat omega_m/ w0

  \param[in] cosmo_fit  The type of fit to do
  \param[out] results_map  A map that will be filled with the
                   results of the fit
  \param[out] maxlikelihood The maximum likelihood

  As a side effect, this will write out a cosgrid as
  a file containing the probability surface.		   
*/
void cosfitter::cosmofit2D_inner( cosfitter::cosmo_fittype cosmo_fit,
				  std::map< param_tags::paramcodes,
				  param_results >& results_map, 
				  double& maxlikelihood ) const {
  if ( sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			  "Can't call this with diagonal errors",1);
  
  if ( fparam.loop_order != fitparam::cosmo_inner )
    throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			  "Can't call this with loop_order set to inner",2);

  std::map< param_tags::paramcodes, param_struct >::const_iterator it1, it2;
  //Point it1, it2 at desired params
  if ( cosmo_fit == cosfitter::omegam_omegade ) {
    it1 = fparam.params.find( param_tags::omegam );
    if (it1 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			    "Can't find Omega_m despite fit type", 4 );
    it2 = fparam.params.find( param_tags::omegade );
    if (it2 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			    "Can't find Omega_DE despite fit type", 4 );
  } else if ( cosmo_fit == cosfitter::flat_omegam_w0 ) {
    it1 = fparam.params.find( param_tags::omegam );
    if (it1 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			    "Can't find Omega_m despite fit type", 4 );
    it2 = fparam.params.find( param_tags::w0 );
    if (it2 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			    "Can't find w0 despite fit type", 4 );
  } else {
    throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			  "Unsupported fit type",8);

  }
  cosgrid2D margprob2D( it1->second, it2->second );
  margprob2D = 0.0;
  
  //Calculate the luminosity distances
  lumdist_array2D lar;
  lar.loadData( margprob2D, sne, fparam, fparam.showprogbar );

  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  param_struct::fittypes alpha_fittype, beta_fittype;
  italpha = fparam.params.find( param_tags::alpha );
  if ( italpha == fparam.params.end() ) 
    alpha_fittype = param_struct::not_set; else
    alpha_fittype = italpha->second.fit;
  itbeta = fparam.params.find( param_tags::beta );
  if ( itbeta == fparam.params.end() ) 
    beta_fittype = param_struct::not_set; else
    beta_fittype = itbeta->second.fit;

  //Allocate working arrays
  unsigned int n1 = margprob2D.getAxisN(1);
  diffarr = new double[ nsn * n1 ];
  prod = new double[ nsn * n1 ];

#ifdef TIMING
  clock_t starttime;
#endif

  double curr_alpha, curr_beta;
  if ( alpha_fittype == param_struct::loop && 
       beta_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on alpha/beta and analytic scriptm
    cosgrid2D working( it1->second, it2->second );
    cosgrid2D nuisance( italpha->second, itbeta->second );
    unsigned int nsteps = nuisance.getAxisN(0);
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.verbose) std::cout << "Calculating Likelihoods" << std::endl;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
#ifdef TIMING
    tot_mult_time = tot_inv_time = tot_nondia_time = 0;
#endif
    for (unsigned int i = 0; i < nsteps; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_alpha = nuisance.getAxisVal(0,i);
      for (unsigned int j = 0; j < nuisance.getAxisN(1); ++j) {
	curr_beta = nuisance.getAxisVal(1,j);
#ifdef TIMING
	starttime = clock();
#endif
	calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
#ifdef TIMING
	tot_nondia_time += clock() - starttime;
#endif
	margprob2D += working;
	nuisance[i][j] = working.getTotal();
      }
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
#ifdef TIMING
    std::cout << "Total time inverting: " << tot_inv_time << " ticks= "
	      << tot_inv_time / static_cast<double>(CLOCKS_PER_SEC) 
	      << " sec" << std::endl;
    std::cout << "Total time doing post-inversion multiply: " 
	      << tot_mult_time << " ticks= "
	      << tot_mult_time / static_cast<double>(CLOCKS_PER_SEC) 
	      << " sec" << std::endl;
    std::cout << "Total time likelihooding: " << tot_nondia_time << " ticks= "
	      << tot_nondia_time / static_cast<double>(CLOCKS_PER_SEC) 
	      << " sec" << std::endl;
#endif
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);

  } else if ( alpha_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on alpha, analytic
    // scriptm and either no beta or a fixed value
    cosgrid2D working( it1->second, it2->second );
    cosgrid1D nuisance( italpha->second );
    if ( beta_fittype == param_struct::not_set ) {
      if ( fparam.colour_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			      "Have Beta cov matrix, but not treating",8);
      curr_beta = 0.0;
    } else if ( beta_fittype == param_struct::fixed ) {
      curr_beta = itbeta->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit2D_inner",
				 "Unknown Beta fit type",16);
    unsigned int nsteps = nuisance.getAxisN();
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i = 0; i < nsteps; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_alpha = nuisance.getAxisVal(i);
      calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
      margprob2D += working;
      nuisance[i] = working.getTotal();
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
  } else if ( beta_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on beta, analytic
    // scriptm and either no alpha or a fixed value
    cosgrid2D working( it1->second, it2->second );
    cosgrid1D nuisance( itbeta->second );
    if ( alpha_fittype == param_struct::not_set ) {
      if ( fparam.width_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			      "Have Alpha cov matrix, but not treating",8);
      curr_alpha = 0.0;
    } else if ( alpha_fittype == param_struct::fixed ) {
      curr_alpha = italpha->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit2D_inner",
				 "Unknown Alpha fit type",16);
    unsigned int nsteps = nuisance.getAxisN();
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i = 0; i < nsteps; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_beta = nuisance.getAxisVal(i);
      calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
      margprob2D += working;
      nuisance[i] = working.getTotal();
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
  } else {
    //Alpha and beta are both fixed or missing
    //This makes things pretty simple
    if ( alpha_fittype == param_struct::not_set ) {
      if ( fparam.width_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			      "Have Alpha cov matrix, but not treating",8);
      curr_alpha = 0.0;
    } else if ( alpha_fittype == param_struct::fixed ) {
      curr_alpha = italpha->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit2D_outer",
				 "Unknown Alpha fit type",16);
    if ( beta_fittype == param_struct::not_set ) {
      if ( fparam.colour_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit2D_inner",
			      "Have Beta cov matrix, but not treating",8);
      curr_beta = 0.0;
    } else if ( beta_fittype == param_struct::fixed ) {
      curr_beta = itbeta->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit2D_inner",
				 "Unknown Beta fit type",16);
    calcLikelihood_nondia( margprob2D, curr_alpha, curr_beta, lar );
  }
  lar.free();

  delete[] prod;
  delete[] diffarr;
  prod = diffarr = NULL;

  unsigned int i1,i2;
  maxlikelihood = margprob2D.getMaximum(i1,i2);

  //Get info about the cosmological parameters
  margprob2D.normalize();
  GetParamMap( margprob2D, results_map );

  //Write probability distribution
  margprob2D.writeFile( fparam.outputfilename, fparam.binaryout );

}

/*!
  Handles 3D cosmological fit in case where the outer loop is
  the cosmoogical parameters (i.e., when the luminosity distance
  is the most expensive thing to calcluate).  Only two types of
  fits are currently supported: omega_m/omega_de/w_0 and 
  flat omega_m/w_0/w_a

  \param[in] cosmo_fit  The type of fit to do
  \param[out] results_map  A map that will be filled with the
                   results of the fit
  \param[out] maxlikelihood The maximum likelihood

  As a side effect, this will write out a cosgrid as
  a file containing the probability surface.		   
*/
void cosfitter::cosmofit3D_outer( cosfitter::cosmo_fittype cosmo_fit,
				  std::map< param_tags::paramcodes,
				  param_results >& results_map, 
				  double& maxlikelihood ) const {
  if ( ! sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			  "Errors must be diagonal",1);
  if ( fparam.loop_order != fitparam::cosmo_outer )
    throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			  "Can't call this with loop_order set to outer",2);

  //Holds current params in order om, ode, w0, wa, with defaults
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

  //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for (it = fparam.params.begin(); it != fparam.params.end(); ++it) 
    if ( it->second.fit == param_struct::fixed )
      utility::updateStandardCosVec( it->first, fparam, current_cosparams, 
				     it->second.fixval );

  //Now point it1, it2, it3 at the vars we will loop over
  std::map< param_tags::paramcodes, param_struct >::const_iterator it1,it2,it3;
  it1 = fparam.params.find( param_tags::omegam );
  if ( it1 == fparam.params.end() )
    throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			  "Can't find Omega_m--all supported 2D fits have it",
			  16);
  if ( it1->second.fit != param_struct::loop )
    throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			  "Omega_m must be loop fit type",32);

  if ( cosmo_fit == cosfitter::omegam_omegade_w0 ) {
    it2 = fparam.params.find( param_tags::omegade );
    if (it2 == fparam.params.end() )
      throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			    "Can't find Omega_de despite specified fit type",
			    64);
    if (it2->second.fit != param_struct::loop )
      throw CosFitterExcept("cosfitter","cosmofit2D_outer",
			    "Omega_DE must have loop fit type for om/ode fit",
			    128);
    it3 = fparam.params.find( param_tags::w0 );
    if (it3 == fparam.params.end() )
      throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			    "Can't find w0 despite specified fit type",
			    64);
    if (it3->second.fit != param_struct::loop )
      throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			    "w0 must have loop fit type for om/ode/w0 fit",
			    128);
  } else if ( cosmo_fit == cosfitter::flat_omegam_w0_wa ) {
    it2 = fparam.params.find( param_tags::w0 );
    if (it2 == fparam.params.end() )
      throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			    "Can't find w0 despite specified fit type",
			    32);
    if (it2->second.fit != param_struct::loop )
      throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			    "w0 must have loop fit type for flat om/w fit",
			    64);
    it3 = fparam.params.find( param_tags::wa );
    if (it3 == fparam.params.end() )
      throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			    "Can't find wa despite specified fit type",
			    32);
    if (it3->second.fit != param_struct::loop )
      throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			    "wa must have loop fit type for flat om/w0/wa fit",
			    64);
  } else {
      throw CosFitterExcept("cosfitter","cosmofit3D_outer",
			    "Unknown fit type",8);
  }

  cosgrid3D margprob3D( it1->second, it2->second, it3->second );
  margprob3D = 0.0;

  param_struct::fittypes alpha_fittype, beta_fittype;
  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  italpha = fparam.params.find( param_tags::alpha );
  if ( italpha == fparam.params.end() ) 
    alpha_fittype = param_struct::not_set; else
    alpha_fittype = italpha->second.fit;
  itbeta = fparam.params.find( param_tags::beta );
  if ( itbeta == fparam.params.end() ) 
    beta_fittype = param_struct::not_set; else
    beta_fittype = itbeta->second.fit;

  int st;
  double cos0, cos1, cos2;
  double ***data, **ptr2, *ptr1;
  param_tags::paramcodes code0 = margprob3D.getAxisSpec(0).first;
  param_tags::paramcodes code1 = margprob3D.getAxisSpec(1).first;
  param_tags::paramcodes code2 = margprob3D.getAxisSpec(2).first;
  unsigned int n0 = margprob3D.getAxisN(0);
  unsigned int n1 = margprob3D.getAxisN(1);
  unsigned int n2 = margprob3D.getAxisN(2);
  data = margprob3D.getData();
  dl.resize(nsn);
  if ( alpha_fittype == param_struct::loop && 
       beta_fittype == param_struct::loop ) {
    cosgrid2D nuisance( italpha->second, itbeta->second );
    nuisance = 0.0;
    cosgrid2D working( italpha->second, itbeta->second );
    unsigned int progbar_stepsize = n0 / 10;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i=0; i < n0; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / n0)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      cos0 = margprob3D.getAxisVal( 0, i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      ptr2 = data[i];
      for (unsigned int j=0; j < n1; ++j) {
	cos1 = margprob3D.getAxisVal( 1, j );
	utility::updateStandardCosVec( code1, fparam, current_cosparams, 
				       cos1 );
	ptr1 = ptr2[j];
	for (unsigned int k=0; k < n2; ++k) {
	  cos2 = margprob3D.getAxisVal( 2, k );
	  utility::updateStandardCosVec( code2, fparam, current_cosparams, 
					 cos2 );
	  st = lm.getLumDist(sne,dl,current_cosparams);
	  if (st == 0) {
	    calcLikelihood_lab( working );
	    nuisance += working; //Update alpha/beta
	    ptr1[k] = working.getTotal();
	  } else ptr1[k] = 0.0;  //Bad lumdist
	}
      }
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);

  } else if ( alpha_fittype == param_struct::loop ) {
    //loop on alpha, analytic scriptm and either no beta or a fixed value
    cosgrid1D nuisance( italpha->second );
    nuisance = 0.0;
    cosgrid1D working( italpha->second );
    unsigned int progbar_stepsize = n0 / 10;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i=0; i < n0; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / n0)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      cos0 = margprob3D.getAxisVal( 0, i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      ptr2 = data[i];
      for (unsigned int j=0; j < n1; ++j) {
	cos1 = margprob3D.getAxisVal( 1, j );
	utility::updateStandardCosVec( code1, fparam, current_cosparams, 
				       cos1 );
	ptr1 = ptr2[j];
	for (unsigned int k=0; k < n2; ++k) {
	  cos2 = margprob3D.getAxisVal( 2, k );
	  utility::updateStandardCosVec( code2, fparam, current_cosparams, 
					 cos2 );
	  st = lm.getLumDist(sne,dl,current_cosparams);
	  if (st == 0) {
	    calcLikelihood_la( working );
	    ptr1[k] = working.getTotal();
	    nuisance += working; //Update alpha
	  } else ptr1[k] = 0.0;  //Bad lumdist
	}
      }
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout) {
      //Well, write alpha, although we won't have beta
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
    }
  } else if ( beta_fittype == param_struct::loop ) {
    //loop on beta, analytic scriptm and either no alpha or a fixed value
    cosgrid1D nuisance( itbeta->second );
    nuisance = 0.0;
    cosgrid1D working( itbeta->second );
    unsigned int progbar_stepsize = n0 / 10;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i=0; i < n0; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / n0)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      cos0 = margprob3D.getAxisVal( 0, i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      ptr2 = data[i];
      for (unsigned int j=0; j < n1; ++j) {
	cos1 = margprob3D.getAxisVal( 1, j );
	utility::updateStandardCosVec( code1, fparam, current_cosparams, 
				       cos1 );
	ptr1 = ptr2[j];
	for (unsigned int k=0; k < n2; ++k) {
	  cos2 = margprob3D.getAxisVal( 2, k );
	  utility::updateStandardCosVec( code2, fparam, current_cosparams, 
					 cos2 );
	  st = lm.getLumDist(sne,dl,current_cosparams);
	  if (st == 0) {
	    calcLikelihood_la( working );
	    ptr1[k] = working.getTotal();
	    nuisance += working; //Update alpha
	  } else ptr1[k] = 0.0;  //Bad lumdist
	}
      }
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout) {
      //Well, write alpha, although we won't have beta
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
    }
  } else {
    //analytic scriptm, and some combination of no alpha or beta or 
    // fixed values
    unsigned int progbar_stepsize = n0 / 10;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i=0; i < n0; ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / n0)+1 << "\%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      cos0 = margprob3D.getAxisVal( 0, i );
      utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
      ptr2 = data[i];
      for (unsigned int j=0; j < n1; ++j) {
	cos1 = margprob3D.getAxisVal( 1, j );
	utility::updateStandardCosVec( code1, fparam, current_cosparams, 
				       cos1 );
	ptr1 = ptr2[j];
	for (unsigned int k=0; k < n2; ++k) {
	  cos2 = margprob3D.getAxisVal( 2, k );
	  utility::updateStandardCosVec( code2, fparam, current_cosparams, 
					 cos2 );
	  st = lm.getLumDist(sne,dl,current_cosparams);
	  if (st == 0) {
	    ptr1[k] = calcLikelihood();
	  } else ptr1[k] = 0.0;  //Bad lumdist
	}
      }
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
  }
  
  unsigned int i1,i2,i3;
  maxlikelihood = margprob3D.getMaximum(i1,i2,i3);

  //Get info about cosmological params
  margprob3D.normalize();
  GetParamMap( margprob3D, results_map );

  //Write probability distribution
  margprob3D.writeFile( fparam.outputfilename, fparam.binaryout );

}

/*!
  Handles 3D cosmological fit in case where the inner loop is
  the cosmoogical parameters (i.e., when the covariance matrix
  inversion is the most expensive thing to calcluate).  Only two types of
  fits are currently supported: omega_m/omega_de/w0 and 
  flat omega_m/w0/wa

  \param[in] cosmo_fit  The type of fit to do
  \param[out] results_map  A map that will be filled with the
                   results of the fit
  \param[out] maxlikelihood The maximum likelihood

  As a side effect, this will write out a cosgrid as
  a file containing the probability surface.		   
*/
//This one is relatively simple because the only fit we support
// for one-D is fixed curvature, omega_m
void cosfitter::cosmofit3D_inner( cosfitter::cosmo_fittype cosmo_fit,
				  std::map< param_tags::paramcodes,
				  param_results >& results_map, 
				  double& maxlikelihood ) const {
  if ( sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			  "Can't call this with diagonal errors",1);
  
  if ( fparam.loop_order != fitparam::cosmo_inner )
    throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			  "Can't call this with loop_order set to inner",2);

  std::map< param_tags::paramcodes, param_struct >::const_iterator it1,it2,it3;
  //Point it1, it2, it3 at desired params
  if ( cosmo_fit == cosfitter::omegam_omegade_w0 ) {
    it1 = fparam.params.find( param_tags::omegam );
    if (it1 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			    "Can't find Omega_m despite fit type", 4 );
    it2 = fparam.params.find( param_tags::omegade );
    if (it2 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			    "Can't find Omega_DE despite fit type", 4 );
    it3 = fparam.params.find( param_tags::w0 );
    if (it3 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			    "Can't find w0 despite fit type", 4 );
  } else if ( cosmo_fit == cosfitter::flat_omegam_w0_wa ) {
    it1 = fparam.params.find( param_tags::omegam );
    if (it1 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			    "Can't find Omega_m despite fit type", 4 );
    it2 = fparam.params.find( param_tags::w0 );
    if (it2 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			    "Can't find w0 despite fit type", 4 );
    it3 = fparam.params.find( param_tags::wa );
    if (it3 == fparam.params.end())
      throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			    "Can't find wa despite fit type", 4 );
  } else {
    throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			  "Unsupported fit type",8);

  }

  cosgrid3D margprob3D( it1->second, it2->second, it3->second );
  margprob3D = 0.0;
  
  //Calculate the luminosity distances
  lumdist_array3D lar;
  lar.loadData( margprob3D, sne, fparam, fparam.showprogbar );

  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  param_struct::fittypes alpha_fittype, beta_fittype;
  italpha = fparam.params.find( param_tags::alpha );
  if ( italpha == fparam.params.end() ) 
    alpha_fittype = param_struct::not_set; else
    alpha_fittype = italpha->second.fit;
  itbeta = fparam.params.find( param_tags::beta );
  if ( itbeta == fparam.params.end() ) 
    beta_fittype = param_struct::not_set; else
    beta_fittype = itbeta->second.fit;

  unsigned int n2 = margprob3D.getAxisN(2);
  diffarr = new double[ nsn * n2 ];
  prod = new double[ nsn * n2 ];

  double curr_alpha, curr_beta;
  if ( alpha_fittype == param_struct::loop && 
       beta_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on alpha/beta and analytic scriptm
    cosgrid3D working( it1->second, it2->second, it3->second );
    cosgrid2D nuisance( italpha->second, itbeta->second );
    unsigned int nsteps = nuisance.getAxisN(0);
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.verbose) std::cout << "Calculating Likelihoods" << std::endl;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i = 0; i < nuisance.getAxisN(0); ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_alpha = nuisance.getAxisVal(0,i);
      for (unsigned int j = 0; j < nuisance.getAxisN(1); ++j) {
	curr_beta = nuisance.getAxisVal(1,j);
	calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
	margprob3D += working;
	nuisance[i][j] = working.getTotal();
      }
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);

  } else if ( alpha_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on alpha, analytic
    // scriptm and either no beta or a fixed value
    cosgrid3D working( it1->second, it2->second, it3->second );
    cosgrid1D nuisance( italpha->second );
    if ( beta_fittype == param_struct::not_set ) {
      if ( fparam.colour_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			      "Have Beta cov matrix, but not treating",8);
      curr_beta = 0.0;
    } else if ( beta_fittype == param_struct::fixed ) {
      curr_beta = itbeta->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit3D_inner",
				 "Unknown Beta fit type",16);
    unsigned int nsteps = nuisance.getAxisN();
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i = 0; i < nuisance.getAxisN(); ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_alpha = nuisance.getAxisVal(i);
      calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
      margprob3D += working;
      nuisance[i] = working.getTotal();
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
  } else if ( beta_fittype == param_struct::loop ) {
    //Omega_m fixed curvature fit with loop on beta, analytic
    // scriptm and either no alpha or a fixed value
    cosgrid3D working( it1->second, it2->second, it3->second );
    cosgrid1D nuisance( itbeta->second );
    if ( alpha_fittype == param_struct::not_set ) {
      if ( fparam.width_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			      "Have Alpha cov matrix, but not treating",8);
      curr_alpha = 0.0;
    } else if ( alpha_fittype == param_struct::fixed ) {
      curr_alpha = italpha->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit2D_outer",
				 "Unknown Alpha fit type",16);
    unsigned int nsteps = nuisance.getAxisN();
    unsigned int progbar_stepsize = nsteps > 10 ? nsteps / 10 : nsteps;
    std::stringstream progbar;
    if (fparam.showprogbar) {
      progbar << "Progress 0%";
      fprintf(stderr,"%s",progbar.str().c_str());
    }
    for (unsigned int i = 0; i < nuisance.getAxisN(); ++i) {
      if ( fparam.showprogbar && i > 0 && i % progbar_stepsize == 0 ) {
	progbar << "..." << (i * 100 / nsteps)+1 << "%";
	fprintf(stderr,"\r%s",progbar.str().c_str());
      }
      curr_beta = nuisance.getAxisVal(i);
      calcLikelihood_nondia( working, curr_alpha, curr_beta, lar );
      margprob3D += working;
      nuisance[i] = working.getTotal();
    }
    if (fparam.showprogbar) fprintf(stderr,"\n");
    working.free();
    nuisance.normalize();
    GetParamMap( nuisance, results_map );
    if (fparam.albetaout)
      nuisance.writeFile(fparam.albetaoutputfilename,
			 fparam.binaryout);
  } else {
    //Alpha and beta are both fixed or missing
    //This makes things pretty simple
    if ( alpha_fittype == param_struct::not_set ) {
      if ( fparam.width_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			      "Have Alpha cov matrix, but not treating",8);
      curr_alpha = 0.0;
    } else if ( alpha_fittype == param_struct::fixed ) {
      curr_alpha = italpha->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit3D_inner",
				 "Unknown Alpha fit type",16);
    if ( beta_fittype == param_struct::not_set ) {
      if ( fparam.colour_covfileset )
	throw CosFitterExcept("cosfitter","cosmofit3D_inner",
			      "Have Beta cov matrix, but not treating",8);
      curr_beta = 0.0;
    } else if ( beta_fittype == param_struct::fixed ) {
      curr_beta = itbeta->second.fixval;
    } else throw CosFitterExcept("cosfitter","cosmofit3D_inner",
				 "Unknown Beta fit type",16);
    calcLikelihood_nondia( margprob3D, curr_alpha, curr_beta, lar );
  }
  lar.free();

  delete[] prod;
  delete[] diffarr;
  prod = diffarr = NULL;

  unsigned int i1,i2, i3;
  maxlikelihood = margprob3D.getMaximum(i1,i2,i3);

  //Get info about the cosmological parameters
  margprob3D.normalize();
  GetParamMap( margprob3D, results_map );

  //Write probability distribution
  margprob3D.writeFile( fparam.outputfilename, fparam.binaryout );

}

/*! Assuming that the luminosity distance has already been
  calculated, loops over the values of alpha and beta in
  the input array and fills it in with the relative probability.
  This version assumes the errors are diagonal.

  \param[out] likelihoodtarg On output, has the probabilities.  On input,
      the axis must be initialized with alpha and beta.
*/
void cosfitter::calcLikelihood_lab( cosgrid2D& likelihoodtarg ) const {

  if ( !sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","calcLikelihood_la",
			  "Can't handle non-diagonal errors",1);
  if ( fparam.errsfixed ) 
    throw CosFitterExcept("cosfitter","calcLikelihood_la",
			  "Errors shouldn't have been pre-computed",2);

  //We expect the chisquare to be vaguely close to nsn, so to make the
  // probabilities closer to 1 and thus slightly easier to manage,
  // we subtract this off before converting to probability
  //The half is because of how this is used
  double halfchi0 = 0.5 * static_cast<double>(nsn);

  //Initialize to zero to avoid contamination from previous calls
  likelihoodtarg = 0.0;
  
  //There is a slight modification here from the Goliath equations.
  //As written in the paper, they rely on subtracting large numbers
  // from each other, which is not something worth avoiding.
  //So the A and B terms are modified by guessing a value for
  // scriptm and evaluating everything around that instead of around zero,
  // which makes for pretty minor modifications to the formulae
  //So -- first estimate scriptm by forming weighted mean without
  // correctly propagated errors
  double wtval,scriptm0;
  wtval = 0.0;
  scriptm0 = 0.0;
  for (unsigned int i = 0; i < nsn; ++i) {
    wtval += 1.0 / pre_vars[i];
    scriptm0 += (sne[i].mag - dl[i])
      / pre_vars[i];
  }
  scriptm0 /= wtval;

  double alpha, beta, alpha_beta, alpha_sq, beta_sq, chisq;
  double amarg_A, amarg_B, diffmag, error_sq, curr_invvar;
  
  for (unsigned int i=0; i < likelihoodtarg.getAxisN(0); ++i) {
    alpha = likelihoodtarg.getAxisVal(0,i);
    alpha_sq = alpha*alpha;

    for (unsigned int j=0; j < likelihoodtarg.getAxisN(1); ++j) {
      beta = likelihoodtarg.getAxisVal(1,j);
      beta_sq = beta*beta;
      alpha_beta = alpha * beta;

      amarg_A = amarg_B = amarg_E = 0.0;
      for (unsigned int k = 0; k < nsn; ++k) {
	diffmag = sne[k].mag - dl[k] + alpha * (sne[k].widthpar-1)
	  - beta * sne[k].colourpar - scriptm0;
	error_sq = pre_vars[k]
          + alpha_sq * sne[k].var_widthpar
	  + beta_sq * sne[k].var_colourpar
          + 2.0 * alpha * sne[k].cov_mag_widthpar
	  - 2.0 * beta * sne[k].cov_mag_colourpar
          - 2.0 * alpha_beta * sne[k].cov_widthpar_colourpar;
	curr_invvar = 1.0 / error_sq;

	amarg_A += diffmag * diffmag * curr_invvar;
	amarg_B += diffmag * curr_invvar;
	amarg_E += curr_invvar;
      }
      chisq = amarg_A + log(amarg_E * inv_twopi) 
	- amarg_B * amarg_B / amarg_E;
      //convert to relative probability
      likelihoodtarg[i][j] = exp( -0.5 * chisq + halfchi0 );
    }
  } 
}


/*! Assuming that the luminosity distance has already been
  calculated, loops over the values of alpha in
  the input array and fills it in with the relative probability.
  This version assumes the errors are diagonal.

  \param[out] likelihoodtarg On output, has the probabilities.  On input,
      the axis must be initialized with the alpha.
*/
void cosfitter::calcLikelihood_la( cosgrid1D& likelihoodtarg ) const {
  //See calcLikelihood_lab for some additional algorithm comments

  if ( !sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","calcLikelihood_la",
			  "Can't handle non-diagonal errors",1);
  if ( fparam.errsfixed ) 
    throw CosFitterExcept("cosfitter","calcLikelihood_la",
			  "Errors shouldn't have been pre-computed",2);
  double halfchi0 = 0.5 * static_cast<double>(nsn);
  likelihoodtarg = 0.0;
  double wtval,scriptm0;
  wtval = 0.0;
  scriptm0 = 0.0;
  for (unsigned int i = 0; i < nsn; ++i) {
    wtval += 1.0 / pre_vars[i];
    scriptm0 += (sne[i].mag - dl[i])
      / pre_vars[i];
  }
  scriptm0 /= wtval;

  //Check to see if we have beta
  //Note that if it is fixed, then the beta part of the error
  // was already included in pre_vars (see calcPreErrs), but
  // the cross term wasn't
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  it = fparam.params.find( param_tags::beta );
  if ( it != fparam.params.end() ) {
    if (it->second.fit != param_struct::fixed)
      throw CosFitterExcept("cosfitter","calcLikelihood_la",
			    "Non fixed beta value",4);
    
    double alpha, beta, alpha_beta, alpha_sq, beta_sq, chisq;
    double amarg_A, amarg_B, diffmag, error_sq, curr_invvar;
    beta = it->second.fixval; 
    beta_sq = beta*beta;

    for (unsigned int i=0; i < likelihoodtarg.getAxisN(); ++i) {
      alpha = likelihoodtarg.getAxisVal(i);
      alpha_sq = alpha*alpha;
      alpha_beta = alpha * beta;
      
      amarg_A = amarg_B = amarg_E = 0.0;
      for (unsigned int k = 0; k < nsn; ++k) {
	diffmag = sne[k].mag - dl[k] + alpha * (sne[k].widthpar-1)
	  - beta * sne[k].colourpar - scriptm0;
	error_sq = pre_vars[k]
          + alpha_sq * sne[k].var_widthpar
          + 2.0 * alpha * sne[k].cov_mag_widthpar
          - 2.0 * alpha_beta * sne[k].cov_widthpar_colourpar;
	curr_invvar = 1.0 / error_sq;

	amarg_A += diffmag * diffmag * curr_invvar;
	amarg_B += diffmag * curr_invvar;
	amarg_E += curr_invvar;
      }
      chisq = amarg_A + log(amarg_E * inv_twopi) 
	- amarg_B * amarg_B / amarg_E;
      likelihoodtarg[i] = exp( -0.5 * chisq + halfchi0 );
    }
  } else {
    //No beta at all
    double alpha, alpha_sq;
    double amarg_A, amarg_B, diffmag, error_sq, curr_invvar, chisq;
    for (unsigned int i=0; i < likelihoodtarg.getAxisN(); ++i) {
      alpha = likelihoodtarg.getAxisVal(i);
      alpha_sq = alpha*alpha;
      amarg_A = amarg_B = amarg_E = 0.0;
      for (unsigned int k = 0; k < nsn; ++k) {
	diffmag = sne[k].mag - dl[k] + alpha * (sne[k].widthpar-1)
	  - scriptm0;
	error_sq = pre_vars[k]
          + alpha_sq * sne[k].var_widthpar
          + 2.0 * alpha * sne[k].cov_mag_widthpar;
	curr_invvar = 1.0 / error_sq;
	amarg_A += diffmag * diffmag * curr_invvar;
	amarg_B += diffmag * curr_invvar;
	amarg_E += curr_invvar;
      }
      chisq = amarg_A + log(amarg_E * inv_twopi) 
	- amarg_B * amarg_B / amarg_E;
      likelihoodtarg[i] = exp( -0.5 * chisq + halfchi0 );
    }
  }
}

/*! Assuming that the luminosity distance has already been
  calculated, loops over the values of beta in
  the input array and fills it in with the relative probability.
  This version assumes the errors are diagonal.

  \param[out] likelihoodtarg On output, has the probabilities.  On input,
      the axis must be initialized with the beta.
*/
void cosfitter::calcLikelihood_lb( cosgrid1D& likelihoodtarg ) const {
  //See calcLikelihood_lab for some additional algorithm comments
  
  if ( !sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","calcLikelihood_lb",
			  "Can't handle non-diagonal errors",1);

  if ( fparam.errsfixed ) 
    throw CosFitterExcept("cosfitter","calcLikelihood_lb",
			  "Errors shouldn't have been pre-computed",2);
  double halfchi0 = 0.5 * static_cast<double>(nsn);
  likelihoodtarg = 0.0;
  double wtval,scriptm0;
  wtval = 0.0;
  scriptm0 = 0.0;
  for (unsigned int i = 0; i < nsn; ++i) {
    wtval += 1.0 / pre_vars[i];
    scriptm0 += (sne[i].mag - dl[i])
      / pre_vars[i];
  }
  scriptm0 /= wtval;

  //Check to see if we have alpha
  //Note that if it is fixed, then the alpha part of the error
  // was already included in pre_vars (see calcPreErrs), but
  // the cross term wasn't
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  it = fparam.params.find( param_tags::alpha );
  if ( it != fparam.params.end() ) {
    if (it->second.fit != param_struct::fixed)
      throw CosFitterExcept("cosfitter","calcLikelihood_lb",
			    "Non fixed alpha value",4);
    
    double alpha, beta, alpha_beta, alpha_sq, beta_sq, chisq;
    double amarg_A, amarg_B, diffmag, error_sq, curr_invvar;
    alpha = it->second.fixval; 
    alpha_sq = alpha*alpha;
    for (unsigned int i=0; i < likelihoodtarg.getAxisN(); ++i) {
      beta = likelihoodtarg.getAxisVal(i);
      beta_sq = beta*beta;
      alpha_beta = alpha * beta;
      amarg_A = amarg_B = amarg_E = 0.0;
      for (unsigned int k = 0; k < nsn; ++k) {
	diffmag = sne[k].mag - dl[k] + alpha * (sne[k].widthpar-1)
	  - beta * sne[k].colourpar - scriptm0;
	error_sq = pre_vars[k]
          + beta_sq * sne[k].var_colourpar
          - 2.0 * beta * sne[k].cov_mag_colourpar
          - 2.0 * alpha_beta * sne[k].cov_widthpar_colourpar;
	curr_invvar = 1.0 / error_sq;

	amarg_A += diffmag * diffmag * curr_invvar;
	amarg_B += diffmag * curr_invvar;
	amarg_E += curr_invvar;
      }
      chisq = amarg_A + log(amarg_E * inv_twopi) 
	- amarg_B * amarg_B / amarg_E;
      likelihoodtarg[i] = exp( -0.5 * chisq + halfchi0 );

    }
  } else {
    //No alpha at all
    double beta, beta_sq;
    double amarg_A, amarg_B, diffmag, error_sq, curr_invvar, chisq;
    for (unsigned int i=0; i < likelihoodtarg.getAxisN(); ++i) {
      beta = likelihoodtarg.getAxisVal(i);
      beta_sq = beta*beta;
      amarg_A = amarg_B = amarg_E = 0.0;
      for (unsigned int k = 0; k < nsn; ++k) {
	diffmag = sne[k].mag - dl[k] - beta * sne[k].colourpar
	  - scriptm0;
	error_sq = pre_vars[k]
          + beta_sq * sne[k].var_colourpar
          - 2.0 * beta * sne[k].cov_mag_colourpar;
	curr_invvar = 1.0 / error_sq;
	amarg_A += diffmag * diffmag * curr_invvar;
	amarg_B += diffmag * curr_invvar;
	amarg_E += curr_invvar;
      }
      chisq = amarg_A + log(amarg_E * inv_twopi) 
	- amarg_B * amarg_B / amarg_E;
      likelihoodtarg[i] = exp( -0.5 * chisq + halfchi0 );
    }
  }
}


/*! Assuming that the luminosity distance has already been
  calculated, loops over the values of alpha/beta in
  the input array and fills it in with the relative probability.
  This version assumes the errors are diagonal.
*/
double cosfitter::calcLikelihood() const {
  //See calcLikelihood_lab for some additional algorithm comments

  if ( !sne.areErrorsDiagonal() ) 
    throw CosFitterExcept("cosfitter","calcLikelihood",
			  "Can't handle non-diagonal errors",1);
  if ( ! fparam.errsfixed ) 
    throw CosFitterExcept("cosfitter","calcLikelihood",
			  "Errors should be fixed",2);
  double halfchi0 = 0.5 * static_cast<double>(nsn);
  double wtval,scriptm0,curr_invvar;
  wtval = 0.0;
  scriptm0 = 0.0;
  for (unsigned int i = 0; i < nsn; ++i) {
    curr_invvar = 1.0/pre_vars[i];
    wtval += curr_invvar;
    scriptm0 += (sne[i].mag - dl[i]) * curr_invvar;
  }
  scriptm0 /= wtval;

  //Check to see if we have alpha or beta
  //If present they must be fixed.  If so, then their contribution
  // to the error has already been calculated and added in calcPreErrs,
  // including the cross-term if present.  amarg_E should also be
  // pre-loaded
  std::map< param_tags::paramcodes, param_struct >::const_iterator italpha;
  std::map< param_tags::paramcodes, param_struct >::const_iterator itbeta;
  italpha = fparam.params.find( param_tags::alpha );
  itbeta = fparam.params.find( param_tags::beta );
  double amarg_A, amarg_B, diffmag;
  amarg_A = amarg_B = 0.0;
  if ( italpha != fparam.params.end() ) {
    if ( italpha->second.fit != param_struct::fixed ) 
      throw CosFitterExcept("cosfitter","calcLikelihood",
			    "Non fixed alpha",4);
    double alpha = italpha->second.fixval;
    if ( itbeta != fparam.params.end() ) {
      //Alpha and beta
      if ( itbeta->second.fit != param_struct::fixed ) 
	throw CosFitterExcept("cosfitter","calcLikelihood",
			      "Non fixed beta",4);
      double beta = itbeta->second.fixval;
      for (unsigned int i=0; i < nsn; ++i) {
	diffmag = sne[i].mag - dl[i] + alpha * (sne[i].widthpar-1)
	  - beta * sne[i].colourpar - scriptm0;
	curr_invvar = 1.0 / pre_vars[i]; //Remember -- already includes al,beta
	amarg_A += diffmag * diffmag * curr_invvar;
	amarg_B += diffmag * curr_invvar;
      }
    } else {
      //Alpha, no beta
       for (unsigned int i=0; i < nsn; ++i) {
	diffmag = sne[i].mag - dl[i] + alpha * (sne[i].widthpar-1) - scriptm0;
	curr_invvar = 1.0 / pre_vars[i]; //Remember -- already includes al,beta
	amarg_A += diffmag * diffmag * curr_invvar;
	amarg_B += diffmag * curr_invvar;
       }
    }
  } else if ( itbeta != fparam.params.end() ) {
    //Beta, no alpha
    if ( itbeta->second.fit != param_struct::fixed ) 
      throw CosFitterExcept("cosfitter","calcLikelihood",
			    "Non fixed beta",4);
    double beta = itbeta->second.fixval;
    for (unsigned int i=0; i < nsn; ++i) {
      diffmag = sne[i].mag - dl[i] - beta * sne[i].colourpar - scriptm0;
	curr_invvar = 1.0 / pre_vars[i]; //Remember -- already includes beta
	amarg_A += diffmag * diffmag * curr_invvar;
	amarg_B += diffmag * curr_invvar;
    }
  } else {
    //Finally, truly no alpha or beta
    for (unsigned int i=0; i < nsn; ++i) {
      diffmag = sne[i].mag - dl[i] - scriptm0;
      curr_invvar = 1.0 / pre_vars[i]; 
      amarg_A += diffmag * diffmag * curr_invvar;
      amarg_B += diffmag * curr_invvar;
    }
  }

  double chi2 = amarg_A + log(amarg_E * inv_twopi) 
    - amarg_B * amarg_B / amarg_E;

  //std::cout << amarg_A << " " << amarg_E << " " << amarg_B << 
  // " chisq: " << chi2 << std::endl;


  return exp( -0.5 * chi2 + halfchi0 );
}

/*!
  This handles the non-diagonal 0D case, and is overloaded with
  the dimensional cases.  This inverts the covariance matrix
  for given values of \f$\alpha, \beta\f$.

  \param[in] alpha The value of \f$\alpha\f$ to use
  \param[in] beta The value of \f$\beta\f$ to use
  \param[in] lar The precomputed array of luminosity distances
  \return The likelihood of the result
*/
double cosfitter::calcLikelihood_nondia( double alpha, double beta, 
					 const lumdist_array0D& lar) const {
  double likelihood;
  if (!fparam.errsfixed) computeInvCovMatrix( alpha, beta );
  
  //Matrix version
  //Get corrected mags
  std::vector< double > corr_mag; //Holds alpha, beta corrected mags
  corr_mag.resize( nsn );
  for (unsigned int i=0; i < nsn; ++i)
    corr_mag[i] = sne[i].mag + alpha * (sne[i].widthpar-1)
      - beta * sne[i].colourpar;

  //Make array of differences, indexing is SN by cosmparam
  // held in Column-major (i.e., Fortran) order
  // diffarr must be pre-allocated
  std::vector< double > invvars(nsn);
  for (unsigned int i=0; i < nsn; ++i) invvars[i] = invcovmatrix[i][i];
  double const* ldptr;
  double wtval,scriptm0;
  if (lar.isGood() == 0) {
    ldptr = lar.getLumdists();
    for (unsigned int i = 0; i < nsn; ++i)
      diffarr[i] = corr_mag[i] - ldptr[i];
    //Now figure out approximate scriptm and adjust
    wtval = 0.0;
    scriptm0 = 0.0;
    for (unsigned int i = 0; i < nsn; ++i) {
      wtval += invvars[i];
      scriptm0 += diffarr[i] * invvars[i];
    }
    scriptm0 /= wtval;
    for (unsigned int i = 0; i < nsn; ++i)
      diffarr[i] -= scriptm0;
  } else for (unsigned int i = 0; i < nsn; ++i)
    diffarr[i] = 0.0;
  
  //Now we do the matrix multiply 
  //This is V^{-1}_{ij} m_{jk} where ij are SN indicies, k is the
  // cosmo index.  The result should be nsn by ncosmo
#if USEMKL || USEATLAS
  //V is symmetric, so we use DSYMM.  Also, for V row vs. column
  // ordering doesn't matter by symmetry
  //The product matrix (prod) must be pre-allocated
  
  int in0 = 1;
  int insn = static_cast<int>( nsn );
  CBLAS_SIDE side = CblasLeft;
  CBLAS_ORDER order = CblasColMajor;
  CBLAS_UPLO uplo = CblasUpper;
  cblas_dsymm( order, side, uplo, insn, in0, 1.0, invcovmatrix.getData(),
	       insn, diffarr, insn, 0.0, prod, insn );
#elif USEACCELERATE
  //There is something horribly wrong with dsymm in accelerate -- 
  // it is very, very slow.  So we use dgemm
  int in0 = 1;
  int insn = static_cast<int>( nsn );
  CBLAS_ORDER order = CblasColMajor;
  CBLAS_TRANSPOSE trans = CblasNoTrans;
  cblas_dgemm( order, trans, trans, insn, in0, insn, 1.0, 
	       invcovmatrix.getData(), insn, diffarr, insn, 0.0, prod, insn );
#else
  double tmpval;
  double *rowptr;
  for (unsigned int i = 0; i < nsn; ++i) {
    rowptr = invcovmatrix[i];
    tmpval = rowptr[0] * diffarr[0];
    for (unsigned int k = 1; k < nsn; ++k)
      tmpval += rowptr[k] * diffarr[k];
    prod[ i ] = tmpval;
  }
#endif
  //Then do the non matrix prod part, prod is nsn by ncosmo (nsn,n0)
  //The result should be a vector which is ncosmo long with the
  // chisqs
  double amarg_A, amarg_B;
  double chisq, halfchi0;
  halfchi0 = 0.5 * static_cast<double>(nsn);
  if (lar.isGood() == 0) {
    amarg_A = diffarr[0] * prod[0];
    for (unsigned int i = 1; i < nsn; ++i )
      amarg_A += diffarr[i] * prod[i];
    amarg_B = diffarr[0] * invcovmatrix_rowsums[0];
    for (unsigned int i = 1; i < nsn; ++i )
      amarg_B += diffarr[i] * invcovmatrix_rowsums[i];
    chisq = amarg_A + log( amarg_E * inv_twopi ) - 
      amarg_B*amarg_B/amarg_E;
    likelihood = exp( -0.5 * chisq + halfchi0 );
  } else likelihood = 0.0;

  return likelihood;

}



/*!
  This handles the non-diagonal 1D case, and is overloaded with
  the multi-dimensional cases.  This inverts the covariance matrix
  for given values of \f$\alpha, \beta\f$ and loops over the cosmological 
  parameters specified by likelihoodtarg.

  \param[out] likelihoodtarg On input must be loaded with the axis
    information for the cosmological parameter, on output is filled with
    the relative probability values.
  \param[in] alpha The value of \f$\alpha\f$ to use
  \param[in] beta The value of \f$\beta\f$ to use
  \param[in] lar The precomputed array of luminosity distances

*/
void cosfitter::calcLikelihood_nondia( cosgrid1D& likelihoodtarg, double alpha,
				       double beta, const lumdist_array1D& lar) const {

#ifdef TIMING
  clock_t starttime = clock();
#endif
  if (!fparam.errsfixed) computeInvCovMatrix( alpha, beta );
#ifdef TIMING
  tot_inv_time += clock() - starttime;  
#endif

  //Matrix version
  unsigned int n0 = likelihoodtarg.getAxisN();
  //Get corrected mags
  std::vector< double > corr_mag; //Holds alpha, beta corrected mags
  corr_mag.resize( nsn );
  for (unsigned int i=0; i < nsn; ++i)
    corr_mag[i] = sne[i].mag + alpha * (sne[i].widthpar-1)
      - beta * sne[i].colourpar;

  //Make array of differences, indexing is SN by cosmparam
  // held in Column-major (i.e., Fortran) order
  // diffarr must be pre-allocated
  std::vector< double > invvars(nsn);
  for (unsigned int i=0; i < nsn; ++i) invvars[i] = invcovmatrix[i][i];
  double const* ldptr;
  double* diffptr;
  double* prodptr;
  double wtval,scriptm0;
  for (unsigned int j = 0; j < n0; ++j) {
    diffptr = diffarr + j * nsn;
    if (lar.isGood(j) == 0) {
      ldptr = lar.getLumdists( j );
      for (unsigned int i = 0; i < nsn; ++i) {
	diffptr[i] = corr_mag[i] - ldptr[i];
      }
      //Now figure out approximate scriptm and adjust
      wtval = 0.0;
      scriptm0 = 0.0;
      for (unsigned int i = 0; i < nsn; ++i) {
	wtval += invvars[i];
	scriptm0 += diffptr[i] * invvars[i];
      }
      scriptm0 /= wtval;
      for (unsigned int i = 0; i < nsn; ++i)
	diffptr[i] -= scriptm0;
    } else for (unsigned int i = 0; i < nsn; ++i)
      diffptr[i] = 0.0;
  }

#ifdef TIMING  
  starttime = clock();
#endif

  //Now we do the matrix multiply 
  //This is V^{-1}_{ij} m_{jk} where ij are SN indicies, k is the
  // cosmo index.  The result should be nsn by ncosmo
#if USEMKL || USEATLAS
  //V is symmetric, so we use DSYMM.  Also, for V row vs. column
  // ordering doesn't matter by symmetry
  //The product matrix (prod) must be pre-allocated
  int in0 = static_cast<int>( n0 );
  int insn = static_cast<int>( nsn );
  CBLAS_SIDE side = CblasLeft;
  CBLAS_ORDER order = CblasColMajor;
  CBLAS_UPLO uplo = CblasUpper;
  //CBLAS_TRANSPOSE trans = CblasNoTrans;
  cblas_dsymm( order, side, uplo, insn, in0, 1.0, invcovmatrix.getData(),
	       insn, diffarr, insn, 0.0, prod, insn );
#elif USEACCELERATE
  //There is something horribly wrong with dsymm in accelerate -- 
  // it is very, very slow.  So we use dgemm
  int in0 = static_cast<int>( n0 );
  int insn = static_cast<int>( nsn );
  CBLAS_ORDER order = CblasColMajor;
  CBLAS_TRANSPOSE trans = CblasNoTrans;
  cblas_dgemm( order, trans, trans, insn, in0, insn, 1.0, 
	       invcovmatrix.getData(), insn, diffarr, insn, 0.0, prod, insn );
#else
  double tmpval;
  double *rowptr;
  for (unsigned int i = 0; i < nsn; ++i) {
    rowptr = invcovmatrix[i];
    for (unsigned int j = 0; j < n0; ++j) {
      diffptr = diffarr + j*nsn;
      tmpval = rowptr[0] * diffptr[0];
      for (unsigned int k = 1; k < nsn; ++k)
        tmpval += rowptr[k] * diffptr[k];
      prod[ j * nsn + i ] = tmpval;
    }
  } 
#endif

  //Then do the non matrix prod part, prod is nsn by ncosmo (nsn,n0)
  //The result should be a vector which is ncosmo long with the
  // chisqs
  double amarg_A, amarg_B;
  double chisq, halfchi0;
  halfchi0 = 0.5 * static_cast<double>(nsn);
  for (unsigned int j = 0; j < n0; ++j) {
    if (lar.isGood(j) == 0) {
      diffptr = diffarr + j * nsn;
      prodptr = prod + j * nsn;
      amarg_A = diffptr[0]*prodptr[0];
      for (unsigned int i = 1; i < nsn; ++i )
	amarg_A += diffptr[i] * prodptr[i];
      amarg_B = diffptr[0] * invcovmatrix_rowsums[0];
      for (unsigned int i = 1; i < nsn; ++i )
	amarg_B += diffptr[i] * invcovmatrix_rowsums[i];
      chisq = amarg_A + log( amarg_E * inv_twopi ) - 
	amarg_B*amarg_B/amarg_E;
      likelihoodtarg[j] = exp( -0.5 * chisq + halfchi0 );
    } else likelihoodtarg[j] = 0.0;
  }

#ifdef TIMING
  tot_mult_time += clock() - starttime;  
#endif
}

/*!
  This handles the non-diagonal 2D case, and is overloaded with
  the other-dimensional cases.  This inverts the covariance matrix
  for given values of \f$\alpha, \beta\f$ and loops over the cosmological 
  parameters specified by likelihoodtarg.

  \param[out] likelihoodtarg On input must be loaded with the axis
    information for the cosmological parameter, on output is filled with
    the relative probability values.
  \param[in] alpha The value of \f$\alpha\f$ to use
  \param[in] beta The value of \f$\beta\f$ to use
  \param[in] lar The precomputed array of luminosity distances
*/
void cosfitter::calcLikelihood_nondia( cosgrid2D& likelihoodtarg, double alpha,
				       double beta, const lumdist_array2D& lar) const {

#ifdef TIMING
  clock_t starttime = clock();
#endif
  if (! fparam.errsfixed ) computeInvCovMatrix( alpha, beta );
#ifdef TIMING
  tot_inv_time += clock() - starttime;  
#endif

  //Matrix version
  std::vector< double > corr_mag; //Holds alpha, beta corrected mags
  corr_mag.resize( nsn );
  for (unsigned int i=0; i < nsn; ++i)
    corr_mag[i] = sne[i].mag + alpha * (sne[i].widthpar-1)
      - beta * sne[i].colourpar;

  unsigned int n0, n1;
  n0 = likelihoodtarg.getAxisN(0);
  n1 = likelihoodtarg.getAxisN(1);

  //Pre-load
  std::vector< double > invvars(nsn);
  for (unsigned int i=0; i < nsn; ++i) invvars[i] = invcovmatrix[i][i];

  //Cosmo loop
  double wtval, scriptm0;
  double const* ldptr;
  double *diffptr, *prodptr;
  likelihoodtarg = 0.0;
  dl.resize( nsn );
#ifdef TIMING  
  starttime = clock();
#endif
  for (unsigned int i=0; i < n0; ++i) {
    //Now, for each of these we make an array of differences on the second
    // cosmological parameter.  diffarr is nsn by ncosparam and in
    // column major (fortran) ordering.  It must be pre-allocated
    for (unsigned int j=0; j < n1; ++j) {
      diffptr = diffarr + j * nsn;
      if (lar.isGood(i,j) == 0) {
	//Lumdist is good, copy it into dl
	ldptr = lar.getLumdists( i, j );
	for (unsigned int k = 0; k < nsn; ++k)
	  diffptr[k] = corr_mag[k] - ldptr[k];
	//Now figure out approximate scriptm and adjust
	wtval = 0.0;
	scriptm0 = 0.0;
	for (unsigned int k = 0; k < nsn; ++k) {
	  wtval += invvars[k];
	  scriptm0 += diffptr[k] * invvars[k];
	}
	scriptm0 /= wtval;
	for (unsigned int k = 0; k < nsn; ++k)
	  diffptr[k] -= scriptm0;
      } else for (unsigned int k = 0; k < nsn; ++k)
	diffptr[k] = 0.0;
    }

    //Now we do the matrix multiply 
    //This is V^{-1}_{ij} m_{jk} where ij are SN indicies, k is the
    // cosmo index.  The result should be nsn by ncosmo
#if USEMKL || USEATLAS 
    //V is symmetric, so we use DSYMM.  Also, for V row vs. column
    // ordering doesn't matter by symmetry
    //The product matrix (prod) must be pre-allocated
    
    int in1 = static_cast<int>( n1 );
    int insn = static_cast<int>( nsn );
    CBLAS_SIDE side = CblasLeft;
    CBLAS_ORDER order = CblasColMajor;
    CBLAS_UPLO uplo = CblasUpper;
    cblas_dsymm( order, side, uplo, insn, in1, 1.0, invcovmatrix.getData(),
		 insn, diffarr, insn, 0.0, prod, insn );
#elif USEACCELERATE
    //There is something horribly wrong with dsymm in accelerate -- 
    // it is very, very slow.  So we use dgemm
    int in1 = static_cast<int>( n1 );
    int insn = static_cast<int>( nsn );
    CBLAS_ORDER order = CblasColMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    cblas_dgemm( order, trans, trans, insn, in1, insn, 1.0, 
		 invcovmatrix.getData(), insn, diffarr, insn, 0.0, prod, insn );
#else
    //Do by hand
    double tmpval;
    double *rowptr;
    for (unsigned int k = 0; k < nsn; ++k) {
      rowptr = invcovmatrix[k];
      for (unsigned int j = 0; j < n1; ++j) {
	diffptr = diffarr + j*nsn;
	tmpval = rowptr[0] * diffptr[0];
	for (unsigned int m = 1; m < nsn; ++m)
	  tmpval += rowptr[m] * diffptr[m];
	prod[ j*nsn + k ] = tmpval;
      }
    }  
#endif
    //Then do the non matrix prod part, prod is nsn by ncosmo (nsn,n1)
    //The result should be a vector which is ncosmo long with the
    // chisqs
    double amarg_A, amarg_B;
    double chisq, halfchi0;
    halfchi0 = 0.5 * static_cast<double>(nsn);
    for (unsigned int j = 0; j < n1; ++j) {
      if (lar.isGood(i,j) == 0) {
	diffptr = diffarr + j * nsn;
	prodptr = prod + j * nsn;
	amarg_A = diffptr[0] * prodptr[0];
	for (unsigned int k = 1; k < nsn; ++k )
	  amarg_A += diffptr[k] * prodptr[k];
	amarg_B = diffptr[0] * invcovmatrix_rowsums[0];
	for (unsigned int k = 1; k < nsn; ++k )
	  amarg_B += diffptr[k] * invcovmatrix_rowsums[k];
	chisq = amarg_A + log( amarg_E * inv_twopi ) - 
	  amarg_B*amarg_B/amarg_E;
	likelihoodtarg[i][j] = exp( -0.5 * chisq + halfchi0 );
      } else likelihoodtarg[i][j] = 0.0;
    }
  }
#ifdef TIMING
  tot_mult_time += clock() - starttime;  
#endif
}

/*!
  This handles the non-diagonal 3D case, and is overloaded with
  the other-dimensional cases.  This inverts the covariance matrix
  for given values of \f$\alpha, \beta\f$ and loops over the cosmological 
  parameters specified by likelihoodtarg.

  \param[out] likelihoodtarg On input must be loaded with the axis
    information for the cosmological parameter, on output is filled with
    the relative probability values.
  \param[in] alpha The value of \f$\alpha\f$ to use
  \param[in] beta The value of \f$\beta\f$ to use
  \param[in] lar The precomputed array of luminosity distances
*/
void cosfitter::calcLikelihood_nondia( cosgrid3D& likelihoodtarg, double alpha,
				       double beta, const lumdist_array3D& lar) const {
  
#ifdef TIMING
  clock_t starttime = clock();
#endif
  if (! fparam.errsfixed ) computeInvCovMatrix( alpha, beta );
#ifdef TIMING
  tot_inv_time += clock() - starttime;  
#endif

  std::vector< double > corr_mag; //Holds alpha, beta corrected mags
  corr_mag.resize( nsn );
  for (unsigned int i=0; i < nsn; ++i)
    corr_mag[i] = sne[i].mag + alpha * (sne[i].widthpar-1)
      - beta * sne[i].colourpar;

  unsigned int n0, n1, n2;
  n0 = likelihoodtarg.getAxisN(0);
  n1 = likelihoodtarg.getAxisN(1);
  n2 = likelihoodtarg.getAxisN(2);

  //Pre-load
  std::vector< double > invvars(nsn);
  for (unsigned int i=0; i < nsn; ++i) invvars[i] = invcovmatrix[i][i];

  //Cosmo loop
  double chisq, halfchi0, wtval, scriptm0;
  double const* ldptr;
  double *diffptr, *prodptr;
  halfchi0 = 0.5 * static_cast<double>(nsn);
  likelihoodtarg = 0.0;
  dl.resize(nsn);

#ifdef TIMING  
  starttime = clock();
#endif
  for (unsigned int i=0; i < n0; ++i) {
    for (unsigned int j=0; j < n1; ++j) {
      //Now, for each of these we make an array of differences on the second
      // cosmological parameter.  diffarr is nsn by ncosparam and in
      // column major (fortran) ordering.  It must be pre-allocated
      for (unsigned int k=0; k < n2; ++k) {
	diffptr = diffarr + k * nsn;
	if (lar.isGood(i,j,k) == 0) {
	  //Lumdist is good, copy it into dl
	  ldptr = lar.getLumdists( i, j, k );
	  for (unsigned int m = 0; m < nsn; ++m)
	    diffptr[m] = corr_mag[m] - ldptr[m];
	  //Now figure out approximate scriptm and adjust
	  wtval = 0.0;
	  scriptm0 = 0.0;
	  for (unsigned int m = 0; m < nsn; ++m) {
	    wtval += invvars[m];
	    scriptm0 += diffptr[m] * invvars[m];
	  }
	  scriptm0 /= wtval;
	  for (unsigned int m = 0; m < nsn; ++m)
	  diffptr[m] -= scriptm0;
	} else for (unsigned int m = 0; m < nsn; ++m)
	  diffptr[m] = 0.0;
      }

      //Now we do the matrix multiply 
      //This is V^{-1}_{ij} m_{jk} where ij are SN indicies, k is the
      // cosmo index.  The result should be nsn by ncosmo
#if USEMKL || USEATLAS
      //V is symmetric, so we use DSYMM.  Also, for V row vs. column
      // ordering doesn't matter by symmetry
      //The product matrix (prod) must be pre-allocated
      int in2 = static_cast<int>( n2 );
      int insn = static_cast<int>( nsn );
      CBLAS_SIDE side = CblasLeft;
      CBLAS_ORDER order = CblasColMajor;
      CBLAS_UPLO uplo = CblasUpper;
      cblas_dsymm( order, side, uplo, insn, in2, 1.0, invcovmatrix.getData(),
		   insn, diffarr, insn, 0.0, prod, insn );
#elif USEACCELERATE
      //There is something horribly wrong with dsymm in accelerate -- 
      // it is very, very slow.  So we use dgemm
      int in2 = static_cast<int>( n2 );
      int insn = static_cast<int>( nsn );
      CBLAS_ORDER order = CblasColMajor;
      CBLAS_TRANSPOSE trans = CblasNoTrans;
      cblas_dgemm( order, trans, trans, insn, in2, insn, 1.0, 
		   invcovmatrix.getData(), insn, diffarr, insn, 0.0, 
		   prod, insn );
#else
      //Do by hand
      double tmpval;
      double *rowptr;
      for (unsigned int m = 0; m < nsn; ++m) {
	rowptr = invcovmatrix[m];
	for (unsigned int k = 0; k < n2; ++k) {
	  diffptr = diffarr + k*nsn;
	  tmpval = rowptr[0] * diffptr[0];
	  for (unsigned int n = 1; n < nsn; ++n)
	    tmpval += rowptr[n] * diffptr[n];
	  prod[ k*nsn + m ] = tmpval;
	}
      }
#endif
      //Then do the non matrix prod part, prod is nsn by ncosmo (nsn,n1)
      //The result should be a vector which is ncosmo long with the
      // chisqs
      double amarg_A, amarg_B;
      for (unsigned int k = 0; k < n2; ++k) {
	if (lar.isGood(i,j,k) == 0) {
	  diffptr = diffarr + k * nsn;
	  prodptr = prod + k * nsn;
	  amarg_A = diffptr[0] * prodptr[0];
	  for (unsigned int m = 1; m < nsn; ++m )
	    amarg_A += diffptr[m] * prodptr[m];
	  amarg_B = diffptr[0] * invcovmatrix_rowsums[0];
	  for (unsigned int m = 1; m < nsn; ++m )
	    amarg_B += diffptr[m] * invcovmatrix_rowsums[m];
	  chisq = amarg_A + log( amarg_E * inv_twopi ) - 
	    amarg_B*amarg_B/amarg_E;
	  likelihoodtarg[i][j][k] = exp( -0.5 * chisq + halfchi0 );
	} else likelihoodtarg[i][j][k] = 0.0;
      }
    }
  }
#ifdef TIMING
  tot_mult_time += clock() - starttime;  
#endif
}


/*!
  \param[in] results_map Map of information about cosmological parameters
  \param[in] intrinsic Assumed intrinsic 'dispersion' in magnitudes
*/
void cosfitter::indivchi( const std::map< param_tags::paramcodes, 
			  param_results >& results_map, 
			  const std::map< unsigned int, double>& 
			  intrinsic) const {
  double fitmag, chisqcurr, varcurr, dmag, totchisq, emptyfac;
  string fmt;
  const double zfacsq = (5.0/log(10.0))*(5.0/log(10.0));
  std::map< param_tags::paramcodes, param_results >::const_iterator it;

  //Grab the params
  double w0, wa, om, ode, alpha, beta, scriptm;
  it = results_map.find( param_tags::w0 );
  if ( it != results_map.end() ) w0 = it->second.getVal(); else w0 = -1.0;
  it = results_map.find( param_tags::wa );
  if ( it != results_map.end() ) wa = it->second.getVal(); else wa = 0.0;
  it = results_map.find( param_tags::omegam );
  if ( it != results_map.end() ) om = it->second.getVal(); else om = 1.0;
  it = results_map.find( param_tags::omegade );
  if ( it != results_map.end() ) ode = it->second.getVal(); else 
    ode = 1.0 - fparam.ocurv - om;
  it = results_map.find( param_tags::alpha );
  if ( it != results_map.end() ) alpha = it->second.getVal(); else 
    alpha = 0.0;
  it = results_map.find( param_tags::beta );
  if ( it != results_map.end() ) beta = it->second.getVal(); else 
    beta = 0.0;
  it = results_map.find( param_tags::scriptm );
  if ( it != results_map.end() ) scriptm = it->second.getVal(); else 
    scriptm = 0.0;

  lm.getLumDist( sne, dl, om, ode, w0, wa );

  totchisq = 0.0;

  printf("For parameters\n");
  printf(" w0: %6.3f wa: %6.3f Om: %6.3f Ol: %6.3f\n", w0, wa, om, ode);
  printf(" alpha: %6.3f beta: %6.3f scriptM: %6.3f\n", alpha, beta, scriptm);
  if (fparam.usekomatsuform) 
    printf("Using Komatsu w(a) form with atrans: %6.4f\n",fparam.atrans);
  if (!sne.areErrorsDiagonal())
    printf("Warning -- individual values don't include covariances\n");

  fmt = "%-11s %-6s %-7s %-6s %-7s %-6s %8s %-8s %-8s\n";
  printf(fmt.c_str(),
	 "SN","zcmb","mag","str","fitmag","errs","diff","chisq","sigma");
  fmt = "%-11s %-6.4f %-7.4f %-6.3f %-7.4f %-6.4f %8.4f %-8.4f %-8.4f\n";
  double dzerrsq,rms;
  rms = 0.0;

  double currintrinsicsq;
  std::map<unsigned int, double>::const_iterator pos;

  for (unsigned int i = 0; i < nsn; ++i) {
    fitmag = dl[i] - ( alpha*(sne[i].widthpar-1.0) ) 
      + ( beta*sne[i].colourpar ) + scriptm;

    pos = intrinsic.find( sne[i].dataset );
    if (pos == intrinsic.end()) {
      std::stringstream errstr;
      errstr << "Unknown dataset number " << sne[i].dataset;
      errstr << " for sn: " << sne[i].name;
      throw CosFitterExcept("cosfitter","indivchi",
			    errstr.str(),1);
    }
    currintrinsicsq = pos->second * pos->second;

    //Error part
    emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
    dzerrsq = fparam.pecz * fparam.pecz + sne[i].var_z;
    dzerrsq *= zfacsq*emptyfac*emptyfac;
    varcurr = sne[i].var_mag + dzerrsq + currintrinsicsq + 
      alpha*alpha*sne[i].var_widthpar +
      beta*beta*sne[i].var_colourpar +
      2*alpha*sne[i].cov_mag_widthpar 
      - 2*beta*sne[i].cov_mag_colourpar
      - 2*alpha*beta*sne[i].cov_widthpar_colourpar;
    dmag = fitmag - sne[i].mag;
    rms += dmag*dmag;
    chisqcurr = dmag*dmag/varcurr;
    totchisq += chisqcurr;

    printf( fmt.c_str() ,
	    sne[i].name.c_str(),sne[i].zcmb,sne[i].mag,sne[i].widthpar,
	    fitmag,sqrt(varcurr), dmag,
	    chisqcurr, sqrt(chisqcurr)); 
  }
  rms = sqrt( rms / static_cast<double>(nsn) );
  
  printf("Total ChiSquare: %7.3f for %u SN\n",totchisq,nsn);
  printf("RMS around best fit: %7.3f\n",rms);

}

/*!
  \param outstream Output stream to write to
  \param[in] fittag Type of fit performed
  \param[in] results_map Map of results of fit
*/
//If you make changes in here, make sure you update convert_to_fits,
// which reads this output file.
void cosfitter::print_results( std::ostream& outstream, 
			       const std::string& fittag,
			       const std::map< param_tags::paramcodes,
			       param_results >& results_map ) const {
  const double zfacsq = 25.0/( log(10.0)*log(10.0) );
  double fitmag, varcurr, totchisq, emptyfac;
  std::map< param_tags::paramcodes, param_results >::const_iterator it;

  bool fixedintrinsic; //!< True if there is only one for all sne
  double currintrinsicsq;
  if ( fparam.nintrinsicdisp == 0 ) {
    fixedintrinsic = 1;
    std::map<unsigned int,double>::const_iterator pos = 
      fparam.intrinsicdisp.find( 0 );
    currintrinsicsq = pos->second * pos->second;
  } else {
    currintrinsicsq = 0.0;
    fixedintrinsic = 0;
  }

  //Grab the params
  double w0, wa, om, ode, alpha, beta, scriptm;
  it = results_map.find( param_tags::w0 );
  if ( it != results_map.end() ) w0 = it->second.getVal(); else w0 = -1.0;
  it = results_map.find( param_tags::wa );
  if ( it != results_map.end() ) wa = it->second.getVal(); else wa = 0.0;
  it = results_map.find( param_tags::omegam );
  if ( it != results_map.end() ) om = it->second.getVal(); else om = 1.0;
  it = results_map.find( param_tags::omegade );
  if ( it != results_map.end() ) ode = it->second.getVal(); else 
    ode = 1.0 - fparam.ocurv - om;
  it = results_map.find( param_tags::alpha );
  if ( it != results_map.end() ) alpha = it->second.getVal(); else 
    alpha = 0.0;
  it = results_map.find( param_tags::beta );
  if ( it != results_map.end() ) beta = it->second.getVal(); else 
    beta = 0.0;
  it = results_map.find( param_tags::scriptm );
  if ( it != results_map.end() ) scriptm = it->second.getVal(); else 
    scriptm = 0.0;

  lm.getLumDist( sne, dl, om, ode, w0, wa );

  totchisq = 0.0;

  ios::fmtflags oldFlags = outstream.flags();

  outstream << "Fit type: " << fittag << std::endl;
  outstream << "Parameters: " << std::endl;
  outstream << "w0: " << std::setprecision(3) << w0 << std::endl;
  outstream << "wa: " << std::setprecision(3) << wa << std::endl;
  outstream << "Om: " << std::setprecision(4) << om << std::endl;
  outstream << "Ode: " << std::setprecision(4) << ode << std::endl;
  outstream << "Sm: " << std::setprecision(5) << scriptm << std::endl;
  outstream << "Alpha: "  << std::setprecision(5) << alpha << std::endl;
  outstream << "Beta: " << std::setprecision(5) << beta << std::endl;
  if (fparam.nintrinsicdisp == 0) {
    outstream << "Intrinsicdisp: " << std::setprecision(3) 
	      << sqrt(currintrinsicsq) << std::endl; 
  } else {
    outstream << "Number of intrinsic dispersion datasets: "
	      << fparam.nintrinsicdisp << std::endl;
    for (std::map<unsigned int,double>::const_iterator pos = 
	   fparam.intrinsicdisp.begin(); pos != fparam.intrinsicdisp.end();
	 ++pos) 
      outstream << " Dataset: " << pos->first << " intrinsicdisp: " 
		<< pos->second << std::endl;
  }
  outstream << "Pecz: " << std::setprecision(4) << fparam.pecz << std::endl;
  outstream << "Number of objects: " << nsn << std::endl;
  if (!sne.areErrorsDiagonal())
    outstream << "Warning -- individual values don't include covariances"
	      << std::endl;
  outstream << left << setw(9) << "SN" << right << setw(7) << "zcmb" << 
    setw(7) << "zhel" << setw(9) << "raw_mag" << setw(7) << "err" << 
    setw(8) << "D_L+sm" << setw(8) << "corr" << setw(8) << "fitmag" << 
    setw(7) << "err" << setw(8) << "diff" << setw(8) << "chisq" << 
    setw(8) << "sigma" << setw(9) << "width" << setw(9) << "dwidth" << 
    setw(8) << "colour" << setw(8) << "dcolour" << setw(4) <<
    "set" << std::endl;
  
  double rms;
  double dzerrsq, dmag, chisqcurr, corr;
  totchisq = 0.0;
  rms = 0.0;
  for (unsigned int i = 0; i < nsn; ++i) {
    corr = ( -alpha * (sne[i].widthpar-1.0) ) + 
      ( beta * sne[i].colourpar );
    fitmag = dl[i] + corr + scriptm;

    emptyfac = (1.0+sne[i].zcmb)/(sne[i].zcmb*(1+0.5*sne[i].zcmb));
    dzerrsq = fparam.pecz * fparam.pecz + sne[i].var_z;
    dzerrsq *= zfacsq*emptyfac*emptyfac;

    varcurr = sne[i].var_mag + dzerrsq 
      + (alpha*alpha*sne[i].var_widthpar)
      + (beta*beta*sne[i].var_colourpar)
      + (2.0*alpha*sne[i].cov_mag_widthpar)
      - (2.0*beta*sne[i].cov_mag_colourpar)
      - (2.0*alpha*beta*sne[i].cov_widthpar_colourpar);

    if (! fixedintrinsic) {
      std::map<unsigned int,double>::const_iterator pos = 
	fparam.intrinsicdisp.find( sne[i].dataset );
      //Assume it's there, or the previous code would have crashed
      currintrinsicsq = pos->second * pos->second;
    }
    varcurr += currintrinsicsq;

    dmag = fitmag - sne[i].mag ;
    rms += dmag*dmag;

    chisqcurr = dmag*dmag/varcurr;
    totchisq += chisqcurr;

    outstream << std::fixed << std::left << std::setw(9) <<
      sne[i].name << std::right << 
      std::setw(7) << std::setprecision(3) << sne[i].zcmb << 
      std::setw(7) << std::setprecision(3) << sne[i].zhel <<
      std::setw(9) << std::setprecision(3) << sne[i].mag << 
      std::setw(7) << std::setprecision(3) << sqrt(sne[i].var_mag) << 
      std::setw(8) << std::setprecision(3) << dl[i]+scriptm << 
      std::setw(8) << std::setprecision(3) << corr <<
      std::setw(8) << std::setprecision(3) << fitmag <<
      std::setw(7) << std::setprecision(3) << sqrt(varcurr) << 
      std::setw(8) << std::setprecision(3) << dmag <<      
      std::setw(8) << std::setprecision(3) << chisqcurr <<
      std::setw(8) << std::setprecision(3) << sqrt(chisqcurr) << 
      std::setw(9) << std::setprecision(3) << sne[i].widthpar << 
      std::setw(9) << std::setprecision(3) << sqrt(sne[i].var_widthpar) <<
      std::setw(8) << std::setprecision(3) << sne[i].colourpar <<
      std::setw(8) << std::setprecision(3) << sqrt(sne[i].var_colourpar) <<
      std::setw(4) << sne[i].dataset << std::endl;
  }
  rms = sqrt( rms / static_cast<double>( nsn ) );

  outstream << "Total chisquare: " << fixed << std::setprecision(3) <<
    totchisq << std::endl;
  outstream << "RMS around fit: " << std::fixed << std::setprecision(3) <<
    rms << std::endl;

  outstream.flags( oldFlags );  //Restore old status

}

/*!
  \param[in] results_map Map of parameter values
*/
void  
cosfitter::write_paramsummary(const std::map< param_tags::paramcodes,
			      param_results >& results_map ) const {
  
  if (! fparam.paramsummary ) {
    string errstr = "Param summary filename must be set";
    throw CosFitterExcept("cosfitter","write_paramsummary",errstr,1);
  }
  
  ofstream ofs( fparam.paramsummaryfilename.c_str() );
  if (! ofs ) {
    std::stringstream errstrng;
    errstrng << "Error opening param summary output file: " << 
      fparam.paramsummaryfilename;
    throw CosFitterExcept("cosfitter","write_paramsummary",
			  errstrng.str(),2);
  }

  //Grab the parameter values
  std::map< param_tags::paramcodes, param_results >::const_iterator it;
  it = results_map.find( param_tags::w0 );
  if ( it != results_map.end() ) ofs << it->second << std::endl;
  it = results_map.find( param_tags::wa );
  if ( it != results_map.end() )  ofs << it->second << std::endl;
  it = results_map.find( param_tags::omegam );
  if ( it != results_map.end() )  ofs << it->second << std::endl;
  it = results_map.find( param_tags::omegade );
  if ( it != results_map.end() )  ofs << it->second << std::endl;
  it = results_map.find( param_tags::alpha );
  if ( it != results_map.end() )  ofs << it->second << std::endl;
  it = results_map.find( param_tags::beta );
  if ( it != results_map.end() )  ofs << it->second << std::endl;
  it = results_map.find( param_tags::scriptm );
  if ( it != results_map.end() )  ofs << it->second << std::endl;

  ofs.close();

}

/*!
  \param[in] results_map Results of cosmological fit
  \return The estimate of \f${\mathcal M}\f$ given the other parameters.
 */
double cosfitter::estimate_scriptm( const std::map< param_tags::paramcodes, 
				    param_results >& results_map ) const {
  //Grab the parameter values
  double w0, wa, om, ode, alpha, beta, scriptm;
  std::map< param_tags::paramcodes, param_results >::const_iterator it;
  it = results_map.find( param_tags::w0 );
  if ( it != results_map.end() ) w0 = it->second.getVal(); else w0 = -1.0;
  it = results_map.find( param_tags::wa );
  if ( it != results_map.end() ) wa = it->second.getVal(); else wa = 0.0;
  it = results_map.find( param_tags::omegam );
  if ( it != results_map.end() ) om = it->second.getVal(); else om = 1.0;
  it = results_map.find( param_tags::omegade );
  if ( it != results_map.end() ) ode = it->second.getVal(); else 
    ode = 1.0 - fparam.ocurv - om;
  it = results_map.find( param_tags::alpha );
  if ( it != results_map.end() ) alpha = it->second.getVal(); else 
    alpha = 0.0;
  it = results_map.find( param_tags::beta );
  if ( it != results_map.end() ) beta = it->second.getVal(); else 
    beta = 0.0;
  it = results_map.find( param_tags::scriptm );
  if ( it != results_map.end() ) scriptm = it->second.getVal(); else 
    scriptm = 0.0;

  lm.getLumDist( sne, dl, om, ode, w0, wa );

  if ( sne.areErrorsDiagonal() ) {
    double fitmag, invvar, wt, wtmean;
    wt = 0.0;
    wtmean = 0.0;
    if ( fparam.errsfixed ) {
      for (unsigned int i = 0; i < nsn; ++i) {
	fitmag = dl[i]
	  - ( alpha*(sne[i].widthpar-1) )
	  + ( beta*sne[i].colourpar );
	invvar = 1.0 / pre_vars[i];
	wt += invvar;
	wtmean += (sne[i].mag - fitmag) * invvar;
      }
    } else {
      for (unsigned int i = 0; i < nsn; ++i) {
	fitmag = dl[i] 
	  - alpha*(sne[i].widthpar-1) + beta*sne[i].colourpar;
	if ( have_fixed_alpha ) {
	  //Alpha term already in there, beta and cross term not
 	  invvar = 1.0 / (pre_vars[i] + (beta*beta*sne[i].var_colourpar)
			  - (2.0*beta*sne[i].cov_mag_colourpar) 
			  - (2.0*alpha*beta*sne[i].cov_widthpar_colourpar) );
	} else if ( have_fixed_beta ) {
	  //Beta in there, alpha and cross term no
 	  invvar = 1.0 / (pre_vars[i] + (alpha*alpha*sne[i].var_widthpar) 
			  + (2.0*alpha*sne[i].cov_mag_widthpar)
			  - (2.0*alpha*beta*sne[i].cov_widthpar_colourpar) );
	} else {
	  //Neither alpha nor beta in there
 	  invvar = 1.0 / (pre_vars[i] + (alpha*alpha*sne[i].var_widthpar) 
			  + (beta*beta*sne[i].var_colourpar)
			  + (2.0*alpha*sne[i].cov_mag_widthpar)
			  - (2.0*beta*sne[i].cov_mag_colourpar) 
			  - (2.0*alpha*beta*sne[i].cov_widthpar_colourpar) );
	}
	wt += invvar;
	wtmean += (sne[i].mag - fitmag) * invvar;
      }
    }
    return wtmean / wt;
  } else {
    if (!fparam.errsfixed) computeInvCovMatrix( alpha, beta );
    double diffmag, amarg_B;  //amarg_E computed by computeInvCovMatrix
    amarg_B = 0.0;
    for (unsigned int i = 0; i < nsn; ++i) {
      diffmag = sne[i].mag - dl[i] + alpha * (sne[i].widthpar-1)
	- beta * sne[i].colourpar;
      amarg_B += diffmag * invcovmatrix_rowsums[i];
    }
    return amarg_B / amarg_E;
  }
} 

/*!
  Read in the parameter file.
*/
void cosfitter::readParamFile() {

  fparam.readFile( paramfilename );

}

/*!
  Reads in parameters from parameterfile, then loads the supernova
  data and calculates the non-parameter dependent variances.
 */
void cosfitter::prep(const std::string& paramfile) {
  paramfilename = paramfile;

  readParamFile();
  
  sne.setName("Supernova");
  sne.readData( fparam.datafilename );

  nsn = sne.size();

  //Handle reading full style cov matrix
  if ( fparam.mag_covfileset || fparam.width_covfileset || 
       fparam.colour_covfileset || fparam.magwidth_covfileset || 
       fparam.magcolour_covfileset || fparam.widthcolour_covfileset ) {
    std::vector< std::string > FileNames(6);
    if (fparam.mag_covfileset) FileNames[0] = fparam.mag_covfilename; else
      FileNames[0] = "";
    if (fparam.mag_covfileset) FileNames[1] = fparam.width_covfilename; else
      FileNames[1] = "";
    if (fparam.mag_covfileset) FileNames[2] = fparam.colour_covfilename; else
      FileNames[2] = "";
    if (fparam.mag_covfileset) FileNames[3] = fparam.magwidth_covfilename; else
      FileNames[3] = "";
    if (fparam.mag_covfileset) 
      FileNames[4] = fparam.magcolour_covfilename; else FileNames[4] = "";
    if (fparam.mag_covfileset) 
      FileNames[5] = fparam.widthcolour_covfilename; else FileNames[5] = "";
    sne.readCovData( FileNames );
  }  


  if ( fparam.woodbury_covfileset ) {
    sne.readWoodburyCovData( fparam.woodbury0_covfilename,
			     fparam.woodburya_covfilename,
			     fparam.woodburyb_covfilename );
    if ( sne.getWoodburyCovMatrixSNeRef().getNsn() != nsn ) 
      throw CosFitterExcept("cosfitter","prep",
			    "Woodbury covmatrix not same size as SN data",8);
  }  

  sne.zcmbsort(); 

  //Calculate variances
  calcPreErrs();

  //Set up Luminosity distance parameters
  if (fparam.usekomatsuform) {
    lm.setUseKomatsuForm();
    lm.setAtrans( fparam.atrans );
  } else lm.unsetUseKomatsuForm();
  lm.setDz( fparam.dzint ); //In case the trapezoidal rule is used

  isprepped = true;
  
}

//Main driver
void cosfitter::dofit() const {

  if (! isprepped )
    throw CosFitterExcept("cosfitter","dofit",
			  "Attempt to run fit without prepping",1);

  //Now figure out what sort of fit we are going to do
  if (fparam.verbose) cout << "Starting fit" << std::endl;

  //Figure out how many cosmological loop params we have
  unsigned int cosmo_fit_dimension = 0;
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for ( it = fparam.params.begin(); it != fparam.params.end(); ++it)
    if (it->second.param_spec.second == param_tags::cosmological && 
	it->second.fit == param_struct::loop ) ++cosmo_fit_dimension;

  //Results variables
  std::string fittag;
  std::map< param_tags::paramcodes, param_results > results_map;
  double maxlikelihood;

  LoadFixedParams( fparam.params, results_map );

  switch (cosmo_fit_dimension) {
  case 0 :
    cosmofit0D( results_map, maxlikelihood, fittag );
    break; 
  case 1 : 
    cosmofit1D( results_map, maxlikelihood, fittag );
    break; 
  case 2: 
    cosmofit2D( results_map, maxlikelihood, fittag );
    break;
  case 3: 
    cosmofit3D( results_map, maxlikelihood, fittag );
    break;
  default :
    std::stringstream errstrng;
    errstrng << "Currently unsupported number of cosmological params: " 
	     << cosmo_fit_dimension;
    throw CosFitterExcept("cosfitter","dofit",
			  errstrng.str(),4);
    break;
  }

  //Add scriptm
  double sm = estimate_scriptm( results_map );
  param_results sm_param;
  sm_param.name = param_tags::scriptmtag;
  sm_param.value = sm;
  sm_param.mostlikelyval = sm;
  sm_param.fit = param_struct::analytic;
  sm_param.spec.first = param_tags::scriptm;
  sm_param.spec.second = param_tags::nuisance;
  results_map[ param_tags::scriptm ] = sm_param;

  //Verbosity stuff
  if (fparam.paramsummary)
    write_paramsummary( results_map );
  if (fparam.verbose || fparam.extendedout) {
    if (fparam.verbose) {
      cout << "Number of SNe: " << nsn << std::endl;
      cout << "Maximum likelihood value: " << maxlikelihood << std::endl;
      cout << "Equivalent chisquare (assuming gaussian): " <<
        -2.0*log( maxlikelihood ) + nsn << std::endl;
      if (fparam.fixcurv) cout << "Value of Omega_curv: " 
			       << fparam.ocurv << std::endl;
      if (fparam.usekomatsuform) 
	std::cout << "Using Komatsu form of w(a) with "
		  << " atrans: " << fparam.atrans << std::endl;
    }

    if (fparam.verbose) {
      for (map<param_tags::paramcodes, param_results>::const_iterator i = 
	     results_map.begin(); i != results_map.end(); ++i) 
	std::cout << i->second << std::endl;

      indivchi( results_map, fparam.intrinsicdisp );
    }

    if (fparam.extendedout) {
      ofstream ofs( fparam.extendedoutfilename.c_str() );
      if (!ofs) {
        stringstream errstrng("");
        errstrng << "Error opening extended output file: " <<
          fparam.extendedoutfilename;
        throw CosFitterExcept("cosfitter","cosmofit1D",
                              errstrng.str(),2);
      }
      print_results( ofs, fittag, results_map );
      ofs.close();
    }

  }
}
