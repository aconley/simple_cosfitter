//lumdist.cc

//Luminosity distance integrator
//Given a list of SNe redshifts, calculates and returns an array
// of D_l evaluated at each position (for a fixed Omega_m and Omega_lambda).
// D_l is the 'hubble constant free' luminosity distance (c/H0 has been
//  removed).

#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>

#include <lumdist.h>
#include <cosfitterexcept.h>

using namespace std;

const double lumdist::prec_lum = 1e-5;

/*! 
  \brief General w=-1 integrand in GSL form
  \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] params Array of parameters 
     (\f$\Omega_m\f$, \f$\Omega_{\Lambda}\f$)
  \returns The value of the integral term 
             \f$ \frac{1}{\sqrt{ \left(1+z\right)^2\left(1+\Omega_m z \right)
                                 - z \left( 2 + z \right) \Omega_{\Lambda}}}\f$

 */
inline double calcVal_gsl(double z, void* params) {
  double zpo,val;
  const double *in_params = static_cast<const double*>(params);
  double om = in_params[0];
  double ol = in_params[1];
  zpo=(1.0+z);
  val=zpo*zpo*(1+om*z)-z*(2+z)*ol;
  return 1./sqrt(val);
}

/*! 
  \brief General w=-1 integrand in GSL ODE form
  \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] y Current value of integral.  Not used
  \param[out] f The value of the integral term
             \f$ \frac{1}{\sqrt{ \left(1+z\right)^2\left(1+\Omega_m z \right)
                  - z \left( 2 + z \right) \Omega_{\Lambda}}}\f$
  \param[in] params Array of parameters 
     (\f$\Omega_m\f$, \f$\Omega_{\Lambda}\f$)
 */
inline int calcVal_gsl_ode(double z, const double y[1], 
			   double f[1], void* params) {
  double zpo,val;
  const double *in_params = static_cast<const double*>(params);
  double om = in_params[0];
  double ol = in_params[1];
  zpo=(1.0+z);
  val=zpo*zpo*(1+om*z)-z*(2+z)*ol;
  f[0] = 1./sqrt(val);
  return GSL_SUCCESS;
}

/*! 
  \brief Integrand for constant w in GSL form
  \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] params Array of parameters 
     (\f$w, \Omega_m, \Omega_{DE}\f$)
  \returns The value of the integral term
     \f$ \frac{1}{\sqrt{ \left(1+z\right)^2 \Omega_m + 
             \left(1+z\right)^2 \Omega_k +
             \left(1+z\right)^{3\left(1+w\right)} \Omega_{DE} } } \f$

 */
inline double calcValW_gsl(double z, void* params) {
  double zpo,val;
  const double *in_params = static_cast<const double*>(params);
  double w = in_params[0];
  double om = in_params[1];
  double ow = in_params[2];
  double ok = 1.0 - om - ow;
  zpo=(1.0+z);
  val = zpo*zpo * ( zpo*om + ok ) + ow*pow(zpo,3.0*(1+w));
  return 1./sqrt(val);
}

/*! 
  \brief General w=const integrand in GSL ODE form
  \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] y Current value of integral.  Not used
  \param[out] f The value of the integral term
     \f$ \frac{1}{\sqrt{ \left(1+z\right)^3 \Omega_m + 
             \left(1+z\right)^2 \Omega_k +
             \left(1+z\right)^{3\left(1+w\right)} \Omega_{DE} } } \f$
  \param[in] params Array of parameters 
     (\f$w\f$, \f$\Omega_m\f$, \f$\Omega_{DE}\f$)
 */
inline int calcValW_gsl_ode(double z, const double y[1], 
			   double f[1], void* params) {
  double zpo,val;
  const double *in_params = static_cast<const double*>(params);
  double w = in_params[0];
  double om = in_params[1];
  double ow = in_params[2];
  double ok = 1.0 - om - ow;
  zpo=(1.0+z);
  val = zpo*zpo * ( zpo*om + ok ) + ow*pow(zpo,3.0*(1+w));
  f[0] = 1./sqrt(val);
  return GSL_SUCCESS;
}

/*! 
  \brief Integrand for w0, wa in GSL form
  \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] params Array of parameters 
     (\f$w0, wa, \Omega_m, \Omega_{DE}\f$)
  \returns The value of the integral term
*/
inline double calcValW0WA_gsl(double z, void* params) {
  double zpo,val;
  const double *in_params = static_cast<const double*>(params);
  double w0 = in_params[0];
  double wa = in_params[1];
  double om = in_params[2];
  double ow = in_params[3];
  double ok = 1.0 - om - ow;
  zpo=(1.0+z);
  val = zpo*zpo * ( zpo*om + ok ) + 
    ow * pow( zpo,3.0*(1+w0+wa) ) * exp( -3.0 * wa * z/zpo );
  return 1./sqrt(val);
}

/*! 
  \brief General w0, wa integrand in GSL ODE form
  \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] y Current value of integral.  Not used
  \param[out] f The value of the integral term
  \param[in] params Array of parameters 
     (\f$w_0, w_a, \Omega_m, \Omega_{DE}\f$)
 */
inline int calcValW0WA_gsl_ode(double z, const double y[1], 
			       double f[1], void* params) {
  double zpo,val;
  const double *in_params = static_cast<const double*>(params);
  double w0 = in_params[0];
  double wa = in_params[1];
  double om = in_params[2];
  double ow = in_params[3];
  double ok = 1.0 - om - ow;
  zpo=(1.0+z);
  val = zpo*zpo * ( zpo*om + ok ) + 
    ow * pow( zpo,3.0*(1+w0+wa) ) * exp( -3.0 * wa * z/zpo );
  f[0] = 1./sqrt(val);
  return GSL_SUCCESS;
}

/*! 
  \brief Integrand for w0, wa in GSL form, using Komatsu parameterization
  \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] params Array of parameters 
     (\f$w0, wa, \Omega_m, \Omega_{DE}, atrans\f$)
  \returns The value of the integral term
*/
inline double calcValW0WAKom_gsl(double z, void* params) {
  double zpo,val;
  const double *in_params = static_cast<const double*>(params);
  double w0 = in_params[0];
  double wa = in_params[1];
  double om = in_params[2];
  double ow = in_params[3];
  double atrans = in_params[4];
  double ok = 1.0 - om - ow;
  double opat = 1.0 + atrans;
  zpo=(1.0+z);
  val = zpo*zpo * ( zpo*om + ok ) + 
    ow * exp( -3.0*( z*wa/zpo + (1+w0+opat*wa)*log( (1/zpo + atrans)/opat ) ) );
  return 1./sqrt(val);
}

/*! 
  \brief General w0, wa integrand in GSL ODE form
  \ingroup lumdist
  Not in lumdist class because the GSL does not accept pointers to member
  functions.
  \param[in] z Redshift
  \param[in] y Current value of integral.  Not used
  \param[out] f The value of the integral term
  \param[in] params Array of parameters 
     (\f$w_0, w_a, \Omega_m, \Omega_{DE}, atrans\f$)
 */
inline int calcValW0WAKom_gsl_ode(double z, const double y[1], 
				  double f[1], void* params) {
  double zpo,val;
  const double *in_params = static_cast<const double*>(params);
  double w0 = in_params[0];
  double wa = in_params[1];
  double om = in_params[2];
  double ow = in_params[3];
  double atrans = in_params[4];
  double ok = 1.0 - om - ow;
  double opat = 1.0 + atrans;
  zpo=(1.0+z);
  val = zpo*zpo * ( zpo*om + ok ) + 
    ow * exp( -3.0*( z*wa/zpo + (1+w0+opat*wa)*log( (1/zpo + atrans)/opat ) ) );
  f[0] = 1./sqrt(val);
  return GSL_SUCCESS;
}


void lumdist::init_gsl() {
  //Avoid aborting on error handler
  gsl_set_error_handler_off();
  
  T45 = const_cast<gsl_odeiv_step_type*>(gsl_odeiv_step_rkf45);
  s_ode = gsl_odeiv_step_alloc(T45,1);
  c_ode = gsl_odeiv_control_y_new(prec_lum, 0.0);
  e_ode = gsl_odeiv_evolve_alloc(1);

  intspace = gsl_integration_workspace_alloc(512);
}

lumdist::lumdist() : param_tol(1e-3), useKomatsuForm(false),
		     atrans(0.1), dz(0.0001) {
  init_gsl();
}

lumdist::lumdist(const std::vector<double>& zval) : 
  param_tol(1e-3), useKomatsuForm(false), atrans(0.1), dz(0.0001) {
  init_gsl();
  zcmb.reserve( zval.size() );
  zcmb.assign( zval.begin(), zval.end() );
}

lumdist::~lumdist() {
  gsl_odeiv_evolve_free(e_ode);
  gsl_odeiv_control_free(c_ode);
  gsl_odeiv_step_free(s_ode);
  gsl_integration_workspace_free(intspace);
}


void lumdist::setZcmb(const std::vector<double>& zvals) const {
  zcmb.reserve( zvals.size() );
  zcmb.assign( zvals.begin(), zvals.end());
}

void lumdist::setZcmb(const SNeData& sne) const {
  zcmb.resize( sne.size() );
  for (unsigned int i = 0; i < zcmb.size(); ++i)
    zcmb[i] = sne[i].zcmb;
}

/*!
  Evaluates only the integral part of the luminosity distance term
  \param[in] om \f$\Omega_m\f$
  \param[in] ol \f$\Omega_{\Lambda}\f$
  \param[out] intval Filled with the integral values
  \returns 0 on success

  Uses a brute-force trapezoidal algorithm.
*/
int lumdist::lumDistSumInt(std::vector<double>& intval, double om, 
			   double ol) const {
  bool canbounce;
  unsigned int nz = zcmb.size();
  double zmax, olbounce, temp, inner=0.0;
  double z, currz, val, currint, running_sum, last_grid;

  intval.resize(nz);

  zmax = zcmb[nz-1]; //Assume sorted
  
  if (om == 0 && ol >= 1.0) return 1;
  if (om > 0.5) {
    inner = cos( acos( (1-om)/om )/3.0 );
    olbounce = 4*om*pow( inner , 3 );
  } else {
    temp = (1-om)/om;
    inner = cosh( log( temp + sqrt(temp*temp - 1) ) / 3.0 );
    olbounce = 4*om*pow( inner , 3 );
  }
  canbounce = (ol >= olbounce);

  //If we can bounce, we have a lower limit on the maximum allowable
  // redshift, so we can toss out at least some cases, although not all
  //if (canbounce) {
  //  if (zmax > (2*inner-1)) return 1;
  //} 
  if (canbounce) return 1;

  //Initialize sum_running to -1/2 first point (first term in trapezoidal
  // rule)
  running_sum = - calcVal(0.0,om,ol) * dz / 2.0;

  last_grid = 0.0;

  //Main sum loop
  val = 0.0;
  for (unsigned int i = 0; i < nz; ++i) {
    currz = zcmb[i];

    for (z = last_grid; z < currz; z += dz) {
      val = calcVal(z,om,ol);
      running_sum += val * dz;
    }

    //z will have overstepped by dz as the exit condition
    last_grid = z - dz;

    //Subtract of the 1/2 end point for the trapezoidal rule
    currint = running_sum -  val * dz / 2.0;

    //Since we will start at the last grid point again, we need to
    // subtract off the last added value.  Otherwise we will add
    // the last grid point twice
    running_sum -= val*dz;

    //Add in the last non-grid point value evaluated at the midpoint
    // between the two
    val = calcVal( (last_grid + currz)/2.0, om, ol);
    currint += val * (currz - last_grid);

    intval[i] = currint;
  }

  return 0;
}

/*!
  Uses a fairly brute-force trapezoidal algorithm.
  \param[in] sne List of SN to find luminosity distance integrals for.  Must
                  be sorted by SNeDataEntry::zcmb
  \param[in] om \f$\Omega_m\f$
  \param[in] ol \f$\Omega_{\Lambda}\f$
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success
*/
int lumdist::lumDistSum(const SNeData& sne, std::vector<double>& dl,
			double om, double ol) const {
  
  const double pi = 3.14159265358979323846264338327950288419716939937510582;

  setZcmb( sne );
  unsigned int nsn = sne.size();

  int st = lumDistSumInt( dl, om, ol );
  if (st) return st;

  //Apply curvature stuff
  double sqrtk, omplusol;
  omplusol = om + ol;
  if (fabs(omplusol - 1.0) < 0.001) {
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * dl[i] );
    }
  } else if (omplusol < 1.0) {
    sqrtk = sqrt( fabs( 1 - omplusol ) );
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10 ((1.0 + sne[i].zhel) * sinh(sqrtk * dl[i])/sqrtk );
    }
  } else {
    sqrtk = sqrt( fabs( 1 - omplusol ) );
    double compval = pi / sqrtk;
    for (unsigned int i = 0; i < nsn; ++i) 
      if (dl[i] >= compval || isnan(dl[i]) ) return 1;
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * sin(sqrtk * dl[i])/sqrtk );
    }
  }

  return 0;

}


/*!
  Evaluates only the integral part of the luminosity distance term
  \param[in] om \f$\Omega_m\f$
  \param[in] ol \f$\Omega_{\Lambda}\f$
  \param[out] intval Filled with the integral values
  \returns 0 on success

  Uses a Runge-Kutta-Felhberg 4-5 method
*/
int lumdist::lumDistSumInt_rk4(std::vector<double>& intval, double om, 
			       double ol) const {
  bool canbounce;
  unsigned int nz;
  double z, zmax, olbounce, temp, inner=0.0;

  nz = zcmb.size();
  intval.resize(nz);

  zmax = zcmb[nz-1]; //Assume sorted

  //Check bounciness
  if (om == 0 && ol >= 1.0) return 1; //Bounces
  if (om > 0.5) {
    inner = cos( acos( (1-om)/om )/3.0 );
    olbounce = 4*om*pow( inner , 3 );
  } else {
    temp = (1-om)/om;
    inner = cosh( log( temp + sqrt(temp*temp - 1) ) / 3.0 );
    olbounce = 4*om*pow( inner , 3 );
  }
  canbounce = (ol >= olbounce);

  //If we can bounce, we have a lower limit on the maximum allowable
  // redshift, so we can toss out at least some cases, although not all
  //if (canbounce) {
  //  if (zmax > (2*inner-1)) return 1;
  //} 
  if (canbounce) return 1;

  double params[2];
  params[0] = om; params[1] = ol;

  gsl_odeiv_system sys = {&calcVal_gsl_ode, 0, 1, &params};

  z = 0;
  double h = 1e-5; //Initial step guess
  if ( h < zcmb[0] ) h = zcmb[0]/2.0;
  double y[1] = { 0.0 };  //Inital value of integral

  int st;
  double zcurrtarg;
  for (unsigned int i = 0; i < nz; ++i) {
    zcurrtarg = zcmb[i];
    while (z < zcurrtarg) {
      st = gsl_odeiv_evolve_apply (e_ode, c_ode, s_ode, &sys, &z, zcurrtarg, 
				   &h, y);
      if (st) return 1; //Failure
    }
    intval[i] = y[0];
  }
  return 0;
}

/*!
  Uses Runge-Kutta-Felhberg 4-5 method.  Faster than trapezoidal algorithm
  \param[in] sne List of SN to find luminosity distance integrals for.  Must
                  be sorted by SNeDataEntry::zcmb
  \param[in] om \f$\Omega_m\f$
  \param[in] ol \f$\Omega_{\Lambda}\f$, the density parameter of the 
                   cosmological constant
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success
*/
int lumdist::lumDistSum_rk4(const SNeData& sne,std::vector<double>& dl,
			    double om, double ol) const {

  const double pi = 3.14159265358979323846264338327950288419716939937510582;

  setZcmb( sne );
  unsigned int nsn = zcmb.size();

  int st = lumDistSumInt_rk4( dl, om, ol );
  if (st) return st;

  //Apply curvature stuff
  double sqrtk, omplusol;
  omplusol = om + ol;
  if (fabs(omplusol - 1.0) < 0.001) {
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * dl[i] );
    }
  } else if (omplusol < 1.0) {
    sqrtk = sqrt( fabs( 1 - omplusol ) );
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10 ((1.0 + sne[i].zhel) * sinh(sqrtk * dl[i])/sqrtk );
    }
  } else {
    sqrtk = sqrt( fabs( 1 - omplusol ) );
    double compval = pi/sqrtk;
    for (unsigned int i = 0; i < nsn; ++i)
      if (dl[i] >= compval || isnan(dl[i]) ) return 1;
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * sin(sqrtk * dl[i])/sqrtk );
    }
  }

  return 0;
 

}

/*!
  Based on calculations presented in Kantowski et al., 2000 ApJ 545, 549
  \param[in] zcmb Redshift to evaluate at
  \param[in] om \f$\Omega_m\f$
  \param[in] ol \f$\Omega_{\Lambda}\f$
  \param[out] integralval The integral value
  \returns 0 on success

  These equations require ol >= 0.
*/
int lumdist::singleLumDistInt(double zcmb, double om, double ol,
			      double& integralval) const {
  const double paramtol = 1e-5;

  double ok; //Omega_k
  ok = 1.0 - om - ol;
  double kappa;
  kappa = - ok/abs(ok);  //-1 for open, 0 flat, 1 close

  if (ol < 0.0) return 1;

  double bparam;  //Equation (8) in Kantowski
  bparam = - 27.0/2.0 * om*om*ol/(ok*ok*ok);

  if (ol == 0) {
    //Mattig distance, corrected for difference between zhel and zcmb,
    // something that is usually ignored
    integralval = 2.0/(om*om) * (om*zcmb + (om-2.0)*(sqrt(1+om*zcmb)-1) );
    integralval /= (1 + zcmb);
    return 0;
  } else if (om == 0) {
    if (ol <= 0.0) return 1; //Can't do this value
    if (ol >= 1.0) return 1; //Nor this
    integralval = (1+zcmb)*(1+zcmb-sqrt(ol + (1+zcmb)*(1+zcmb)*(1-ol)))/ol;
    return 0;
  } else if ( abs(ok) < paramtol ) {
    //Assume flat universe, but ol = 0 already taken care of
    double k = sqrt( 0.5 + sqrt(3.0)/4.0 );
    double deltaphi, top, bottom, tmpfac;
    tmpfac = pow( 1.0/om - 1, 1.0/3.0 );
    top = -zcmb * sqrt( sqrt(3.0)*om*tmpfac ) *
      sqrt( 1.0 + zcmb / (1.0 + tmpfac ) );
    bottom = 1.0 + zcmb / (1.0 + tmpfac) + 
      sqrt( 1.0 + om*zcmb*(3 + 3*zcmb + zcmb*zcmb) );
    deltaphi = 2 * atan( top/bottom );
    
    integralval = 1.0 / sqrt( sqrt(3.0) * om*tmpfac ) *
      ( - gsl_sf_ellint_F( deltaphi, k, GSL_PREC_DOUBLE ) );
    return 0;
  } else if ( (bparam < 0) || (bparam > 2) ) {
    //Equations (9-15) and (7) in Kantowski
    double vk;
    vk = pow( kappa*(bparam-1) + sqrt( bparam*(bparam-2) ), 1.0/3.0 );
    double y1;
    y1 = 1.0/3.0 * ( -1 + kappa*(vk + 1.0/vk) );
    double A;
    A = sqrt( y1 * (3*y1 + 2 ) );
    double g = 1.0/sqrt(A);
    double k = g * sqrt( (2*A + kappa*(1+3*y1)) / 4.0 );
    double deltaphiz, top, bottom, opz;
    opz = 1.0 + zcmb;
    top = -zcmb * sqrt( A * abs(ok ) * (1 + zcmb / (1 - ok * y1 / om ) ) );
    bottom = 1 + zcmb / (1 - ok * y1 / om ) + 
      sqrt( opz*opz*(1+om*zcmb) - zcmb*(zcmb+2)*ol );
    deltaphiz = 2 * atan( top/bottom );
    
    integralval = -g*gsl_sf_ellint_F( deltaphiz, k, GSL_PREC_DOUBLE );
    return 0;
    
  } else if (bparam > 0 && bparam < 2) {

    //Check against bouncing universes
    double inner, olbounce;
    if (om > 0.5) {
      inner = cos( acos( (1-om)/om )/3.0 );
      olbounce = 4*om*pow( inner , 3 );
    } else {
      double temp = (1-om)/om;
      inner = cosh( log( temp + sqrt(temp*temp - 1) ) / 3.0 );
      olbounce = 4*om*pow( inner , 3 );
    }
    if (ol >= olbounce) return 1; //Bouncing universe
    
    //Otherwise we are in lower right hand sliver
    double y1,y2,y3;

    y1 = 1.0/3.0 * ( -1.0 + cos( acos( 1.0 - bparam )/3.0 ) +
		     sqrt(3.0) * sin( acos( 1.0 - bparam )/3.0 ) );
    y2 = 1.0/3.0 * ( -1.0 - 2 * cos( acos( 1.0 - bparam )/3.0 ) );
    y3 = 1.0/3.0 * ( -1.0 + cos( acos( 1.0 - bparam ) / 3.0 ) -
		     sqrt(3.0) * sin( acos( 1.0 - bparam ) / 3.0 ) );
    double g = 2 / sqrt( y1 - y2 );
    double k = sqrt( (y1-y3)/(y1-y2) );
    double deltaphiz, top, bottom, opz;
    opz = 1.0 + zcmb;
    top = sqrt( y1 - y2 ) * ( sqrt(y3 - om/ok) - sqrt(y3 - opz)*om/ok );
    bottom = sqrt( (y1 - om/ok) * (y2 - opz*om/ok) ) +
      sqrt( (y2 - om/ok) * (y1 - opz*om/ok) );
    deltaphiz = 2.0 * atan( top/bottom );
    
    integralval = -g*gsl_sf_ellint_F( deltaphiz, k, GSL_PREC_DOUBLE );
    return 0;
    
  } else if (bparam == 2) {
    //Check against bouncing universes
    double inner, olbounce;
    if (om > 0.5) {
      inner = cos( acos( (1-om)/om )/3.0 );
      olbounce = 4*om*pow( inner , 3 );
    } else {
      double temp = (1-om)/om;
      inner = cosh( log( temp + sqrt(temp*temp - 1) ) / 3.0 );
      olbounce = 4*om*pow( inner , 3 );
    }
    if (ol >= olbounce) return 1; //Bouncing universe
    
    double top, bottom, opz;
    opz = 1.0 + zcmb;
    top = ( sqrt( 1.0/3.0 - om/ok ) + 1.0 ) * 
      ( sqrt(1.0/3.0 - opz * om / ok ) - 1.0 );
    bottom = ( sqrt( 1.0/3.0 - om/ok) - 1.0 ) *
      ( sqrt(1.0/3.0 - opz * om/ok) + 1.0 );
    
    integralval = log( top / bottom );
    return 0;
  }
  
  //Hopefully we won't get here, but if we do throw the
  // full integral at the problem
  int gsl_retcode;
  size_t neval;
  double error;
  gsl_function F;
  double params[2];
  params[0] = om; params[1] = ol;
  F.function = &calcVal_gsl;
  F.params = &params;
  
  gsl_retcode = gsl_integration_qng(&F,0,zcmb,0,prec_lum,&integralval,
				    &error,&neval);
  
  if (gsl_retcode) {
    //Simple integration didn't work -- call adaptive routine
      gsl_retcode = 
	gsl_integration_qags(&F,0,zcmb,0,prec_lum,512,intspace,&integralval,
			     &error);
      if (gsl_retcode) {
	std::cerr << "Error -- integral not evaluated" << std::endl;
	return 2;
      }
  }
  return 0;

}

/*!
  Based on calculations presented in Kantowski et al., 2000 ApJ 545, 549
  \param[in] SN Object with SN information
  \param[in] om \f$\Omega_m\f$
  \param[in] ol \f$\Omega_{\Lambda}\f$
  \param[out] dl luminosity distance (c/H0 removed).
            Actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success

  This requires ol >= 0.0.
*/
int lumdist::singleLumDist(const SNeDataEntry& SN, double om, double ol,
			   double &dl) const {

  const double pi = 3.14159265358979323846264338327950288419716939937510582;

  double zcmb = SN.zcmb;
  double zhel = SN.zhel;
  double dlval;

  double ok = 1.0 - om - ol;

  int st = singleLumDistInt(zcmb,om,ol,dlval);
  if (st) return st;

  //Note that dlval already includes one sqrtk factor in this case
  if (fabs(ok) < 0.001) {
    dl = 5.0 * log10( (1.0 + zhel) * dlval );
  } else if (ok > 0) {
    double sqrtk = sqrt( fabs( ok ) );
    dl = 5*log10( (1.0 + zhel) * sinh(dlval)/sqrtk );
  } else {
    double sqrtk = sqrt( fabs( ok ) );
    double compval = pi / sqrtk;
    if (dlval >= compval || isnan(dlval) ) return 1;
    dl = 5*log10( (1.0 + zhel) * sin(dlval)/sqrtk );
  }
  return 0;

}

/*!
  Evaluates only the integral part of the luminosity distance term
  \param[in] w \f$w\f$
  \param[in] om \f$\Omega_m\f$
  \param[in] ol \f$\Omega_{\Lambda}\f$
  \param[out] intval Filled with the integral values
  \returns 0 on success

  Uses a brute-force trapezoidal algorithm.
*/
int lumdist::lumDistSumWInt(std::vector<double>& intval, double w,
			    double om, double ol) const {
  //We can no longer test for bouncing cases because w complicates
  // this significantly
  
  //Test for cosmological constant case
  if (w == -1.0) 
    return lumDistSumInt(intval,om,ol);

  unsigned int nz = zcmb.size();
  intval.resize(nz);

  double ok, zmax;
  double z, wfac, currz, val, currint, running_sum, last_grid;

  intval.resize(nz);
  
  zmax = zcmb[nz-1]; //Assume sorted
  
  //Prepare some values that aren't worth recalculating in the inner loop
  ok = 1.0 - om - ol;
  wfac = 3.0*(1.0+w);

  //Initialize sum_running to -1/2 first point (first term in trapezoidal
  // rule)
  running_sum = - calcValW(0.0,om,ol,ok,wfac) * dz / 2.0;

  last_grid = 0.0;

  //Main sum loop
  val = 0.0;
  for (unsigned int i = 0; i < nz; ++i) {
    currz = zcmb[i];

    for (z = last_grid; z < currz; z += dz) {
      val = calcValW(z,om,ol,ok,wfac);
      running_sum += val * dz;
    }

    //z will have overstepped by dz as the exit condition
    last_grid = z - dz;

    //Subtract of the 1/2 end point for the trapezoidal rule
    currint = running_sum -  val * dz / 2.0;

    //Since we will start at the last grid point again, we need to
    // subtract off the last added value.  Otherwise we will add
    // the last grid point twice
    running_sum -= val*dz;

    //Add in the last non-grid point value evaluated at the midpoint
    // between the two
    val = calcValW( (last_grid + currz)/2.0, om, ol, ok, wfac);
    currint += val * (currz - last_grid);

    intval[i] = currint;
  }

  return 0;

}


/*!
  Uses a fairly brute-force trapezoidal algorithm.
  \param[in] sne List of SN to find luminosity distance integrals for.  Must
                  be sorted by SNeDataEntry::zcmb
  \param[in] w \f$w\f$, the equation of state parameter of the dark energy
  \param[in] om \f$\Omega_m\f$
  \param[in] ow \f$\Omega_{DE}\f$, density parameter for dark energy
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success
*/
int lumdist::lumDistSumW(const SNeData& sne, std::vector<double>& dl, 
			 double w, double om, double ow) const {

  const double pi = 3.14159265358979323846264338327950288419716939937510582;

  setZcmb( sne );
  unsigned int nsn = zcmb.size();

  int st = lumDistSumWInt( dl, w, om, ow );
  if (st) return st;

  //Apply curvature stuff
  double ok = 1.0 - om - ow;
  if (fabs(ok) < 0.001) {
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * dl[i] );
    }
  } else if (ok > 0) {
    double sqrtk = sqrt( fabs( ok ) );
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10 ((1.0 + sne[i].zhel) * sinh(sqrtk * dl[i])/sqrtk );
    }
  } else {
    double sqrtk = sqrt( fabs( ok ) );
    for (unsigned int i = 0; i < nsn; ++i)
      if (dl[i]*sqrtk >= pi) return 1;
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * sin(sqrtk * dl[i])/sqrtk );
    }
  }

  return 0;

}

/*!
  Uses GSL numeric integration to find luminosity distance.
  \param[in] SN Object with info about SN
  \param[in] w \f$w\f$, the equation of state parameter of the dark energy
  \param[in] om \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$, the density parameter of the dark energy
  \param[out] dl luminosity distance (c/H0 removed).
            Actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success
*/
int lumdist::singleLumDistW(const SNeDataEntry& SN, double w, double om, 
			    double ode, double &dl) const {

  double zhel = SN.zhel;
  double zcmb = SN.zcmb;

  //w is not -1, so we have to throw the full integral at the problem
  const double paramtol = 1e-5;

  int gsl_retcode;
  size_t neval;
  double currint, error, sqrtk;
  gsl_function F;
  double params[3];
  params[0] = w; params[1] = om; params[2] = ode;
  F.function = &calcValW_gsl;
  F.params = &params;
  
  double ok; //Omega_k
  ok = 1.0 - om - ode;

  gsl_retcode = gsl_integration_qng(&F,0,zcmb,0,prec_lum,&currint,
				    &error,&neval);
  
  if (gsl_retcode) {
    //Simple integration didn't work -- call adaptive routine
    gsl_retcode = 
      gsl_integration_qags(&F,0,zcmb,0,prec_lum,512,intspace,&currint,&error);
    if (gsl_retcode) {
      //Well, that didn't work either.  Give up.
      return 1;
    }
  }

  if (fabs(ok) < paramtol) {
    dl = 5 * log10( (1.0 + zhel)*currint );
  } else if (ok > 0) {
    sqrtk = sqrt( fabs( ok ) );
    dl = 5*log10( (1.0 + zhel) * sinh(sqrtk * currint) / sqrtk );
  } else {
    sqrtk = sqrt( fabs( ok ) );
    dl = 5*log10( (1.0 + zhel) * sin(sqrtk * currint) / sqrtk );
  }
  
  return 0;
}

/*!
  Evaluates only the integral part of the luminosity distance term
  \param[in] w \f$w\f$
  \param[in] om \f$\Omega_m\f$
  \param[in] ol \f$\Omega_{\Lambda}\f$
  \param[out] intval Filled with the integral values
  \returns 0 on success

  Uses Runge-Kutta-Felhberg 4-5 method.  Faster than trapezoidal algorithm.
*/
int lumdist::lumDistSumWInt_rk4(std::vector<double>& intval, double w,
				double om, double ol) const {
  //We can no longer test for bouncing cases because w complicates
  // this significantly
  
  //Test for cosmological constant case
  if (w == -1.0) 
    return lumDistSumInt_rk4(intval,om,ol);

  unsigned int nz = zcmb.size();
  intval.resize(nz);

  double params[3];
  params[0] = w; params[1] = om; params[2] = ol;

  gsl_odeiv_system sys = {&calcValW_gsl_ode, 0, 1, &params};

  double z = 0;
  double h = 1e-5; //Initial step guess
  if ( h < zcmb[0] ) h = zcmb[0]/2.0;
  double y[1] = { 0.0 };  //Inital value of integral

  int st;
  double zcurrtarg;
  for (unsigned int i = 0; i < nz; ++i) {
    zcurrtarg = zcmb[i];
    while (z < zcurrtarg) {
      st = gsl_odeiv_evolve_apply(e_ode, c_ode, s_ode, &sys, &z, zcurrtarg, 
				  &h, y);
      if (st != 0) return 1; //Failure
    }
    intval[i] = y[0];
  }

  return 0;
}


/*!
  Uses Runge-Kutta-Felhberg 4-5 method.  Faster than trapezoidal algorithm
  \param[in] sne List of SN to find luminosity distance integrals for.  Must
                  be sorted by SNeDataEntry::zcmb
  \param[in] w \f$w\f$, the equation of state parameter for the dark energy
  \param[in] om \f$\Omega_m\f$
  \param[in] ow \f$\Omega_{DE}\f$, the density parameter of the 
                  dark energy.
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success
*/
int lumdist::lumDistSumW_rk4(const SNeData& sne, std::vector<double>& dl, 
			     double w, double om, double ow) const {

  const double pi = 3.14159265358979323846264338327950288419716939937510582;

  setZcmb( sne );
  unsigned int nsn = zcmb.size();

  int st = lumDistSumWInt_rk4( dl, w, om, ow );
  if (st) return st;

  //Apply curvature stuff
  double ok = 1.0 - om - ow;
  if (fabs(ok) < 0.001) {
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * dl[i] );
    }
  } else if (ok > 0) {
    double sqrtk = sqrt( fabs( ok ) );
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10 ((1.0 + sne[i].zhel) * sinh(sqrtk * dl[i])/sqrtk );
    }
  } else {
    double sqrtk = sqrt( fabs( ok ) );
    double compval = pi / sqrtk;
    for (unsigned int i = 0; i < nsn; ++i) 
      if (dl[i] >= compval || isnan(dl[i]) ) return 1;
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * sin(sqrtk * dl[i])/sqrtk );
    }
  }

  return 0;
}

/*!
  Evaluates only the integral part of the luminosity distance term
  \param[in] w0 \f$w_0\f$
  \param[in] wa \f$w_a\f$
  \param[in] om \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \param[out] intval Filled with the integral values
  \returns 0 on success

  Uses a brute-force trapezoidal algorithm.
*/
int lumdist::lumDistSumW0WAInt(std::vector<double>& intval, double w0,
			       double wa, double om, double ode) const {
  //We can no longer test for bouncing cases because w complicates
  // this significantly
  
  //Test for constant w case
  if ( abs(wa) < param_tol ) return lumDistSumWInt( intval, w0, om, ode );

  unsigned int nz = zcmb.size();
  intval.resize(nz);

  double ok, zmax;
  double z, wfac, currz, val, currint, running_sum, last_grid;

  intval.resize(nz);
  
  zmax = zcmb[nz-1]; //Assume sorted
  
  //Prepare some values that aren't worth recalculating in the inner loop
  ok = 1.0 - om - ode;
  wfac = 3.0*(1.0+w0+wa);

  //Initialize sum_running to -1/2 first point (first term in trapezoidal
  // rule)
  if (useKomatsuForm) 
    running_sum = - calcValW0WAKom(0.0,om,ode,ok,w0,wa) * dz / 2.0; else
    running_sum = - calcValW0WA(0.0,om,ode,ok,wfac,wa) * dz / 2.0;

  last_grid = 0.0;

  //Main sum loop
  val = 0.0;
  for (unsigned int i = 0; i < nz; ++i) {
    currz = zcmb[i];

    if (useKomatsuForm) {
      for (z = last_grid; z < currz; z += dz) {
	val = calcValW0WAKom(z,om,ode,ok,w0,wa);
	running_sum += val * dz;
      }
    } else {
      for (z = last_grid; z < currz; z += dz) {
	val = calcValW0WA(z,om,ode,ok,wfac,wa);
	running_sum += val * dz;
      }
    }

    //z will have overstepped by dz as the exit condition
    last_grid = z - dz;

    //Subtract of the 1/2 end point for the trapezoidal rule
    currint = running_sum -  val * dz / 2.0;

    //Since we will start at the last grid point again, we need to
    // subtract off the last added value.  Otherwise we will add
    // the last grid point twice
    running_sum -= val*dz;

    //Add in the last non-grid point value evaluated at the midpoint
    // between the two
    if (useKomatsuForm)
      val = calcValW0WAKom( (last_grid + currz)/2.0, om, ode, ok, 
			    w0, wa); else
      val = calcValW0WA( (last_grid + currz)/2.0, om, ode, ok, wfac, wa);
    currint += val * (currz - last_grid);

    intval[i] = currint;
  }

  return 0;

}


/*!
  Uses a fairly brute-force trapezoidal algorithm.
  \param[in] sne List of SN to find luminosity distance integrals for.  Must
                  be sorted by SNeDataEntry::zcmb
  \param[in] w0 \f$w_0\f$
  \param[in] wa \f$w_a\f$
  \param[in] om \f$\Omega_m\f$
  \param[in] ow \f$\Omega_{DE}\f$, density parameter for dark energy
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success
*/
int lumdist::lumDistSumW0WA(const SNeData& sne, std::vector<double>& dl, 
			    double w0, double wa, double om, double ow) const {

  const double pi = 3.14159265358979323846264338327950288419716939937510582;

  setZcmb( sne );
  unsigned int nsn = zcmb.size();

  int st = lumDistSumW0WAInt( dl, w0, wa, om, ow );
  if (st) return st;

  //Apply curvature stuff
  double ok = 1.0 - om - ow;
  if (fabs(ok) < 0.001) {
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * dl[i] );
    }
  } else if (ok > 0) {
    double sqrtk = sqrt( fabs( ok ) );
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10 ((1.0 + sne[i].zhel) * sinh(sqrtk * dl[i])/sqrtk );
    }
  } else {
    double sqrtk = sqrt( fabs( ok ) );
    for (unsigned int i = 0; i < nsn; ++i)
      if (dl[i]*sqrtk >= pi) return 1;
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * sin(sqrtk * dl[i])/sqrtk );
    }
  }

  return 0;

}

/*!
  Uses GSL numeric integration to find luminosity distance.
  \param[in] SN Object with info about SN
  \param[in] w0 \f$w_0\f$
  \param[in] wa \f$w_a\f$
  \param[in] om \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$, the density parameter of the dark energy
  \param[out] dl luminosity distance (c/H0 removed).
            Actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success
*/
int lumdist::singleLumDistW0WA(const SNeDataEntry& SN, double w0, double wa,
			       double om, double ode, double &dl) const {

  double zhel = SN.zhel;
  double zcmb = SN.zcmb;

  //w is not -1, so we have to throw the full integral at the problem
  const double paramtol = 1e-5;

  int gsl_retcode;
  size_t neval;
  double currint, error, sqrtk;
  gsl_function F;
  double params[5];
  params[0] = w0; params[1] = wa; params[2] = om; params[3] = ode;
  if (useKomatsuForm) { 
    params[4] = atrans;
    F.function = &calcValW0WAKom_gsl;
  } else F.function = &calcValW0WA_gsl;
  F.params = &params;
  
  double ok; //Omega_k
  ok = 1.0 - om - ode;

  gsl_retcode = gsl_integration_qng(&F,0,zcmb,0,prec_lum,&currint,
				    &error,&neval);
  
  if (gsl_retcode) {
    //Simple integration didn't work -- call adaptive routine
    gsl_retcode = 
      gsl_integration_qags(&F,0,zcmb,0,prec_lum,512,intspace,&currint,&error);
    if (gsl_retcode) {
      //Well, that didn't work either.  Give up.
      return 1;
    }
  }

  if (fabs(ok) < paramtol) {
    dl = 5 * log10( (1.0 + zhel)*currint );
  } else if (ok > 0) {
    sqrtk = sqrt( fabs( ok ) );
    dl = 5*log10( (1.0 + zhel) * sinh(sqrtk * currint) / sqrtk );
  } else {
    sqrtk = sqrt( fabs( ok ) );
    dl = 5*log10( (1.0 + zhel) * sin(sqrtk * currint) / sqrtk );
  }
  
  return 0;
}

/*!
  Evaluates only the integral part of the luminosity distance term
  \param[in] w0 \f$w_0\f$
  \param[in] wa \f$w_a\f$
  \param[in] om \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \param[out] intval Filled with the integral values
  \returns 0 on success

  Uses Runge-Kutta-Felhberg 4-5 method.  Faster than trapezoidal algorithm,
  except for large numbers of SN (where the overheads kill you).
*/
int lumdist::lumDistSumW0WAInt_rk4(std::vector<double>& intval, double w0,
				   double wa, double om, double ode) const {
  //We can no longer test for bouncing cases because w complicates
  // this significantly
  
  //Test for constant w case
  if ( abs(wa) < param_tol ) return lumDistSumWInt_rk4( intval, w0, om, ode );

  unsigned int nz = zcmb.size();
  intval.resize(nz);

  double params[5];
  params[0] = w0; params[1] = wa; params[2] = om; params[3] = ode;

  gsl_odeiv_system sys;
  if (useKomatsuForm) {
    params[4] = atrans;
    sys.function = &calcValW0WAKom_gsl_ode;
  } else sys.function = calcValW0WA_gsl_ode;
  sys.jacobian = 0;
  sys.dimension = 1;
  sys.params = &params;

  double z = 0;
  double h = 1e-5; //Initial step guess
  if ( h < zcmb[0] ) h = zcmb[0]/2.0;
  double y[1] = { 0.0 };  //Inital value of integral

  int st;
  double zcurrtarg;
  for (unsigned int i = 0; i < nz; ++i) {
    zcurrtarg = zcmb[i];
    while (z < zcurrtarg) {
      st = gsl_odeiv_evolve_apply(e_ode, c_ode, s_ode, &sys, &z, zcurrtarg, 
				  &h, y);
      if (st != 0) return 1; //Failure
    }
    intval[i] = y[0];
  }

  return 0;
}


/*!
  Uses Runge-Kutta-Felhberg 4-5 method.  Usually faster than trapezoidal 
  algorithm
  \param[in] sne List of SN to find luminosity distance integrals for.  Must
                  be sorted by SNeDataEntry::zcmb
  \param[in] w0 \f$w_0\f$
  \param[in] wa \f$w_a\f$
  \param[in] om \f$\Omega_m\f$
  \param[in] ow \f$\Omega_{DE}\f$, the density parameter of the 
                  dark energy.
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success
*/
int lumdist::lumDistSumW0WA_rk4(const SNeData& sne, std::vector<double>& dl, 
				double w0, double wa, double om, 
				double ow) const {

  const double pi = 3.14159265358979323846264338327950288419716939937510582;

  setZcmb( sne );
  unsigned int nsn = zcmb.size();

  int st = lumDistSumW0WAInt_rk4( dl, w0, wa, om, ow );
  if (st) return st;

  //Apply curvature stuff
  double ok = 1.0 - om - ow;
  if (fabs(ok) < 0.001) {
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * dl[i] );
    }
  } else if (ok > 0) {
    double sqrtk = sqrt( fabs( ok ) );
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10 ((1.0 + sne[i].zhel) * sinh(sqrtk * dl[i])/sqrtk );
    }
  } else {
    double sqrtk = sqrt( fabs( ok ) );
    double compval = pi / sqrtk;
    for (unsigned int i = 0; i < nsn; ++i) 
      if (dl[i] >= compval || isnan(dl[i]) ) return 1;
    for (unsigned int i = 0; i < nsn; ++i) {
      dl[i] = 5*log10( (1.0 + sne[i].zhel) * sin(sqrtk * dl[i])/sqrtk );
    }
  }

  return 0;
}


/*!
  \param[in] sne List of SN to find luminosity distance integrals for.  Must
                  be sorted by SNeDataEntry::zcmb
  \param[in] om \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$, the density parameter of the 
                  dark energy.
  \param[in] w0 \f$w_0\f$, the equation of state parameter for the dark energy
                  at the current epoch
  \param[in] wa \f$w_a\f$ Minus the derivative of \f$w\f$ with respect to 
                 \f$a\f$
  \param[in] method Which integration method to use.  If the specified 
                      method is not possible, an error is returned.
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success

  The model here is \f$w \left(a\right) = w_0 + w_a \left( 1 - a \right)\f$
  (i.e., Linder (2003)).
*/
int lumdist::getLumDist( const SNeData& sne, std::vector<double>& dl,
			 double om, double ode, double w0, double wa,
			 lumdist::int_method method) const {

  //Single is most efficient for small numbers of SN, then RK4,
  // and for large enough sets the trapezoidal rule is better
  //These numbers are approximately true on my machine, but might
  // not be for others.
  const unsigned int nmaxsingle = 50;
  const unsigned int nmaxrk4 = 600; 

  //Figure out which integrals we can use
  bool islambda, isconstw;
  if ( (!useKomatsuForm) && abs(wa) < param_tol ) {
    isconstw = true;
    if ( abs( w0 + 1.0 ) < param_tol ) islambda = true; else islambda=false;
  } else {
    isconstw = false;
    islambda = false;
  }

  unsigned int nsn = sne.size();
  
  if ( method == lumdist::AUTO ) {
    //Decide what method to use based on number of SN, etc.
    if (islambda && ode >= 0.0) {
      //Single is viable
      if (nsn <nmaxsingle) {
	method = lumdist::SINGLE; 
      } else if (nsn < nmaxrk4) {
	method = lumdist::RK4;
      } else method = lumdist::TRAPEZOIDAL;
    } else {
      //Single not viable
      if (nsn < nmaxrk4) method = lumdist::RK4; else 
	method=lumdist::TRAPEZOIDAL;
    }
  } else {
    //User specified.  If SINGLE, make sure we can
    if ( method == lumdist::SINGLE && ! islambda ) return 1;
  }

  //Now do the actual calculation
  int retval;
  if (method == lumdist::SINGLE) {
    retval = 0;
    for (unsigned int i=0; i < nsn; ++i) {
      retval |= singleLumDist( sne[i], om, ode, dl[i] );
    }
    return retval;
  } else if (method == lumdist::RK4) {
    if (islambda) {
      return lumDistSum_rk4( sne, dl, om, ode );
    } else if (isconstw) {
      return lumDistSumW_rk4( sne, dl, w0, om, ode );
    } else {
      return lumDistSumW0WA_rk4( sne, dl, w0, wa, om, ode );
    }
  } else if (method == lumdist::TRAPEZOIDAL) {
    if (islambda) {
      return lumDistSum( sne, dl, om, ode );
    } else if (isconstw) {
      return lumDistSumW( sne, dl, w0, om, ode );
    } else {
      return lumDistSumW0WA( sne, dl, w0, wa, om, ode );
    }
  } else return 1; //Shouldn't happen


}

/*!
  \param[in] sne List of SN to find luminosity distance integrals for.  Must
                  be sorted by SNeDataEntry::zcmb
  \param[in] params 4 element vector of parameter values in the order
              \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \param[in] method Which integration method to use.  If the specified 
                      method is not possible, an error is returned.
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned for each SN
  \returns 0 on success

  The model here is \f$w \left(a\right) = w_0 + w_a \left( 1 - a \right)\f$
  (i.e., Linder (2003)) or the Komatsu et al. (2008) form.
*/
int lumdist::getLumDist( const SNeData& sne, std::vector<double>& dl, 
			 const std::vector<double>& params,
			 int_method method ) const {
  if ( params.size() != 4 )
    throw CosFitterExcept("lumdist","getLumDist",
			  "params vec must have 4 elements",1);
  return getLumDist( sne, dl, params[0], params[1], params[2], params[3],
		     method );
}
/*!
  \param[in] sn  Single SN to find luminosity distance for. 
  \param[in] om \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$, the density parameter of the 
                  dark energy.
  \param[in] w0 \f$w_0\f$, the equation of state parameter for the dark energy
                  at the current epoch
  \param[in] wa \f$w_a\f$ Minus the derivative of \f$w\f$ with 
                 respect to \f$a\f$
  \param[out] dl luminosity distance (c/H0 removed).  Assumed pre-allocated.
            actually, 5 * log10 (dl) is returned
  \returns 0 on success

  The model here is \f$w \left(a\right) = w_0 + w_a \left( 1 - a \right)\f$
  (i.e., Linder (2003)) or the Komatsu et al. (2008) form.
*/
int lumdist::getLumDist( const SNeDataEntry& sn, double& dl,
			 double om, double ode, double w0, double wa) const {
  //Figure out which function to call
  bool islambda, isconstw;
  //Note that the Komatsu form doesn't allow a truly constant w
  if ( abs(wa) < param_tol && (! useKomatsuForm) ) {
    isconstw = true;
    if ( abs( w0 + 1.0 ) < param_tol ) islambda = true; else islambda=false;
  } else {
    isconstw = false;
    islambda = false;
  }

  if (islambda) {
    //The elliptic integral method doesn't do well for Omega_m=0
    if ( abs(om) < param_tol ) return singleLumDistW( sn, w0, om, ode, dl );
    else return singleLumDist( sn, om, ode, dl );
  }
  if (isconstw) return singleLumDistW( sn, w0, om, ode, dl );
  return singleLumDistW0WA( sn, w0, wa, om, ode, dl );
}
