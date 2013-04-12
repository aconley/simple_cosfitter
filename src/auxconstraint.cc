
#include<string>
#include<sstream>

#ifdef USEMKL
#include<math.h>
#else
#include<cmath>
#endif

#include "param_tags.h"
#include "auxconstraint.h"
#include "cosfitterexcept.h"

#include <gsl/gsl_errno.h>

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
  (\f$\Omega_m\f$,\f$\Omega_{\Lambda}\f$)
  \returns 1/E term used in many cosmological expressions
*/
double auxconstraint::one_over_E( double z, void *params ) {
  double opz, om, ol, ok;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + static_cast<double>(z);
  om = in_params[0];
  ol = in_params[1];
  ok = 1.0 - om - ol;
  return 1.0/sqrt( (om*opz + ok)*opz*opz + ol );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
  (\f$w\f$, \f$\Omega_m\f$,\f$\Omega_{DE}\f$)
  \returns 1/E term used in many cosmological expressions
*/
double auxconstraint::one_over_E_w( double z, void *params ) {
  double opz, w, om, ode, ok;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + static_cast<double>(z);
  w = in_params[0];
  om = in_params[1];
  ode = in_params[2];
  ok = 1.0 - om - ode;
  return 1.0/sqrt( (om*opz + ok)*opz*opz + ode*pow(opz,3*(1+w)) );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
  (\f$w_0\f$, \f$w_a\f$ \f$\Omega_m\f$,\f$\Omega_{DE}\f$)
  \returns 1/E term used in many cosmological expressions
*/
double auxconstraint::one_over_E_w0wa( double z, void *params ) {
  double opz, w0, wa, om, ode, ok;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + static_cast<double>(z);
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  ok = 1.0 - om - ode;
  return 1.0/sqrt( (om*opz + ok)*opz*opz + 
		   ode*pow(opz,3*(1+w0+wa))*exp(-3 * wa * z/opz) );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
  (\f$w_0\f$, \f$w_a\f$ \f$\Omega_m\f$,\f$\Omega_{DE}, atrans\f$)
  \returns 1/E term used in many cosmological expressions

  This uses the Komatsu et al. (2008) form for the expansion history
*/
double auxconstraint::one_over_E_w0waKom( double z, void *params ) {
  double opz, w0, wa, om, ode, ok, atrans, val;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + static_cast<double>(z);
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  atrans = in_params[4];
  ok = 1.0 - om - ode;
  double opat = 1.0 + atrans;

  val = opz*opz * ( opz*om + ok ) + 
    ode * exp( -3.0*( z*wa/opz + (1+w0+opat*wa)*
		      log( (1/opz + atrans)/opat ) ) );
  return 1./sqrt(val);
}


/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
  (\f$\Omega_m\f$,\f$\Omega_{\Lambda}, \Omega_r\f$)
  \returns 1/E term used in many cosmological expressions
  This version includes \f$\Omega_r\f$
*/
double auxconstraint::one_over_E_rad( double z, void *params ) {
  double opz, om, ol, orad, ok;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + static_cast<double>(z);
  om = in_params[0];
  ol = in_params[1];
  orad = in_params[2];
  ok = 1.0 - om - ol - orad;
  return 1.0/sqrt( ((orad*opz+om)*opz + ok)*opz*opz + ol );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
  (\f$w\f$, \f$\Omega_m\f$,\f$\Omega_{DE}\f$,\f$\Omega_r\f$)
  \returns 1/E term used in many cosmological expressions
  This version includes \f$\Omega_r\f$
*/
double auxconstraint::one_over_E_w_rad( double z, void *params ) {
  double opz, w, om, ode, ok, orad;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + static_cast<double>(z);
  w = in_params[0];
  om = in_params[1];
  ode = in_params[2];
  orad = in_params[3];
  ok = 1.0 - om - ode - orad;
  return 1.0/sqrt( ((orad*opz+om)*opz + ok)*opz*opz + ode*pow(opz,3*(1+w)) );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
  (\f$w_0\f$, \f$w_a\f$ \f$\Omega_m\f$,\f$\Omega_{DE}\f$, \f$\Omega_r\f$)
  \returns 1/E term used in many cosmological expressions
  This version includes \f$\Omega_r\f$
*/
double auxconstraint::one_over_E_w0wa_rad( double z, void *params ) {
  double opz, w0, wa, om, ode, ok, orad;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + static_cast<double>(z);
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  orad = in_params[4];
  ok = 1.0 - om - ode - orad;
  return 1.0/sqrt( ((orad*opz+om)*opz + ok)*opz*opz + 
		   ode*pow(opz,3*(1+w0+wa))*exp(-3 * wa * z/opz) );
}

/*!
  This is in the form for the GSL

  \param[in] z Redshift
  \param[in] params Cosmological parameters 
  (\f$w_0\f$, \f$w_a\f$ \f$\Omega_m\f$,\f$\Omega_{DE}\f$, \f$\Omega_r, 
  atrans\f$)
  \returns 1/E term used in many cosmological expressions
  This version includes \f$\Omega_r\f$, and is in the Komatsu form
  for w(a)
*/
double auxconstraint::one_over_E_w0waKom_rad( double z, void *params ) {
  double opz, w0, wa, om, ode, ok, orad, atrans,val;
  double *in_params = static_cast<double*>(params);
  opz = 1.0 + static_cast<double>(z);
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  orad = in_params[4];
  atrans = in_params[5];
  ok = 1.0 - om - ode - orad;
  double opat = 1.0+atrans;
  val = ( (orad*opz + om)*opz + ok) * opz*opz + 
    ode*exp( -3.0*( z*wa/opz + (1+w0+opat*wa)*log( (1/opz + atrans)/opat ) ) );
  return 1./sqrt(val);
}


/*!
  This is in the form for the GSL

  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
  (\f$\Omega_m , \Omega_{\Lambda} , \Omega_b h^2 , h\f$)
  \returns The integral term for the \f$r_s \left( z_{cmb} \right)\f$ equation
*/
double auxconstraint::rsint( double a, void *params ) {
  double om, ol, ok, obh2, h;
  double *in_params = static_cast<double*>(params);

  const double zfac = 26053.28054114; //2.5e4 * (tcmb/2.7)^-4, tcmb=2.728
  const double rfac = 32827.1334818364; //31500 * (tcmb/2.7)^-4

  om = in_params[0];
  ol = in_params[1];
  ok = 1.0 - om - ol;
  obh2 = in_params[2];
  h = in_params[3];

  double aeq = 1.0 / (1 + zfac * om * h * h); //Matter, radiation equality
  double avRb = rfac * obh2;

  double a4e2 = om * (aeq + a) + (ol * a * a + ok) * a * a;

  return 1.0/sqrt( 3.0 * (1 + avRb*a) * a4e2 );
}

/*!
  This is in the form for the GSL

  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
  (\f$w, \Omega_m , \Omega_{DE} , \Omega_b h^2 , h\f$)
  \returns The integral term for the \f$r_s \left( z_{cmb} \right)\f$ equation
*/
double auxconstraint::rsint_w( double a, void *params ) {
  double w, om, ode, ok, obh2, h;
  double *in_params = static_cast<double*>(params);

  const double zfac = 26053.28054114; //2.5e4 * (tcmb/2.7)^-4, tcmb=2.728
  const double rfac = 32827.1334818364; //31500 * (tcmb/2.7)^-4
  
  w  = in_params[0];
  om = in_params[1];
  ode = in_params[2];
  ok = 1.0 - om - ode;
  obh2 = in_params[3];
  h = in_params[4];

  double aeq = 1.0 / (1 + zfac * om * h * h); //Matter, radiation equality
  double avRb = rfac * obh2;

  double a4e2 = om * (aeq + a) + ok * a * a + ode * pow( a, 1 - 3 * w );

  return 1.0/sqrt( 3.0 * (1 + avRb*a) * a4e2 );
}

/*!
  This is in the form for the GSL

  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
  (\f$w_0 , w_a, \Omega_m , \Omega_{DE} , \Omega_b h^2 , h\f$)
  \returns The integral term for the \f$r_s \left( z_{cmb} \right)\f$ equation
*/
double auxconstraint::rsint_w0wa( double a, void *params ) {
  double w0, wa, om, ode, ok, obh2, h;
  double *in_params = static_cast<double*>(params);

  const double zfac = 26053.28054114; //2.5e4 * (tcmb/2.7)^-4, tcmb=2.728
  const double rfac = 32827.1334818364; //31500 * (tcmb/2.7)^-4
  
  w0  = in_params[0];
  wa  = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  ok = 1.0 - om - ode;
  obh2 = in_params[4];
  h = in_params[5];

  double aeq = 1.0 / (1 + zfac * om * h * h); //Matter, radiation equality
  double avRb = rfac * obh2;

  double a4e2 = om * (aeq + a) + ok * a * a + 
    ode * pow( a, 1 - 3 * (w0 + wa) ) * exp( -3 * wa * ( 1.0 - a ) );

  return 1.0/sqrt( 3.0 * (1 + avRb*a) * a4e2 );
}

/*!
  This is in the form for the GSL

  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
  (\f$w_0 , w_a, \Omega_m , \Omega_{DE} , \Omega_b h^2 , h, atrans\f$)
  \returns The integral term for the \f$r_s \left( z_{cmb} \right)\f$ equation

  Uses the Komatsu et al. (2008) form for w(a)
*/
double auxconstraint::rsint_w0waKom( double a, void *params ) {
  double w0, wa, om, ode, ok, obh2, h, atrans;
  double *in_params = static_cast<double*>(params);

  const double zfac = 26053.28054114; //2.5e4 * (tcmb/2.7)^-4, tcmb=2.728
  const double rfac = 32827.1334818364; //31500 * (tcmb/2.7)^-4
  
  w0  = in_params[0];
  wa  = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  ok = 1.0 - om - ode;
  obh2 = in_params[4];
  h = in_params[5];
  atrans = in_params[6];
  double opat = 1.0+atrans;
  double aeq = 1.0 / (1 + zfac * om * h * h); //Matter, radiation equality
  double avRb = rfac * obh2;
  double a4e2 = om * (aeq + a) + ok * a * a + 
    ode * exp( -3.0*( (1-a)*wa + (1+w0+opat*wa)*log( (a + atrans)/opat )));

  return 1.0/sqrt( 3.0 * (1 + avRb*a) * a4e2 );
}

/*!
  This is in the form for the GSL
  
  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
  (\f$\Omega_m , \Omega_{\Lambda} , \Omega_b h^2, h \f$)
  \returns The integral term for \f$r_s H_0 / c\f$  Does not include the
  \f$\sqrt(3)\f$ term, and is done in terms of \f$a\f$
*/
double auxconstraint::rsHovercint( double a, void *params ) {
  double om, ode, obh2, h, orad, ok;
  double *in_params = static_cast<double*>(params);

  om = in_params[0];
  ode = in_params[1];
  obh2 = in_params[2];
  h = in_params[3];

  orad = oradh2 / (h*h); 

  ok = 1.0 - om - ode - orad;

  //a^2 * H/H0
  double a2E = sqrt( orad + a * (om + a*(ok + a*a*ode) ) );
  
  return 1.0 / ( a2E*sqrt(1.0 + 0.75*a*(obh2/ogammah2)) );
}

/*!
  This is in the form for the GSL
  
  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
  (\f$w, Omega_m , \Omega_{\Lambda} , \Omega_b h^2, h \f$)
  \returns The integral term for \f$r_s H_0 / c\f$  Does not include the
  \f$\sqrt(3)\f$ term, and is done in terms of \f$a\f$
*/
double auxconstraint::rsHovercint_w( double a, void *params ) {
  double w, om, ode, obh2, h, orad, ok;
  double *in_params = static_cast<double*>(params);
  
  w = in_params[0];
  om = in_params[1];
  ode = in_params[2];
  obh2 = in_params[3];
  h = in_params[4];

  orad = oradh2 / (h*h); 

  ok = 1.0 - om - ode - orad;

  //a^2 * H/H0
  double a2E = sqrt( orad + a * (om + a*ok + ode*pow(a,-3.0*w)) );
  
  return 1.0 / ( a2E*sqrt(1.0 + 0.75*a*(obh2/ogammah2)) );
}


/*!
  This is in the form for the GSL
  
  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
  (\f$w_0, w_a, Omega_m , \Omega_{\Lambda} , \Omega_b h^2, h \f$)
  \returns The integral term for \f$r_s H_0 / c\f$  Does not include the
  \f$\sqrt(3)\f$ term, and is done in terms of \f$a\f$
*/
double auxconstraint::rsHovercint_w0wa( double a, void *params ) {
  double w0, wa, om, ode, obh2, h, orad, ok;
  double *in_params = static_cast<double*>(params);
  
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  obh2 = in_params[4];
  h = in_params[5];

  orad = oradh2 / (h*h); 

  ok = 1.0 - om - ode - orad;

  //a^2 * H/H0
  double a2E = sqrt( orad + a * (om + a*ok + 
				ode*pow(a,-3.0*(w0+wa))*exp(-3.0*wa*(1-a)) ) );
  
  return 1.0 / ( a2E*sqrt(1.0 + 0.75*a*(obh2/ogammah2)) );
}


/*!
  This is in the form for the GSL
  
  \param[in] a Scale factor
  \param[in] params Cosmological parameters 
  (\f$w_0, w_a, Omega_m , \Omega_{\Lambda} , \Omega_b h^2, h, atrans \f$)
  \returns The integral term for \f$r_s H_0 / c\f$  Does not include the
  \f$\sqrt(3)\f$ term, and is done in terms of \f$a\f$

  Uses the Komatsu et al. (2008) form for w(a)
*/
double auxconstraint::rsHovercint_w0waKom( double a, void *params ) {
  double w0, wa, om, ode, obh2, h, orad, ok, atrans, opat;
  double *in_params = static_cast<double*>(params);
  
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  obh2 = in_params[4];
  h = in_params[5];
  atrans = in_params[6];

  opat = 1.0 + atrans;

  orad = oradh2 / (h*h); 

  ok = 1.0 - om - ode - orad;

  //a^2 * H/H0
  double val = ode*exp( -3.0*( (1.0-a)*wa + (1+w0+opat*wa) * 
			      log( (a + atrans)/opat ) ) );
  double a2E = sqrt( orad + a*(om + a*( ok + a*a*val ) ) );

  return 1.0 / ( a2E*sqrt(1.0 + 0.75*a*(obh2/ogammah2)) );
}


/*!
  Based on formulae from
  <A HREF="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1992ARA%26A..30..499C&amp;db_key=AST&amp;data_type=HTML&amp;format=&amp;high=4270f9326909263">Carroll and Press, AARA 1993, 30: 499.</A>
  Ignores \f$\Omega_{rad}\f$.  This is only valid for
  the cosmological constant case (\f$w=-1\f$).
*/
bool auxconstraint::isRebounding(double om, double ol) {

  const double minom = 1e-4; //Smaller than this, use the small \om expansion
  if (om < minom) {
    return ol > 1.0 - om;
  }

  //Otherwise do the full calculation
  double arg = (1.0 - om)/om;
  if ( om >= 0.5 ) {
    return ol >= 4*om*pow( cos( acos( arg ) / 3.0 ), 3);
  } else {
    // cmath doesn't include inverse hyperbolic functions, so this
    // looks funny, but cosh^{-1}(z) = ln( z + sqrt(z+1)*sqrt(z-1) )
    return ol >= 4*om*pow( cosh( log( arg + sqrt(arg*arg-1) )/ 3.0),3);
  }
}

/*!
  \param[in] v Vector of \f$h, \Omega_b h^2\f$
  \param[in] params Other parameters, very indirectly
  \returns \f$\chi^2\f$
*/
double auxconstraint::percival_chisq(const gsl_vector *v, void *params) {

  const unsigned int nz = 2;
  const double percival_z[nz] = { 0.2, 0.35 };
  const double percival_params[nz] = { 0.1905, 0.1097 };
  const double badval = 1e50;

  double h, h2, obh2;
  h = gsl_vector_get(v,0);
  obh2 = gsl_vector_get(v,1);
  h2 = h*h;
  
  void **vptr = static_cast<void**>(params);

  distance_helper *const dh = static_cast<distance_helper *const>(vptr[0]);
  
  double* in_params = static_cast<double*>(vptr[1]);
  double om, ode, w0, wa;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  
  double omh2;
  omh2 = om * h2;

  double orad;
  orad = oradh2 / h2;

  //Bad value test
  if ( om < 0 || obh2 < 0 || h < 0 || obh2 / h2 > om) 
    return badval;

  double zdrag;
  zdrag = dh->GetZdrag( obh2, omh2 );

  //Get H0 rs(zd)/c
  double K1, K2_0, K2_1;
  bool success_K1, success_K2;
  K1=dh->GetK1( zdrag, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1) return badval; //Failure

  //We need two values of K2, one for each z
  K2_0=dh->GetK2( percival_z[0], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_0 == 0) return badval; //Failure
  K2_1=dh->GetK2( percival_z[1], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_1 == 0) return badval; //Failure


  //Really, we hold H0/c CV = [K2^2 z/E]^1/3
  double DV_0, DV_1;
  bool success_DV;
  DV_0 = pow( K2_0 * K2_0 * percival_z[0] * 
	      dh->GetOneOverE( percival_z[0], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_0 == 0 ) return badval;  //Failure
  DV_1 = pow( K2_1 * K2_1 * percival_z[1] * 
	      dh->GetOneOverE( percival_z[1], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_1 == 0 ) return badval;  //Failure
  
  double delta_params[2]; //Difference between values and Percival fits
  delta_params[0] = K1/DV_0 - percival_params[0];
  delta_params[1] = K1/DV_1 - percival_params[1];

  double invcov[2][2];
  invcov[0][0] = 30124.0;
  invcov[0][1] = invcov[1][0] = -17227.0;
  invcov[1][1] = 86977.0;

  double chisq;
  chisq = 0.0;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      chisq += delta_params[i] * invcov[i][j] * delta_params[j];

  return chisq;
}


/*!
  \param[in] v Vector of \f$h, \Omega_b h^2\f$
  \param[in] params Other parameters, very indirectly
  \returns \f$\chi^2\f$
*/
double auxconstraint::wmap7_shift_chisq(const gsl_vector *v, void *params) {

  const double pi = 3.14159265358979323846264338327950288419716939937510582;
  const double wmap7_params[3] = { 302.09, 1.725, 1091.3 };
  const double badval = 1e50;

  double h, h2, obh2;
  h = gsl_vector_get(v,0);
  obh2 = gsl_vector_get(v,1);
  h2 = h*h;
  
  void **vptr = static_cast<void**>(params);

  distance_helper *const dh = static_cast<distance_helper *const>(vptr[0]);
  
  double* in_params = static_cast<double*>(vptr[1]);
  double om, ode, w0, wa;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  
  double omh2;
  omh2 = om * h2;

  double orad;
  orad = oradh2 / h2;

  //Bad value test
  if ( om < 0 || obh2 < 0 || h < 0 || obh2 / h2 > om) 
    return badval;

  double zs;
  zs = dh->GetZstar( obh2, omh2 );

  double K1, K2;
  bool success_K1, success_K2;
  K1=dh->GetK1( zs, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1) return badval; //Failure
  if (K1 == 0) return badval; //Failure
  K2=dh->GetK2( zs, om, ode, orad, w0, wa, success_K2 );
  if (! success_K2) return badval; //Failure

  double la, R;
  la = pi * K2 / K1;
  R  = sqrt( om )*K2;

  double delta_params[3]; //Between params and WMAP fits
  delta_params[0] = la - wmap7_params[0];
  delta_params[1] = R - wmap7_params[1];
  delta_params[2] = zs - wmap7_params[2];

  //Table 10 of Komatsu et al.
  double invcov[3][3];
  invcov[0][0] = 2.305;
  invcov[0][1] = invcov[1][0] = 29.698;
  invcov[0][2] = invcov[2][0] = -1.333;
  invcov[1][1] = 6825.270;
  invcov[1][2] = invcov[2][1] = -113.180;
  invcov[2][2] = 3.414;

  double chisq;
  chisq = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      chisq += delta_params[i] * invcov[i][j] * delta_params[j];

  return chisq;
}


/*!
  \param[in] v Vector of \f$h, \Omega_b h^2\f$
  \param[in] params Other parameters, very indirectly
  \returns \f$\chi^2\f$ of Percival BAO and WMAP7
*/
double auxconstraint::percival_wmap7_chisq(const gsl_vector *v, void *params) {

  const unsigned int nz = 2;
  const double percival_z[nz] = { 0.2, 0.35 };
  const double percival_params[nz] = { 0.1905, 0.1097 };
  const double pi = 3.14159265358979323846264338327950288419716939937510582;
  const double wmap7_params[3] = { 302.09, 1.725, 1091.3 };
  const double badval = 1e50;

  double h, h2, obh2;
  h = gsl_vector_get(v,0);
  obh2 = gsl_vector_get(v,1);
  h2 = h*h;
  
  void **vptr = static_cast<void**>(params);

  distance_helper *const dh = static_cast<distance_helper *const>(vptr[0]);
  
  double* in_params = static_cast<double*>(vptr[1]);
  double om, ode, w0, wa;
  w0 = in_params[0];
  wa = in_params[1];
  om = in_params[2];
  ode = in_params[3];
  
  double omh2;
  omh2 = om * h2;

  double orad;
  orad = oradh2 / h2;

  //Bad value test
  if ( om < 0 || obh2 < 0 || h < 0 || obh2 / h2 > om) 
    return badval;

  double zdrag, zstar;
  zdrag = dh->GetZdrag( obh2, omh2 );
  zstar = dh->GetZstar( obh2, omh2 );

  //Get H0 rs(zd)/c
  double K1_wmap, K1_bao, K2_wmap, K2_bao_0, K2_bao_1;
  bool success_K1, success_K2;
  K1_bao=dh->GetK1( zdrag, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1) return badval; //Failure
  K1_wmap=dh->GetK1( zstar, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1 || K1_wmap == 0.0) return badval; //Failure

  K2_wmap=dh->GetK2( zstar, om, ode, orad, w0, wa, success_K2 );
  if (! success_K2) return badval; //Failure
  K2_bao_0=dh->GetK2( percival_z[0], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_bao_0 == 0) return badval; //Failure
  K2_bao_1=dh->GetK2( percival_z[1], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_bao_1 == 0) return badval; //Failure

  //WMAP shift params
  double la, R;
  la = pi * K2_wmap / K1_wmap;
  R  = sqrt( om )*K2_wmap;

  //Really, we hold H0/c CV = [K2^2 z/E]^1/3
  double DV_0, DV_1;
  bool success_DV;
  DV_0 = pow( K2_bao_0 * K2_bao_0 * percival_z[0] * 
	      dh->GetOneOverE( percival_z[0], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_0 == 0 ) return badval;  //Failure
  DV_1 = pow( K2_bao_1 * K2_bao_1 * percival_z[1] * 
	      dh->GetOneOverE( percival_z[1], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_1 == 0 ) return badval;  //Failure

  double delta_params_bao[2]; //Difference between values and Percival fits
  delta_params_bao[0] = K1_bao/DV_0 - percival_params[0];
  delta_params_bao[1] = K1_bao/DV_1 - percival_params[1];

  double invcov_bao[2][2];
  invcov_bao[0][0] = 30124.0;
  invcov_bao[0][1] = invcov_bao[1][0] = -17227.0;
  invcov_bao[1][1] = 86977.0;

  double chisq;
  chisq = 0.0;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 2; ++j)
      chisq += delta_params_bao[i] * invcov_bao[i][j] * delta_params_bao[j];

  double delta_params_wmap[3]; //Between params and WMAP fits
  delta_params_wmap[0] = la - wmap7_params[0];
  delta_params_wmap[1] = R - wmap7_params[1];
  delta_params_wmap[2] = zstar - wmap7_params[2];

  //Table 10 of comatsu
  double invcov_wmap[3][3];
  invcov_wmap[0][0] = 2.305;
  invcov_wmap[0][1] = invcov_wmap[1][0] = 29.698;
  invcov_wmap[0][2] = invcov_wmap[2][0] = -1.333;
  invcov_wmap[1][1] = 6825.270;
  invcov_wmap[1][2] = invcov_wmap[2][1] = -113.180;
  invcov_wmap[2][2] = 3.414;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      chisq += delta_params_wmap[i] * invcov_wmap[i][j] * 
	delta_params_wmap[j];

  return chisq;
}



////////////////////////////////////////////////////////////////
//                   Distance Helper
////////////////////////////////////////////////////////////////
auxconstraint::distance_helper::distance_helper(bool fcurv, double ok) :
  useKomatsuForm(false), atrans(0.1) {
  work = gsl_integration_workspace_alloc(1000); 
}

auxconstraint::distance_helper::~distance_helper() {
  gsl_integration_workspace_free(work);
}

/*!
  \param[in] z Redshift
  \param[in] obh2 \f$ \Omega_b h^2 \f$
  \param[in] om \f$ \Omega_m \f$
  \param[in] ode \f$ \Omega_{DE}\f$
  \param[in] w0 \f$w_{0}\f$
  \param[in] wa \f$w_{a}\f$
  \param[in] h \f$H_{0} = 100 h\f$ km/sec
  \param[out] success 1 if computation was successful, 0 if failed
  \returns \f$K_1 = H_{0} r_{s} \left( z \right) / c\f$

*/
double
auxconstraint::distance_helper::GetK1( double z, double obh2, double om,
				       double ode, double w0, double wa, 
				       double h, bool& success ) const {
  const double doubletol = 1e-5; //Used for various comso parameter tolerances
  const double prefac = 1.0 / sqrt(3.0);
  const double badval = 1e30;

  success = 0;

  double a = 1.0 / (1.0 + z);

  //Quick exit tests for bad inputs (bouncing, etc.)
  bool constw = fabs(wa) < doubletol;
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return badval;
  if (om < 0.0) return badval;

  //Set up integral
  gsl_function F;
  if ( cosmoconst ) {
    double params[4];
    params[0] = om;
    params[1] = ode;
    params[2] = obh2;
    params[3] = h;
    F.function = rsHovercint;
    F.params = &params;
  } else if (constw) {
    double params[5];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    params[3] = obh2;
    params[4] = h;
    F.function = rsHovercint_w;
    F.params = &params;
  } else if (useKomatsuForm) {
    //w(a), Komatsu form
    double params[7];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = obh2;
    params[5] = h;
    params[6] = atrans;
    F.function = rsHovercint_w0waKom;
    F.params = &params;
  } else {
    //w(a), Linder form
    double params[6];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = obh2;
    params[5] = h;
    F.function = rsHovercint_w0wa;
    F.params = &params;
  }
  
  //This integral is more difficult than the BAO one, 
  // and requires adaptive integration
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  double intval, error;
  int st = gsl_integration_qag(&F,0,a,0,1e-7,1000,6,work,&intval,&error);
  gsl_set_error_handler( old_handler );

  //Failure if st != 0
  if (st) return badval;

  success = 1;
  return prefac * intval;
}


/*!
  \param[in] z Redshift
  \param[in] om \f$ \Omega_m \f$
  \param[in] ode \f$ \Omega_{DE}\f$
  \param[in] orad \f$ \Omega_r\f$, the density of radiation (not just photons)
  \param[in] w0 \f$w_{0}\f$
  \param[in] wa \f$w_{a}\f$
  \param[out] success 1 if computation was successful, 0 if failed
  \returns \f$K_2 = \left(1 + z \right) H_{0} 
  D_{A} \left( z_{\star} \right) / c\f$, where \f$D_{A}\f$ is
  the comoving angular diameter distance
*/
double
auxconstraint::distance_helper::GetK2( double z, double om,
				       double ode, double orad, double w0, 
				       double wa, bool& success ) const {
  const double doubletol = 1e-5; //Used for various comso parameter tolerance
  const double pi = 3.14159265358979323846264338327950288419716939937510582;
  const double badval = 1e30;

  success = 0;
  //Quick exit tests for bad inputs (bouncing, etc.)
  bool constw = fabs(wa) < doubletol;
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return badval;
  if (om < 0.0) return badval;

  //Set up integral
  gsl_function F;
  if ( cosmoconst ) {
    double params[3];
    params[0] = om;
    params[1] = ode;
    params[2] = orad;
    F.function = one_over_E_rad;
    F.params = &params;
  } else if (constw) {
    double params[4];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    params[3] = orad;
    F.function = one_over_E_w_rad;
    F.params = &params;
  } else if (useKomatsuForm) {
    //w(a), Komatsu Form
    double params[6];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = orad;
    params[5] = atrans;
    F.function = one_over_E_w0waKom_rad;
    F.params = &params;
  } else {
    //w(a), Linder form
    double params[5];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = orad;
    F.function = one_over_E_w0wa_rad;
    F.params = &params;
  }
  
  //This integral is more difficult than the BAO one, 
  // and requires adaptive integration
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  double intval, error;
  int st = gsl_integration_qag(&F,0,z,0,1e-7,1000,6,work,&intval,&error);
  gsl_set_error_handler( old_handler );

  //Failure if st != 0
  if (st) return badval;
#if USEMKL
  if ( isnan(intval) ) return badval;
#else
  if (std::isnan(intval) ) return badval;
#endif

  double sqrtk, ok;
  double finalval;
  ok = 1.0 - om - ode - orad;
  if (fabs(ok) < doubletol) {
    finalval = intval;
  } else if (ok > 0.0) {
    //Open
    sqrtk = sqrt( fabs( ok ) );
    finalval = sinh( sqrtk * intval )/sqrtk;
  } else {
    //Closed
    sqrtk = sqrt( fabs( ok ) );

    //proper distance is zero or negative at this redshift
    double compval = pi / sqrtk;
    if (intval >= compval) return badval;

    finalval = sin( sqrtk * intval )/sqrtk;
  }

  success = 1;
  return finalval;
}

/*!
  \param[in] z Redshift
  \param[in] om \f$ \Omega_m \f$
  \param[in] ode \f$ \Omega_{DE}\f$
  \param[in] orad \f$ \Omega_r\f$, the density of radiation (not just photons)
  \param[in] w0 \f$w_{0}\f$
  \param[in] wa \f$w_{a}\f$
  \param[out] success 1 if computation was successful, 0 if failed
  \returns \f$1/E = H_{0} / H\left( z \right)\f$
*/
double
auxconstraint::distance_helper::GetOneOverE( double z, double om,
					     double ode, double orad, 
					     double w0, double wa, 
					     bool& success ) const {
  const double doubletol = 1e-5; //Used for various comso parameter tolerance
  const double badval = 1e30;

  success = 0;
  //Quick exit tests for bad inputs (bouncing, etc.)
  bool constw = fabs(wa) < doubletol;
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return badval;
  if (om < 0.0) return badval;

  //Set up integral
  if ( cosmoconst ) {
    double params[3];
    params[0] = om;
    params[1] = ode;
    params[2] = orad;
    success = 1;
    return one_over_E_rad( z, params );
  } else if (constw) {
    double params[4];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    params[3] = orad;
    success = 1;
    return one_over_E_w_rad( z, params );
  } else if (useKomatsuForm) {
    //w(a), Komatsu Form
    double params[6];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = orad;
    params[5] = atrans;
    success = 1;
    return one_over_E_w0waKom_rad( z, params );
  } else {
    //w(a), Linder form
    double params[5];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = orad;
    success = 1;
    return one_over_E_w0wa_rad( z, params );
  }

  return badval;
}



/*!
  \param[in] obh2  \f$\Omega_b h^2\f$
  \param[in] omh2  \f$\Omega_m h^2\f$
  \returns The value of \f$z_{\star}\f$, the decopuling redshift

  Fitting forumula from Hu & Sugiyama (1996)
*/
double
auxconstraint::distance_helper::GetZstar( double obh2, double omh2 ) const {
  double g1 = 0.0783*pow(obh2,-0.238) / ( 1.0 + 39.5*pow(obh2,0.763) );
  double g2 = 0.560 / ( 1.0 + 21.1 * pow(obh2,1.81)  );
  return 1048.0*(1 + 0.00124*pow(obh2, -0.738))*(1+g1*pow(omh2,g2));
}

/*!
  \param[in] obh2  \f$\Omega_b h^2\f$
  \param[in] omh2  \f$\Omega_m h^2\f$
  \returns The value of \f$z_{d}\f$, the drag redshift

  Fitting formula from Eisenstein & Hu (1998)
*/
double
auxconstraint::distance_helper::GetZdrag( double obh2, double omh2 ) const {
  double b1 = 0.313*pow(omh2,-0.419)*(1.0 + 0.607*pow(omh2,0.674));
  double b2 = 0.238*pow(omh2,0.223);
  return 1291.0*pow(omh2,0.251) * (1.0 + b1*pow(obh2,b2)) / 
    (1.0+0.659*pow(omh2,0.828));
}

////////////////////////////////////////////////////////////////
//                    base_auxconstraint
////////////////////////////////////////////////////////////////
auxconstraint::base_auxconstraint::base_auxconstraint(bool fcurv, double ok) :
  fixcurv(fcurv), ocurv(ok) { }

auxconstraint::base_auxconstraint::~base_auxconstraint() {}

/*!
  Little utility function to update a vector with values
  of the cosmological parameters in the order that getLumDist
  understands.

  \param[in] code The parameter code from param_tags
  \param[out] vec The vector to update
  \param[in] val The value to insert

  vec is in the order \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$.
*/
void 
auxconstraint::base_auxconstraint::updateStandardCosVec( param_tags::paramcodes code, 
							 std::vector<double>& vec,
							 double val) const {
  if ( vec.size() != 4 ) vec.resize(4);
  switch ( code ) {
  case param_tags::omegam :
    vec[0] = val;
    if (fixcurv) vec[1] = 1.0 - ocurv - val;
    break;
  case param_tags::omegade :
    vec[1] = val;
    break;
  case param_tags::w0 :
    vec[2] = val;
    break;
  case param_tags::wa :
    vec[3] = val;
    break;
  default :
    //Do nothing -- don't know parameter
    break;
  }
}


void
auxconstraint::base_auxconstraint::makeProbSurface( const std::map< param_tags::paramcodes,
						    param_struct >& params,
						    const std::string& fname) {
  
  //Figure out how many cosmological params we are dealing with
  unsigned int cosmo_fit_dimension = 0;
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for ( it = params.begin(); it != params.end(); ++it)
    if (it->second.param_spec.second == param_tags::cosmological && 
	it->second.fit == param_struct::loop ) ++cosmo_fit_dimension;

  switch ( cosmo_fit_dimension ) {
  case 1 :
    makeProbSurface1D( params, fname );
    break;
  case 2 :
    makeProbSurface2D( params, fname );
    break;
  case 3 :
    makeProbSurface3D( params, fname );
    break;
  default :
    std::stringstream errstrng;
    errstrng << "Currently unsupported number of cosmological params: " 
	     << cosmo_fit_dimension;
    throw CosFitterExcept("auxconstraing::base_auxconstraint",
			  "MakeProbSurface",
			  errstrng.str(),1);
    break;
  }
}

void 
auxconstraint::base_auxconstraint::makeProbSurface1D(const std::map< param_tags::paramcodes,
						     param_struct >& params,
						     const std::string& fname,
						     bool binaryout ) {
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

  //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for (it = params.begin(); it != params.end(); ++it) 
    if ( it->second.fit == param_struct::fixed )
      updateStandardCosVec( it->first, current_cosparams, it->second.fixval );


  //Now point it at the var we will loop over, in this case only
  // omega_m
  it = params.find( param_tags::omegam );
  if ( it == params.end() )
    throw CosFitterExcept("auxconstraint::base_auxconstraint",
			  "makeProbSurface1D",
			  "Only supported 1D is Omega_m, flat", 1);
  if ( !fixcurv )
    throw CosFitterExcept("auxconstraint::base_auxconstraint",
			  "makeProbSurface1D",
			  "Only supported 1D is flat", 1);

  cosgrid1D prob1D( it->second );

  param_tags::paramcodes code0 = prob1D.getAxisSpec().first;
  unsigned int n0 = prob1D.getAxisN();
  double cos0;
  for (unsigned int i = 0; i < n0; ++i) {
    cos0 = prob1D.getAxisVal(i);
    updateStandardCosVec( code0, current_cosparams, cos0 );
    prob1D[i] = exp( -0.5 * GetChiSq( current_cosparams ) );
  }
  prob1D.normalize();
  prob1D.writeFile( fname, binaryout );
}

void 
auxconstraint::base_auxconstraint::makeProbSurface2D(const std::map< param_tags::paramcodes,
						     param_struct >& params,
						     const std::string& fname,
						     bool binaryout ) {

  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

  std::map< param_tags::paramcodes, param_struct >::const_iterator it1, it2;
  //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  for (it1 = params.begin(); it1 != params.end(); ++it1) 
    if ( it1->second.fit == param_struct::fixed )
      updateStandardCosVec( it1->first, current_cosparams, 
			    it1->second.fixval );

  //Point it1 and it2 at the vars we will loop over
  it1 = params.find( param_tags::omegam );
  if (it1 == params.end() || it1->second.fit != param_struct::loop ) {
    std::string errstr = "All supported 2D fit types use Omega_m ";
    throw CosFitterExcept("auxconstraint::base_auxconstraint",
			  "makeProbSurface2D",
			  errstr,1);
  }

  it2 = params.find( param_tags::omegade );
  if (it2 != params.end() ) {
    //We know there were 2 cosmo params with loops, and we already
    // made sure omega_m was present
    if (fixcurv) 
      throw CosFitterExcept("auxconstraint::base_auxconstraint",
			    "makeProbSurface2D",
			    "Have Omega_m and Omega_de in 2D, but fixed curvature",2);
  } else {
    it2 = params.find( param_tags::w0 );
    if (it2 == params.end() ) {
      std::string errstr = "Can't determine fit type: have omega_m, ";
      errstr += "no omega_de, no w0";
      throw CosFitterExcept("auxconstraint::base_auxconstraint",
			    "makeProbSurface2D",
			    errstr,4);
    }
    if (!fixcurv)
      throw CosFitterExcept("auxconstraint::base_auxconstraint",
			    "makeProbSurface2D",
			    "Need fixed curvature to 2D fit omega_m, w0",16);
  }

  cosgrid2D prob2D( it1->second, it2->second );

  param_tags::paramcodes code0 = prob2D.getAxisSpec(0).first;
  param_tags::paramcodes code1 = prob2D.getAxisSpec(1).first;
  unsigned int n0 = prob2D.getAxisN(0);
  unsigned int n1 = prob2D.getAxisN(1);
  double cos0, cos1;
  for (unsigned int i = 0; i < n0; ++i) {
    cos0 = prob2D.getAxisVal(0,i);
    updateStandardCosVec( code0, current_cosparams, cos0 );
    for (unsigned int j = 0; j < n1; ++j) {
      cos1 = prob2D.getAxisVal(1,j);
      updateStandardCosVec( code1, current_cosparams, cos1 );
      prob2D[i][j] = exp( -0.5 * GetChiSq( current_cosparams ) );
    }
  }
  prob2D.normalize();
  prob2D.writeFile( fname, binaryout );
}

void 
auxconstraint::base_auxconstraint::makeProbSurface3D(const std::map< param_tags::paramcodes,
						     param_struct >& params,
						     const std::string& fname,
						     bool binaryout ) {

  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

  std::map< param_tags::paramcodes, param_struct >::const_iterator it1,it2,it3;
  //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  for (it1 = params.begin(); it1 != params.end(); ++it1) 
    if ( it1->second.fit == param_struct::fixed )
      updateStandardCosVec( it1->first, current_cosparams, 
			    it1->second.fixval );

  //Point it1/2/3 at the vars we will loop over
  it1 = params.find( param_tags::omegam );
  if (it1 == params.end() || it1->second.fit != param_struct::loop ) {
    std::string errstr = "All supported 2D fit types use Omega_m ";
    throw CosFitterExcept("auxconstraint::base_auxconstraint",
			  "makeProbSurface2D",errstr,1);
  }

  it2 = params.find( param_tags::omegade );
  if (it2 != params.end() ) {
    if (fixcurv) 
      throw CosFitterExcept("auxconstraint::base_auxconstraint",
			    "makeProbSurface3D",
			    "Have Omega_m and Omega_de in 3D, but fixed curvature",2);
    it3 = params.find( param_tags::w0 );
    if (it3 == params.end() ) 
      throw CosFitterExcept("auxconstraint::base_auxconstraint",
			    "makeProbSurface3D",
			    "Can't decide fit type: have om, ode, no w0",4);
  } else {
    it2 = params.find( param_tags::w0 );
    if (it2 == params.end() ) {
      std::string errstr = "Can't determine fit type: have omega_m, ";
      errstr += "no omega_de, no w0";
      throw CosFitterExcept("auxconstraint::base_auxconstraint",
			    "makeProbSurface3D",
			    errstr,8);
    }
    it3 = params.find( param_tags::wa );
    if (it2 == params.end() ) {
      std::string errstr = "Can't determine fit type: have omega_m, ";
      errstr += "no w0, no wa";
      throw CosFitterExcept("auxconstraint::base_auxconstraint",
			    "makeProbSurface3D",
			    errstr,16);
    }
    if (!fixcurv)
      throw CosFitterExcept("auxconstraint::base_auxconstraint",
			    "makeProbSurface3D",
			    "Need fixed curvature to 3D fit omega_m, w0, wa",
			    32);
  }
 

  cosgrid3D prob3D( it1->second, it2->second, it3->second );

  param_tags::paramcodes code0 = prob3D.getAxisSpec(0).first;
  param_tags::paramcodes code1 = prob3D.getAxisSpec(1).first;
  param_tags::paramcodes code2 = prob3D.getAxisSpec(2).first;
  unsigned int n0 = prob3D.getAxisN(0);
  unsigned int n1 = prob3D.getAxisN(1);
  unsigned int n2 = prob3D.getAxisN(2);
  double cos0, cos1, cos2;
  for (unsigned int i = 0; i < n0; ++i) {
    cos0 = prob3D.getAxisVal(0,i);
    updateStandardCosVec( code0, current_cosparams, cos0 );
    for (unsigned int j = 0; j < n1; ++j) {
      cos1 = prob3D.getAxisVal(1,j);
      updateStandardCosVec( code1, current_cosparams, cos1 );
      for (unsigned int k = 0; k < n2; ++k) {
	cos2 = prob3D.getAxisVal(2,k);
	updateStandardCosVec( code2, current_cosparams, cos2 );
	prob3D[i][j][k] = exp( -0.5 * GetChiSq( current_cosparams ) );
      }
    }
  }
  prob3D.normalize();
  prob3D.writeFile( fname, binaryout );
}

/*!
  \param[in] params 4 element vector of parameter values in the order
              \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \returns \f$\chi^2\f$ 
 */
double
auxconstraint::base_auxconstraint::GetChiSq(const std::vector<double>& params) {
  if ( params.size() != 4 )
    throw CosFitterExcept("auxconstraint::base_auxconstraint","GetChiSq",
                          "Param vector wrong size",1);
  return GetChiSq( params[2], params[3], params[0], params[1] );
}


////////////////////////////////////////////////////////////////
//                         baoE05
////////////////////////////////////////////////////////////////

auxconstraint::baoE05::baoE05(double z, double aval, double aerr, 
			      double n, double f, bool fcurv,
			      double ok) :
  base_auxconstraint(fcurv,ok), qag_alloc(false), z(z), a_meas(aval), 
  a_err(aerr), ns(n), fnu(f), useKomatsuForm( false ), atrans( 0.1 ) {}

auxconstraint::baoE05::~baoE05() {
  if (qag_alloc) gsl_integration_workspace_free(work);
}

//This only addressed mutable variables, so it's const
void auxconstraint::baoE05::allocateQAGWork() const {
  work = gsl_integration_workspace_alloc(1000); 
  qag_alloc = true;
}


/*!
  \param[in] w0   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns Eisenstein A parameter
*/
double auxconstraint::baoE05::aval(double w0, double wa,
				   double om, double ode) const {

  const double doubletol = 1e-5; //Used for various comso parameter tolerances
  const double badval = 1e30;

  double ok = 1.0 - om - ode;

  if (om < 0.0) return badval;  //Doesn't work

  //Komatsu form doesn't allow truly constant w
  bool constw = (! useKomatsuForm) && (fabs(wa) < doubletol);
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);

  //Bouncing universe test -- hard to have the BAO constraint,
  // which assumes things about the CMB
  if (cosmoconst && isRebounding( om, ode ) ) return badval;

  //Set up and do integration
  int int_fail;
  size_t neval;
  double intval, error, ooverE;
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  gsl_function F;
  if ( cosmoconst ) {
    double params[2];
    params[0] = om;
    params[1] = ode;
    F.function = one_over_E;
    F.params = &params;
    int_fail = gsl_integration_qng(&F,0,z,0,1e-7,&intval,&error,&neval);
    if (int_fail) {
      //Try adaptive routine
      if (!qag_alloc) allocateQAGWork();
      int_fail = gsl_integration_qag(&F,0,z,0,1e-7,1000,4,work,
				     &intval,&error);
    }
    ooverE = one_over_E( z, params );
  } else if (constw) {
    double params[3];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    F.function = one_over_E_w;
    F.params = &params;
    int_fail = gsl_integration_qng(&F,0,z,0,1e-7,&intval,&error,&neval);
    if (int_fail) {
      //Try adaptive routine
      if (!qag_alloc) allocateQAGWork();
      int_fail = gsl_integration_qag(&F,0,z,0,1e-7,1000,4,work,
				     &intval,&error);
    }
    ooverE = one_over_E_w( z, params );
  } else if (useKomatsuForm) {
    //w(a), Komatsu form
    double params[5];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = atrans;
    F.function = one_over_E_w0waKom;
    F.params = &params;
    int_fail = gsl_integration_qng(&F,0,z,0,1e-7,&intval,&error,&neval);
    if (int_fail) {
      //Try adaptive routine
      if (!qag_alloc) allocateQAGWork();
      int_fail = gsl_integration_qag(&F,0,z,0,1e-7,1000,4,work,
				     &intval,&error);
    }
    ooverE = one_over_E_w0waKom( z, params );
  } else {
    //w(a), Linder form
    double params[4];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    F.function = one_over_E_w0wa;
    F.params = &params;
    int_fail = gsl_integration_qng(&F,0,z,0,1e-7,&intval,&error,&neval);
    if (int_fail) {
      //Try adaptive routine
      if (!qag_alloc) allocateQAGWork();
      int_fail = gsl_integration_qag(&F,0,z,0,1e-7,1000,4,work,
				     &intval,&error);
    }
    ooverE = one_over_E_w0wa( z, params );
  }

  gsl_set_error_handler( old_handler );

  double curra;
  if (int_fail) {
    //Integral didn't evaluate
    return badval;
  } else {
    //Recall that f1 is 1/E, not E from Eisenstein
    if ( fabs(ok) < doubletol ) {
      //Treat as flat
      curra = sqrt( om ) * pow( ooverE, 1.0/3.0 ) *
        pow( intval / z, 2.0/3.0 );
    } else if ( ok < 0.0 ) {
      //Closed Universe
      curra = sqrt( om ) *
        pow( ooverE / fabs( ok ), 1.0/3.0 ) *
        pow( 1.0/z * sin( sqrt(fabs(ok)) * intval ), 2.0/3.0 );
    } else {
      //Open Universe
      curra = sqrt( om ) *
        pow( ooverE / ok, 1.0/3.0 ) *
        pow( 1.0/z * sinh( sqrt(ok) * intval ), 2.0/3.0 );
    }
  }
  return curra;
}

/*!
  \param[in] params 4 element vector of parameter values in the order
  \f$\Omega_m, \Omega_{DE}, w_{0}, w_{a}\f$
  \returns Eisenstein A parameter
*/
double auxconstraint::baoE05::aval( const std::vector<double>& params ) const {
  if ( params.size() != 4 )
    throw CosFitterExcept("auxconstraint::baoE05","aval",
			  "Param vector wrong size",1);
  return aval( params[2], params[3], params[0], params[1] );
}

/*!
  \param[in] w0   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns \f$\chi^2\f$ at parameter values

  If you specify a neutrino fraction, this will be modified
  following Goobar et al. (2006).
  Note that an expansion is used for the neutrino fraction, so if
  the value is large the expansion is probably invalid
*/
double
auxconstraint::baoE05::GetChiSq(double w0, double wa,
				double om, double ode) {

  const double doubletol = 1e-5; //Used for various comso parameter tolerances

  double a_int; //!< Internal value of a

  //Modify for different scalar index
  a_int = a_meas;
  if (fabs(ns - 0.98) > doubletol) a_int *= pow( ns/0.98, -0.35 );

  //And neutrino mass
  if (fabs(fnu) > doubletol) a_int *= (1.0 + 0.94*fnu);

  double chisq = ( a_meas-aval(w0,wa,om,ode) )/a_err;
  return chisq*chisq;
}

////////////////////////////////////////////////////////////////////
//               WMAP3 Shift parameter
////////////////////////////////////////////////////////////////////


auxconstraint::wmap3yr_dls::wmap3yr_dls(bool fcurv, double ok) :
  base_auxconstraint(fcurv,ok), useKomatsuForm(false),  atrans(0.1) {
  work = gsl_integration_workspace_alloc(1000); 
  h = 0.72; h_error = 0.16; //HST Key project, errors doubled
  nhsteps = 50;
  nombh2steps = 50;
  usela = false;
}

auxconstraint::wmap3yr_dls::~wmap3yr_dls() {
  gsl_integration_workspace_free(work);
}

std::pair<double,double> auxconstraint::wmap3yr_dls::GetHubblePrior() const {
  return std::pair<double,double>(h,h_error);
}

/*!
  \param[in] w0  \f$w_0\f$
  \param[in] wa  \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns Value of \f${\mathcal R}\f$, the shift parameter, at the 
  given values of the cosmological parameters

*/
double
auxconstraint::wmap3yr_dls::rval(double w0, double wa, double om, double ode) const
{

  const double doubletol = 1e-5; //Used for various comso parameter tolerances

  double zval = 1089; //z of last scattering surface

  double ok = 1.0 - om - ode;

  //Rebounding test
  bool constw = (! useKomatsuForm) && (fabs(wa) < doubletol);
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return 10000;

  if (om < 0.0) return 10000;

  //Set up and do integration
  gsl_function F;
  if ( cosmoconst ) {
    double params[2];
    params[0] = om;
    params[1] = ode;
    F.function = one_over_E;
    F.params = &params;
  } else if (constw) {
    double params[3];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    F.function = one_over_E_w;
    F.params = &params;
  } else if (useKomatsuForm) {
    double params[5];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = atrans;
    F.function = one_over_E_w0waKom;
    F.params = &params;
  } else {
    //w(a), Linder form
    double params[4];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    F.function = one_over_E_w0wa;
    F.params = &params;
  }

  //This integral is more difficult than the BAO one, 
  // and requires adaptive integration
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  double intval, error;
  int st = gsl_integration_qag(&F,0,zval,0,1e-7,1000,6,work,&intval,&error);
  gsl_set_error_handler( old_handler );

  double curr_r;
  if (st) {
    //Integral didn't evaluate
    return 11000;
  } else {
    if ( fabs(ok) < doubletol ) {
      //Treat as flat
      curr_r = sqrt( om ) * intval;
    } else if ( ok < 0.0 ) {
      //Closed Universe
      curr_r = sqrt( om/fabs(ok) ) * sin( sqrt(fabs(ok)) * intval );
    } else {
      //Open Universe
      curr_r = sqrt( om/ok ) * sin( sqrt(ok) * intval );
    }
  }
  return curr_r;
}

/*!
  \param[in] params 4 element vector of parameter values in the order
  \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \returns Shift parameter \f$R\f$
*/
double 
auxconstraint::wmap3yr_dls::rval( const std::vector<double>& params ) const {
  if ( params.size() != 4 )
    throw CosFitterExcept("auxconstraint::wmap3yr_dls","rval",
			  "Param vector wrong size",1);
  return rval( params[2], params[3], params[0], params[1] );
}

/*!
  \param[in] w0  \f$w_0\f$
  \param[in] wa  \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \param[in] obh2 \f$\Omega_b h^2\f$
  \param[in] h \f$h\f$
  \returns Value of \f$l_a\f$, the second shift parameter, at the 
  given values of the cosmological parameters

*/
double
auxconstraint::wmap3yr_dls::laval(double w0, double wa, double om, 
				  double ode, double obh2, double h) const
{

  const double doubletol = 1e-5; //Used for various comso parameter tolerances


  //Quick exit tests for bad inputs (bouncing, etc.)
  bool constw = (! useKomatsuForm) && (fabs(wa) < doubletol);
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return 10000;
  if (om < 0.0) return 10000;

  double r;
  r = rval( w0, wa, om, ode );
  
  return laval( r, w0, wa, om, ode, obh2, h );

}

/*!
  \param[in] params 4 element vector of parameter values in the order
  \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \param[in] obh2 \f$\Omega_b h^2\f$
  \param[in] h \f$h\f$	      
  \returns Shift parameter \f$l_a\f$, the second shift parameter
*/
double 
auxconstraint::wmap3yr_dls::laval( const std::vector<double>& params,
				   double obh2, double h ) const {
  if ( params.size() != 4 )
    throw CosFitterExcept("auxconstraint::wmap3yr_dls","laval",
			  "Param vector wrong size",1);
  return laval( params[2], params[3], params[0], params[1], obh2, h );
}

/*!
  \param[in] r \f$R\f$ shift parameter from rval
  \param[in] w0  \f$w_0\f$
  \param[in] wa  \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \param[in] obh2 \f$\Omega_b h^2\f$
  \param[in] h \f$h\f$
  \returns Value of \f$l_a\f$, the second shift parameter, at the 
  given values of the cosmological parameters

*/
double
auxconstraint::wmap3yr_dls::laval(double r, double w0, double wa, double om, 
				  double ode, double obh2, double h) const
{

  const double doubletol = 1e-5; //Used for various comso parameter tolerances

  const double zval = 1089; //z of last scattering surface
  const double aval = 1.0 / ( 1.0 + zval ); //a of last scatterin surfact

  //Quick exit tests for bad inputs (bouncing, etc.)
  bool constw = fabs(wa) < doubletol;
  bool cosmoconst = constw && (fabs(w0 + 1) < doubletol);
  if (cosmoconst && isRebounding( om, ode ) ) return 10000;
  if (om < 0.0) return 10000;

  //Set up r_s integral
  gsl_function F;
  if ( cosmoconst ) {
    double params[4];
    params[0] = om;
    params[1] = ode;
    params[2] = obh2;
    params[3] = h;
    F.function = rsint;
    F.params = &params;
  } else if (constw) {
    double params[5];
    params[0] = w0;
    params[1] = om;
    params[2] = ode;
    params[3] = obh2;
    params[4] = h;
    F.function = rsint_w;
    F.params = &params;
  } else if (useKomatsuForm) {
    double params[7];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = obh2;
    params[5] = h;
    params[6] = atrans;
    F.function = rsint_w0waKom;
    F.params = &params;
  } else {
    //w(a), Linder form
    double params[6];
    params[0] = w0;
    params[1] = wa;
    params[2] = om;
    params[3] = ode;
    params[4] = obh2;
    params[5] = h;
    F.function = rsint_w0wa;
    F.params = &params;
  }
  
  //This integral is more difficult than the BAO one, 
  // and requires adaptive integration
  gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
  double intval, error;
  int st = gsl_integration_qag(&F,0,aval,0,1e-7,1000,6,work,&intval,&error);
  gsl_set_error_handler( old_handler );

  if (st) {
    //Integral didn't evaluate
    return 11000;
  } 
  
  return 3.14159 * r / ( sqrt( om ) * intval );
}

/*!
  \param[in] r \f$R\f$ shift parameter from rval
  \param[in] params 4 element vector of parameter values in the order
  \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$
  \param[in] obh2 \f$\Omega_b h^2\f$
  \param[in] h \f$h\f$	      
  \returns Shift parameter \f$l_a\f$, the second shift parameter
*/
double 
auxconstraint::wmap3yr_dls::laval( double r, const std::vector<double>& params,
				   double obh2, double h ) const {
  if ( params.size() != 4 )
    throw CosFitterExcept("auxconstraint::wmap3yr_dls","laval",
			  "Param vector wrong size",1);
  return laval( r, params[2], params[3], params[0], params[1], obh2, h );
}

/*!
  \param[in] w0   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns \f$\chi^2\f$ at parameter values

  This version only does the first shift parameter, which does not
  require an external prior on \f$h\f$.
*/
double
auxconstraint::wmap3yr_dls::GetChiSq_Ronly(double w0, double wa, 
					   double om, double ode) const
{

  //Updated for 2007 W&M values
  //The discontinuity is a little unnerving
  double ival;
  if (fabs(om + ode)-1.0 < 0.001) {
    ival = 1.70;
  } else {
    ival = 1.71;
  }
  double ival_error = 0.03;

  double chisq = ( ival - rval( w0, wa, om, ode ) ) / ival_error;
  return chisq*chisq;
}

/*!
  \param[in] w0   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns \f$\chi^2\f$ at parameter values

  This version incorporates the second shift paramter \f$l_a\f$ of
  Wang and Mukherjee 2007, astro-ph/0703780, and numerically marginalizes
  over \f$\Omega_b h^2\f$ and \f$h\f$ using the WMAP3 value for the first
  and some sort of external prior on the second.

*/
double
auxconstraint::wmap3yr_dls::GetChiSq_Rla(double w0, double wa, 
					 double om, double ode) const
{
  double R, la, ombh2;
  const double ombh2_error = 0.00082;
  double invcov[4][4];

  if (fabs(om + ode)-1 < 0.001) {
    //Assume flat
    R = 1.70;
    la = 302.2;
    ombh2 = 0.022;

    //This is the inverse covariance matrix for the parameters
    // R, l_a, ombh2, h.  We can pre-compute this because h is assumed
    // to be independent of the others.  An alternative would be to use
    // the WMAP constraints on h, but these are quite degenerate with
    // the other params and hence much of the benefit of the second
    // shift parameter will be lost.
    invcov[0][0] = 1128.97;
    invcov[0][1] = invcov[1][0] = 4.07564;
    invcov[0][2] = invcov[2][0] = 1272.01;
    invcov[0][3] = invcov[3][0] = 0.0;
    invcov[1][1] = 1.28161;
    invcov[1][2] = invcov[2][1] = 1250.85;
    invcov[1][3] = invcov[3][1] = 0.0;
    invcov[2][2] = 2.71459e06;
    invcov[2][3] = invcov[3][2] = 0.0;
    invcov[3][3] = 1.0 / (h_error*h_error);

  } else {
    R = 1.71;
    la = 302.5;
    ombh2 = 0.02173;

    invcov[0][0] = 1131.3220;
    invcov[0][1] = invcov[1][0] = 4.8061031;
    invcov[0][2] = invcov[2][0] = 5234.4160;
    invcov[0][3] = invcov[3][0] = 0.0;
    invcov[1][1] = 1.1678060;
    invcov[1][2] = invcov[2][1] = 1077.2189;
    invcov[1][3] = invcov[3][1] = 0.0;
    invcov[2][2] = 2481446.0;
    invcov[2][3] = invcov[3][2] = 0.0;
    invcov[3][3] = 1.0 / (h_error*h_error);
  }
  
  //And now we numerically marginalize over ombh2 and h, which is
  // slow, but such is life.  We go from -3 to +3 sigma in each.
  double minombh2, deltaombh2;
  double minh, deltah;
  minh = h - 3.0 * h_error;
  deltah = 6.0 * h_error / static_cast<double>( nhsteps - 1 );
  minombh2 = ombh2 - 3.0 * ombh2_error;
  deltaombh2 = 6.0 * ombh2_error / static_cast<double>( nombh2steps - 1 );

  double currh, currombh2;
  double delta_params[4];
  double currR = rval( w0, wa, om, ode );
  double chisq, prob, pullh, pullombh2;
  delta_params[0] = R - currR;
  prob = 0.0;
  for (unsigned int i = 0; i < nhsteps; ++i) {
    currh = minh + deltah * static_cast<double>(i);
    delta_params[3] = h - currh;
    for (unsigned int j = 0; j < nombh2steps; ++j) {
      currombh2 = minombh2 + deltaombh2 * static_cast<double>(j);
      delta_params[1] = la - laval( currR, w0, wa, om, ode, currombh2,
				    currh );
      delta_params[2] = ombh2 - currombh2;

      //And do the matrix multiply to get the chisq.  
      // This is a small matrix, so just do it by hand.
      chisq = 0.0;
      for (unsigned int k = 0; k < 4; ++k)
	for (unsigned int m = 0; m < 4; ++m)
	  chisq += delta_params[k] * invcov[k][m] * delta_params[m];

      pullh = delta_params[3]/h_error;
      pullombh2 = delta_params[2]/ombh2_error;

      //Now we have to do the prior marginalization bit
      //We assume our priors in ombh2 and h are uncorrelated,
      // or this would be much messier.
      prob += exp( -0.5 * ( chisq + pullh*pullh + pullombh2*pullombh2 ) );
    }
  }

  //Do the normalization bit.  This is the Gaussian prefactor and the
  // size of the boxes
  prob *= deltah * deltaombh2 / (2.0*3.14159274101257*h_error*ombh2_error);

  //And back to a chisq
  return -2.0 * log( prob );
}
/*!
  \param[in] w0   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns \f$\chi^2\f$ at parameter values

*/
double
auxconstraint::wmap3yr_dls::GetChiSq(double w0, double wa, double om,
				     double ode) {
  if (usela) 
    return GetChiSq_Rla(w0,wa,om,ode);
  else 
    return GetChiSq_Ronly(w0, wa, om, ode);
}

////////////////////////////////////////////////////////////////////
//               WMAP7 Shift parameters
////////////////////////////////////////////////////////////////////

const gsl_multimin_fminimizer_type* auxconstraint::wmap7yr_dls::T =
  gsl_multimin_fminimizer_nmsimplex;

auxconstraint::wmap7yr_dls::wmap7yr_dls(bool fcurv, double ok) :
  base_auxconstraint(fcurv,ok) {
  v = gsl_vector_alloc(2);
  ss = gsl_vector_alloc(2);
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc(T, 2);
}
auxconstraint::wmap7yr_dls::~wmap7yr_dls() {
  gsl_vector_free(v);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
}

/*!
  \param[in] w0   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns \f$\chi^2\f$ at parameter values for WMAP7 shift parameters

  This follows the prescription of Komatsu et al. (2008).  The
  nuisance parameters \f$\Omega_b h^2\f$ and \f$h\f$ are marginalized
  over by minimizing the \f$\chi^2\f$ for these parameters at each
  value of the cosmological parameters (see Lampton et al. (1976) for
  discussion).
*/
double
auxconstraint::wmap7yr_dls::GetChiSq(double w0, double wa,
				     double om, double ode) {
  
  const size_t maxiter = 500;
  const double badval = 1e30;
  const double sizetol1 = 1e-4, sizetol2 = 1e-6;

  if (om <= 0) return badval;  //Can't handle no matter universe
  

  //Do an initial evalution of K2 to make sure we aren't in a bad region
  // using default values.  Basically, this is to test for bad closed
  // universes with 0 or negative D_A
  bool K2_test_success;
  dhelp.GetK2( 1090.0, om, ode, 5e-5, w0, wa, K2_test_success );
  if (!K2_test_success) return badval;

  //Set initial conditions for h, obh2 to WMAP7 values
  gsl_vector_set(v, 0, 0.724);
  gsl_vector_set(v, 1, 0.02268);

  //Set initial step sizes
  gsl_vector_set(ss, 0, 0.05);
  gsl_vector_set(ss, 1, 0.004);
  
  //Make par vector
  double params[4];
  params[0] = w0;
  params[1] = wa;
  params[2] = om;
  params[3] = ode;
  
  //We have to jump through some serious hoops to pass
  // in an instance of the class as a parameter
  void **varr;
  varr = new void*[2];
  varr[0] = static_cast<void*>(&dhelp); //This is evil, but such is c++->c
  varr[1] = static_cast<void*>(params);
  void *vptr;
  vptr = static_cast<void*>(varr);

  //Setup function
  minex_func.n = 2;
  minex_func.f = &wmap7_shift_chisq;
  minex_func.params = vptr;

  gsl_multimin_fminimizer_set(s, &minex_func, v, ss);

  int status;
  size_t iter = 0;
  double size;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
           
    if (status) 
      break;
     
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, sizetol1);
     
    /*
      printf ("%5d %10.3e %10.3e f() = %10.3e size = %.7f\n", 
      iter,
      gsl_vector_get(s->x, 0), 
      gsl_vector_get(s->x, 1), 
      s->fval, size);
    */
  } while (status == GSL_CONTINUE && iter < maxiter);

  //We use non-converged minima anyways, but warn
  if (status != GSL_SUCCESS) {
    printf("Warning -- simplex did not converge\n");
    printf(" at w0: %6.3f wa: %6.3f om: %5.3f ode: %6.3f\n",
	   w0,wa,om,ode);
  } else {
    //Following the recommendation of NumRec, if we succeeded
    // we restart the minimization at the current guess

    //reset initial step sizes to half old values
    gsl_vector_set(ss, 0, 0.025);
    gsl_vector_set(ss, 1, 0.002);

    gsl_multimin_fminimizer_set(s, &minex_func, s->x, ss);
    iter = static_cast<size_t>(0);
    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status) 
	break;
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, sizetol2);
    } while (status == GSL_CONTINUE && iter < maxiter);
    if (status != GSL_SUCCESS) {
      printf("Warning -- second simplex did not converge\n");
      printf(" at w0: %6.3f wa: %6.3f om: %5.3f ode: %6.3f\n",
	     w0,wa,om,ode);
    }
  }

  //Get chisq at minimum
  double chisq;
  chisq = wmap7_shift_chisq( s->x, varr );

  //Clean up
  delete[] varr;

  return chisq;

}

///////////////////////////////////////////////////////////////
//                         omprior                           //
///////////////////////////////////////////////////////////////
auxconstraint::omprior::omprior(double om, double omerr,
				bool fixcurv, double ok) :
  base_auxconstraint(fixcurv,ok) {
  ommean = om;
  omerror = omerr;
}

auxconstraint::omprior::omprior(std::pair<double,double> vals,
				bool fixcurv, double ok) :
  base_auxconstraint(fixcurv,ok) {
  ommean = vals.first;
  omerror = vals.second;
}

double auxconstraint::omprior::GetChiSq(double omval, double a, 
					double b, double c) {
  //Doesn't get much simpler than this.  We ignore the other arguments
  double pull = (omval - ommean) / omerror;
  return pull*pull;
}


////////////////////////////////////////////////////////////////////
//               Percival (09) BAO constraints
////////////////////////////////////////////////////////////////////

const gsl_multimin_fminimizer_type* auxconstraint::baoP09::T =
  gsl_multimin_fminimizer_nmsimplex;

auxconstraint::baoP09::baoP09(bool fcurv, double ok) :
  base_auxconstraint(fcurv,ok) {
  v = gsl_vector_alloc(2);
  ss = gsl_vector_alloc(2);
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc(T, 2);
}
auxconstraint::baoP09::~baoP09() {
  gsl_vector_free(v);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
}

/*!
  \param[in] w0   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns \f$\chi^2\f$ at parameter values for Percival BAO parameters

  The nuisance parameters \f$\Omega_b h^2\f$ and \f$h\f$ are marginalized
  over by minimizing the \f$\chi^2\f$ for these parameters at each
  value of the cosmological parameters (see Lampton et al. (1976) for
  discussion).
*/
double
auxconstraint::baoP09::GetChiSq(double w0, double wa,
				double om, double ode) {
  
  const size_t maxiter = 500;
  const double badval = 1e30;
  const double sizetol1 = 1e-4, sizetol2 = 1e-6;

  if (om <= 0) return badval;  //Can't handle no matter universe
  

  //Do an initial evalution of K2 to make sure we aren't in a bad region
  // using default values.  Basically, this is to test for bad closed
  // universes with 0 or negative D_A
  bool K2_test_success;
  dhelp.GetK2( 1090.0, om, ode, 5e-5, w0, wa, K2_test_success );
  if (!K2_test_success) return badval;

  //Set initial conditions for h, obh2 to WMAP7 values
  gsl_vector_set(v, 0, 0.724);
  gsl_vector_set(v, 1, 0.02268);

  //Set initial step sizes
  gsl_vector_set(ss, 0, 0.05);
  gsl_vector_set(ss, 1, 0.004);
  
  //Make par vector
  double params[4];
  params[0] = w0;
  params[1] = wa;
  params[2] = om;
  params[3] = ode;
  
  //We have to jump through some serious hoops to pass
  // in an instance of the class as a parameter
  void **varr;
  varr = new void*[2];
  varr[0] = static_cast<void*>(&dhelp); //This is evil, but such is c++->c
  varr[1] = static_cast<void*>(params);
  void *vptr;
  vptr = static_cast<void*>(varr);

  //Setup function
  minex_func.n = 2;
  minex_func.f = &percival_chisq;
  minex_func.params = vptr;

  gsl_multimin_fminimizer_set(s, &minex_func, v, ss);

  int status;
  size_t iter = 0;
  double size;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
           
    if (status) 
      break;
     
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, sizetol1);
     
    /*
      printf ("%5d %10.3e %10.3e f() = %10.3e size = %.7f\n", 
      iter,
      gsl_vector_get(s->x, 0), 
      gsl_vector_get(s->x, 1), 
      s->fval, size);
    */
  } while (status == GSL_CONTINUE && iter < maxiter);

  //We use non-converged minima anyways, but warn
  if (status != GSL_SUCCESS) {
    printf("Warning -- simplex did not converge\n");
    printf(" at w0: %6.3f wa: %6.3f om: %5.3f ode: %6.3f\n",
	   w0,wa,om,ode);
  } else {
    //Following the recommendation of NumRec, if we succeeded
    // we restart the minimization at the current guess

    //reset initial step sizes to half old values
    gsl_vector_set(ss, 0, 0.025);
    gsl_vector_set(ss, 1, 0.002);

    gsl_multimin_fminimizer_set(s, &minex_func, s->x, ss);
    iter = static_cast<size_t>(0);
    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status) 
	break;
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, sizetol2);
    } while (status == GSL_CONTINUE && iter < maxiter);
    if (status != GSL_SUCCESS) {
      printf("Warning -- second simplex did not converge\n");
      printf(" at w0: %6.3f wa: %6.3f om: %5.3f ode: %6.3f\n",
	     w0,wa,om,ode);
    }
  }

  //Get chisq at minimum
  double chisq;
  chisq = percival_chisq( s->x, varr );

  //Clean up
  delete[] varr;

  return chisq;

}

////////////////////////////////////////////////////////////////////
//               Percival (09) + WMAP7 BAO constraints
////////////////////////////////////////////////////////////////////

const gsl_multimin_fminimizer_type* auxconstraint::baoP09_wmap7yr_dls::T =
  gsl_multimin_fminimizer_nmsimplex;

auxconstraint::baoP09_wmap7yr_dls::baoP09_wmap7yr_dls(bool fcurv, double ok) :
  base_auxconstraint(fcurv,ok) {
  v = gsl_vector_alloc(2);
  ss = gsl_vector_alloc(2);
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc(T, 2);
}
auxconstraint::baoP09_wmap7yr_dls::~baoP09_wmap7yr_dls() {
  gsl_vector_free(v);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
}

/*!
  \param[in] w0   \f$w_0\f$
  \param[in] wa   \f$w_a\f$
  \param[in] om  \f$\Omega_m\f$
  \param[in] ode \f$\Omega_{DE}\f$
  \returns \f$\chi^2\f$ at parameter values for Percival BAO + WMAP7 parameters

  The nuisance parameters \f$\Omega_b h^2\f$ and \f$h\f$ are marginalized
  over by minimizing the \f$\chi^2\f$ for these parameters at each
  value of the cosmological parameters (see Lampton et al. (1976) for
  discussion).
*/
double
auxconstraint::baoP09_wmap7yr_dls::GetChiSq(double w0, double wa,
					    double om, double ode) {
  
  const size_t maxiter = 500;
  const double badval = 1e30;
  const double sizetol1 = 1e-4, sizetol2 = 1e-6;

  if (om <= 0) return badval;  //Can't handle no matter universe

  //Do an initial evalution of K2 to make sure we aren't in a bad region
  // using default values.  Basically, this is to test for bad closed
  // universes with 0 or negative D_A
  bool K2_test_success;
  dhelp.GetK2( 1090.0, om, ode, 5e-5, w0, wa, K2_test_success );
  if (!K2_test_success) return badval;

  //Set initial conditions for h, obh2 to WMAP7 values
  gsl_vector_set(v, 0, 0.724);
  gsl_vector_set(v, 1, 0.02268);

  //Set initial step sizes
  gsl_vector_set(ss, 0, 0.05);
  gsl_vector_set(ss, 1, 0.004);
  
  //Make par vector
  double params[4];
  params[0] = w0;
  params[1] = wa;
  params[2] = om;
  params[3] = ode;
  
  //We have to jump through some serious hoops to pass
  // in an instance of the class as a parameter
  void **varr;
  varr = new void*[2];
  varr[0] = static_cast<void*>(&dhelp); //This is evil, but such is c++->c
  varr[1] = static_cast<void*>(params);
  void *vptr;
  vptr = static_cast<void*>(varr);

  //Setup function
  minex_func.n = 2;
  minex_func.f = &percival_wmap7_chisq;
  minex_func.params = vptr;

  gsl_multimin_fminimizer_set(s, &minex_func, v, ss);

  int status;
  size_t iter = 0;
  double size;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
           
    if (status) 
      break;
     
    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, sizetol1);
     
    /*
      printf ("%5d %10.3e %10.3e f() = %10.3e size = %.7f\n", 
      iter,
      gsl_vector_get(s->x, 0), 
      gsl_vector_get(s->x, 1), 
      s->fval, size);
    */
  } while (status == GSL_CONTINUE && iter < maxiter);

  //We use non-converged minima anyways, but warn
  if (status != GSL_SUCCESS) {
    printf("Warning -- simplex did not converge\n");
    printf(" at w0: %6.3f wa: %6.3f om: %5.3f ode: %6.3f\n",
	   w0,wa,om,ode);
  } else {
    //Following the recommendation of NumRec, if we succeeded
    // we restart the minimization at the current guess

    //reset initial step sizes to half old values
    gsl_vector_set(ss, 0, 0.025);
    gsl_vector_set(ss, 1, 0.002);

    gsl_multimin_fminimizer_set(s, &minex_func, s->x, ss);
    iter = static_cast<size_t>(0);
    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status) 
	break;
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, sizetol2);
    } while (status == GSL_CONTINUE && iter < maxiter);
    if (status != GSL_SUCCESS) {
      printf("Warning -- second simplex did not converge\n");
      printf(" at w0: %6.3f wa: %6.3f om: %5.3f ode: %6.3f\n",
	     w0,wa,om,ode);
    }
  }

  //Get chisq at minimum
  double chisq;
  chisq = percival_wmap7_chisq( s->x, varr );

  //Clean up
  delete[] varr;

  return chisq;

}

