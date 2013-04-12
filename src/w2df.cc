//Creates a 2dF probability surface in w/om space, assuming flatness


//Uses code for calculating f due to Eric Linder

#include <cmath>
#include "w2df.h"
#include "cosfitterexcept.h"
#include "param_tags.h"

using namespace std;

/*!
  \param[in] zfin Redshift to evaluate growth factor at
  \param[in] cosparams Cosmological parameters to use
  \param[in] nstep Number of steps to use.  10000 is reasonable
  \returns The growth factor at redshift zfin
*/
double w2dFGRS::fcalc(double zfin, const w2df_params& cosparams, 
		      int nstep) {
  
  //Starting and ending value of integration over scale factor a
  const double starta = 1e-4;
  
  //Variables
  int nstep1 = nstep + 1;
  double a;                 //Integration variable
  double enda;              //Final value of a
  double v[2], dv[2];
  double stepsize;
  double *aa;               //Array of a values
  double *y0;               //Array of growth factors
  double *y1;               //Array of derivatives of growth factors
  double *f;                //Array of f

  enda = 1.0/(1.0 + zfin);

  //Initialize arrays
  //These will be of nstep+1 in size, as the final entry
  // will hold the results of the integration
  aa = new double[nstep1];
  y0 = new double[nstep1];
  y1 = new double[nstep1];
  f  = new double[nstep1];
  
  //Initial conditions
  y0[0] = v[0] = 1.0; //growth
  y1[0] = v[1] = 0.0; //derivative
  aa[0] = a = starta;
  stepsize = (enda-starta)/nstep; 

  //Integration loop
  for (int i = 0; i < nstep; ++i) {
    derivs(a,v,dv,cosparams);
    rk4(v,dv,a,stepsize,v,cosparams);
    if ( (a + stepsize) == a )
      throw CosFitterExcept("w2df","fcalc","stepsize not large enough",1);

    a += stepsize; 
    aa[i+1] = a;
    y0[i+1] = v[0];
    y1[i+1] = v[1];
  }
      
  //Evaluate velocity, or log, growth factor 1+d\ln G/d\ln a
  for (int i = 0; i < nstep1; ++i) {
    f[i] = 1.0 + aa[i]*y1[i]/y0[i];
  }

  //Clean up
  double retval = f[nstep];
  delete[] aa; delete[] y0; delete[] y1; delete[] f;

  return retval;
}

/*!
  4th order Runge-Kutta integrator
  From numerical recipies, specialized to case of 2 variables
  and to this particular problem

  \param[in] y Initial value for variables at x
  \param[in] dydx Derivatives of variables at x
  \param[in] x Starting position
  \param[in] h Amount to advance integral
  \param[in] params Cosmological parameters for this problem, used to get
                     derivatives
  \param[out] yout New values of variables
*/
void w2dFGRS::rk4( double y[2], double dydx[2], double x,
		   double h, double yout[2], const w2df_params& params) {

  double yt[2],dyt[2],dym[2]; //Working version of input and 2 x derivs

  double hh;  //half step size
  double h6;  //Sixth step size
  double xh;  //x + h/2
  
  hh = h*0.5;
  h6 = h/6.0;
  xh = x + hh;

  //First step
  yt[0] = y[0] + hh*dydx[0];
  yt[1] = y[1] + hh*dydx[1];

  //Second step
  derivs(xh,yt,dyt,params);
  yt[0] = y[0] + hh*dyt[0];
  yt[1] = y[1] + hh*dyt[1];

  //Third step
  derivs(xh,yt,dym,params);
  yt[0] = y[0] + h*dym[0];
  yt[1] = y[1] + h*dym[1];
  dym[0] += dym[0];
  dym[1] += dym[1];

  //Fourth and final step, accumulating previous increments
  derivs(x+h,yt,dyt,params);
  yout[0] = y[0] + h6*(dydx[0]+dyt[0]+2*dym[0]);
  yout[1] = y[1] + h6*(dydx[1]+dyt[1]+2*dym[1]);
  
}


/*!
  \param[in] x Place to evaluate derivative
  \param[in] y Input values
  \param[out] dydx Derivatives at x
  \param[in] params Cosmological parameters

  For this particular problem, x is a, y[0] is the growth factor,
  y[1] is d growth_factor/da.
*/
void w2dFGRS::derivs( double x, double y[2], double dydx[2],
		      const w2df_params& params) {
  double wofa, powarg, exparg, xxx;

  wofa = params.w0 + params.wa*(1.0 - x); //w at current a
  powarg = 3.0 * ( params.w0 + params.wa );
  exparg = 3.0 * params.wa * (1.0 - x);
  xxx = params.om/(1.0-params.om) * pow(x,powarg) * exp( exparg );
  
  dydx[0] = y[1];
  dydx[1] = -(7.0/2.0 - 1.5 * wofa /(1+xxx)) * y[1]/x;
  dydx[1] -= 3.0/2.0*(1-wofa)/(1+xxx)*y[0]/(x*x);
}

/*!
  Interface to find growth factor at a given redshift
  \param[in] zfin Redshift to get growth factor at
  \param[in] om \f$\Omega_m\f$
  \param[in] w0 w at z=0
  \param[in] wa Derivative of w with a.  Specifically, 
     \f$w\left(a\right) = w0 + wa \left(1 - a\right)\f$.  a is the scale
     factor: \f$ a = \frac{1}{1+z}\f$.
  \returns The growth factor at zfin
*/
double w2dFGRS::get_growth_factor(double zfin, double om, double w0, 
				  double wa) {
  w2df_params params(om,w0,wa);
  return fcalc(zfin,params,10000);
}


/*!
  Fill in grid of probability values by comparing growth factor
  measurement with the theoretical value.  This is specialized
  to the flat \f$ w \f$, \f$ \Omega_m \f$ case.

  \param[in] fval The measured growth of structure parameter
  \param[in] dfval The error in fval
  \param[in] z The redshift fval was measured at
  \param[in] nw The number of \f$ w \f$ points to calculate
  \param[in] minw The minimum \f$ w \f$ value to consider
  \param[in] maxw The maximum \f$ w \f$ value to consider
  \param[in] nom The number of \f$ \Omega_m \f$ points to calculate
  \param[in] minom The minimum \f$ \Omega_m \f$ value to consider
  \param[in] maxom The maximum \f$ \Omega_m \f$ value to consider
  \returns A cosgrid2D loaded up with the relative probability levels
    for various values of the cosmological parameters.
*/
cosgrid2D w2dFGRS::make_prob_grid(double fval, double dfval, double z,
				  int nw, double minw, double maxw,
				  int nom, double minom, double maxom) {
  

  //Caluclate step sizes
  double dw, dom;
  dw = dom = 0.0;
  if (nw > 1) dw = (maxw - minw)/(nw - 1);
  if (nom > 1) dom = (maxom - minom)/(nom - 1);

  param_tags::paramspec spec1, spec2;
  spec1 = std::make_pair( param_tags::w0, param_tags::cosmological );
  spec1 = std::make_pair( param_tags::omegam, param_tags::cosmological );
  cosgrid2D outgrid(nw,minw,maxw,dw,param_tags::wtag, spec1,
		    nom,minom,maxom,dom,param_tags::omtag,spec2);
  
  //Calculation loop
  double currom, currw, ftheory, fvar, diff_f;
  fvar = 1.0 / (dfval*dfval);
  for (int i = 0; i < nw; ++i) {
    currw = minw + i * dw;
    for (int j = 0; j < nom; ++j) {
      currom = minom + j*dom;
      ftheory = get_growth_factor(z,currom,currw,0.0);
      diff_f = fval - ftheory;
      //The probability is proportional to -1/2 chisquare
      outgrid[i][j] = exp( - 0.5 * diff_f*diff_f *fvar );
    }
  }

  //Normalize grid
  outgrid.normalize();

  return outgrid;

}


/*!
  Fill in grid of probability values by comparing growth factor
  measurement with the theoretical value.  This is specialized
  to the flat \f$w\f$, \f$\Omega_m\f$ case.

  \param[in] fval The measured growth of structure parameter
  \param[in] dfval The error in fval
  \param[in] z The redshift fval was measured at
  \param[in] ingrid cosgrid2D used to set the values the grid is calculated at
  \returns A cosgrid2D loaded up with the relative probability levels
    for various values of the cosmological parameters.
*/
cosgrid2D w2dFGRS::make_prob_grid(double fval, double dfval, double z,
				  const cosgrid2D& ingrid) {

  //Check to make sure that grid corresponds to w/om
  if (ingrid.getAxisLabel(0) != param_tags::wtag || 
      ingrid.getAxisLabel(1)!= param_tags::omtag ) 
    throw CosFitterExcept("w2df","make_prob_grid",
			  "Input cosgrid not of correct type",1);

  cosgrid2D outgrid( ingrid );
  
  //Calculation loop
  double currom, currw, ftheory, fvar, diff_f;
  fvar = 1.0 / (dfval*dfval);
  for (unsigned int i = 0; i < outgrid.getAxisN(0); ++i) {
    currw = outgrid.getAxisVal(0,i);
    for (unsigned int j = 0; j < outgrid.getAxisN(1); ++j) {
      currom = outgrid.getAxisVal(1,j);
      ftheory = get_growth_factor(z,currom,currw,0.0);
      diff_f = fval - ftheory;
      //The probability is proportional to -1/2 chisquare
      outgrid[i][j] = exp( - 0.5 * diff_f*diff_f *fvar );
    }
  }

  //Normalize grid
  outgrid.normalize();

  return outgrid;
}


