//User interface driver for single luminosity
// distances

#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include <string>

#include "auxconstraint.h"

using namespace std;

int main( int argc, char** argv ) {

  int c;
  
  bool userode, flatflag, komatsuform;
  double om, ode, w0, wa, h, atrans, obh2;
  ode = 0.75;
  w0 = -1.0;
  wa = 0.0;
  h = 0.65;
  atrans = 0.1;
  userode = false;
  flatflag = false;
  komatsuform = false;

  const unsigned int nz = 2;
  const double percival_z[nz] = { 0.2, 0.35 };
  const double pi = 3.14159265358979323846264338327950288419716939937510582;


  int option_index = 0;
  static struct option long_options[] = { 
    {"help", no_argument, 0, 'h'},
    {"komatsu",no_argument,0,'k'},
    {"tol",required_argument,0,'t'},
    {"omega_m",required_argument, 0, 'm'},
    {"omega_de",required_argument, 0, 'd'},
    {"w0",required_argument,0, '2'},
    {"wa",required_argument,0, '3'},
    {"atrans",required_argument,0, 'a'},
    { 0 , 0 , 0 , 0 }
  };
  
  while ( ( c = getopt_long( argc, argv, "hkm:d:2:3:a:",
			     long_options, &option_index ) ) != -1 )
    switch (c) {
    case 'h' :
      cerr << "NAME\n";
      cerr << "\tget_acoustic -- Returns acoustic distances."<<endl;
      cerr << "SYNOPSIS\n";
      cerr << "\tget_acoustic om h obh2 [ -d | --omega_de=OMEGA_DE]" << endl;
      cerr << "\t [ --w0=W0 ] [ --wa=WA ] [ --h=H ]" << std::endl;
      cerr << "\t [ -k ] [ -a | --atrans=ATRANS ]" << endl;
      cerr << "DESCRIPTION\n";
      cerr << "\tThis calculates various BAO and CMB parameters given the cosmological" << endl;
      cerr << "\tparameters." << std::endl;
      cerr << "OPTIONS" << endl;
      cerr << "\t-d, --omega_de OMEGA_DE\n\t\tOmega_DarkEnergy (def: 0.75)" << endl;
      cerr << "\t--w0 W0\n\t\tThe constant equation of state parameter for the"
	   << std::endl;
      cerr << "\t\t dark energy. (def: -1)" << std::endl;
      cerr << "\t--wa W0A\n\t\tThe a derivative of w (def: 0)." << endl;
      cerr << "\t--h H\n\t\tThe reduced Hubble constant (H0/100 km s^-1 Mpc^-1) (def: 0.65)" << endl;
      cerr << "\t--t, --tol TOL\n\t\tTolerance for parameters.  For example "
	   << "if Omega_k < TOL" << std::endl;
      cerr << "\t\tthen the fit is treated as flat" << std::endl;
      cerr << "\t-f\n\t\tAssume a flat universe for the fit" << endl;
      cerr << "\t-k\n\t\tUse Komatsu et al. form for w(a)" << endl;
      cerr << "\t--a, --atrans ATRANS\n\t\tTransition scale factor for"
	   << " Komatsu w(a) form (def: 0.1)" << std::endl;
      return 0;
      break;
    case 'd' :
      ode = atof( optarg );
      userode = true;
      break;
    case '2' :
      w0 = atof( optarg );
      break;
    case '3' :
      wa = atof( optarg );
      break;
    case 'k' :
      komatsuform = true;
      break;
    case 'a' :
      atrans = atof(optarg);
      break;
    }

  if (optind > argc - 3) {
    cerr << "get_acoustic :: Required arguments not provided" << endl;
    return 1;
  }

  om   = atof(argv[optind]);
  h    = atof(argv[optind+1]);
  obh2 = atof(argv[optind+2]);
  
  double h2 = h * h;
  double omh2 = om * h2;

  if (userode) {
    if ( fabs( om + ode - 1.0 ) > 0.001 ) flatflag=true;
  } else ode = 1.0 - om;    

  auxconstraint::distance_helper dhelp;
  
  if (komatsuform) {
    dhelp.SetUseKomatsuForm();
    dhelp.SetAtrans(atrans);
  }

  double zstar, zdrag;

  zstar = dhelp.GetZstar( obh2, omh2 );
  zdrag = dhelp.GetZdrag( obh2, omh2 );

  double orad = auxconstraint::oradh2 / h2;

  double K1_wmap, K1_bao, K2_wmap, K2_bao_0, K2_bao_1;
  bool success_K1, success_K2;
  K1_bao=dhelp.GetK1( zdrag, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1) exit(1);
  K1_wmap=dhelp.GetK1( zstar, obh2, om, ode, w0, wa, h, success_K1 );
  if (! success_K1 || K1_wmap == 0.0) exit(1); //Failure
  K2_wmap=dhelp.GetK2( zstar, om, ode, orad, w0, wa, success_K2 );
  if (! success_K2) exit(1); //Failure
  K2_bao_0=dhelp.GetK2( percival_z[0], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_bao_0 == 0) exit(1); //Failure
  K2_bao_1=dhelp.GetK2( percival_z[1], om, ode, orad, w0, wa, success_K2 );
  if (! success_K2 || K2_bao_1 == 0) exit(1); //Failure

  double la, R;
  la = pi * K2_wmap / K1_wmap;
  R  = sqrt( om )*K2_wmap;

  //Really, we hold H0/c CV = [K2^2 z/E]^1/3
  double DV_0, DV_1;
  bool success_DV;
  DV_0 = pow( K2_bao_0 * K2_bao_0 * percival_z[0] * 
	      dhelp.GetOneOverE( percival_z[0], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_0 == 0 ) exit(1);  //Failure
  DV_1 = pow( K2_bao_1 * K2_bao_1 * percival_z[1] * 
	      dhelp.GetOneOverE( percival_z[1], om, ode, orad, w0, wa, 
			       success_DV ), 1.0/3.0 );
  if ( ! success_DV || DV_1 == 0 ) exit(1);  //Failure

  double rsoverDV_0, rsoverDV_1;
  rsoverDV_0 = K1_bao/DV_0;
  rsoverDV_1 = K1_bao/DV_1;

  double cinvH = 2997.92458 / h; 
  printf("For w0: %8.5f wa: %8.5f om: %8.5f ode: %8.5f\n",
	 w0,wa,om,ode );
  printf("    omh2: %8.5f obh2: %9.6f h: %5.3f\n",omh2,obh2,h);
  printf("zstar: %9.3f l_a: %8.5f R: %8.5f\n", zstar, la, R);
  printf("zdrag: %9.3f rs(zd): %7.3f Mpc\n", zdrag, cinvH * K1_bao);
  printf("rs/Dv (0.2): %7.5f rs/Dv (0.35): %7.5f\n",rsoverDV_0,
	 rsoverDV_1);
  return 0;
}
