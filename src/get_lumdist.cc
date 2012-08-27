//User interface driver for single luminosity
// distances

#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include <string>

#include <lumdist.h>
#include <snedata.h>

using namespace std;

int main( int argc, char** argv ) {

  int c;
  
  bool userom, userode, flatflag, logflag, komatsuform;
  float redshift;
  double om, ode, w0, wa, h, atrans;
  om = 0.25;
  ode = 0.75;
  w0 = -1.0;
  wa = 0.0;
  h = 0.65;
  atrans = 0.1;
  logflag = false;
  userom = false;
  userode = false;
  flatflag = false;
  komatsuform = false;

  lumdist lm;

  int option_index = 0;
  static struct option long_options[] = { 
    {"help", no_argument, 0, 'h'},
    {"flat", no_argument, 0, 'f'},
    {"log", no_argument, 0, 'l'},
    {"komatsu",no_argument,0,'k'},
    {"tol",required_argument,0,'t'},
    {"h", required_argument, 0, '1' },
    {"omega_m",required_argument, 0, 'm'},
    {"omega_de",required_argument, 0, 'd'},
    {"w0",required_argument,0, '2'},
    {"wa",required_argument,0, '3'},
    {"atrans",required_argument,0, 'a'},
    { 0 , 0 , 0 , 0 }
  };
  
  while ( ( c = getopt_long( argc, argv, "hflkm:d:1:2:3:t:a:",
			     long_options, &option_index ) ) != -1 )
    switch (c) {
    case 'h' :
      cerr << "NAME\n";
      cerr << "\tget_lumdist -- Returns luminosity distance."<<endl;
      cerr << "SYNOPSIS\n";
      cerr << "\tget_lumdist [ -m | --omega_m=OMEGA_MASS ] [ -d | --omega_de=OMEGA_DE]" << endl;
      cerr << "\t [ --w0=W0 ] [ --wa=WA ] [ --h=H ] [ -t | --tol=TOL ] [ -f ] [ -l ]" << std::endl;
      cerr << "\t [ -k ] [ -a | --atrans=ATRANS ] redshift [ redshift2 ... ]" << endl;
      cerr << "DESCRIPTION\n";
      cerr << "\tThis calculates the luminosity distance(s) given the cosmological" << endl;
      cerr << "\tparameters. The only required input is the redshift, since the" << endl; 
      cerr << "\tother parameters have defaults.  Relativistic particles" << endl;
      cerr << "\t(neutrinos, etc.) are not supported.  Multiple redshifts are possible." << endl;
      cerr << "OPTIONS" << endl;
      cerr << "\t-m, --omega_m OMEGA_MASS\n\t\tOmega_mass (more precisely, the density parameter" << endl;
      cerr << "\t\tof non-relativistic matter (def: 0.25)\n";
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
      cerr << "\t-l\n\t\tReturn 5 log_10 (dl) instead of dl" << endl;
      cerr << "\t-k\n\t\tUse Komatsu et al. form for w(a)" << endl;
      cerr << "\t--a, --atrans ATRANS\n\t\tTransition scale factor for"
	   << " Komatsu w(a) form (def: 0.1)" << std::endl;
      return 0;
      break;
    case '1' : 
      h = atof( optarg );
      break;
    case 'm' :
      om = atof( optarg );
      userom = true;
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
    case 't' :
      lm.setParamTol( atof(optarg) );
      break;
    case 'f' :
      flatflag = true;
      break;
    case 'l' :
      logflag = true;
      break;
    case 'k' :
      komatsuform = true;
      break;
    case 'a' :
      atrans = atof(optarg);
      break;
    }

  if (optind > argc - 1) {
    cerr << "get_lumdist :: Required argument redshift not provided" << endl;
    return 1;
  }

  if (flatflag) {
    if (userode && userom && ( fabs(om + ode - 1.0) > 0.0001) ) {
      std::cerr << "get_lumdist :: Flat flag set but user specified " 
		<< "values do not allow" << std::endl;
      return 2;
    }
    if (userode) { om = 1.0 - ode; } else if (userom) { ode = 1.0 - om; }
  }

  int nsn = argc - optind;
  SNeData sne;
  sne.resize(nsn);
  std::vector<double> distances(nsn);
  for (int i = 0; i < nsn; ++i) {
    redshift = atof( argv[optind+i] );
    sne[i].zhel = redshift;
    sne[i].zcmb = redshift;
  }
  if (komatsuform) lm.setUseKomatsuForm();
  lm.setAtrans(atrans);
  lm.getLumDist( sne, distances, om, ode, w0, wa, lumdist::AUTO  );

  double cinvH = 2997.92458 / h; 
  printf("For w0: %8.5f wa: %8.5f om: %8.5f ode: %8.5f\n",
	 w0,wa,om,ode );
  for (int i = 0; i < nsn; ++i) {
    //Convert this to Mpc
    double dmpc = cinvH * pow( 10.0, 0.2 * distances[i] );
    if (logflag) {
      dmpc = 5 * log10( dmpc );
      printf(" z: %8.5f 5 log(dl): %9.6f\n",sne[i].zhel,dmpc );
    } else {
      printf(" z: %8.5f dl: %12.6f Mpc\n",sne[i].zhel,dmpc );
    }
  }
  return 0;
}
