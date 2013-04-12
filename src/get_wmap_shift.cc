//User interface driver for giving the WMAP shift params

#include <cstdlib>
#include <iostream>
#include <getopt.h>
#include <string>

#include "auxconstraint.h"

using namespace std;

int main( int argc, char** argv ) {

  int c;
  
  bool userom, userode, flatflag, getla;
  double om, ode, w0, wa, h, obh2;
  om = 0.25;
  ode = 0.75;
  w0 = -1.0;
  wa = 0.0;
  obh2 = 0.022;
  h = 0.72;
  userom = false;
  userode = false;
  flatflag = false;
  getla = false;

  auxconstraint::wmap3yr_dls wmap;

  int option_index = 0;
  static struct option long_options[] = { 
    {"help", no_argument, 0, 'h'},
    {"flat", no_argument, 0, 'f'},
    {"h", required_argument, 0, '1' },
    {"omega_b_h2", required_argument, 0, '2'},
    {"omega_m",required_argument, 0, 'm'},
    {"omega_de",required_argument, 0, 'd'},
    {"w0",required_argument,0, '3'},
    {"wa",required_argument,0, '4'},
    {"getla", no_argument, 0, '5'},
    { 0 , 0 , 0 , 0 }
  };
  
  while ( ( c = getopt_long( argc, argv, "hfl5m:d:1:2:3:4:t:",
			     long_options, &option_index ) ) != -1 )
    switch (c) {
    case 'h' :
      cerr << "NAME\n";
      cerr << "\tget_wmap_shift -- Returns WMAP shift parameter(s)."<<endl;
      cerr << "SYNOPSIS\n";
      cerr << "\twmap_shift [ -m | --omega_m=OMEGA_MASS ] [ -d | --omega_de=OMEGA_DE]" << endl;
      cerr << "\t [ --w0=W0 ] [ --wa=WA ] [ --h=H ] [ --omega_b_h2=OMEGA_B_H2 ]" << std::endl;
      cerr << "\t [ -f ] [ --getla ] " << std::endl;
      cerr << "DESCRIPTION\n";
      cerr << "\tThis calculates the wmap shift parameter(s) given the cosmological" << endl;
      cerr << "\tparameters. All of the input parameter have defaults.  Relativistic " << std::endl; 
      cerr << " massive particles (neutrinos, etc.) are not supported." << std::endl;
      cerr << "OPTIONS" << endl;
      cerr << "\t-m, --omega_m OMEGA_MASS\n\t\tOmega_mass (more precisely, the density parameter" << endl;
      cerr << "\t\tof non-relativistic matter (def: 0.25)\n";
      cerr << "\t-d, --omega_de OMEGA_DE\n\t\tOmega_DarkEnergy (def: 0.75)" << endl;
      cerr << "\t--w0 W0\n\t\tThe constant equation of state parameter for the dark energy. (def: -1)" << endl;
      cerr << "\t--wa W0A\n\t\tThe a derivative of w (def: 0)." << endl;
      cerr << "\t--h H\n\t\tThe reduced Hubble constant (H0/100 km s^-1 Mpc^-1) (def: 0.72)" << endl;
      cerr << "\t --omega_b_h2 OMEGA_B_H2\n\t\tThe density of baryons (def: 0.022)" << std::endl;
      cerr << "\t\tthen the fit is treated as flat" << std::endl;
      cerr << "\t-f\n\t\tAssume a flat universe for the fit" << endl;
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
      obh2 = atof( optarg );
      break;
    case '3' :
      w0 = atof( optarg );
      break;
    case '4' :
      wa = atof( optarg );
      break;
    case '5' :
      getla = true;
      break;
    case 'f' :
      flatflag = true;
      break;
    }

  if (flatflag) {
    if (userode && userom && ( fabs(om + ode - 1.0) > 0.0001) ) {
      std::cerr << "get_wmap_shift :: Flat flag set but user specified " 
		<< "values do not allow" << std::endl;
      return 2;
    }
    if (userode) { om = 1.0 - ode; } else if (userom) { ode = 1.0 - om; }
  }

  double rval = wmap.rval(w0,wa,om,ode);
  printf("For w0: %8.5f wa: %8.5f om: %8.5f ode: %8.5f\n obh2: %8.5f h: %8.5f:\n",
	 w0,wa,om,ode,obh2,h );
  printf("R: %6.3f\n",rval);
  if (getla) {
    printf("l_a: %6.2f\n", wmap.laval(rval,w0,wa,om,ode,obh2,h) );
  }
  return 0;
}
