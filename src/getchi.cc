//getchi.cc

//Gives the chisquare at a specific value of om/ol/alpha/beta,
//Unlike indivchi, this 1) Doesn't take scriptm as an argument,
// instead analytically marginalizing over it's value using
// the Goliath et al. 2001 trick 2) Doesn't output individual
// contributions to the chisq, since this is a slightly tricky
// issue using this method

#include <string>
#include <iostream>
#include <cstdlib>
#include <getopt.h>

#include "fitter.h"
#include "lumdist.h"
#include "cosfitterexcept.h"

using namespace std;

int main(int argc, char** argv) {

  string paramfile;
  double w = -1, om=1.0, ol=0.0, alpha=1.45, beta=4.1;
  int c;
  bool flatflag = false;
  
  //Option parsing.  Use getopt_long because simple getopt has
  // trouble with negative inputs
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"w", required_argument,0, '0'},
    {"om",required_argument,0, 'm'},
    {"ol",required_argument,0, 'l'},
    {"alpha", required_argument,0, 'a'},
    {"beta", required_argument, 0, 'b'},
    {"flat", no_argument, 0, 'f'},
    { 0 , 0 , 0 , 0 } 
  };

  while ( (c = getopt_long(argc, argv, "h0:m:l:a:b:f",
			   long_options,&option_index )) != -1 )
    switch(c) {
    case 'h' :
      cerr << "NAME" << endl;
      cerr << "\tgetchi" << endl;
      cerr << "SYNOPSIS" << endl;
      cerr << "\tgetchi [ --w=W ] [ -m | --om=OM ] [ -l | --ol=OL ]" <<
	endl;
      cerr << "\t [ -a | --alpha=ALPHA ] [ -b | --beta=BETA ]"
	   << " [ -f ] PARAMFILE" << endl;
      cerr << "DESCRIPTION" << endl;
      cerr << "\tThis program calculates the chisquare of a SN data set " <<
	"relative" << endl;
      cerr << "\tto the specified parameters, using analytical marginalization"
	   << endl;
      cerr << "\tover the nuisance parameter scriptM.  The PARAMFILE is a"
	   << endl;
      cerr << "\tsimple_cosfitter style parameter file used to determine " << 
	" where" << endl;
      cerr << "\tthe data is read from and the intrinsic dispersion.  Note " <<
	" that" << endl;
      cerr << "\tit in no way determines the values of the parameters used." <<
	endl;
      cerr << "OPTIONS:\n";
      cerr << "\t--w W\n\t\tDark energy equation of state (def: -1.0)" <<
	endl;
      cerr << "\t-m, --om OM\n\t\tDensity parameter of non-relativistic " <<
	"matter (OMEGA_M) (def: 1)" << endl;
      cerr << "\t-l, --ol OL\n\t\tDensity parameter of dark energy " <<
	"(OMEGA_LAMBDA) (def: 0)" << endl;
      cerr << "\t-a, --alpha ALPHA\n\t\tSlope of the stretch-luminosity " <<
	"relation (def: 1.45)" << endl;
      cerr << "\t-b, --beta BETA\n\t\tSlope of the colour-luminosity " <<
	"relation (def: 4.1)" << endl;
      cerr << "\t-f\n\t\tFlat universe flag, forcing OL = 1.0 - OM" << endl;
      cerr << endl; 
      exit(1);
      break;
    case '0' :
      w = atof(optarg);
      break;
    case 'm' :
      om = atof( optarg );
      break;
    case 'l' :
      ol = atof( optarg );
      break;
    case 'a' :
      alpha = atof( optarg );
      break; 
    case 'b' :
      beta = atof( optarg );
      break;
    case 'f' :
      flatflag = true;
      break;
    }

  if (optind > argc-1) {
    cerr << "ERROR: Param file not provided\n";
    exit(2);
  }
  paramfile = string( argv[optind] );
  
  if (flatflag) ol = 1.0 - om;

  cosfitter sf;
  try {
    sf.prep( paramfile );
    printf("For parameters\n");
    printf("w: %6.3f Om: %6.3f Ol: %6.3f Al: %7.4f Beta: %7.4f\n",
	   w,om,ol,alpha,beta);
    printf("Number of Objects: %u\n", sf.getNsn());
    printf("ScriptM Marginalized ChiSquare: %10.5f\n",
	   sf.getchi(w,om,ol,alpha,beta));
  } catch (const CosFitterExcept& ex) {
    cerr << "Error encountered in indivchi\n";
    cerr << ex << endl;
    exit(1);
  }
 
}
