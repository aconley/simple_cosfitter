//User interface driver for WMAP shift parameter

#include<cstdlib>
#include<cmath>
#include<iostream>
#include<getopt.h>
#include<string>

#include "param_tags.h"
#include "paramfile.h"
#include "auxconstraint.h"
#include "cosfitterexcept.h"

using namespace std;

int main(int argc, char **argv) {

  int c;

  std::string paramfile, outfile;
  bool usela;
  bool userh, userherror;
  double hval, hval_error;

  usela = userh = userherror = false;
  hval = 0.72; hval_error = 0.16;  //Note -- not used.

  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"usela",no_argument, 0, '1'},
    {"hmean",required_argument, 0, '2'},
    {"herror",required_argument, 0, '3'},
    { 0 , 0 , 0 , 0 } 
  };

  while ( ( c = getopt_long( argc, argv, "h12:3:",
			     long_options, &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      cerr << "NAME\n";
      cerr << "\twmap_grid -- Produces WMAP 3yr probability surface using"
	   << std::endl;
      cerr << "\tthe CMB shift parameter." << std::endl;
      cerr << "SYNOPSIS\n";
      cerr << "\twmap_grid param_file [ --usela ] [ --hmean=HMEAN ] "
	   << "[ --herror=HERROR ]" << std::endl;
      cerr << "\t output_file"<<endl;
      cerr << "DESCRIPTION\n";
      cerr << "\tThis program is designed to calculate the probability surface\n";
      cerr << "\tof the WMAP 3rd year results in the form of the the CMB shift"
	   << std::endl;
      cerr << "\tshift parameter as defined in Bond et al. (1997), with the"
	   << std::endl;
      cerr << "\tvalue given by Wang and Mukherjee (2006). The first argument"
	   << std::endl; 
      cerr << "\tshould be the name of a simple_cosfitter style parameter file, which"
	   << std::endl;
      cerr << "\tis used to set the parameters of the fit.  The last argument should be" << std::endl;
      cerr << "\tThe desired name of the output file." << std::endl;
      cerr << "OPTIONS" << std::endl;
      cerr << "\t--usela" << std::endl;
      cerr << "\t\tUse the second shift parameter l_a as well as R." 
	   << std::endl;
      cerr << "\t--hmean=HMEAN" << std::endl;
      cerr << "\t\tMean value for h prior, if using l_a. (Def: 0.72)"
	   << std::endl;
      cerr << "\t--herror=HERROR" << std::endl;
      cerr << "\t\tSigma value for h prior, if using l_a. (Def: 0.16)"
	   << std::endl;
      cerr << "NOTES" << std::endl;
      cerr << "\tOne should be careful using the shift parameter that ones use is" << std::endl;
      cerr << "\tconfined to the parameters that were considered when the shift" << std::endl;
      cerr << "\tvalue was determined.  For example, trying to determine w, Omega_m," << std::endl;
      cerr << "\tand Omega_DE at the same time using this gives a very unrealistic" << std::endl;
      cerr << "\tconstraint because they weren't included in the derivation."
	   << std::endl;
      return 0;
      break;
    case '1' :
      usela = true;
      break;
    case '2' :
      userh = true;
      hval = atof( optarg );
      break;
    case '3' :
      userherror = true;
      hval_error = atof( optarg );
      break;
    }
  
  if (optind > argc - 2) {
    cerr << "wmap_grid :: Required arguments not provided\n";
    return 1;
  }

  paramfile = std::string( argv[optind] );
  outfile = std::string( argv[optind+1] ); 

  fitparam fparam;
  auxconstraint::wmap3yr_dls wmap;

  if (usela) {
    wmap.SetUseLa();
    if (userh) wmap.SetHubblePriorMean(hval);
    if (userherror) wmap.SetHubblePriorSigma(hval_error);
  }

  try {
    fparam.readFile( paramfile, true );
    if (fparam.fixcurv)
      wmap.SetFixCurv( fparam.ocurv );
    if (fparam.usekomatsuform) {
      wmap.SetUseKomatsuForm();
      wmap.SetAtrans( fparam.atrans );
    }
    wmap.makeProbSurface( fparam.params, outfile );
  } catch (const CosFitterExcept& ex) {
    cerr << "Error encountered making WMAP constraints\n";
    cerr << ex << endl;
    cerr << "Aborting" << endl;
    exit(1);
  }
  return 0;
}
