//User interface driver for WMAP9 shift parameters

#include<cstdlib>
#include<cmath>
#include<iostream>
#include<getopt.h>
#include<string>

#include <param_tags.h>
#include <paramfile.h>
#include <auxconstraint.h>
#include <cosfitterexcept.h>

using namespace std;

int main(int argc, char **argv) {

  int c;

  std::string paramfile, outfile;

  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    { 0 , 0 , 0 , 0 } 
  };

  while ( ( c = getopt_long( argc, argv, "h",
			     long_options, &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      cerr << "NAME\n";
      cerr << "\twmap9_grid -- Produces WMAP 9yr probability surface using"
	   << std::endl;
      cerr << "\tthe CMB shift parameters l_a and R, and the redshift of" 
	   << std::endl;
      cerr << "\tdecopuling, z_star." << std::endl;
      cerr << "SYNOPSIS\n";
      cerr << "\twmap9_grid param_file output_file" << std::endl;
      cerr << "DESCRIPTION\n";
      cerr << "\tThis program is designed to calculate the probability surface\n";
      cerr << "\tof the WMAP 9th year results in the form of the the CMB shift"
	   << std::endl;
      cerr << "\tshift parameters as defined in Hinshaw et al. (2012). The first"
	   << std::endl;
      cerr << "\targument should be the name of a simple_cosfitter style parameter file,"
	   << std::endl;
      cerr << "\twhich is used to set the parameters of the fit.  The last argument" << std::endl;
      cerr << "\tshould be the desired name of the output file." << std::endl;
      cerr << "NOTES" << std::endl;
      cerr << "\tOne should be careful using the shift parameter that ones use is" << std::endl;
      cerr << "\tconfined to the parameters that were considered when the shift" << std::endl;
      cerr << "\tvalue was determined.  Right now, nonflat, w0, wa are supported"
	   << std::endl;
      return 0;
      break;
    }
  
  if (optind > argc - 2) {
    cerr << "wmap9_grid :: Required arguments not provided\n";
    return 1;
  }

  paramfile = std::string( argv[optind] );
  outfile = std::string( argv[optind+1] ); 

  fitparam fparam;
  auxconstraint::wmap9yr_dls wmap;

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
    cerr << "Error encountered making WMAP9 constraints\n";
    cerr << ex << endl;
    cerr << "Aborting" << endl;
    exit(1);
  }
  return 0;
}
