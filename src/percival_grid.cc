//User interface driver for percival BAO parameters
// possibly also including WMAP5

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
  bool includewmap;

  includewmap = false;

  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"includewmap", no_argument, 0, 'i'},
    { 0 , 0 , 0 , 0 } 
  };

  while ( ( c = getopt_long( argc, argv, "hi",
			     long_options, &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      cerr << "NAME\n";
      cerr << "\tpercival_grid -- Produces a Percival BAO probability surface"
	   << std::endl;
      cerr << "\tbased on Percival et al. (2009)." << std::endl;
      cerr << "SYNOPSIS\n";
      cerr << "\tpercival_grid [ -i, --includewmap ] param_file output_file" 
	   << std::endl;
      cerr << "DESCRIPTION\n";
      cerr << "\tThis program is designed to calculate the probability surface\n";
      cerr << "\tof the Percival et al. (2009) BAO constraints.  The first"
	   << std::endl;
      cerr << "\targument should be the name of a simple_cosfitter style parameter file,"
	   << std::endl;
      cerr << "\twhich is used to set the parameters of the fit.  The last argument" << std::endl;
      cerr << "\tshould be the desired name of the output file." << std::endl;
      cerr << "OPTIONS" << std::endl;
      cerr << "\t-i, --includewmap" << std::endl;
      cerr << "\t\tAlso uses the WMAP9 constraints" << std::endl;
      return 0;
      break;
    case 'i' :
      includewmap = true;
      break;
    }
  
  if (optind > argc - 2) {
    cerr << "percival_grid :: Required arguments not provided\n";
    return 1;
  }

  paramfile = std::string( argv[optind] );
  outfile = std::string( argv[optind+1] ); 

  fitparam fparam;
  try {
    fparam.readFile( paramfile, true );
    if (includewmap) {
      auxconstraint::baoP09_wmap9yr_dls bao;
      if (fparam.fixcurv)
	bao.SetFixCurv( fparam.ocurv );
      if (fparam.usekomatsuform) {
	bao.SetUseKomatsuForm();
	bao.SetAtrans( fparam.atrans );
      }
      bao.makeProbSurface( fparam.params, outfile );
    } else {
      auxconstraint::baoP09 bao;
      if (fparam.fixcurv)
	bao.SetFixCurv( fparam.ocurv );
      if (fparam.usekomatsuform) {
	bao.SetUseKomatsuForm();
	bao.SetAtrans( fparam.atrans );
      }
      bao.makeProbSurface( fparam.params, outfile );
    }
  } catch (const CosFitterExcept& ex) {
    cerr << "Error encountered making Percival BAO constraints\n";
    cerr << ex << endl;
    cerr << "Aborting" << endl;
    exit(1);
  }
  return 0;
}
