#include <string>
#include <iostream>

#include <getopt.h>

#include "fitter.h"
#include "cosfitterexcept.h"

int main(int argc, char **argv) {
  
  int c;
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {0,0,0,0}
  };

  if ( argc == 1 ) {
    std::cerr << "Required argument paramfile not provided" << std::endl;
    return 1;
  }

  while ( ( c = getopt_long( argc, argv, "h",long_options,
			     &option_index ) ) != -1 )
    switch(c) {
    case 'h' :
      std::cerr << "NAME" << std::endl;
      std::cerr << "\trunfitter_twoscriptm -- Runs simple_cosfitter on a data set" 
		<< std::endl;
      std::cerr << "\twith two scriptm values." << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\trunfitter_twoscriptm PARAMFILE [PARAMFILE2 ...]" << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tCommand line routine for doing cosmology fits." 
		<< std::endl;
      std::cerr << "\tFits are completely controlled through the parameter "
		<< "file(s) PARAMFILE." << std::endl;
      std::cerr << "\tFor details of what should be in the parameter file, "
		<< "please" << std::endl;
      std::cerr << "\tconsult the Doxygen documentation (by running make docs)."
	   << std::endl; 
      std::cerr << "\tExamples can be found in the samples subdirectory."
	   << std::endl;
      return 0;
      break;
    }

  int nparamfiles = argc - 1;

  for(int i = 0; i < nparamfiles; ++i) {

    std::string paramfile = std::string( argv[optind+i] );

    cosfitter fitter;
    try {
      fitter.prep( paramfile );
    } catch (const CosFitterExcept& ex) {
      std::cerr << "Error initializing simple_cosfitter\n";
      std::cerr << ex << std::endl;
      std::cerr << "Aborting fit\n";
      return 1;
    }
    
    if (nparamfiles > 1 && (fitter.fparam.verbose || 
			    fitter.fparam.showprogbar))
      std::cout << "Doing fit to param file: " << paramfile << std::endl;

    try {
      fitter.dofit();
    } catch (const CosFitterExcept& ex) {
      std::cerr << "Error performing fit\n";
      std::cerr << ex << std::endl;
      std::cerr << "Aborting\n";
      return 2;
    }
  }

  return 0;

}
