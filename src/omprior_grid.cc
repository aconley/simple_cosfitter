//User interface driver for omprior

#include<cstdlib>
#include<cmath>
#include<iostream>
#include<getopt.h>
#include<string>

#include <param_tags.h>
#include <paramfile.h>
#include <auxconstraint.h>
#include <cosfitterexcept.h>

int main(int argc, char **argv) {

  int c;

  double omval, omerr;
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
      std::cerr << "NAME\n";
      std::cerr << "\tomprior_grid -- Produces Omega_m prior grid"<<std::endl;
      std::cerr << "SYNOPSIS\n";
      std::cerr << "\tomprior_grid omval omerr paramfile output_file" << 
	std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tThis program is designed to calculate the probability "
		<< "surface" << std::endl;
      std::cerr << "\tfor an Omega_m prior.  The third argument should "
		<< "be a simple_cosfitter" << std::endl;
      std::cerr << "\tstyle parameter file which will be used to determine " 
		<< "the parameters." << std::endl;
      std::cerr << "\tof the grid.  The final argument should be the "
		<< "desired name of the" << std::endl;
      std::cerr << "\toutput file.  The first two arguments describe "
		<< "the Gaussian prior." << std::endl;
      return 0;
      break;
    }
  
  if (optind > argc - 4) {
    std::cerr << "baryonpeak_grid :: Required arguments not provided\n";
    return 1;
  }

  omval = atof( argv[optind] );
  omerr = atof( argv[optind+1] );
  paramfile = std::string( argv[optind+2] );
  outfile = std::string( argv[optind+3] ); 

  if (omval < 0.0) {
    std::cerr << "Error in omprior_grid: om value less than 0" 
	      << std::endl;
    exit(1);
  }
  if (omerr <= 0.0) {
    std::cerr << "Error in omprior_grid: om error less than or equal "
	      << "to zero" << std::endl;
    exit(2);
  }

  fitparam fparam;
  auxconstraint::omprior omp(omval, omerr,fparam.fixcurv,fparam.ocurv);

  try {
    fparam.readFile( paramfile, true );
    omp.makeProbSurface( fparam.params, outfile );
  } catch (const CosFitterExcept& ex) {
    std::cerr << "Error encountered making omprior constraints\n";
    std::cerr << ex << std::endl;
    std::cerr << "Aborting" << std::endl;
    exit(4);
  }
  return 0;
}
