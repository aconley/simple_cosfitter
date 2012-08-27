//User interface driver for baryonpeak

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

  double z, aval, aerr;
  double fnu, ns;
  std::string paramfile, outfile;
  bool userns, userfnu, userz, useraval, useraerr;

  //Defaults
  z = 0.35; //For Eisenstein
  userz = false;
  aval = 0.469; //For Eisenstein
  useraval = false;
  aerr = 0.017; //For Eisenstein
  useraerr = false;
  ns = 0.98; //Power index for eisenstein measurement
  userns = false;
  fnu = 0.0;
  userfnu = false;

  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"usekomatsu",no_argument,0,'k'},
    {"ns", required_argument, 0, 'n'},
    {"fnu", required_argument, 0, 'f'},
    {"z", required_argument, 0, 'z'},
    {"aval", required_argument, 0, 'a'},
    {"aerr", required_argument, 0, 'e'},
    {"atrans", required_argument, 0, 't'},
    { 0 , 0 , 0 , 0 } 
  };

  while ( ( c = getopt_long( argc, argv, "hn:f:z:a:e:",
			     long_options, &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      cerr << "NAME\n";
      cerr << "\tbaryonpeak_grid -- Produces SDSS baryon acoustic peak"<<endl;
      cerr << "\tprobability surface, based on Eisenstein et al. (2005)" << endl;
      cerr << "SYNOPSIS\n";
      cerr << "\tbaryonpeak_grid [ -f | --fnu=FNU ] [ -n | --ns=NS ] " 
	   << std::endl;
      cerr << "\t [ -z, --z=Z ] [ -a, --aval=AVAL ] [ --aerr=AERR ] " 
	   << std::endl;
      cerr << "\t paramfile output_file" << std::endl;
      cerr << "DESCRIPTION\n";
      cerr << "\tThis program is designed to calculate the probability surface\n";
      cerr << "\tof the SDSS Baryon Peak results in Omega_m, Omega_Lambda space.\n";
      cerr << "\tIt is based on the A parameter described in Eisenstein et al.\n";
      cerr << "\t2005, ApJ 633, 560.  The first argument should be a simple_cosfitter\n";
      cerr << "\tstyle parameter file which will be used to determine the parameters." << std::endl;
      cerr << "\tThe final argument should be the desired name of the output file.\n";
      cerr << "OPTIONS:\n";
      cerr << "\t-f, --fnu FNU" << std::endl;
      cerr << "\t\tNeutrinto mass energy fraction, using Goobar et al. (2006)."
	   << std::endl;
      cerr << "\t\tThis is only valid for small values. (Def: 0.0)" 
	   << std::endl;
      cerr << "\t-n, --ns NS\n\t\tAssumed power spectrum index for Eistenstein measurement"
	   << endl;
      cerr << "\t\tThe default value of A is actually 0.469 (ns/0.98)^{-0.35}."
	   << endl;
      cerr << "\t\t(def: 0.98)" << endl;
      cerr << "\t-z, --z Z\n\t\tRedshift of measurement (def: 0.35)" << endl;
      cerr << "\t-a, --aval AVAL" << std::endl;
      cerr << "\t\tA measurement to compare against (def: 0.469)" << endl;
      cerr << "\t--aerr AERR" << std::endl;
      cerr << "\t\tA measurement error (def: 0.017)" << endl;
      return 0;
      break;
    case 'f' :
      fnu = atof( optarg );
      userfnu = true;
      break;
    case 'n' : 
      ns = atof( optarg );
      userns = true;
      break;
    case 'z' :
      userz = true;
      z = atof( optarg );
      break;
    case 'a' :
      useraval = true;
      aval = atof( optarg );
      break;
    case 'e' :
      useraerr = true;
      aerr = atof( optarg );
      break;
    }
  
  if (optind > argc - 2) {
    cerr << "baryonpeak_grid :: Required arguments not provided\n";
    return 1;
  }

  paramfile = std::string( argv[optind] );
  outfile = std::string( argv[optind+1] ); 

  fitparam fparam;

  auxconstraint::baoE05 sdss_bao;
  if (userns) sdss_bao.SetNs( ns );
  if (userfnu) sdss_bao.SetFnu( fnu );
  if (userz) sdss_bao.SetZ( z );
  if (useraval) sdss_bao.SetA( aval );
  if (useraerr) sdss_bao.SetAErr( aerr );
  try {
    fparam.readFile( paramfile, true );
    if (fparam.fixcurv)
      sdss_bao.SetFixCurv( fparam.ocurv );
    if (fparam.usekomatsuform) {
      sdss_bao.SetUseKomatsuForm();
      sdss_bao.SetAtrans( fparam.atrans );
    }
    sdss_bao.makeProbSurface( fparam.params, outfile );
  } catch (const CosFitterExcept& ex) {
    cerr << "Error encountered making baryon constraints\n";
    cerr << ex << endl;
    cerr << "Aborting" << endl;
    exit(1);
  }
  return 0;
}
