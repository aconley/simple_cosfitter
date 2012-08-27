//Command line program for calculating w/om probability
// grid from growth of structure parameter

#include <iostream>
#include <string>
#include <getopt.h>
#include <stdlib.h>

#include <w2df.h>
#include <cosgrids.h>
#include <cosfitterexcept.h>

using namespace std;

int main(int argc, char **argv) {

  int c;

  bool ingridset, binaryfile;
  int nw, nom;
  double minw, maxw, minom, maxom;
  double fval, dfval, zval;
  string ingridfile, outfile;

  //Defaults (from Hawkins et al.)
  zval = 0.15;
  fval = 0.5096;
  dfval = 0.1080;

  //Other defaults
  ingridset = false;
  binaryfile = false;
  nw = 50; minw = -3; maxw = 0;
  nom = 50; minom = 0; maxom = 1;

  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"zval", required_argument, 0, 'z'},
    {"fval", required_argument, 0, 'f'},
    {"dfval", required_argument, 0, 'd'},
    {"nw", required_argument, 0, '1'},
    {"minw",required_argument, 0, '2'},
    {"maxw",required_argument, 0, '3'},
    {"nom", required_argument, 0, '4'},
    {"minom", required_argument, 0, '5'},
    {"maxom", required_argument, 0, '6'},
    {"ingrid", required_argument, 0, '7'},
    {"binaryfile",no_argument,0,'b'},
    {0,0,0,0}
  };

  while ( ( c = getopt_long( argc, argv, "hz:f:d:1:2:3:4:5:6:7:b",
			     long_options, &option_index ) ) != -1 ) 
    switch (c) {
    case 'h' :
      cerr << "NAME\n";
      cerr << "\tw2df_grid -- Produces 2dFGRS Omega_M, w contours assuming\n";
      cerr << "\t a flat Universe\n";
      cerr << "SYNOPSIS\n";
      cerr << "\tw2df_grid [ -z | --zval=ZVAL ] [ -f | --fval=FVAL ]" << endl;
      cerr << "\t [ --d | --dfval=DFVAL ] [ --nw=NW ] [ --minw=MINW ]" << endl;
      cerr << "\t [ --maxw=MAXW ] [ --nom=NOM ] [ --minom=MINOM ]\n";
      cerr << "\t [ --maxom=MAXOM ] [ --ingrid=INGRID ] [ -b, --binaryfile ] " 
	   << endl;
      cerr << "\t output_file\n";
      cerr << "DESCRIPTION\n";
      cerr << "\tThis program is designed to calculate the probability surface\n";
      cerr << "\tof the 2dFGRS results in Omega_m,.w space, assuming a flat\n";
      cerr << "\tUniverse and a constant value for the dark energy equation of\n";
      cerr << "\tstate parameter w.  The final argument is the name of an output\n";
      cerr << "\tfile, which is written using the cosgrid2D file format.\n";
      cerr << "\tIf binaryfile is set, the output file will be in binary\n";
      cerr << "\tformat -- and the input file must also be binary if ingrid\n";
      cerr << "\tis used in conjunction with binaryfile.\n";
      cerr << "OPTIONS:\n";
      cerr << "\t-z, --zval ZVAL\n\t\tRedshift of measurment (def: 0.15)\n";
      cerr << "\t-f, --fval FVAL\n\t\tGrowth parameter f (d\\ln D/d\\ln a) " << 
	"(def: 0.5096)\n";
      cerr << "\t-d, --dval=DFVAL\n\t\tError in fval (def: 0.108)\n";
      cerr << "\t--nw NW\n\t\tNumber of w bins (def: 50)\n";
      cerr << "\t--minw MINW\n\t\tMinimum w value (def: -3)\n";
      cerr << "\t--maxw MAXW\n\t\tMaximum w value (def: 0)\n";
      cerr << "\t--nom NOM\n\t\tNumber of Omega_m bins (def: 50)\n";
      cerr << "\t--minom MINOM\n\t\tMinimum Omega_m (def: 0)\n";
      cerr << "\t--maxom MAXOM\n\t\tMaximum Omega_m (def: 1)\n";
      cerr << "\t--ingrid INGRID\n\t\tName of file containing cosgrid2D" <<
	"which can be used to set" << endl ;
      cerr << "\t\t parameter spacing. If set, nw,minw,maxw,etc. " << 
	"are ignored" << endl;
      cerr << "\t--b, --binaryfile\n\t\tOutput as binary file.  If ingrid "
	   << "is set, assume it is binary." << endl;
      return 0;
      break;
    case 'z' :
      zval = atof( optarg );
    case 'f' : 
      fval = atof( optarg );
      break;
    case 'd' :
      dfval = atof( optarg );
      break;
    case '1' :
      nw = atoi( optarg );
      break;
    case '2' :
      minw = atof( optarg );
      break;
    case '3' : 
      maxw = atof( optarg );
      break;
    case '4' :
      nom = atoi( optarg );
      break;
    case '5' :
      minom = atof( optarg );
      break;
    case '6' :
      maxom = atof( optarg );
      break;
    case '7' :
      ingridfile = string( optarg );
      ingridset = true;
      break;
    case 'b' :
      binaryfile = true;
      break;
    }
  
  if (optind > argc - 1) {
    cerr << "w2df_grid :: Required arguments not provided\n";
    return 1;
  }
  
  outfile = string(argv[optind]); 

  cosgrid2D outgrid;
  try {
    if (ingridset) {
      cosgrid2D ingrid;
      ingrid.readFile( ingridfile, binaryfile );
      outgrid = w2dFGRS::make_prob_grid(fval,dfval,zval,ingrid);
    } else {
      outgrid = w2dFGRS::make_prob_grid(fval,dfval,zval,nw,minw,maxw,
					nom,minom,maxom);
    }
  } catch (const CosFitterExcept& ex) {
    cerr << "Error encountered in w2df_grid" << endl;
    cerr << ex;
    cerr << "Aborting\n";
    exit(1);
  }

  outgrid.writeFile( outfile, binaryfile );

  return 0;
}
