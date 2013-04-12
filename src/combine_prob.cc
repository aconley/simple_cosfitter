#include<iostream>
#include<vector>
#include<string>
#include<getopt.h>

#include "cosgrids.h"
#include "cosfitterexcept.h"

using namespace std;

int main(int argc, char **argv) {
  int c;
  bool normalize, binaryfile, silent;

  //Defaults
  normalize = true;
  silent = false;
  binaryfile = false;

  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"binaryfile", no_argument, 0, 'b' },
    {"nonormalize", no_argument, 0, 'n'},
    {"silent", no_argument, 0, 's'},
    {0,0,0,0}
  };

  while ( ( c = getopt_long( argc, argv, "hnb",long_options,
			     &option_index ) ) != -1 )
    switch (c) {
    case 'h' :
      cerr << "NAME" << endl;
      cerr << "\tcombine_prob -- Combined probability surfaces\n";
      cerr << "SYNOPSIS" << endl;
      cerr << "\tcombine_prob [ -n | --nonormalize ] [ -b | --binaryfile ]"
	   << endl; 
      cerr << "\t [ --silent ] INFILE1 [ INFILE2 ... ] OUTFILE" << endl;
      cerr << "DESCRIPTION" << endl;
      cerr << "\tTakes serialized contour1/2/3D probability grids and combines" <<
	endl;
      cerr << "\tthem via multiplication.  INFILES are a list of previously" <<
	endl;
      cerr << "\tprepared probability surfaces in the cosgrid file format" <<
	endl;
      cerr << "\twhich are multiplied and then output at OUTFILE." << endl;
      cerr << "\tBy default these are all assumed to be text files.  If they"
	   << endl;
      cerr << "\tare binary, then set the --binaryfile argument.  The output"
	   << endl;
      cerr << "\twill also then be in cosgrid binary format.\n";
      cerr << "OPTIONS" << endl;
      cerr << "\t-n, --nonormalize" << endl;
      cerr << "\t\tDon't renormalize the probabilities after multiplication" <<
	endl;
      cerr << "\t-b, --binaryfile" << endl;
      cerr << "\t\tInputs and output as binary files" << endl;
      cerr << "\t--silent" << endl;
      cerr << "\t\tRun in silent mode" << endl;
      return 0;
      break;
    case 'n' :
      normalize = false;
      break;
    case 's' :
      silent = true;
      break;
    case 'b' :
      binaryfile = true;
      break;
    }

  if (optind > argc - 2) {
    cerr << "ERROR in combine_prob -- not enough arguments provided\n";
    return 1;
  }

  unsigned int nin = static_cast<unsigned int>(argc - optind - 1);
  vector<std::string> infiles(nin);

  for (unsigned int i = 0; i < nin; ++i)
    infiles[i] = string( static_cast<char*>(argv[optind+i]));

  string outfile = string( static_cast<char*>(argv[optind+nin]) );

  unsigned int naxes0, cnaxes;
  try {
    naxes0 = NAxesInFile( infiles[0] );
    for (unsigned int i = 1; i < nin; ++i) {
      cnaxes = NAxesInFile( infiles[i] );
      if ( cnaxes != naxes0 ) {
	std::cerr << "Error -- file " << infiles[i] << " has " << cnaxes
		  << " axes" << std::endl;
	std::cerr << " not the same as " << infiles[0] << " which has "
		  << naxes0 << std::endl;
	return 2;
      }
    }
  } catch (CosFitterExcept& ex) {
    cerr << "ERROR in combine_prob: number of axes problem" << std::endl;
    cerr << ex << endl;
    return 1;
  }

  try {
    if (naxes0 == 1) {
      cosgrid1D ingrid, outgrid;
      if (!silent) cout << "Processing: " << infiles[0] << endl;
      outgrid.readFile( infiles[0], binaryfile );
      for (unsigned int i = 1; i < nin; ++i ) {
	if (!silent) cout << "Processing: " << infiles[i] << endl;
	ingrid.readFile( infiles[i], binaryfile );
	outgrid *= ingrid;
      }
      if (normalize) outgrid.normalize();
      outgrid.writeFile( outfile, binaryfile );      
    } else if (naxes0 == 2) {
      cosgrid2D ingrid, outgrid;
      if (!silent) cout << "Processing: " << infiles[0] << endl;
      outgrid.readFile( infiles[0], binaryfile );
      for (unsigned int i = 1; i < nin; ++i ) {
	if (!silent) cout << "Processing: " << infiles[i] << endl;
	ingrid.readFile( infiles[i], binaryfile );
	outgrid *= ingrid;
      }
      if (normalize) outgrid.normalize();
      outgrid.writeFile( outfile, binaryfile );
    } else if (naxes0 == 3) {
      cosgrid3D ingrid, outgrid;
      if (!silent) cout << "Processing: " << infiles[0] << endl;
      outgrid.readFile( infiles[0], binaryfile );
      for (unsigned int i = 1; i < nin; ++i ) {
	if (!silent) cout << "Processing: " << infiles[i] << endl;
	ingrid.readFile( infiles[i], binaryfile );
	outgrid *= ingrid;
      }
      if (normalize) outgrid.normalize();
      outgrid.writeFile( outfile, binaryfile );
    } else {
      std::cerr << "Unsupported number of axes: " << naxes0 << std::endl;
      return 4;
    }
  } catch (CosFitterExcept& ex) {
    cerr << "ERROR in combine_prob\n";
    cerr << ex << endl;
    return 1;
  }

  return 0;
}
