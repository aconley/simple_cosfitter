//indivchi.cc

//Gives the chisquare at a specific value of om/ol/sm/alpha/beta,
// also providing the chisquare contribution of each object

#include <string>
#include <iostream>
#include <cstdlib>
#include <map>
#include <set>
#include <getopt.h>

#include <param_tags.h>
#include <param_results.h>
#include <fitter.h>
#include <lumdist.h>
#include <cosfitterexcept.h>

using namespace std;

int main(int argc, char** argv) {

  string paramfile;
  double intrinsicval = 0.0;
  double w0 = -1, wa = 0.0, om=1.0, ol=0.0, sm=23.9, al=1.45, beta=4.1;
  int c;
  bool flatflag = false, intrinsicflag = false;
  
  //Option parsing.  Use getopt_long because simple getopt has
  // trouble with negative inputs
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"w0", required_argument,0, '0'},
    {"wa", required_argument,0, '1'},
    {"om",required_argument,0, 'm'},
    {"ol",required_argument,0, 'l'},
    {"sm", required_argument, 0, 's'},
    {"alpha", required_argument,0, 'a'},
    {"beta", required_argument, 0, 'b'},
    {"flat", no_argument, 0, 'f'},
    {"intrinsic", required_argument, 0, 'i'},
    { 0 , 0 , 0 , 0 } 
  };

  while ( (c = getopt_long(argc, argv, "h0:1:m:l:s:a:b:fi:",
			   long_options,&option_index )) != -1 )
    switch(c) {
    case 'h' :
      cerr << "NAME" << endl;
      cerr << "\t indivchi" << endl;
      cerr << "SYNOPSIS" << endl;
      cerr << "\tindivchi [ --w0=W0 ] [ --wa=WA ] [ -m | --om=OM ] "
	   << "[ -l | --ol=OL ]" << std::endl;
      cerr << "\t [ -s | --sm=SM ] [ -a | --alpha=ALPHA ] [ -b | --beta=BETA ]"
	   << endl;
      cerr << "\t [ -f ] [ -i | --intrinsic=INTRINSIC ] PARAMFILE" << endl;
      cerr << "DESCRIPTION" << endl;
      cerr << "\tThis program calculates the chisquare of a SN data set " <<
	"relative" << endl;
      cerr << "\tto the specified parameters, also providing the individual"
	   << endl;
      cerr << "\tchisquare contributions of each event.  The PARAMFILE is a"
	   << endl;
      cerr << "\tsimple_cosfitter style parameter file used to determine " <<
	" where" << endl;
      cerr << "\tthe data is read from and the intrinsic dispersion.  Note " <<
	" that" << endl;
      cerr << "\tit in no way determines the values of the parameters used." <<
	endl;
      cerr << "OPTIONS:\n";
      cerr << "\t--w0 W0\n\t\tConstant dark energy equation of "
	   << "state (def: -1.0)" << std::endl;
      cerr << "\t--wa WA\n\t\tA derivative of dark energy equation of "
	   << "state (def: 0.0)" << std::endl;
      cerr << "\t-m, --om OM\n\t\tDensity parameter of non-relativistic " <<
	"matter (OMEGA_M) (def: 1)" << endl;
      cerr << "\t-l, --ol OL\n\t\tDensity parameter of dark energy " <<
	"(OMEGA_LAMBDA) (def: 0)" << endl;
      cerr << "\t-s, --sm SM\n\t\tScript-M, luminosity/Hubble constant " <<
	"nusiance param (SCRIPTM)\n\t\t(def: 23.9)" << endl;
      cerr << "\t-a, --alpha ALPHA\n\t\tSlope of the stretch-luminosity " <<
	"relation (def: 1.45)" << endl;
      cerr << "\t-b, --beta BETA\n\t\tSlope of the colour-luminosity " <<
	"relation (def: 4.1)" << endl;
      cerr << "\t-f\n\t\tFlat universe flag, forcing OL = 1.0 - OM" << endl;
      cerr << "\t-i, --intrinsic INTRINSIC\n\t\tAssumed intrinsic dispersion." 
	   << " Otherwise set to value" << std::endl;
      cerr << "\t\tfrom PARAMFILE" << endl;
      cerr << endl; 
      exit(1);
      break;
    case '0' :
      w0 = atof(optarg);
      break;
    case '1' :
      wa = atof(optarg);
      break;
    case 'm' :
      om = atof( optarg );
      break;
    case 'l' :
      ol = atof( optarg );
      break;
    case 's' :
      sm = atof( optarg );
      break;
    case 'a' :
      al = atof( optarg );
      break; 
    case 'b' :
      beta = atof( optarg );
      break;
    case 'f' :
      flatflag = true;
      break;
    case 'i' :
      intrinsicval = atof( optarg );
      intrinsicflag = true;
      break;
    }

  if (optind > argc-1) {
    cerr << "ERROR: Param file not provided\n";
    exit(2);
  }
  paramfile = string( argv[optind] );
  
  if (flatflag) ol = 1.0 - om;

  cosfitter sf;

  //We have to build a map of values
  std::map< param_tags:: paramcodes, param_results > results_map;
  param_results pres;
  pres.fit = param_struct::fixed;

  pres.name = param_tags::w0tag; pres.value = w0;
  results_map[ param_tags::w0 ] = pres;
  pres.name = param_tags::watag; pres.value = wa;
  results_map[ param_tags::wa ] = pres;
  pres.name = param_tags::omtag; pres.value = om;
  results_map[ param_tags::omegam ] = pres;
  pres.name = param_tags::odetag; pres.value = ol;
  results_map[ param_tags::omegade ] = pres;
  pres.name = param_tags::alphatag; pres.value = al;
  results_map[ param_tags::alpha ] = pres;
  pres.name = param_tags::betatag; pres.value = beta;
  results_map[ param_tags::beta ] = pres;
  pres.name = param_tags::scriptmtag; pres.value = sm;
  results_map[ param_tags::scriptm ] = pres;

  std::map<unsigned int, double> intrinsic;
  try {
    sf.prep( paramfile );
    if (! intrinsicflag ) intrinsic = sf.fparam.intrinsicdisp; else {
      //We need to know which data set entries are present
      std::set<unsigned int> dataset = sf.getDataSetList();
      for (std::set<unsigned int>::iterator pos = dataset.begin();
	   pos != dataset.end(); ++pos)
	intrinsic[ *pos ] = intrinsicval;
    }
    sf.indivchi( results_map, intrinsic );
  } catch (const CosFitterExcept& ex) {
    cerr << "Error encountered in indivchi\n";
    cerr << ex << endl;
    exit(1);
  }
 
}
