//paramfile.cc

//Stuff for reading parameter files
#include "paramfile.h"
#include "utility.h"
#include "cosfitterexcept.h"

#include <fstream>
#include <sstream>
#include <algorithm> /* For swap */

using namespace std;


/////////////////////////////////////////////////
//              param_struct                   //
/////////////////////////////////////////////////
param_struct::param_struct() : name(""), fit( param_struct::not_set ),
			       fixval(0.0), min(0.0), max(0.0), n(0), 
			       dval(0.0) {
  param_spec = std::pair< param_tags::paramcodes, param_tags::paramtype>
    ( param_tags::unknown, param_tags::unset );
}


/////////////////////////////////////
//            fitparam             //
/////////////////////////////////////
fitparam::fitparam() : ocurv(0.0), loop_order( fitparam::unknown ),
		       strategy( fitparam::not_set ), errsfixed( false ) {
  param_struct pstruct;

  //Only one we always have
  pstruct.name=param_tags::scriptmtag;
  pstruct.param_spec.first=param_tags::scriptm;
  pstruct.param_spec.second = param_tags::nuisance;
  pstruct.fit = param_struct::analytic; 
  params[ param_tags::scriptm ] = pstruct;

  ocurv = 0.0; //i.e., flat
  dzint = 0.0005;
  pecz = 0.001; //300 km/sec
  atrans = 0.1; //transition scale factor for Komatsu form

  nintrinsicdisp = 0;

  usekomatsuform = false;
  flatonlyfit = false;
  fixcurv = false;
  albetaout = false;
  extendedout = false;
  paramsummary = false;
  mag_covfileset = false;
  width_covfileset = false;
  colour_covfileset = false;
  magwidth_covfileset = false;
  magcolour_covfileset = false;
  widthcolour_covfileset = false;
  woodbury_covfileset = false;
  binaryout = false;
  verbose = false;
  showprogbar = false;

  outputfilename = std::string("");
}

/*!
  See main page of documentation for description of file format
  \param[in] paramfilename Name of file to read from
  \param[in] silent Don't print status text if if verbose is set in
    param file
*/
void fitparam::readFile(const std::string& paramfilename, bool silent) {
  ifstream ifp;
  stringstream str("");
  
  ifp.open(paramfilename.c_str());
  if (ifp.fail()) {
    stringstream errstrng("");
    errstrng << "Error opening parameter file: " << paramfilename;
    throw CosFitterExcept("cosfitter","readParamFile",
			  errstrng.str(),1);
  }

  //Input flags
  bool dataset, outfileset, outdataset;
  dataset = outfileset = outdataset = false;
  bool alpha_present, beta_present;

  //Read in and process lines
  std::string line, word0;
  vector<string> words;
  param_struct pstruct;
  while (!ifp.eof()) {
    getline(ifp,line);
    utility::stringwords(line,words);

    if (words.size() == 0) continue; //Nothing on line
    if (words[0][0] == '#') continue;  //Comment line

    word0 = utility::uppercase(words[0]);

    //This is a big if statement which handles each possible command
    // line option
    if (word0=="DATAFILE") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : DATAFILE filepath",2);
      }
      str.str(words[1]); str.clear(); str >> datafilename;
      dataset = true;
    } 
    else if (word0 == "MAGCOVFILE" || word0=="COVFILE" || 
	     word0 == "MAGCOVMATRIX" || word0=="COVMATRIX") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : MAGCOVFILE filepath",2);
      }
      str.str(words[1]); str.clear(); str >> mag_covfilename;
      mag_covfileset = true;
    }
    else if (word0 == "WIDTHCOVFILE" || word0 == "WIDTHCOVMATRIX" ||
	     word0 == "ALPHACOVFILE" || word0 == "ALPHA_COVFILE" ||
	     word0 == "ALCOVFILE" || word0 == "AL_COVFILE" ||
	     word0 == "ALPHACOVMATRIX" ) {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : WIDTHCOVFILE filepath",2);
      }
      str.str(words[1]); str.clear(); str >> width_covfilename;
      width_covfileset = true;
    }
    else if (word0 == "COLOURCOVFILE" || word0 == "COLOURCOVMATRIX" ||
	     word0 == "COLORCOVFILE" || word0 == "COLORCOVMATRIX" ||
	     word0 == "BETACOVFILE" || word0 == "BETA_COVFILE" ||
	     word0 == "BECOVFILE" || word0 == "BE_COVFILE" ||
	     word0 == "BETACOVMATRIX" ) {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : COLOURCOVFILE filepath",2);
      }
      str.str(words[1]); str.clear(); str >> colour_covfilename;
      colour_covfileset = true;
    }
    else if (word0 == "MAGWIDTHCOVFILE" || word0 == "MAGWIDTHCOVMATRIX" ) {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : MAGWIDTHCOVFILE filepath",2);
      }
      str.str(words[1]); str.clear(); str >> magwidth_covfilename;
      magwidth_covfileset = true;
    }
    else if (word0 == "MAGCOLOURCOVFILE" || word0 == "MAGCOLOURCOVMATRIX" ||
	     word0 == "MAGCOLORCOVFILE" || word0 == "MAGCOLORCOVMATRIX" ) {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : MAGCOLOURCOVFILE filepath",2);
      }
      str.str(words[1]); str.clear(); str >> magcolour_covfilename;
      magcolour_covfileset = true;
    }
    else if (word0 == "WIDTHCOLOURCOVFILE" || 
	     word0 == "WIDTHCOLOURCOVMATRIX" ||
	     word0 == "WIDTHCOLORCOVFILE" || 
	     word0 == "WIDTHCOLORCOVMATRIX" ) {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : WIDTHCOLOURCOVFILE filepath",2);
      }
      str.str(words[1]); str.clear(); str >> widthcolour_covfilename;
      widthcolour_covfileset = true;
    }
    else if (word0 == "WOODBURYCOVFILES" || word0 == "WOODBURY_COVFILES" ||
	     word0 == "WOODCOVFILES" || word0 == "WOOD_COVFILES") {
      if (words.size()!=4) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : WOODBURYCOVFILE file0 filea fileb",2);
      }
      str.str(words[1]); str.clear(); str >> woodbury0_covfilename;
      str.str(words[2]); str.clear(); str >> woodburya_covfilename; 
      str.str(words[3]); str.clear(); str >> woodburyb_covfilename;
      woodbury_covfileset = true;
    }
    else if (word0=="OUTFILE") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : OUTFILE filepath",2);
      }
      str.str(words[1]); str.clear(); str >> outputfilename;
      outfileset = true;
    }
    else if (word0=="W" || word0=="W0" || word0 == "W_0") {
      pstruct.name = param_tags::w0tag; 
      pstruct.param_spec.first=param_tags::w0;
      pstruct.param_spec.second = param_tags::cosmological;
      if (words.size() == 2) {
	pstruct.fit = param_struct::fixed;
	str.str(words[1]); str.clear(); str >> pstruct.fixval;
      } else if (words.size() == 4) {
	pstruct.fit = param_struct::loop;
	str.str(words[1]); str.clear(); str >> pstruct.min;
	str.str(words[2]); str.clear(); str >> pstruct.max;
	str.str(words[3]); str.clear(); str >> pstruct.n;
	if (pstruct.min > pstruct.max) swap(pstruct.min,pstruct.max);
      } else {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE: W ( min max num | val )",2);
      }
      params[ pstruct.param_spec.first ] = pstruct;
    }
    else if (word0=="WA" || word0 == "W_A") {
      pstruct.name = param_tags::watag; 
      pstruct.param_spec.first=param_tags::wa;
      pstruct.param_spec.second = param_tags::cosmological;
      if (words.size() == 2) {
	pstruct.fit = param_struct::fixed;
	str.str(words[1]); str.clear(); str >> pstruct.fixval;
      } else if (words.size() == 4) {
	pstruct.fit = param_struct::loop;
	str.str(words[1]); str.clear(); str >> pstruct.min;
	str.str(words[2]); str.clear(); str >> pstruct.max;
	str.str(words[3]); str.clear(); str >> pstruct.n;
	if (pstruct.min > pstruct.max) swap(pstruct.min,pstruct.max);
      } else {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE: WA ( min max num | val )",2);
      }
      params[ pstruct.param_spec.first ] = pstruct;
    }
    else if (word0=="OM" || word0=="OMEGAM" || word0 == "OMEGA_M") {
      pstruct.name = param_tags::omtag; 
      pstruct.param_spec.first=param_tags::omegam;
      pstruct.param_spec.second = param_tags::cosmological;
      if (words.size() == 2) {
	pstruct.fit = param_struct::fixed;
	str.str(words[1]); str.clear(); str >> pstruct.fixval;
      } else if (words.size() == 4) {
	pstruct.fit = param_struct::loop;
	str.str(words[1]); str.clear(); str >> pstruct.min;
	str.str(words[2]); str.clear(); str >> pstruct.max;
	str.str(words[3]); str.clear(); str >> pstruct.n;
	if (pstruct.min > pstruct.max) swap(pstruct.min,pstruct.max);
      } else {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE: OM ( min max num | val )",2);
      }
      params[ pstruct.param_spec.first ] = pstruct;
    }
    else if (word0 == "OL" || word0 == "OMEGALAMDA" || word0 == "OMEGA_LAMBDA" 
	     || word0 == "ODE" || word0 =="OMEGADE" || word0 == "OMEGA_DE" ) {
      pstruct.name = param_tags::odetag; 
      pstruct.param_spec.first=param_tags::omegade;
      pstruct.param_spec.second = param_tags::cosmological;
      if (words.size() == 2) {
	pstruct.fit = param_struct::fixed;
	str.str(words[1]); str.clear(); str >> pstruct.fixval;
      } else if (words.size() == 4) {
	pstruct.fit = param_struct::loop;
	str.str(words[1]); str.clear(); str >> pstruct.min;
	str.str(words[2]); str.clear(); str >> pstruct.max;
	str.str(words[3]); str.clear(); str >> pstruct.n;
	if (pstruct.min > pstruct.max) swap(pstruct.min,pstruct.max);
      } else {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE: ODE ( min max num | val )",2);
      }
      params[ pstruct.param_spec.first ] = pstruct;
    }
    else if (word0=="ALPHA" || word0 == "AL") {
      pstruct.name = param_tags::alphatag; 
      pstruct.param_spec.first=param_tags::alpha;
      pstruct.param_spec.second = param_tags::nuisance;
      if (words.size() == 2) {
	pstruct.fit = param_struct::fixed;
	str.str(words[1]); str.clear(); str >> pstruct.fixval;
      } else if (words.size() == 4) {
	pstruct.fit = param_struct::loop;
	str.str(words[1]); str.clear(); str >> pstruct.min;
	str.str(words[2]); str.clear(); str >> pstruct.max;
	str.str(words[3]); str.clear(); str >> pstruct.n;
	if (pstruct.min > pstruct.max) swap(pstruct.min,pstruct.max);
      } else {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE: ALPHA ( min max num | val )",2);
      }
      params[ pstruct.param_spec.first ] = pstruct;
    }
    else if (word0=="BETA" || word0 == "BE") {
      pstruct.name = param_tags::betatag; 
      pstruct.param_spec.first=param_tags::beta;
      pstruct.param_spec.second = param_tags::nuisance;
      if (words.size() == 2) {
	pstruct.fit = param_struct::fixed;
	str.str(words[1]); str.clear(); str >> pstruct.fixval;
      } else if (words.size() == 4) {
	pstruct.fit = param_struct::loop;
	str.str(words[1]); str.clear(); str >> pstruct.min;
	str.str(words[2]); str.clear(); str >> pstruct.max;
	str.str(words[3]); str.clear(); str >> pstruct.n;
	if (pstruct.min > pstruct.max) swap(pstruct.min,pstruct.max);
      } else {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE: BETA ( min max num | val )",2);
      }
      params[ pstruct.param_spec.first ] = pstruct;
    }
    else if (word0=="SNDISP") {
      if (words.size() < 2 || words.size() > 4) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : SNDISP [dataset] intrinsicsndispersion",2);
      }
      double idisp; 
      if (words.size() == 2) {
	str.str(words[1]); str.clear(); str >> idisp;
	intrinsicdisp[ 0 ] = idisp;
      } else {
	unsigned int dset;
	str.str(words[1]); str.clear(); str >> dset;
	str.str(words[2]); str.clear(); str >> idisp;
	intrinsicdisp[ dset ] = idisp;
      }
    }
    else if (word0=="PECZ") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : PECZ peculiar_z",2);
      }
      str.str(words[1]); str.clear(); str >> pecz;
    }
    else if (word0=="DZINT") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : DZINT z_integral_step",2);
      }
      str.str(words[1]); str.clear(); str >> dzint;
    }
    else if (word0=="FLATONLYFIT") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : FLATONLYFIT yes/no",2);
      }
      if (words[1] == "YES" || words[1] == "yes") {
	flatonlyfit = true;
	fixcurv = true;
	ocurv = 0.0;
      } else if (words[1] == "NO" || words[1] == "no") {
	flatonlyfit = false;
      } else {
	stringstream errstrng("");
	errstrng << "Unknown FLATONLYFIT specification: " << words[1];
	throw CosFitterExcept("cosfitter","readParamFile",
			      errstrng.str(),4);
      }
    }
    else if (word0=="FIXCURV" || word0 == "FIXCURVE") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : FIXCURV Omega_curv",2);
      }
      str.str(words[1]); str.clear(); str >> ocurv;
      fixcurv = true;
    } 
    else if (word0=="KOMATSUFORM" || word0=="KOMATSU" || 
	     word0=="USEKOMATSU" || word0=="USEKOMATSUFORM" ) {
      if (words.size() != 2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : KOMATSUFORM yes/no",2);
      }
      if (words[1] == "YES" || words[1] == "yes") {
	usekomatsuform = true;
      } else if (words[1] == "NO" || words[1] == "no") {
	usekomatsuform = false;
      } else {
	stringstream errstrng("");
	errstrng << "Unknown KOMATSUFORM specification: " << words[1];
	throw CosFitterExcept("cosfitter","readParamFile",
			      errstrng.str(),4);
      }     
    }
    else if (word0=="ATRANS") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : ATRANS atrans",2);
      }
      str.str(words[1]); str.clear(); str >> atrans;
    }
    else if (word0=="VERBOSE") {
      if (words.size() != 2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : VERBOSE yes/no",2);
      }
      if (words[1] == "YES" || words[1] == "yes") {
	verbose = true;
      } else if (words[1] == "NO" || words[1] == "no") {
	verbose = false;
      } else {
	stringstream errstrng("");
	errstrng << "Unknown VERBOSE specification: " << words[1];
	throw CosFitterExcept("cosfitter","readParamFile",
			      errstrng.str(),4);
      }     
    } 
    else if (word0=="SHOWPROGBAR") {
      if (words.size() != 2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : SHOWPROGBAR yes/no",2);
      }
      if (words[1] == "YES" || words[1] == "yes") {
	showprogbar = true;
      } else if (words[1] == "NO" || words[1] == "no") {
	showprogbar = false;
      } else {
	stringstream errstrng("");
	errstrng << "Unknown SHOWPROGBAR specification: " << words[1];
	throw CosFitterExcept("cosfitter","readParamFile",
			      errstrng.str(),4);
      }     
    } 
    else if (word0=="ALBETAOUTFILE" || word0 == "ALPHABETAOUTFILE" ||
	     word0=="ALPHABETA" ) {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : ALBETAOUTFILE filename",2);
      }
      str.str(words[1]); str.clear(); str >> albetaoutputfilename;
      albetaout = true;
    } 
    else if (word0 == "EXTENDEDOUTFILE" || word0 == "EXTENDEDOUT" ||
	     word0 == "EXTOUTFILE" || word0 == "EXTOUT") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : EXTENDEDOUTFILE filename",2);
      }
      str.str(words[1]); str.clear(); str >> extendedoutfilename;
      extendedout = true;
    }
    else if (word0 == "PARAMSUMMARYFILE" || word0 == "PARAMSUMMARY") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : PARAMSUMMARYFILE filename",2);
      }
      str.str(words[1]); str.clear(); str >> paramsummaryfilename;
      paramsummary = true;
    }
    else if (word0=="BINARYOUT") {
      if (words.size()!=2) {
	throw CosFitterExcept("cosfitter","readParamFile",
			      "PARAM USAGE : BINARYOUT yes/no",2);
      }
      if (words[1] == "YES" || words[1] == "yes") {
	binaryout = true;
      } else if (words[1] == "NO" || words[1] == "no") {
	binaryout = false;
      } else {
	stringstream errstrng("");
	errstrng << "Unknown BINARYOUT specification: " << words[1];
	throw CosFitterExcept("cosfitter","readParamFile",
			      errstrng.str(),4);
      }
    } 
    else {
      cerr << "WARNING : unparsed parameter " << words[0] << endl;
    }
  }
  ifp.close();

  any_covfile = mag_covfileset || width_covfileset || colour_covfileset ||
    magwidth_covfileset || magcolour_covfileset || widthcolour_covfileset ||
    woodbury_covfileset;

  //Some error checking
  if (! dataset) {
    string errstrng = "A data file must be specified";
    throw CosFitterExcept("cosfitter","readParamFile",errstrng,8);
  }

  //Verbosity/progbar interaction
  if (verbose) showprogbar = true;

  //If somebody only provided one value for alpha or beta, it's the
  // same as having it fixed
  std::map< param_tags::paramcodes, param_struct >::iterator italpha, itbeta;
  italpha = params.find( param_tags::alpha );
  if (italpha == params.end() ) alpha_present=false; else 
    alpha_present=true;
  if ( alpha_present && italpha->second.fit == param_struct::loop &&
       italpha->second.n < 2 ) {
    italpha->second.fit = param_struct::fixed;
    italpha->second.fixval = italpha->second.min;
  }
  itbeta = params.find( param_tags::beta );
  if (itbeta == params.end() ) beta_present=false; else 
    beta_present=true;
  if ( beta_present && itbeta->second.fit == param_struct::loop &&
       itbeta->second.n < 2 ) {
    itbeta->second.fit = param_struct::fixed;
    itbeta->second.fixval = itbeta->second.min;
  }

  //Check to see if albetaout makes sense
  if ( albetaout ) {
    if ( alpha_present && beta_present ) {
      if ( italpha->second.fit == param_struct::fixed &&
	   itbeta->second.fit == param_struct::fixed )
	throw CosFitterExcept("paramfile","readParamFile",
			      "Asking for alpha/beta output but both fixed",
			      128);
    } else if ( alpha_present ) {
      if ( italpha->second.fit == param_struct::fixed )
	throw CosFitterExcept("paramfile","readParamFile",
			      "Asking for alpha/beta output but no beta and alpha fixed",
			      128);
    } else if ( beta_present ) {
      if ( itbeta->second.fit == param_struct::fixed )
	throw CosFitterExcept("paramfile","readParamFile",
			      "Asking for alpha/beta output but no alpha and beta fixed",
			      128);
    } else {
      throw CosFitterExcept("paramfile","readParamFile",
			    "Asking for alpha/beta output but using neither",
			    128);
    }
  } 

  //Print a warning about using a full + woodbury cov file
  if ( (mag_covfileset || width_covfileset || colour_covfileset ||
	magwidth_covfileset || magcolour_covfileset ||
	widthcolour_covfileset ) && woodbury_covfileset ) {
    std::cout << "Warning: you have provided both a standard covariance"
	      << std::endl << " file and a woodbury one.  There will be"
	      << std::endl << " no performance improvement for the woodbury"
	      << std::endl;
  }

  //Decide if our errors are fixed ahead of time.  This is true
  // if we don't have alpha and beta, or if we do have them and
  // they are fixed.  Note we already tried to find them above
  if ( alpha_present ) {
    if ( beta_present ) {
      //Alpha and beta
      if (italpha->second.fit == param_struct::fixed &&
	  itbeta->second.fit == param_struct::fixed ) errsfixed=true; else
	errsfixed=false;
    } else {
      //No beta
      if (italpha->second.fit == param_struct::fixed) errsfixed=true; else
	errsfixed=false;
    }
  } else if (beta_present) {
    //We already know there is no alpha, so we don't need to check that
    if (itbeta->second.fit == param_struct::fixed) errsfixed=true; else
      errsfixed=false;
  } else errsfixed=true; //Have neither alpha nor beta

  //We will make the assumption later that having an alpha/beta covariance
  // matrix means you actually want to do something with them
  if ( width_covfileset && ! alpha_present )
    throw CosFitterExcept("cosfitter","readParamFile",
			  "Have width_covfile but not using alpha",32);
  if ( magwidth_covfileset && ! alpha_present )
    throw CosFitterExcept("cosfitter","readParamFile",
			  "Have magwidth_covfile but not using alpha",32);
  if ( colour_covfileset && ! beta_present )
    throw CosFitterExcept("cosfitter","readParamFile",
			  "Have colour_covfile but not using beta",32);
  if ( magcolour_covfileset && ! beta_present )
    throw CosFitterExcept("cosfitter","readParamFile",
			  "Have magcolour_covfile but not using beta",32);
  if ( widthcolour_covfileset && ! alpha_present )
    throw CosFitterExcept("cosfitter","readParamFile",
			  "Have widthcolour_covfile but not using alpha",32);
  if ( widthcolour_covfileset && ! beta_present )
    throw CosFitterExcept("cosfitter","readParamFile",
			  "Have widthcolour_covfile but not using beta",32);

  //Figure out which order the loops will be in -- if we have
  // covariance matricies then the outer loop should be on
  // alpha/beta since matrix inverson is the most expensive thing
  // to do (for reasonable numbers of SN).
  if ( any_covfile) loop_order = fitparam::cosmo_inner;
  else loop_order = cosmo_outer;

  //Now determine nuisance loop strategy, with a big if statment
  if (alpha_present) {
    if (italpha->second.fit == param_struct::fixed) {
      if (beta_present) {
	if (itbeta->second.fit == param_struct::fixed) 
	  strategy = fitparam::m_analytic;
	else strategy = fitparam::m_analytic_b_loop;
      } else strategy = fitparam::m_analytic;
    } else {
      if (beta_present) {
	if (itbeta->second.fit == param_struct::fixed) 
	  strategy = fitparam::m_analytic_a_loop;
	else strategy = fitparam::m_analytic_ab_loop;
      } else strategy = fitparam::m_analytic_a_loop;
    }
  } else if (beta_present) {
    //We know alpha isn't present
    if (itbeta->second.fit == param_struct::fixed)
      strategy = fitparam::m_analytic;
    else strategy = fitparam::m_analytic_b_loop;
  } else strategy = fitparam::m_analytic; //No alpha or beta
    
  //Now attach a string description
  switch (strategy) {
  case fitparam::not_set :
    throw CosFitterExcept("cosfitter","readParamFile",
			  "Didn't determine strategy for some reason",64);
    break;
  case fitparam::m_analytic : 
    strategy_string = "Analytic on script-M";
    break;
  case fitparam::m_analytic_a_loop :
    strategy_string = "Analytic on script-M, loop on alpha";
    break;
  case fitparam::m_analytic_b_loop :
    strategy_string = "Analytic on script-M, loop on beta";
    break;
  case fitparam::m_analytic_ab_loop :
    strategy_string = "Analytic on script-M, loop on alpha and beta";
    break;
  }

  //Figure out the step sizes
  std::map< param_tags::paramcodes, param_struct >::iterator it;
  for ( it = params.begin(); it != params.end(); ++it) 
    if (it->second.fit == param_struct::loop) {
      if (it->second.n == 1) {
	it->second.dval = 0.0;
      } else it->second.dval = (it->second.max - it->second.min)/
	       static_cast<double>(it->second.n - 1);
    }
  
  //Decide if we are calling ode omega_lambda or omega_de, and w or w0
  //Note that they will still be stored by odetag and w0tag in params,
  // but they will be written differently when output
  it = params.find( param_tags::w0 );
  if (it == params.end()) { //No w0
    it = params.find( param_tags::omegade );
    if (it != params.end() ) it->second.name = param_tags::oltag;
  } 
  //Now w or w0
  it = params.find( param_tags::wa );
  if (it == params.end()) { //No wa
    it = params.find( param_tags::w0 );
    if (it != params.end() ) it->second.name = param_tags::wtag;
  } 

  //Handle the intrinsic dispersions
  nintrinsicdisp = intrinsicdisp.size();
  if (nintrinsicdisp == 0) intrinsicdisp[0] = 0.0; //!< Default value

  //Now print out summary info
  if (verbose && !silent) {
    cout << "From file: " << paramfilename << endl; 
    cout << "PARAMETERS:\n------------------\n";
    for (it = params.begin(); it != params.end(); ++it) {
      if (it->second.fit == param_struct::loop) {
	std::cout << "  " << it->second.name << " = " << it->second.min 
		  << ".." << it->second.max << " , " << it->second.n 
		  << " steps" << std::endl;
      } else if (it->second.fit == param_struct::fixed) {
	std::cout << "  " << it->second.name << " = " << it->second.fixval 
		  << std::endl;
      } else if (it->second.fit == param_struct::analytic) {
	std::cout << "  " << it->second.name << " Analytic marginalization" 
		  << std::endl;
      }

    }
    std::cout << "  Nuisance parameter strategy: " << strategy_string 
	      << std::endl;
    if (nintrinsicdisp == 0) {
      std::cout << "  Intrinsic SN dispersion = "
		<< intrinsicdisp[0] << std::endl;
    } else {
      std::cout << "  Number of intrinsic dispersion datasets: "
		<< nintrinsicdisp << std::endl;
      for (std::map<unsigned int,double>::iterator pos = intrinsicdisp.begin();
	   pos != intrinsicdisp.end(); ++pos) {
	std::cout << "   Dataset: " << pos->first << " disp: "
		  << pos->second << std::endl;
      }
    }

    std::cout << "  Peculiar velocity z = " << pecz << std::endl;
    if (fixcurv) {
      if (flatonlyfit) {
	std::cout << "  Flat universe only" << std::endl;;
      } else {
	std::cout << "  Doing fit with Omega_curvature fixed at: " << ocurv 
	     << std::endl;
      }
    }
    if (usekomatsuform)
      std::cout << "  Using Komatsu form for w(a) with atrans: " << atrans
		<< std::endl;
    if (dataset)
      std::cout << "  Data file = " << datafilename << std::endl;
    if (mag_covfileset)
      std::cout << "  Mag covariance matrix file = " << mag_covfilename 
		<< std::endl;
    if (width_covfileset)
      std::cout << "  Width covariance matrix file = " << width_covfilename 
		<< std::endl;
    if (colour_covfileset)
      std::cout << "  Colour covariance matrix file = " << colour_covfilename 
		<< std::endl;
    if (magwidth_covfileset)
      std::cout << "  Mag-Width covariance matrix file = " 
		<< magwidth_covfilename << std::endl;
    if (magcolour_covfileset)
      std::cout << "  Mag-Colour covariance matrix file = " 
		<< magcolour_covfilename << std::endl;
    if (widthcolour_covfileset)
      std::cout << "  Width-Colour covariance matrix file = " 
		<< widthcolour_covfilename << std::endl;
    if (woodbury_covfileset) {
      std::cout << "  Woodbury m_B covariance matrix file = "
		<< woodbury0_covfilename << std::endl;
      std::cout << "  Woodbury s covariance matrix file = "
		<< woodburya_covfilename << std::endl;
      std::cout << "  Woodbury c covariance matrix file = "
		<< woodburyb_covfilename << std::endl;
    }
    if (outfileset)
      std::cout << "  Output file = " << outputfilename << std::endl;
    if (extendedout)
      std::cout << "  Extended output file = " << extendedoutfilename 
		<< std::endl;
    if (paramsummary)
      std::cout << "  Short parameter summary file = " << paramsummaryfilename 
		<< std::endl;
    if (albetaout) 
      std::cout << "  Output alpha/beta file = " << albetaoutputfilename 
		<< std::endl;
    if (binaryout)
      std::cout << "  Output probabilities will be in binary format" 
		<< std::endl;
    std::cout << std::endl;
  }
}
