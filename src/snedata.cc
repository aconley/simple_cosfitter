#include <string>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <stdexcept>

#include <utility.h>
#include <cosfitterexcept.h>
#include <snedata.h>

using namespace std;

/*!
  If it is necessary to calculate a large number of these, it is more
  efficient to do it by hand than to use this function.

  \param[in] alpha Width (stretch) correction coefficient
  \param[in] beta Colour correction coefficient
  \returns Corrected mangitude value
*/
double SNeDataEntry::getCorrectedMag(double alpha, double beta) const {
  return static_cast<double>(mag + alpha * (widthpar - 1.0) - beta * colourpar);
}

/*!
  If it is necessary to calculate a large number of these, it is more
  efficient to do it by hand than to use this function.

  \param[in] alpha Width (stretch) correction coefficient
  \param[in] beta Colour correction coefficient
  \returns Error in corrected magnitude, not including peculiar velocity
    or intrinsic dispersion
*/
double SNeDataEntry::getCorrectedMagError(double alpha, 
                                         double beta) const {
  double var = var_mag 
    + (alpha*alpha*var_widthpar) 
    + (beta*beta*var_colourpar) 
    + (2.0*alpha*cov_mag_widthpar) 
    - (2.0*beta*cov_mag_colourpar) 
    - (2.0*alpha*beta*cov_widthpar_colourpar);
  return sqrt(var);
}

/*!
  If it is necessary to calculate a large number of these, it is more
  efficient to do it by hand than to use this function

  \param[in] alpha Width (stretch) correction coefficient
  \param[in] beta Colour correction coefficient
  \returns Variance corrected magnitude, not including peculiar velocity
    or intrinsic dispersion
*/
double SNeDataEntry::getCorrectedMagVar(double alpha, 
                                       double beta) const {
  return var_mag + alpha*alpha*var_widthpar + 
    beta*beta*var_colourpar + 2*alpha*cov_mag_widthpar -
    2*beta*cov_mag_colourpar - 2*alpha*beta*cov_widthpar_colourpar;
}

/*!
  Makes sure that all the defined correlation coefficients are
  less than 1.
 */
bool SNeDataEntry::isCovValid() const {
  if ( (cov_mag_widthpar > 0) && (var_mag > 0) && (var_widthpar > 0) &&
       ( fabs(cov_mag_widthpar/sqrt(var_mag*var_widthpar)) > 1 ) ) return false;
  if ( (cov_mag_colourpar > 0) && (var_mag > 0) && (var_colourpar > 0) &&
       ( fabs(cov_mag_colourpar/sqrt(var_mag*var_colourpar)) > 1) ) 
    return false;
  if ( (cov_widthpar_colourpar > 0) && (var_widthpar > 0) && 
       (var_colourpar > 0) &&
       ( fabs(cov_widthpar_colourpar/sqrt(var_widthpar*var_colourpar)) > 1) )
			      return false;
  return true;
}

SNeData::SNeData() : Name(""), isdiagonal( true ), havefullcov( false ),
		     havewoodburycov(false) { }


SNeData::SNeData(const std::string& name) : Name(name), isdiagonal( true ),
					    havefullcov( false ), 
					    havewoodburycov( false ) { }

SNeData::SNeData(const std::string& name, 
		 const std::string &fileName) : Name(name), 
						isdiagonal( true ),
						havefullcov( false ), 
						havewoodburycov(false) {
  readData(fileName,false);
}

void SNeData::zcmbsort() {
  if (points.size() == 0) return; //Nothing to do
  if ( isdiagonal ) {
    sort ( points.begin(), points.end(), SNeSortByZcmb() );
  } else {
    //More complex case -- we have a covariance matrix
    // to drag around
    unsigned int nsn = points.size();

    SNeSortByZcmbIndex srt( *this );
    std::vector<unsigned int> indxarr( nsn );
    for (unsigned int i = 0; i < nsn; ++i) indxarr[i] = i;
    sort( indxarr.begin(), indxarr.end(), srt );

    //Sort data points
    std::vector<SNeDataEntry> newpnts( nsn );
    for (unsigned int i = 0; i < nsn; ++i)
      newpnts[i] = points[ indxarr[i] ];
    points = newpnts;

    //Reorder the covariance matricies
    if ( havefullcov ) covmatrix.reorder( indxarr );
    if ( havewoodburycov ) woodbury_covmatrix.reorder( indxarr );
  }
}

/*!
Each line in the data file should have the format:
snname zcmb zhel dz widthpar dwidthpar colourpar dcolourpar cov_mag_width cov_mag_colourpar cov_widthpar_colourpar [dataset]

snname is a string (quotes are not necessary, and in fact probably not
a good idea), but the others are all
doubleing point numbers.  There can't be any spaces in snname, or bad
things will happen.  The lines must have an entry for each of the
above parameters, although they can be set to zero.

dataset is optional, and is an integer.  This is to support different
 dispersions for different datasets, which is controlled via the
SNDISP option.  It is up to you to ensure that, if you specify datasets,
you specify a dispersion for each.

Lines that begin with # are ignored.

Any old data is eliminated, as is the covariance matrix if present.

Some examples can be found in the test subdirectory.

*/
void SNeData::readData(const std::string& FileName, bool verbose) { 
  ifstream fl(FileName.c_str());

  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Data file " << FileName << " not found";
    throw CosFitterExcept("SNeData","readData",errstrng.str(),2);
  }

  std::string line;
  std::vector<std::string> words;
  std::stringstream str("");
  SNeDataEntry sndata;

  double dz, dmag, dwidthpar, dcolourpar;
  unsigned int wsize;

  points.resize(0);

  //File reading loop
  while(fl) {
    getline(fl,line);
    if (line[0]=='#' || line[0]==';') continue; //Comment line
    utility::stringwords(line,words);
    wsize = words.size();
    if (wsize < 13) continue; //Not enough entries on line
    
    if (wsize > 14) {
      std::stringstream errstrng("");
      errstrng << "Too many entries on line: " << line;
      throw CosFitterExcept("SNeData","readData",
			    errstrng.str(),4);
    }

    sndata.name = words[0];

    //Yes, I'm sure there is a more elegant way to do this.
    str.str(words[1]); str.clear(); str >> sndata.zcmb;
    str.str(words[2]); str.clear(); str >> sndata.zhel;
    str.str(words[3]); str.clear(); str >> dz;
    str.str(words[4]); str.clear(); str >> sndata.mag;
    str.str(words[5]); str.clear(); str >> dmag;
    str.str(words[6]); str.clear(); str >> sndata.widthpar;
    str.str(words[7]); str.clear(); str >> dwidthpar;
    str.str(words[8]); str.clear(); str >> sndata.colourpar;
    str.str(words[9]); str.clear(); str >> dcolourpar;
    str.str(words[10]); str.clear(); str >> sndata.cov_mag_widthpar;
    str.str(words[11]); str.clear(); str >> sndata.cov_mag_colourpar;
    str.str(words[12]); str.clear(); str >> sndata.cov_widthpar_colourpar;

    if (wsize == 14) {
      str.str(words[13]); str.clear(); str >> sndata.dataset;
    }

    sndata.var_z = dz*dz;
    sndata.var_mag = dmag*dmag;
    sndata.var_widthpar = dwidthpar*dwidthpar;
    sndata.var_colourpar = dcolourpar*dcolourpar;

    points.push_back(sndata);
  }
  
  if (points.size() == 0) throw CosFitterExcept("SNeData","readData",
						"No data read",8);

  if (verbose) cout << "Read: " << points.size() << " lines from " << 
		 FileName << endl;

  fl.close();

  //Clean out old covmatricies if we are rereading
  if ( !isdiagonal ) {
    isdiagonal = true;
    if (havefullcov) covmatrix.resize(0);
    if (havewoodburycov) woodbury_covmatrix.resize(0,0);
  }

  //Check covariance validity
  unsigned int badidx;
  if (!isCovValid(badidx)) {
    std::stringstream errstrng("");
    errstrng << "Found invalid correlation coefficient for "
	     << points[badidx].name;
    throw CosFitterExcept("SNeData","readData",
			  errstrng.str(),16);
  }

}

/*!
  \param[in] FileNames Vector of string filenames, matched by position
    to the six cov matricies.  An empty string will result in a skipped read.
*/
void SNeData::readCovData( const std::vector<std::string>& FileNames) {
  covmatrix.readCovData( FileNames );
  unsigned int n = covmatrix.getNsn();
  if ( n != 0 ) {
    if ( n != points.size() )
      throw CosFitterExcept("SNeData","readCovData","Cov matrix wrong size",1);
    isdiagonal = false;
    havefullcov = true;
  }
}

/*!
  Read in a Woodbury style covariance matrix.  Replaces the existing one.
*/
void SNeData::readWoodburyCovData(const std::string& FileName0,
				  const std::string& FileNamea,
				  const std::string& FileNameb) { 
  woodbury_covmatrix.readFile( FileName0, FileNamea, FileNameb );
  unsigned int n = woodbury_covmatrix.getNsn();
  if ( n != 0 ) {
    if ( n != points.size() )
      throw CosFitterExcept("SNeData","readWoodburyCovData",
			    "Woodbury cov matrix wrong size",1);
    isdiagonal = false;
    havewoodburycov = true;
  }
}

/*!
  Uses the same format as SNeData::readData.
  \param[in] outfile File to write to
*/
void SNeData::writeData(const std::string& outfile) const {  
  //Since we are doing formatted output, and the c++ style formatted
  // output sucks out loud, use C style

  FILE *fp;
  fp = fopen(outfile.c_str(),"w");
  if (!fp) {
    std::stringstream errstrng("");
    errstrng <<  "Couldn't open file " << outfile << " for writing";
    throw CosFitterExcept("SNeData","writeFile",errstrng.str(),1);
  }

  std::string hdfmt="%-10s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s ";
  hdfmt += "%-11s %-11s %-11s";
  std::string fmt = "%-10s %8.6f %8.6f %8.6f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f ";
  fmt += "%11.8f %11.8f %11.8f";

  //Write the header
  fprintf(fp,hdfmt.c_str(),"name","zcmb","zhel","dz","mag","dmag",
	  "width","dwidth","colour","dcolour","cov_m_w","cov_m_c",
	  "cov_w_c");

  //And the data
  for (std::vector<SNeDataEntry>::const_iterator sn = points.begin(); 
       sn != points.end(); ++sn)
    fprintf(fp,fmt.c_str(),sn->name.c_str(),sn->zcmb,sn->zhel,sqrt(sn->var_z),
	    sn->mag,sqrt(sn->var_mag),sn->widthpar,sqrt(sn->var_widthpar),
	    sn->colourpar,sqrt(sn->var_colourpar),sn->cov_mag_widthpar,
	    sn->cov_mag_colourpar,sn->cov_widthpar_colourpar);

  fclose(fp);
}


/*!
  Use getCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrixSNe SNeData::getCovMatrixSNe() const {
  if (!havefullcov) return covMatrixSNe();
  return covmatrix;
}

/*!
  Use getWoodburyCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of woodbury_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrixWoodburySNe SNeData::getWoodburyCovMatrixSNe() const {
  if (!havewoodburycov) return covMatrixWoodburySNe();
  return woodbury_covmatrix;
}

/*!
  Note that if a covariance matrix (of either form)
  is present the cross terms between the old SN and the 
  new ones will be left blank.

  \param[in] a List to append to current list
  \returns Reference to current list
*/
SNeData& SNeData::operator+=(const SNeData& a) {

  if (&a == this) {
    //Self addition.  It's an interesting question if
    // we should respect the const reference to a or not --
    // here we just ignore the request
    return *this;
  }
   
  unsigned int oldsize1 = points.size();
  unsigned int oldsize2 = a.points.size();
  unsigned int newsize = oldsize1 + oldsize2;
  points.reserve( newsize );
  for (SNeData::const_iterator i = a.begin(); i != a.end(); ++i)
    points.push_back( *i );


  //Handle cov matricies
  if ( havefullcov || a.havefullcov ) {
    if ( ! havefullcov ) covmatrix.matchCovmats( a.covmatrix, oldsize1 );
    covmatrix.append( a.covmatrix );
  }
  if ( havewoodburycov || a.havewoodburycov ) {
    int status;
    status = 0;
    if ( havewoodburycov ) {
      if ( a.havewoodburycov ) {
	woodbury_covmatrix.append(a.woodbury_covmatrix, status);
      } else {
	//Only the current one has one, so append a zero one
	covMatrixWoodburySNe tmp( oldsize2, woodbury_covmatrix.getNu());
	woodbury_covmatrix.append(tmp, status);
      }
    } else {
      //Only the thing we are adding has one
      woodbury_covmatrix.clear();
      woodbury_covmatrix.resize( oldsize1, a.woodbury_covmatrix.getNu() );
      woodbury_covmatrix.append(a.woodbury_covmatrix, status);
    }
  }

  return *this;
}


/*!
  \param[in] indxarr Array of indicies into old array.  
  \returns A new SNeData with only the specified objects present

  No bounds checking is performed, so it is up to the caller
  to ensure that the indicies in indxarr are valid.
*/
SNeData SNeData::operator[](const std::vector<unsigned int>& indxarr) const {
  unsigned int nelements = indxarr.size();

  if (nelements == 0) return SNeData(Name);

  SNeData retval(Name);

  retval.points.resize( nelements );
  
  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    retval.points[i++] = points[indxarr[j]];

  if (havefullcov)
    covmatrix.select( retval.covmatrix, indxarr );
  if ( havewoodburycov ) 
    woodbury_covmatrix.select( retval.woodbury_covmatrix, indxarr );

  return retval;
}

/*!
  \param[in] namearr Vector of strings corresponding to names
    of SN
  \returns A new SNeData with only the specified objects present

  Error checking is done, and an out_of_range exception is thrown
  if any of the passed in names are not found.
*/
SNeData SNeData::operator[](const std::vector<std::string>& namearr) const {

  unsigned int nelements = namearr.size();

  if (nelements == 0) return SNeData(Name);  

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested
  std::map<std::string,unsigned int> name_index_map;
  unsigned int ncurrent = points.size();
  for (unsigned int i = 0; i < ncurrent; ++i)
    name_index_map[ points[i].name ] = i;

  std::map<std::string,unsigned int>::iterator endpoint = name_index_map.end();
  std::map<std::string,unsigned int>::iterator location;
  
  std::vector< unsigned int > indxarr;
  for (unsigned int i = 0, j = 0; i < nelements; ++i) {
    location = name_index_map.find( namearr[i] );
    if (location == endpoint) {
      //Not found
      std::stringstream s("");
      s << "Name " << namearr[i] << " not found in " << Name;
      throw out_of_range(s.str());
    }
    indxarr[j++] = location->second;
  }
  return (*this)[ indxarr ];
}

/*!
  \param[in] indxarr List of indicies to remove 

  Duplicate indicies are allowed, but no range checking is performed,
  so it is up to the caller to ensure that only valid indicies are passed.
*/
void SNeData::remove( const std::vector<int>& indxarr ) {

  unsigned int nelements = points.size();
  unsigned int nremove = indxarr.size();

  if (nremove == 0) return;
  
  std::vector<SNeDataEntry> newpoints(nelements-nremove);

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  std::vector<bool> included(nelements, true);

  for (unsigned int i =0; i < nremove; ++i)
    included[ indxarr[i] ] = false;
  
  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) newpoints[i++] = points[j];
  
  points = newpoints;
  
  if ( havefullcov || havewoodburycov ) {
    std::vector<unsigned int> indxarr_select( nelements-nremove );
    for (unsigned int i = 0, j = 0; i < nelements; ++i)
      if (included[i]) indxarr_select[j++] = i;
    if (havefullcov) {
      covMatrixSNe tmp;
      tmp.matchCovmats( covmatrix, nelements-nremove );
      covmatrix.select( tmp, indxarr_select );
      covmatrix = tmp;
    }
    if ( havewoodburycov ) {
      covMatrixWoodburySNe tmp( nelements-nremove, 
				woodbury_covmatrix.getNu() );
      woodbury_covmatrix.select( tmp, indxarr_select );
      woodbury_covmatrix = tmp;
    }
  }
}

/*!
  \param[in] namearr List of names to remove.
  \param[in] strict Complain if one of the passed in names is
     not found in the list.

  Duplicate indicies are allowed.  Invalid entries in namearr
  result in a out_of_range error being thrown unless strict
  is set to false
*/
void SNeData::remove( const std::vector<std::string>& namearr,
		      bool strict ) {

  unsigned int nelements = points.size();
  unsigned int nremove = namearr.size();

  if (nremove == 0) return; //Nothing to remove

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested, and because
  // we allow duplicates

  std::map<std::string,unsigned int> name_index_map;
  for (unsigned int i = 0; i < nelements; ++i)
    name_index_map[ points[i].name ] = i;

  std::map<std::string,unsigned int>::iterator endpoint = name_index_map.end();
  std::map<std::string,unsigned int>::iterator location;

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  std::vector<bool> included(nelements, true);

  for (std::vector<std::string>::const_iterator i = namearr.begin();
       i != namearr.end(); ++i) {
    location = name_index_map.find( *i );
    if (location == endpoint) {
      //Not found
      if (strict) {
	std::stringstream s("");
	s << "Name " << *i << " not found in " << Name;
	throw out_of_range(s.str());
      }
      //Otherwise just ignore
    } else included[ location->second ] = false;
  }

  unsigned int nremaining = count( included.begin(), included.end(), true );

  std::vector<SNeDataEntry> newpoints(nremaining);

  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) newpoints[i++] = points[j];

  points = newpoints;

  if ( havefullcov || havewoodburycov ) {
    std::vector<unsigned int> indxarr_select( nelements-nremove );
    for (unsigned int i = 0, j = 0; i < nelements; ++i)
      if (included[i]) indxarr_select[j++] = i;
    if (havefullcov) {
      covMatrixSNe tmp;
      tmp.matchCovmats( covmatrix, nelements-nremove );
      covmatrix.select( tmp, indxarr_select );
      covmatrix = tmp;
    }
    if ( havewoodburycov ) {
      covMatrixWoodburySNe tmp( nelements-nremove, 
				woodbury_covmatrix.getNu() );
      woodbury_covmatrix.select( tmp, indxarr_select );
      woodbury_covmatrix = tmp;
    }
  }
}


/*!
  \param[in] indxarr List of indicies to remove from returned copy

  Duplicate indicies are allowed, but no range checking is performed,
  so it is up to the caller to ensure that only valid indicies are passed.

*/
SNeData SNeData::copy_remove( const std::vector<int>& indxarr ) const {

  SNeData retval(Name);

  unsigned int nelements = points.size();
  unsigned int nremove = indxarr.size();

  if (nremove == 0) return retval = *this; //Nothing removed

  retval.points.resize(nelements-nremove);

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  std::vector<bool> included(nelements, true);

  for (unsigned int i =0; i < nremove; ++i)
    included[ indxarr[i] ] = false;
  
  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) retval.points[i++] = points[j];


  if ( havefullcov || havewoodburycov ) {
    std::vector<unsigned int> indxarr_select( nelements-nremove );
    for (unsigned int i = 0, j = 0; i < nelements; ++i)
      if (included[i]) indxarr_select[j++] = i;
    if (havefullcov) {
      retval.covmatrix.matchCovmats( covmatrix, nelements-nremove );
      covmatrix.select( retval.covmatrix, indxarr_select );
    }
    if ( havewoodburycov ) {
      covMatrixWoodburySNe tmp( nelements-nremove, 
				woodbury_covmatrix.getNu() );
      woodbury_covmatrix.select( tmp, indxarr_select );
      retval.woodbury_covmatrix = tmp;
    }
  }
  return retval;
}

/*!
  \param[in] namearr List of names to remove from returned copy
  \param[in] strict Complain if one of the passed in names is
     not found in the list.
  \returns Copy of SNeData with specified entries removed

  Duplicate indicies are allowed.  Invalid entries in namearr
  result in a out_of_range error being thrown unless strict
  is set to false.
*/
SNeData SNeData::copy_remove( const std::vector<std::string>& namearr,
			      bool strict) const {

  SNeData retval(Name);

  unsigned int nelements = points.size();
  unsigned int nremove = namearr.size();

  if (nremove == 0) return retval = *this; //Nothing removed

  //The indexing is a lot more complex in this case, especially
  // because we want to retain the order requested, and because
  // we allow duplicates

  std::map<std::string,unsigned int> name_index_map;
  for (unsigned int i = 0; i < nelements; ++i)
    name_index_map[ points[i].name ] = i;

  std::map<std::string,unsigned int>::iterator endpoint = name_index_map.end();
  std::map<std::string,unsigned int>::iterator location;

  //Make a vector to hold which elements to include
  //A bitset can't be used because we don't know the number of
  // elements ahead of time
  std::vector<bool> included(nelements, true);

  for (std::vector<std::string>::const_iterator i = namearr.begin();
       i != namearr.end(); ++i) {
    location = name_index_map.find( *i );
    if (location == endpoint) {
      //Not found
      if (strict) {
	std::stringstream s("");
	s << "Name " << *i << " not found in " << Name;
	throw out_of_range(s.str());
      }
      //Otherwise just ignore this one
    } else included[ location->second ] = false;
  }

  unsigned int nremaining = count( included.begin(), included.end(), true );

  retval.points.resize(nremaining);

  for (unsigned int i = 0, j = 0; j < nelements; ++j)
    if ( included[j] ) retval.points[i++] = points[j];

  if ( havefullcov || havewoodburycov ) {
    std::vector<unsigned int> indxarr_select( nelements-nremove );
    for (unsigned int i = 0, j = 0; i < nelements; ++i)
      if (included[i]) indxarr_select[j++] = i;
    if (havefullcov) {
      retval.covmatrix.matchCovmats( covmatrix, nelements-nremove );
      covmatrix.select( retval.covmatrix, indxarr_select );
    }
    if ( havewoodburycov ) {
      covMatrixWoodburySNe tmp( nelements-nremove, 
				woodbury_covmatrix.getNu() );
      woodbury_covmatrix.select( tmp, indxarr_select );
      retval.woodbury_covmatrix = tmp;
    }
  }

  return retval;
}

std::set<unsigned int> SNeData::getDataSetList() const {
  std::set<unsigned int> retset;
  for (std::vector<SNeDataEntry>::const_iterator pos = points.begin();
       pos != points.end(); ++pos) 
    retset.insert( pos->dataset );
  return retset;
}

bool SNeData::isCovValid() const {
  if (points.size() == 0) return true;
  bool is_valid = true;
  for (std::vector<SNeDataEntry>::const_iterator it = points.begin();
       it != points.end(); ++it)
    is_valid &= it->isCovValid();
  return is_valid;
}


bool SNeData::isCovValid(unsigned int& badidx) const {
  if (points.size() == 0) return true;
  bool is_valid = true;
  bool curr_valid;
  badidx = points.size();
  for (unsigned int i = 0; i < points.size(); ++i) {
    curr_valid = points[i].isCovValid();
    is_valid &= curr_valid;
    if (!curr_valid) badidx = i;
  }
  return is_valid;
}
