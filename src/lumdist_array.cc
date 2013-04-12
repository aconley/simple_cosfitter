#include <sstream>
#include <string>
#include <cstdio>
#include <algorithm>

#include "param_tags.h"
#include "lumdist_array.h"
#include "cosfitterexcept.h"
#include "utility.h"

/////////////////////////////////////////////////////////////
//                lumdist_array0D
 
lumdist_array0D::lumdist_array0D() {
  nz = 0;
  data = NULL;
  good = false;
}

/*!
  \param[in] sne      List of SN
  \param[in] fparam   parameters of fit
*/
lumdist_array0D::lumdist_array0D( const SNeData& sne,
				  const fitparam& fparam) {
  nz = 0;
  data = NULL;
  good = false;
  loadData(sne, fparam, false);
}

lumdist_array0D::lumdist_array0D( double om, double ode, double w0, double wa,
				  const SNeData& sne ) {
				  
  nz = 0;
  data = NULL;
  good = false;
  loadData(om, ode, w0, wa, sne, false);
}

lumdist_array0D::~lumdist_array0D() {
  free();
}

void lumdist_array0D::free() {
  if (data != NULL) {
    delete[] data;
  }
  data = NULL;
  good = false;
}

/*!
  \param[in] sn       List of SN
  \param[in] fparam   Parameters of fit
  \param[in] verbose  Show progress meter
*/
int lumdist_array0D::loadData( const SNeData& sn, const fitparam& fparam,
			       bool verbose ) {

  //Holds current params in order om, ode, w0, wa, with defaults
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

    //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for (it = fparam.params.begin(); it != fparam.params.end(); ++it) 
    if ( it->second.fit == param_struct::fixed )
      utility::updateStandardCosVec( it->first, fparam, current_cosparams, 
				     it->second.fixval );
  return loadData( current_cosparams[0], current_cosparams[1], 
		   current_cosparams[2], current_cosparams[3],
		   sn, verbose );
}

/*!
  \param[in] om       \f$\Omega_m\f$
  \param[in] ode      \f$\Omega_{DE}\f$
  \param[in] w0       \f$w_0\f$
  \param[in] wa       \f$w_a\f$
  \param[in] sn       List of SN
  \param[in] verbose  Show progress meter
*/
int lumdist_array0D::loadData( double om, double ode, double w0, double wa,
			       const SNeData& sn, bool verbose ) {
  free();

  nz = sn.size();

  //Fill in redshift arrays
  z_hel.resize(nz);
  for (unsigned int i = 0; i < nz; ++i)
    z_hel[i] = sn[i].zhel;
  z_cmb.resize(nz);
  for (unsigned int i = 0; i < nz; ++i)
    z_cmb[i] = sn[i].zcmb;
  
  //Allocate storage for lumdists
  data = new double[nz];
  
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = om; current_cosparams[1] = ode;
  current_cosparams[2] = w0; current_cosparams[3] = wa;

  //Loop
  lumdist lm;
  int st;
  std::vector<double> dl(nz);
  if (verbose) fprintf(stderr,"Filling in luminosity distances\n");
  st = lm.getLumDist( sn, dl, current_cosparams );
  good = st;
  if (st == 0)
    for (unsigned int j = 0; j < nz; ++j)
      data[j] = dl[j];
  return 0;
}




/////////////////////////////////////////////////////////////
//                lumdist_array1D
 
lumdist_array1D::lumdist_array1D() {
  n0 = nz = 0;
  data = NULL;
  good = NULL;
}

/*!
  \param[in] cosgrid  Grid of cosmological parameters
  \param[in] sne      List of SN
  \param[in] fparam   parameters of fit
*/
lumdist_array1D::lumdist_array1D( const cosgrid1D& cosgrid,
				  const SNeData& sne,
				  const fitparam& fparam) {
  n0 = nz = 0;
  data = NULL;
  good = NULL;
  loadData(cosgrid, sne, fparam, false);
}

lumdist_array1D::~lumdist_array1D() {
  free();
}

void lumdist_array1D::free() {
  if (data != NULL) {
    for (unsigned int i = 0; i < n0; ++i)
      delete[] data[i];
    delete[] data;
  }
  if (good != NULL)
    delete[] good;
  data = NULL;
  good = NULL;
}

unsigned int lumdist_array1D::getNAxis(unsigned int i) const {
  if (i==0) return n0;
  if (i==1) return nz;
  return 0;
}

/*!
  \param[in] cosgrid  Grid of cosmological parameters
  \param[in] sn       List of SN
  \param[in] fparam   Parameters of fit
  \param[in] verbose  Show progress meter
*/
int lumdist_array1D::loadData( const cosgrid1D& cosgrid,
			       const SNeData& sn, const fitparam& fparam,
			       bool verbose ) {
  free();

  unsigned int n0;
  n0 = cosgrid.getAxisN();
  nz = sn.size();

  //Fill in redshift arrays
  z_hel.resize(nz);
  for (unsigned int i = 0; i < nz; ++i)
    z_hel[i] = sn[i].zhel;
  z_cmb.resize(nz);
  for (unsigned int i = 0; i < nz; ++i)
    z_cmb[i] = sn[i].zcmb;
  
  //Allocate storage for lumdists
  data = new double*[n0];
  for (unsigned int i = 0; i < n0; ++i)
    data[i] = new double[nz];
  
  //And for good array
  good = new bool[n0];

  //We have to figure out which cosmological parameter
  // we have, as well as what the other, possibly fixed,
  // params are
  param_tags::paramcodes code0 = cosgrid.getAxisSpec().first;

  //Sanity checks
  if ( code0 != param_tags::omegam && code0 != param_tags::w0)
    throw CosFitterExcept("lumdist_array1D","loadData",
			  "Only Omega-m or w0 supported as param for 1D",2);
  if ( ! fparam.fixcurv )
    throw CosFitterExcept("lumdist_array1D","loadData",
			  "Curvature must be fixed in 1D fit",2);

  //Holds current params in order om, ode, w0, wa, with defaults
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

    //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for (it = fparam.params.begin(); it != fparam.params.end(); ++it) 
    if ( it->second.fit == param_struct::fixed )
      utility::updateStandardCosVec( it->first, fparam, current_cosparams, 
				     it->second.fixval );

  //Progress bar stuff
  unsigned int progbar_stepsize = n0 / 10;
  std::stringstream progbar;
  if (verbose) progbar << "Progress 0%";

  if (progbar_stepsize == 0) progbar_stepsize = 1;

  //Loop
  lumdist lm;
  int st;
  std::vector<double> dl(nz);
  double cos0;
  double *arrptr;
  if (verbose) fprintf(stderr,"Filling in luminosity distances\n");
  for (unsigned int i=0; i < n0; ++i) {
    if ( verbose && i > 0 && i % progbar_stepsize == 0 ) {
      progbar << "..." << (i * 100 / n0)+1 << "\%";
      fprintf(stderr,"\r%s",progbar.str().c_str());
    }
    arrptr = data[i];
    cos0 = cosgrid.getAxisVal(i);
    utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
    st = lm.getLumDist( sn, dl, current_cosparams );
    good[i] = st;
    if (st == 0)
      for (unsigned int j = 0; j < nz; ++j)
	arrptr[j] = dl[j];
  }
  if (verbose) fprintf(stderr,"\n");
  return 0;
}


/////////////////////////////////////////////////////////////
//                lumdist_array2D
 
lumdist_array2D::lumdist_array2D() {
  n0 = n1 = nz = 0;
  data = NULL;
  good = NULL;
}

/*!
  \param[in] cosgrid  Grid of cosmological parameters
  \param[in] sne       List of SN
  \param[in] fparam    Parameters of fit
*/
lumdist_array2D::lumdist_array2D( const cosgrid2D& cosgrid,
				  const SNeData& sne,
				  const fitparam& fparam) {
  n0 = n1 = nz = 0;
  data = NULL;
  good = NULL;
  loadData(cosgrid, sne, fparam, false);
}

lumdist_array2D::~lumdist_array2D() {
  free();
}

void lumdist_array2D::free() {
  if (data != NULL) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j=0; j < n1; ++j)
	delete[] data[i][j];
      delete[] data[i];
    }
    delete[] data;
  }
  if (good != NULL) {
    for (unsigned int i = 0; i < n0; ++i)
      delete[] good[i];
    delete[] good;
  }
  data = NULL;
  good = NULL;
}

/*!
  \param[in] cosgrid  Grid of cosmological parameters
  \param[in] sn       List of SN
  \param[in] fparam   Parameters of fit
  \param[in] verbose  Print progress meter
*/
int lumdist_array2D::loadData( const cosgrid2D& cosgrid,
			       const SNeData& sn, const fitparam& fparam,
			       bool verbose ) {
  free();

  unsigned int n0, n1;
  n0 = cosgrid.getAxisN(0);
  n1 = cosgrid.getAxisN(1);
  nz = sn.size();

  //Fill in redshift arrays
  z_hel.resize(nz);
  for (unsigned int i = 0; i < nz; ++i)
    z_hel[i] = sn[i].zhel;
  z_cmb.resize(nz);
  for (unsigned int i = 0; i < nz; ++i)
    z_cmb[i] = sn[i].zcmb;
  
  //Allocate storage for lumdists
  data = new double**[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double*[n1];
    for (unsigned int j = 0; j < n1; ++j)
      data[i][j] = new double[nz];
  }
  
  //And for good array
  good = new bool*[n0];
  for (unsigned int i = 0; i < n0; ++i)
    good[i] = new bool[n1];

  //We have to figure out which cosmological parameters
  // we have, as well as what the other, possibly fixed,
  // params are
  param_tags::paramcodes code0 = cosgrid.getAxisSpec(0).first;
  param_tags::paramcodes code1 = cosgrid.getAxisSpec(1).first;

  //Sanity checks
  if ( ( ( code0 == param_tags::w0 && code1 == param_tags::omegam ) ||
	 ( code1 == param_tags::w0 && code0 == param_tags::omegam ) )
       && ! fparam.fixcurv ) 
    throw CosFitterExcept("lumdist_array2D","loadData",
			  "Must have fixed curvature for om/w0 fit",1);
  if ( ( ( code0 == param_tags::omegam && code1 == param_tags::omegade ) ||
	 ( code1 == param_tags::omegam && code0 == param_tags::omegade ) )
       && fparam.fixcurv ) 
    throw CosFitterExcept("lumdist_array2D","loadData",
			  "Can't have fixed curvature for om/ode fit",2);

  //Holds current params in order om, ode, w0, wa, with defaults
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

  //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for (it = fparam.params.begin(); it != fparam.params.end(); ++it) 
    if ( it->second.fit == param_struct::fixed )
      utility::updateStandardCosVec( it->first, fparam, current_cosparams, 
				     it->second.fixval );

  //Progress bar stuff
  unsigned int progbar_stepsize = n0 / 10;
  std::stringstream progbar;
  if (verbose) progbar << "Progress 0%";

  if (progbar_stepsize == 0) progbar_stepsize = 1;

  //Loop
  lumdist lm;
  std::vector<double> dl(nz);
  double cos0, cos1;
  int st;
  bool *bptr;
  double **arrptr1, *arrptr2;
  if (verbose) fprintf(stderr,"Filling in luminosity distances\n");
  for (unsigned int i=0; i < n0; ++i) {
    if ( verbose && i > 0 && i % progbar_stepsize == 0 ) {
      progbar << "..." << (i * 100 / n0)+1 << "\%";
      fprintf(stderr,"\r%s",progbar.str().c_str());
    }
    cos0 = cosgrid.getAxisVal(0,i);
    utility::updateStandardCosVec( code0, fparam, current_cosparams, cos0 );
    bptr = good[i];
    arrptr1 = data[i];
    for (unsigned int j=0; j < n1; ++j) {
      arrptr2 = arrptr1[j];
      cos1 = cosgrid.getAxisVal(1,j);
      utility::updateStandardCosVec( code1, fparam, current_cosparams, cos1 );
      st = lm.getLumDist( sn, dl, current_cosparams );
      bptr[j] = st;
      if (st == 0)
	for (unsigned int k = 0; k < nz; ++k)
	  arrptr2[k] = dl[k];
    }
  }
  if (verbose) fprintf(stderr,"\n");
  return 0;

}

//////////////////////////////////////////////////////////
//                lumdist_array3D
lumdist_array3D::lumdist_array3D() {
  n0 = n1 = n2 = nz = 0;
  data = NULL;
  good = NULL;
}

/*!
  \param[in] cosgrid  Grid of cosmological parameters
  \param[in] sne      List of SN
  \param[in] fparam   Parameters of fit
*/
lumdist_array3D::lumdist_array3D( const cosgrid3D& cosgrid,
				  const SNeData& sne,
				  const fitparam& fparam) {
  n0 = n1 = n2 = nz = 0;
  data = NULL;
  good = NULL;
  loadData(cosgrid, sne, fparam, false);
}

lumdist_array3D::~lumdist_array3D() {
  free();
}

void lumdist_array3D::free() {
  if (data != NULL) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j=0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k) {
	  delete[] data[i][j][k];
	}
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }
  if (good != NULL) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j)
	delete[] good[i][j];
      delete[] good[i];
    }
    delete[] good;
  }
  data = NULL;
  good = NULL;
}

/*!
  \param[in] cosgrid  Grid of cosmological parameters
  \param[in] sn       List of SN
  \param[in] fparam   Parameters of fit
  \param[in] verbose  Print progress meter
*/
int lumdist_array3D::loadData( const cosgrid3D& cosgrid,
			       const SNeData& sn, const fitparam& fparam,
			       bool verbose ) {
  free();

  unsigned int n0, n1, n2;
  n0 = cosgrid.getAxisN(0);
  n1 = cosgrid.getAxisN(1);
  n2 = cosgrid.getAxisN(2);
  nz = sn.size();

  //Fill in redshift arrays
  z_hel.resize(nz);
  for (unsigned int i = 0; i < nz; ++i)
    z_hel[i] = sn[i].zhel;
  z_cmb.resize(nz);
  for (unsigned int i = 0; i < nz; ++i)
    z_cmb[i] = sn[i].zcmb;
  
  //Allocate storage for lumdists
  data = new double***[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double**[n1];
    for (unsigned int j = 0; j < n1; ++j) {
      data[i][j] = new double*[n2];
      for (unsigned int k = 0; k < n2; ++k)
	data[i][j][k] = new double[nz];
    }
  }
  
  //And for good array
  good = new bool**[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    good[i] = new bool*[n1];
    for (unsigned int j = 0; j < n1; ++j)
      good[i][j] = new bool[n2];
  }

  //We have to figure out which cosmological parameters
  // we have, as well as what the other, possibly fixed,
  // params are
  std::vector< param_tags::paramcodes > codes(3);
  codes[0] = cosgrid.getAxisSpec(0).first;
  codes[1] = cosgrid.getAxisSpec(1).first;
  codes[2] = cosgrid.getAxisSpec(2).first;

  //Sanity checks
  std::vector< param_tags::paramcodes >::const_iterator itcode;
  if ( find( codes.begin(), codes.end(), param_tags::omegam ) != codes.end() &&
       find( codes.begin(), codes.end(), param_tags::w0 ) != codes.end() &&
       find( codes.begin(), codes.end(), param_tags::wa ) != codes.end() &&
       !fparam.fixcurv )
    throw CosFitterExcept("lumdist_array3D","loadData",
			  "Must have fixed curvature for om/w0/wa fit",2);
  if ( find( codes.begin(), codes.end(), param_tags::omegam ) != codes.end() &&
       find( codes.begin(), codes.end(), param_tags::omegade) != codes.end() &&
       find( codes.begin(), codes.end(), param_tags::w0 ) != codes.end() &&
       fparam.fixcurv )
    throw CosFitterExcept("lumdist_array3D","loadData",
			  "Can't have fixed curvature for om/ode/w0/ fit",4);

  //Holds current params in order om, ode, w0, wa, with defaults
  std::vector<double> current_cosparams(4);
  current_cosparams[0] = 1.0; current_cosparams[1] = 0.0;
  current_cosparams[2] = -1.0; current_cosparams[3] = 0.0;

    //Loop over the possible cosmo params, setting the fixed values
  // if we have some
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for (it = fparam.params.begin(); it != fparam.params.end(); ++it) 
    if ( it->second.fit == param_struct::fixed )
      utility::updateStandardCosVec( it->first, fparam, current_cosparams, 
				     it->second.fixval );

  //Progress bar stuff
  unsigned int progbar_stepsize = n0 / 10;
  std::stringstream progbar;
  if (verbose) progbar << "Progress 0%";

  if (progbar_stepsize == 0) progbar_stepsize = 1;

  //Loop
  lumdist lm;
  std::vector<double> dl(nz);
  double cos0, cos1, cos2;
  int st;
  bool **bptr1, *bptr2;
  double ***arrptr1, **arrptr2, *arrptr3;
  if (verbose) fprintf(stderr,"Filling in luminosity distances\n");
  for (unsigned int i=0; i < n0; ++i) {
    if ( verbose && i > 0 && i % progbar_stepsize == 0 ) {
      progbar << "..." << (i * 100 / n0)+1 << "\%";
      fprintf(stderr,"\r%s",progbar.str().c_str());
    }
    cos0 = cosgrid.getAxisVal(0,i);
    utility::updateStandardCosVec( codes[0], fparam, current_cosparams, cos0 );
    bptr1 = good[i];
    arrptr1 = data[i];
    for (unsigned int j=0; j < n1; ++j) {
      arrptr2 = arrptr1[j];
      bptr2 = bptr1[j];
      cos1 = cosgrid.getAxisVal(1,j);
      utility::updateStandardCosVec( codes[1], fparam, current_cosparams, 
				     cos1 );
      for (unsigned int k = 0; k < n2; ++k ) {
	arrptr3 = arrptr2[k];
	cos2 = cosgrid.getAxisVal(1,j);
	utility::updateStandardCosVec( codes[2], fparam, current_cosparams, 
				       cos2 );
	st = lm.getLumDist( sn, dl, current_cosparams );
	bptr2[j] = st;
	if (st == 0) 
	  for (unsigned int m = 0; m < nz; ++m)
	    arrptr3[m] = dl[m];
      }
    }
  }
  if (verbose) fprintf(stderr,"\n");
  return 0;

}
