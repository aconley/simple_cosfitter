//Routine for converting probability grids into FITS files

#include <fitsio.h>
#include <getopt.h>	
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>

#include <snedata.h>
#include <paramfile.h>
#include <param_tags.h>
#include <cosgrids.h>
#include <cosfitterexcept.h>
#include <utility.h>

using namespace std;

//Adds header keywords to fp reflecting fit parameters
int AddFitParams( fitsfile* fp, fitparam& fparam,
		  int &status ) {

  if (status > 0) return status;

  char keyval[FLEN_VALUE]; //Value of key if length is in doubt.
  std::map< param_tags::paramcodes, param_struct >::iterator it;

  int tmp;
  
  //Files
  strncpy( keyval, fparam.datafilename.c_str(), FLEN_VALUE - 9 );
  fits_write_key( fp, TSTRING,const_cast<char*>("DATAFILE"),keyval,
		  const_cast<char*>("Data file"),&status );

  if (fparam.mag_covfileset) {
    strncpy( keyval, fparam.mag_covfilename.c_str(), FLEN_VALUE - 19 );
    fits_write_key( fp, TSTRING,const_cast<char*>("MCOVFILE"),keyval,
		    const_cast<char*>("Mag Cov Matrix file"),&status );
  }
  if (fparam.width_covfileset) {
    strncpy( keyval, fparam.width_covfilename.c_str(), FLEN_VALUE - 22 );
    fits_write_key( fp, TSTRING,const_cast<char*>("WCOVFILE"),keyval,
		    const_cast<char*>("Width Cov Matrix file"),&status );
  }
  if (fparam.colour_covfileset) {
    strncpy( keyval, fparam.colour_covfilename.c_str(), FLEN_VALUE - 22 );
    fits_write_key( fp, TSTRING,const_cast<char*>("CCOVFILE"),keyval,
		    const_cast<char*>("Colour Cov Matrix file"),&status );
  }
  if (fparam.magwidth_covfileset) {
    strncpy( keyval, fparam.magwidth_covfilename.c_str(), FLEN_VALUE - 28 );
    fits_write_key( fp, TSTRING,const_cast<char*>("MWCOVFIL"),keyval,
		    const_cast<char*>("Mag-Width Cov Matrix file"),&status );
  }
  if (fparam.magcolour_covfileset) {
    strncpy( keyval, fparam.magcolour_covfilename.c_str(), FLEN_VALUE - 28 );
    fits_write_key( fp, TSTRING,const_cast<char*>("MCCOVFIL"),keyval,
		    const_cast<char*>("Mag-Colour Cov Matrix file"),&status );
  }
  if (fparam.widthcolour_covfileset) {
    strncpy( keyval, fparam.widthcolour_covfilename.c_str(), FLEN_VALUE - 28 );
    fits_write_key( fp, TSTRING,const_cast<char*>("WCCOVFIL"),keyval,
		    const_cast<char*>("Width-Colour Cov Matrix file"),&status );
  }


  strncpy( keyval, fparam.outputfilename.c_str(), FLEN_VALUE - 11 );
  fits_write_key( fp, TSTRING,const_cast<char*>("OUTFILE"),keyval,
		  const_cast<char*>("Output file"),&status );

  tmp = static_cast<int>(fparam.binaryout);
  fits_write_key( fp, TLOGICAL, const_cast<char*>("BINOUT"), &tmp,
		  const_cast<char*>("Binary output format"), &status );

  tmp = static_cast<int>(fparam.albetaout);
  fits_write_key( fp, TLOGICAL, const_cast<char*>("ALBETOUT"), &tmp,
		  const_cast<char*>("Alpha/beta output"), &status );
  if (fparam.albetaout) {
    strncpy( keyval, fparam.albetaoutputfilename.c_str(), FLEN_VALUE - 12 );
    fits_write_key( fp, TSTRING,const_cast<char*>("ALBTFILE"),keyval,
		    const_cast<char*>("al/beta file"),&status );
  }


  tmp = static_cast<int>(fparam.extendedout);
  fits_write_key( fp, TLOGICAL, const_cast<char*>("EXTOUT"), &tmp,
		  const_cast<char*>("Extended fit info output"), &status );
  if (fparam.extendedout) {
    strncpy( keyval, fparam.extendedoutfilename.c_str(), FLEN_VALUE - 20 );
    fits_write_key( fp, TSTRING,const_cast<char*>("EXOUTFL"),keyval,
		    const_cast<char*>("extended output file"),
		    &status );
  }

  //Type of fit
  it = fparam.params.find( param_tags::w0 );
  if ( it != fparam.params.end() ) tmp = 1; else tmp = 0;
  fits_write_key( fp, TLOGICAL, const_cast<char*>("WFIT"), &tmp,
		  const_cast<char*>("Was w_0 used"), &status );
  it = fparam.params.find( param_tags::wa );
  if ( it != fparam.params.end() ) tmp = 1; else tmp = 0;
  fits_write_key( fp, TLOGICAL, const_cast<char*>("WAFIT"), &tmp,
		  const_cast<char*>("Was w_a used"), &status );

  tmp = static_cast<int>(fparam.flatonlyfit);
  fits_write_key( fp, TLOGICAL, const_cast<char*>("FLATFIT"), &tmp,
		  const_cast<char*>("Flat universe fit"), &status );
  tmp = static_cast<int>(fparam.fixcurv);
  fits_write_key( fp, TLOGICAL, const_cast<char*>("FIXCURV"), &tmp,
		  const_cast<char*>("Fixed curvature universe fit"), &status );

  it = fparam.params.find( param_tags::alpha );
  if ( it != fparam.params.end() ) tmp = 1; else tmp = 0;
  fits_write_key( fp, TLOGICAL, const_cast<char*>("WIDTHCOR"), &tmp,
		  const_cast<char*>("Width-luminosity relation used"), 
		  &status );
  it = fparam.params.find( param_tags::beta );
  if ( it != fparam.params.end() ) tmp = 1; else tmp = 0;
  fits_write_key( fp, TLOGICAL, const_cast<char*>("CLRCORR"), &tmp,
		  const_cast<char*>("Colour-luminosity relation used"), 
		  &status );

  
  //Axes
  int itmp; 
  float ftmp;
  it = fparam.params.find( param_tags::w0 );
  if (it == fparam.params.end() ) {
    strncpy( keyval, "UNUSED", FLEN_VALUE - 20 );
    fits_write_key( fp, TSTRING,const_cast<char*>("W0TYPE"),keyval,
		    const_cast<char*>("Type of w_0 fit"),
		    &status );
  } else {
    switch (it->second.fit) {
    case param_struct::loop :
      strncpy( keyval, "LOOP", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("W0TYPE"),keyval,
		      const_cast<char*>("Type of w_0 fit"),
		      &status );
      ftmp = it->second.min;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MINW0"), &ftmp, 
		      const_cast<char*>("Min w_0"), 
		      &status );
      ftmp = it->second.max;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MAXW0"), &ftmp, 
		      const_cast<char*>("Max w_0"), 
		      &status );
      ftmp = it->second.dval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("DW0"), &ftmp, 
		      const_cast<char*>("Delta w_0"), 
		      &status );
      itmp = static_cast<int>(it->second.n);
      fits_write_key( fp, TINT, const_cast<char*>("NW0"), &itmp, 
		      const_cast<char*>("Number w_0"), 
		      &status );
      break;
    case param_struct::fixed :
      strncpy( keyval, "FIXED", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("W0TYPE"),keyval,
		      const_cast<char*>("Type of w_0 fit"),
		      &status );
      fits_write_key( fp, TFLOAT, const_cast<char*>("W0VAL"), 
		      &it->second.fixval,
		      const_cast<char*>("Fixed w_0 val"), &status );
      break;  
    default :
      break;
    }
  }

  it = fparam.params.find( param_tags::wa );
  if (it == fparam.params.end() ) {
    strncpy( keyval, "UNUSED", FLEN_VALUE - 20 );
    fits_write_key( fp, TSTRING,const_cast<char*>("WATYPE"),keyval,
		    const_cast<char*>("Type of w_a fit"),
		    &status );
  } else {
    switch (it->second.fit) {
    case param_struct::loop :
      strncpy( keyval, "LOOP", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("WATYPE"),keyval,
		      const_cast<char*>("Type of w_a fit"),
		      &status );
      ftmp = it->second.min;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MINWA"), &ftmp, 
		      const_cast<char*>("Min w_a"), 
		      &status );
      ftmp = it->second.max;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MAXWA"), &ftmp, 
		      const_cast<char*>("Max w_a"), 
		      &status );
      ftmp = it->second.dval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("DWA"), &ftmp, 
		      const_cast<char*>("Delta w_a"), 
		      &status );
      itmp = static_cast<int>(it->second.n);
      fits_write_key( fp, TINT, const_cast<char*>("NWA"), &itmp, 
		      const_cast<char*>("Number w_a"), 
		      &status );
      break;
    case param_struct::fixed :
      strncpy( keyval, "FIXED", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING, const_cast<char*>("WATYPE"),keyval,
		      const_cast<char*>("Type of w_a fit"),
		      &status );
      fits_write_key( fp, TFLOAT, const_cast<char*>("WAVAL"), 
		      &it->second.fixval,
		      const_cast<char*>("Fixed w_a val"), &status );
      break;  
    default :
      break;
    }
  }

  it = fparam.params.find( param_tags::omegam );
  if (it == fparam.params.end() ) {
    strncpy( keyval, "UNUSED", FLEN_VALUE - 20 );
    fits_write_key( fp, TSTRING,const_cast<char*>("OMTYPE"),keyval,
		    const_cast<char*>("Type of Omega_m fit"),
		    &status );
  } else 
    switch( it->second.fit ) {
    case param_struct::loop :
      strncpy( keyval, "LOOP", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("OMTYPE"),keyval,
		      const_cast<char*>("Type of Omega_m fit"),
		      &status );
      ftmp = it->second.min;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MINOM"), 
		      &ftmp, const_cast<char*>("Min Omega_m"), 
		      &status );
      ftmp = it->second.max;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MAXOM"), 
		      &ftmp, const_cast<char*>("Max Omega_m"), 
		      &status );
      ftmp = it->second.dval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("DOM"), 
		      &ftmp, const_cast<char*>("Delta Omega_m"), 
		      &status );
      itmp = static_cast<int>(it->second.n);
      fits_write_key( fp, TINT, const_cast<char*>("NOM"), &itmp, 
		      const_cast<char*>("Number Omega_m"), 
		      &status );
      break;
    case param_struct::fixed :
      strncpy( keyval, "FIXED", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("OMTYPE"),keyval,
		      const_cast<char*>("Type of Omega_m fit"),
		      &status );
      ftmp=it->second.fixval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("OMVAL"), 
		      &ftmp,
		      const_cast<char*>("Fixed Omega_m"), &status );
      break;
    default :
      break;
    }

  it = fparam.params.find( param_tags::omegade );
  if (it == fparam.params.end() ) {
    if (fparam.fixcurv) {
      strncpy( keyval, "DEPENDENT", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("ODETYPE"),keyval,
		      const_cast<char*>("Type of Omega_DE fit"),
		      &status );
    } else {
      strncpy( keyval, "UNUSED", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("ODETYPE"),keyval,
		      const_cast<char*>("Type of Omega_DE fit"),
		      &status );
    }
  } else 
    switch( it->second.fit ) {
    case param_struct::loop :
      strncpy( keyval, "LOOP", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("ODETYPE"),keyval,
		      const_cast<char*>("Type of Omega_DE fit"),
		      &status );
      ftmp = it->second.min;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MINODE"), 
		      &ftmp, const_cast<char*>("Min Omega_DE"), 
		      &status );
      ftmp = it->second.max;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MAXODE"), &ftmp,
		      const_cast<char*>("Max Omega_DE"), 
		      &status );
      ftmp = it->second.dval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("DODE"), &ftmp, 
		      const_cast<char*>("Delta Omega_DE"), 
		      &status );
      itmp = static_cast<int>(it->second.n);
      fits_write_key( fp, TINT, const_cast<char*>("NODE"), &itmp, 
		      const_cast<char*>("Number Omega_DE"), 
		      &status );
      break;
    case param_struct::fixed :
      strncpy( keyval, "FIXED", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("ODETYPE"),keyval,
		      const_cast<char*>("Type of Omega_DE fit"),
		      &status );
      ftmp = it->second.fixval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("ODEVAL"), 
		      &ftmp,
		      const_cast<char*>("Fixed Omega_m"), &status );
      break;
    default :
      break;
    }

  it = fparam.params.find( param_tags::alpha );
  if (it == fparam.params.end() ) {
    strncpy( keyval, "UNUSED", FLEN_VALUE - 20 );
    fits_write_key( fp, TSTRING,const_cast<char*>("ALPHTYPE"),keyval,
		    const_cast<char*>("Type of alpha fit"),
		    &status );
  } else 
    switch( it->second.fit ) {
    case param_struct::loop :
      strncpy( keyval, "LOOP", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("ALPHTYPE"),keyval,
		      const_cast<char*>("Type of alpha fit"),
		      &status );
      ftmp = it->second.min;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MINALPHA"), 
		      &ftmp, const_cast<char*>("Min alpha"), 
		      &status );
      ftmp = it->second.max;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MAXALPHA"), 
		      &ftmp, const_cast<char*>("Max alpha"), 
		      &status );
      ftmp = it->second.dval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("DALPHA"),
		      &ftmp, const_cast<char*>("Delta alpha"), 
		      &status );
      itmp = static_cast<int>(it->second.n);
      fits_write_key( fp, TINT, const_cast<char*>("NALPHA"), &itmp, 
		      const_cast<char*>("Number alpha"), 
		      &status );
      break;
    case param_struct::fixed :
      strncpy( keyval, "FIXED", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("ALPHTYPE"),keyval,
		      const_cast<char*>("Type of alpha fit"),
		      &status );
      ftmp = it->second.fixval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("ALPHAVAL"), 
		      &ftmp,
		      const_cast<char*>("Fixed alpha"), &status );
      break;
    default :
      break;
    }

  it = fparam.params.find( param_tags::beta );
  if (it == fparam.params.end() ) {
    strncpy( keyval, "UNUSED", FLEN_VALUE - 20 );
    fits_write_key( fp, TSTRING,const_cast<char*>("BETATYPE"),keyval,
		    const_cast<char*>("Type of beta fit"),
		    &status );
  } else 
    switch( it->second.fit ) {
    case param_struct::loop :
      strncpy( keyval, "LOOP", FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("BETATYPE"),keyval,
		      const_cast<char*>("Type of beta fit"),
		      &status );
      ftmp = it->second.min;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MINBETA"), 
		      &ftmp, const_cast<char*>("Min beta"), 
		      &status );
      ftmp = it->second.max;
      fits_write_key( fp, TFLOAT, const_cast<char*>("MAXBETA"), 
		      &ftmp, const_cast<char*>("Max beta"), 
		      &status );
      ftmp = it->second.dval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("DBETA"), 
		      &ftmp, const_cast<char*>("Delta beta"), 
		      &status );
      itmp = static_cast<int>(it->second.n);
      fits_write_key( fp, TINT, const_cast<char*>("NBETA"), 
		      &itmp, const_cast<char*>("Number beta"), 
		      &status );
      break;
    case param_struct::fixed :
      strncpy( keyval, const_cast<char*>("FIXED"), FLEN_VALUE - 20 );
      fits_write_key( fp, TSTRING,const_cast<char*>("BETATYPE"),
		      keyval,const_cast<char*>("Type of beta fit"),
		      &status );
      ftmp = it->second.fixval;
      fits_write_key( fp, TFLOAT, const_cast<char*>("BETAVAL"), 
		      &ftmp,
		      const_cast<char*>("Fixed beta"), &status );
      break;
    default :
      break;
    }

  if (fparam.fixcurv) 
    fits_write_key( fp, TFLOAT, const_cast<char*>("OCURV"), 
		    &fparam.ocurv, const_cast<char*>("Omega_curv"), 
		    &status );

  //Other stuff
  ftmp = fparam.pecz;
  fits_write_key( fp, TFLOAT, const_cast<char*>("PECZ"), &ftmp,
		  const_cast<char*>("peculiar velocity in z"), &status );

  return status;

}

//1D cosgrid case
int convert_to_fits1D( fitsfile* fp, const string& infile, bool binaryfile,
		       const string& extname, int& status ) {
  if (status > 0) return status;

  cosgrid1D cg;
  cg.readFile( infile, binaryfile );
  long axissize[1];
  axissize[0] = cg.getAxisN();

  fits_create_img( fp, DOUBLE_IMG, 1, axissize, &status );

  char keyval[FLEN_VALUE]; //Value of key if length is in doubt
    
  //Extension name
  strncpy( keyval, extname.c_str(), FLEN_VALUE );
  fits_write_key( fp, TSTRING, const_cast<char*>("EXTNAME"), keyval,
		  const_cast<char*>("Name of extension"), &status );

  //Add WCS information to header
  float crpix = 1, ftmp;
  strncpy( keyval, cg.getAxisLabel().c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		  const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		  const_cast<char*>("Ref pix of axis 1"), &status );
  ftmp = cg.getAxisMin();
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &ftmp, 
		  const_cast<char*>("val at ref pix"), &status );
  ftmp = cg.getAxisD();
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &ftmp, 
		  const_cast<char*>("delta along axis 1"), 
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &ftmp, 
		  const_cast<char*>("WCS matrix element 1 1"),
		  &status );
  
  //Do data writing
  long fpixel[1] = { 1 };
  fits_write_pix( fp, TDOUBLE, fpixel, cg.getAxisN(), cg.getData(), &status );

  return status;

}

int convert_to_fits2D( fitsfile* fp, const string& infile, bool binaryfile,
		       const string& extname, int& status ) {
  if (status > 0) return status;

  cosgrid2D cg;
  cg.readFile( infile, binaryfile );
  long axissize[2];
  axissize[0] = cg.getAxisN(0); axissize[1] = cg.getAxisN(1);

  fits_create_img( fp, DOUBLE_IMG, 2, axissize, &status );

  char keyval[FLEN_VALUE]; //Value of key if length is in doubt.

  strncpy( keyval, extname.c_str(), FLEN_VALUE );
  fits_write_key( fp, TSTRING, const_cast<char*>("EXTNAME"), keyval,
		  const_cast<char*>("Name of extension"), &status );
  fits_write_key( fp, TSTRING, const_cast<char*>("HDUNAME"), keyval,
		  const_cast<char*>("Name of extension"), &status );

  //Add WCS information to header
  float crpix = 1, zero = 0.0, ftmp;
  strncpy( keyval, cg.getAxisLabel(0).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		  const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		  const_cast<char*>("Ref pix of axis 1"), &status );
  ftmp = cg.getAxisMin(0);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRVAL1"), &ftmp, 
		  const_cast<char*>("val at ref pix"), &status );
  ftmp = cg.getAxisD(0);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRDELT1"), &ftmp, 
		  const_cast<char*>("delta along axis 1"), 
		  &status );
  strncpy( keyval, cg.getAxisLabel(1).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE2"), keyval,
		  const_cast<char*>("Type of Data axis 2"),&status);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX2"), &crpix, 
		  const_cast<char*>("Ref pix of axis 2"), &status );
  ftmp = cg.getAxisMin(1);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRVAL2"), &ftmp, 
		  const_cast<char*>("val at ref pix"), &status );
  ftmp = cg.getAxisD(1);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRDELT2"), &ftmp, 
		  const_cast<char*>("delta along axis 2"), 
		  &status );
  ftmp = cg.getAxisD(0);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CD1_1"), &ftmp, 
		  const_cast<char*>("WCS matrix element 1 1"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_2"), &zero, 
		  const_cast<char*>("WCS matrix element 1 2"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD2_1"), &zero, 
		  const_cast<char*>("WCS matrix element 2 1"),
		  &status );
  ftmp = cg.getAxisMin(1);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CD2_2"), &ftmp, 
		  const_cast<char*>("WCS matrix element 2 2"),
		  &status );

  //Do data writing.  We have to make a transposed copy of the
  // data to do this, which is irritating as hell
  double *tmpdata = new double[ cg.getAxisN(0) ];
  long fpixel[2] = { 1, 1 };
  for ( unsigned int i = 0; i < cg.getAxisN(1); ++i ) {
    for (unsigned int j = 0; j < cg.getAxisN(0); ++j) 
      tmpdata[j] = cg[j][i];
    fpixel[1] = static_cast<long>(i+1);
    fits_write_pix( fp, TDOUBLE, fpixel, cg.getAxisN(0), tmpdata, &status );
  }
  delete[] tmpdata;

  return status;
}

int convert_to_fits3D( fitsfile* fp, const string& infile, bool binaryfile,
		       const string& extname, int& status ) {
  if (status > 0) return status;

  cosgrid3D cg;
  cg.readFile( infile, binaryfile );
  long axissize[3];
  axissize[0] = cg.getAxisN(0); axissize[1] = cg.getAxisN(1); 
  axissize[2] = cg.getAxisN(2);

  fits_create_img( fp, DOUBLE_IMG, 3, axissize, &status );

  char keyval[FLEN_VALUE]; //Value of key if length is in doubt.
  strncpy( keyval, extname.c_str(), FLEN_VALUE );
  fits_write_key( fp, TSTRING, const_cast<char*>("EXTNAME"), keyval,
		  const_cast<char*>("Name of extension"), &status );

  float crpix = 1, zero = 0.0, ftmp;
  strncpy( keyval, cg.getAxisLabel(0).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		  const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key( fp, TFLOAT,const_cast<char*>("CRPIX1"), &crpix, 
		  const_cast<char*>("Ref pix of axis 1"), &status );
  ftmp = cg.getAxisMin(0);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRVAL1"), &ftmp, 
		  const_cast<char*>("val at ref pix"), &status );
  ftmp = cg.getAxisD(0);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRDELT1"), &ftmp, 
		  const_cast<char*>("delta along axis 1"), 
		  &status );
  strncpy( keyval, cg.getAxisLabel(1).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE2"), keyval,
		  const_cast<char*>("Type of Data axis 2"),&status);
  fits_write_key( fp, TDOUBLE,const_cast<char*>("CRPIX2"), &crpix, 
		  const_cast<char*>("Ref pix of axis 2"), &status );
  ftmp = cg.getAxisMin(1);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRVAL2"), &ftmp, 
		  const_cast<char*>("val at ref pix"), &status );
  ftmp = cg.getAxisD(1);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRDELT2"), &ftmp, 
		  const_cast<char*>("delta along axis 2"), 
		  &status );
  strncpy( keyval, cg.getAxisLabel(2).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE3"), keyval,
		  const_cast<char*>("Type of Data axis 3"),&status);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX3"), &crpix, 
		  const_cast<char*>("Ref pix of axis 3"), &status );
  ftmp = cg.getAxisMin(2);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRVAL3"), &ftmp, 
		  const_cast<char*>("val at ref pix"), &status );
  ftmp = cg.getAxisD(2);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CRDELT3"), &ftmp, 
		  const_cast<char*>("delta along axis 3"), 
		  &status );
  ftmp = cg.getAxisD(0);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CD1_1"), &ftmp, 
		  const_cast<char*>("WCS matrix element 1 1"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_2"), &zero, 
		  const_cast<char*>("WCS matrix element 1 2"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_3"), &zero, 
		  const_cast<char*>("WCS matrix element 1 3"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD2_1"), &zero, 
		  const_cast<char*>("WCS matrix element 2 1"),
		  &status );
  ftmp = cg.getAxisD(1);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CD2_2"), &ftmp, 
		  const_cast<char*>("WCS matrix element 2 2"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD2_3"), &zero, 
		  const_cast<char*>("WCS matrix element 2 3"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD3_1"), &zero, 
		  const_cast<char*>("WCS matrix element 3 1"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD3_2"), &zero, 
		  const_cast<char*>("WCS matrix element 3 2"),
		  &status );
  ftmp = cg.getAxisD(2);
  fits_write_key( fp, TDOUBLE, const_cast<char*>("CD3_3"), &ftmp, 
		  const_cast<char*>("WCS matrix element 3 3"),
		  &status );

  //Do data writing.  We have to make a transposed copy of the
  // data to do this, which is irritating as hell. cfitsio handles
  // arrays in FORTRAN order, which is not how I store them internally
  long fpixel[3] = { 1, 1, 1 };
  double *tmpdata = new double[ cg.getAxisN(0) ];
  for (unsigned int i = 0; i < cg.getAxisN(2); ++i) {
    fpixel[2] = static_cast<long>(i+1);
    for (unsigned int j = 0; j < cg.getAxisN(1); ++j) {
      fpixel[1] = static_cast<long>(j+1);
      for (unsigned int k = 0; k < cg.getAxisN(0); ++k) 
	tmpdata[k] = cg[k][j][i];
      fits_write_pix( fp, TDOUBLE, fpixel, cg.getAxisN(0), tmpdata, &status );
    }
  }
  delete[] tmpdata;

  return status;
}

//Add in supernova data
int  AddSNeData( fitsfile* fp, const fitparam& fparam, int& status ) {

  if (status > 0) return status;

  SNeData sne;
  sne.setName("Supernova");
  sne.readData( fparam.datafilename );
  sne.zcmbsort();
  long nsn = static_cast<long>(sne.size());
  
  //Make binary table to hold data
  int tfields = 20;
  const char* ttype[] = { "name", "zcmb", "zhel", "dz", "var_z", "mag", "dmag", 
			  "var_mag", "width", "dwidth", "var_width", "colour",
			  "dcolour","var_colour", "thirdvar", "dthirdvar",
			  "var_thirdvar", "cov_m_w", "cov_m_c",
			  "cov_w_c" };
  const char* tform[] = { "12A", "E", "E", "E", "E", "E", "E", "E", "E", "E", 
			  "E", "E", "E", "E", "E", "E", "E", "E", "E", "E" };

  fits_create_tbl( fp, BINARY_TBL, nsn, tfields, const_cast<char**>(ttype), 
		   const_cast<char**>(tform),
		   NULL, const_cast<char*>("SNEDATA"), &status);

  //Build an array of SN names
  char** snnames;
  snnames = new char*[nsn];
  for (long i = 0; i < nsn; ++i) {
    snnames[i] = new char[FLEN_VALUE];
    strncpy( snnames[i], sne[i].name.c_str(), FLEN_VALUE-1 );
  }    
  fits_write_col( fp, TSTRING, 1, 1, 1, nsn, snnames, &status );
  for (long i = 0; i < nsn; ++i) delete[] snnames[i];
  delete[] snnames;

  //Now lots of other things.  Sadly, we have to do lots of copying
  float *vals;
  vals = new float[ nsn ];

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].zcmb;
  fits_write_col( fp, TFLOAT, 2, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].zhel;
  fits_write_col( fp, TFLOAT, 3, 1, 1, nsn, vals, &status );
  
  for (long i = 0; i < nsn; ++i) vals[i] = sqrt( sne[i].var_z );
  fits_write_col( fp, TFLOAT, 4, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].var_z;
  fits_write_col( fp, TFLOAT, 5, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].mag;
  fits_write_col( fp, TFLOAT, 6, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sqrt( sne[i].var_mag );
  fits_write_col( fp, TFLOAT, 7, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].var_mag;
  fits_write_col( fp, TFLOAT, 8, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].widthpar;
  fits_write_col( fp, TFLOAT, 9, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sqrt( sne[i].var_widthpar );
  fits_write_col( fp, TFLOAT, 10, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].var_widthpar;
  fits_write_col( fp, TFLOAT, 11, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].colourpar;
  fits_write_col( fp, TFLOAT, 12, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sqrt( sne[i].var_colourpar );
  fits_write_col( fp, TFLOAT, 13, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].var_colourpar;
  fits_write_col( fp, TFLOAT, 14, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].thirdpar;
  fits_write_col( fp, TFLOAT, 15, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sqrt( sne[i].var_thirdpar );
  fits_write_col( fp, TFLOAT, 16, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].var_thirdpar;
  fits_write_col( fp, TFLOAT, 17, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].cov_mag_widthpar;
  fits_write_col( fp, TFLOAT, 18, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].cov_mag_colourpar;
  fits_write_col( fp, TFLOAT, 19, 1, 1, nsn, vals, &status );

  for (long i = 0; i < nsn; ++i) vals[i] = sne[i].cov_widthpar_colourpar;
  fits_write_col( fp, TFLOAT, 20, 1, 1, nsn, vals, &status );

  delete[] vals;

  return status;
}

//Reads in and parses the extended output file
int AddExtendedData( fitsfile* fp, const fitparam& fparam, int& status ) {

  if (status > 0) return status;

  if ( !fparam.extendedout )
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "No extended output file set",1);
  
  ifstream fl;
  fl.open( fparam.extendedoutfilename.c_str() );
  if (!fl) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << fparam.extendedoutfilename << 
      " for reading";
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  errstrng.str(),1);
  }

  string line;
  vector<string> words;
  stringstream str;
  float w0, wa, om, ode, sm1, sm2, alpha, beta;

  //The first line is the fit type
  getline(fl,line);

  int first = line.find_first_of(":");
  int last=line.find_last_not_of(" ");
  line = line.substr( first, last-first+1 );
  first = line.find_first_of(":");
  string fittype = line.substr( first, line.size()-first );
  
  //Next line is ignored
  getline(fl,line);
  
  //Then the best fit parameters w0, wa, ,om,ode,sm1, sm2, alpha, beta
  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 2)
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with w0: line",2);
  str.str( words[1] );
  str.clear();
  str >> w0;

  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 2)
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with wa: line",2);
  str.str( words[1] );
  str.clear();
  str >> wa;

  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 2)
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with om: line",4);
  str.str( words[1] );
  str.clear();
  str >> om;

  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 2)
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with ode: line",5);
  str.str( words[1] );
  str.clear();
  str >> ode;
  
  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 2)
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with sm1: line",6);
  str.str( words[1] );
  str.clear();
  str >> sm1;

  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 2)
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with sm2: line",6);
  str.str( words[1] );
  str.clear();
  str >> sm2;

  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 2)
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with alpha: line",7);
  str.str( words[1] );
  str.clear();
  str >> alpha;

  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 2)
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with beta: line",8);
  str.str( words[1] );
  str.clear();
  str >> beta;

  //Determine how many intrinsic disp datasets
  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() == 6) {
    str.str(words[5]);
    str.clear();
    unsigned int ndisp;
    str >> ndisp;
    if (ndisp > 0) for (unsigned int i = 0; i < ndisp; ++i)
		     getline(fl,line);
  }

  //Skip the next two (intrinsicdisp, pecz)
  getline(fl,line);

  //Then grab the number of SN
  //The line is Number of objects: nsn
  unsigned int nsn;
  getline(fl,line);
  utility::stringwords(line,words);
  if (words.size() != 4) {
    throw CosFitterExcept("convert_to_fits","AddExtendedData",
			  "Trouble with nsn line",2);
  }
  str.str( words[3] );
  str.clear();
  str >> nsn;


  //Skip the warning and header line
  getline(fl,line);
  getline(fl,line);

  //Allocate the arrays to hold the results
  char **snnames;
  snnames = new char*[ nsn ];
  for (unsigned int i = 0; i < nsn; ++i ) snnames[i] = new char[FLEN_VALUE];

  const int ntags = 16;
  float *flttags[ntags-1]; //-1 because the first tag is a string
  for (unsigned int i = 0; i < ntags-1; ++i) flttags[i] = new float[nsn];

  //Start reading
  for (unsigned int i = 0; i < nsn; ++i) {
    getline(fl,line);
    utility::stringwords(line,words);

    //Get the snname
    str.str(words[0]); str.clear();
    strncpy( snnames[i], str.str().c_str(), FLEN_VALUE-1 );

    //Now get floating tags
    for (unsigned int j = 0; j < ntags-1; ++j) {
      str.str( words[j+1] ); str.clear(); str >> flttags[j][i];
    }
  }

  fl.close();

  //Now start buidling binary table
  const char* ttype[ntags] = { "Name", "zcmb", "zhel", "raw_mag", "err_raw",
			       "D_L+sm", "corr", "fitmag", "err_fitmag", "diff",
			       "chisq", "sigma", "width", "dwidth", "colour",
			       "dcolour" };
  const char* tform[ntags] = { "12A", "E", "E", "E", "E", "E", "E", "E", "E",
			       "E", "E", "E", "E", "E", "E", "E" };

  fits_create_tbl( fp, BINARY_TBL, nsn, ntags, const_cast<char**>(ttype), 
		   const_cast<char**>(tform), NULL,
		   const_cast<char*>("BESTFIT"), &status );

  //Add the best fit values as header keywords
  fits_write_key( fp, TFLOAT, const_cast<char*>("W0BEST"), &w0, 
		  const_cast<char*>("w_0 best fit"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("WABEST"), &wa, 
		  const_cast<char*>("w_a best fit"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("OMBEST"), &om, 
		  const_cast<char*>("Omega_m best fit"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("ODEBEST"), &ode, 
		  const_cast<char*>("Omega_de best fit"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("SM1BEST"), &sm1, 
		  const_cast<char*>("scriptm best fit"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("SM2BEST"), &sm2, 
		  const_cast<char*>("scriptm best fit"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("ALBEST"), &alpha, 
		  const_cast<char*>("alpha best fit"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("BETABEST"), &beta, 
		  const_cast<char*>("beta best fit"),
		  &status );
  fits_write_key( fp, TINT, const_cast<char*>("NSN"), &nsn, 
		  const_cast<char*>("Number of SNe"), &status);

  //And write the columns
  fits_write_col( fp, TSTRING, 1, 1, 1, nsn, snnames, &status );

  for (unsigned int i = 0; i < ntags-1; ++i) 
    fits_write_col( fp, TFLOAT, i+2, 1, 1, nsn, flttags[i], &status );

  //Clean up
  for (unsigned int i = 0; i < nsn; ++i ) delete[] snnames[i];
  delete[] snnames;
  for (unsigned int i = 0; i < ntags-1; ++i ) delete[] flttags[i];

  return status;
}

int main( int argc, char **argv ) {

  int c;
  string inputfile, outputfile;
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {0,0,0,0}
  };

  while ( ( c = getopt_long( argc, argv, "h",long_options,
			     &option_index ) ) != -1 )
    switch (c) {
    case 'h' :
      cerr << "NAME" << endl;
      cerr << "\tconvert_to_fits -- Converts output from fit to FITS file"
	   << endl;
      cerr << "SYNOPSIS" << endl;
      cerr << "\tconvert_to_fits INFILE OUTFILE\n";
      cerr << "DESCRIPTION" << endl;
      cerr << "\tConverts the results of a cosmology fit to a FITS file"
	   << " OUTFILE" << endl;
      cerr << "\tOnly 1,2, and 3D surfaces currently supported." << endl;
      cerr << "\tINFILE is the parameter file of the fit." << endl;
      cerr << "\tAs usual with cfitsio, will fail if file already exists," 
	   << endl;
      cerr << "\tunless a ! is prepended to the filename.  Also supports"
	   << endl;
      cerr << "\tcompression, etc. through the normal cfitsio extended" 
	   << endl;
      cerr << "\tfilename syntax for the output file." << endl;
      return 0;
      break;
    }

  if (optind > argc - 2) {
    cerr << "ERROR in convert_to_fits -- not enough arguments provided\n";
    return 1;
  }

  inputfile = string( argv[optind] );
  outputfile = string( argv[optind+1] );

  //Process the parameter file
  fitparam fparam;
  try {
    fparam.readFile( inputfile, true );
  } catch (const CosFitterExcept& ex) {
    cerr << "ERROR reading parameter file in convert_to_fits\n";
    cerr << ex << endl;
    return 2;
  }

  //Open the file
  int status = 0;
  fitsfile *fp;
  fits_create_file( &fp, outputfile.c_str(), &status );
  
  //Write the prob surface
  unsigned int naxes;
  bool binaryfile = fparam.binaryout;
  try {
    naxes = NAxesInFile( fparam.outputfilename, binaryfile );
  } catch (const CosFitterExcept& ex) {
    cerr << "ERROR in convert_to_fits\n";
    cerr << ex << endl;
    return 4;
  }

  //Now read it in
  try {
    switch (naxes) {
    case 1 :
      convert_to_fits1D( fp, fparam.outputfilename, binaryfile, 
			 "PROB", status );
      break;
    case 2 :
      convert_to_fits2D( fp, fparam.outputfilename, binaryfile, 
			 "PROB", status );
      break;
    case 3 :
      convert_to_fits3D( fp, fparam.outputfilename, binaryfile, 
			 "PROB", status );
      break;
    default :
      cerr << "Number of axes in " << fparam.outputfilename 
	   << " (" << naxes << ") not"
	   << " currently supported" << endl;
      return 8;
    }

    AddFitParams( fp, fparam, status );
    
    //Do SNe data table
    AddSNeData( fp, fparam, status );
    AddFitParams( fp, fparam, status );
    
    if (fparam.extendedout) {
      AddExtendedData( fp, fparam, status );
      AddFitParams( fp, fparam, status );
    }
    
    //Do marginalized alpha/beta stuff
    if (fparam.albetaout) {
      //We have the possibility of this being 1D now
      it = fparam.params.find( param_tags::alpha );
      if ( it == fparam.params.end() ) {
	//No alpha, must have beta
	convert_to_fits1D( fp, fparam.albetaoutputfilename,
			   binaryfile, "BETAPROB", status );
	AddFitParams( fp, fparam, status );
      } else {
	//Have alpha
	it = fparam.params.find( param_tags::beta );
	if ( it == fparam.params.end() ) {
	  //But no beta
	  convert_to_fits1D( fp, fparam.albetaoutputfilename,
			     binaryfile, "ALPHPROB", status );
	  AddFitParams( fp, fparam, status );
	} else {
	  //Have both
	  convert_to_fits2D( fp, fparam.albetaoutputfilename, 
			     binaryfile, "ALPHBETA", status );
	  AddFitParams( fp, fparam, status );
	}
      }
    }
  } catch (const CosFitterExcept& ex) {
    cerr << "ERROR in convert_to_fits\n";
    cerr << ex << endl;
    return 8;
  }

  fits_close_file( fp, &status );

  if (status) {
    fits_report_error(stderr, status); /* print any error message */
    return 1024;
  }

  return 0;
}
