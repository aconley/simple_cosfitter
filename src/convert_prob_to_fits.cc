//Routine for converting probability grids into FITS files

#include <fitsio.h>
#include <getopt.h>
#include <string>
#include <string.h>
#include <iostream>

#include "cosgrids.h"
#include "cosfitterexcept.h"

using namespace std;

//1D cosgrid case
int convert_to_fits1D( fitsfile* fp, const string& inputfile, 
		       bool binaryfile, bool marg, int& status ) {
  if (status > 0) return status;  //No need to do anything, since FITS bad

  cosgrid1D cg;
  cg.readFile( inputfile, binaryfile );
  long axissize[1];
  axissize[0] = static_cast<long>(cg.getAxisN());

  fits_create_img( fp, DOUBLE_IMG, 1, axissize, &status );

  //Add WCS information to header
  float crpix = 1, tmpval;
  char keyval[FLEN_VALUE]; //Value of key if length is in doubt.
  strncpy( keyval, cg.getAxisLabel().c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		  const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		  const_cast<char*>("Ref pix of axis 1"), &status );
  tmpval = cg.getAxisMin();
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		  const_cast<char*>("val at ref pix"), &status );
  tmpval = cg.getAxisD();
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &tmpval,
		  const_cast<char*>("delta along axis 1"), &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &tmpval, 
		  const_cast<char*>("WCS matrix element 1 1"),&status );
  
  //Do data writing
  long fpixel[1] = { 1 };
  fits_write_pix( fp, TDOUBLE, fpixel, cg.getAxisN(), cg.getData(), 
		  &status );

  return status;

}

int convert_to_fits2D( fitsfile *fp, const string& inputfile,
		       bool binaryfile, bool marg, int& status ) {
  if (status > 0) return status;  //No need to do anything, since FITS bad

  cosgrid2D cg;
  cg.readFile( inputfile, binaryfile );
  long axissize[2];
  axissize[0] = static_cast<long>(cg.getAxisN(0)); 
  axissize[1] = static_cast<long>(cg.getAxisN(1));

  fits_create_img( fp, DOUBLE_IMG, 2, axissize, &status );

  //Add WCS information to header
  float crpix = 1, zero = 0.0, tmpval;
  char keyval[FLEN_VALUE]; //Value of key if length is in doubt.
  strncpy( keyval, cg.getAxisLabel(0).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		  const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		  const_cast<char*>("Ref pix of axis 1"), &status );
  tmpval = cg.getAxisMin(0);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		  const_cast<char*>("val at ref pix"), &status );
  tmpval = cg.getAxisD(0);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &tmpval,
		  const_cast<char*>("delta along axis 1"), &status );
  strncpy( keyval, cg.getAxisLabel(1).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE2"), keyval,
		  const_cast<char*>("Type of Data axis 2"),&status);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX2"), &crpix, 
		  const_cast<char*>("Ref pix of axis 2"), &status );
  tmpval = cg.getAxisMin(1);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL2"), &tmpval, 
		  const_cast<char*>("val at ref pix"), &status );
  tmpval = cg.getAxisD(1);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT2"), &tmpval, 
		  const_cast<char*>("delta along axis 2"), &status );
  tmpval = cg.getAxisD(0);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &tmpval,
		  const_cast<char*>("WCS matrix element 1 1"),&status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_2"), &zero, 
		  const_cast<char*>("WCS matrix element 1 2"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD2_1"), &zero, 
		  const_cast<char*>("WCS matrix element 2 1"),
		  &status );
  tmpval = cg.getAxisD(1);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD2_2"), &tmpval, 
		  const_cast<char*>("WCS matrix element 2 2"),&status );

  //Do data writing.  We have to make a transposed copy of the
  // data to do this, which is irritating as hell
  double *tmpdata = new double[ cg.getAxisN(0) ];
  long fpixel[2] = { 1, 1 };
  for ( unsigned int i = 0; i < cg.getAxisN(1); ++i ) {
    for (unsigned int j = 0; j < cg.getAxisN(0); ++j) tmpdata[j] = cg[j][i];
    fpixel[1] = static_cast<long>(i+1);
    fits_write_pix( fp, TDOUBLE, fpixel, cg.getAxisN(0), tmpdata, &status );
  }
  delete[] tmpdata;

  //Do marginalized 1D surfaces
  if (marg) {
    long iaxissize[1];
    long ifpixel[1] = { 1 };
    cosgrid1D cg_col;
    cg_col = cg.collapseAlongAxis(1);
    iaxissize[0] = static_cast<long>(cg_col.getAxisN() );
    fits_create_img( fp, DOUBLE_IMG, 1, iaxissize, &status );
    strncpy( keyval, cg.getAxisLabel(0).c_str(), FLEN_VALUE - 19 );
    fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		    const_cast<char*>("Type of Data axis 1"),
		    &status);
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		    const_cast<char*>("Ref pix of axis 1"), 
		    &status );
    tmpval = cg_col.getAxisMin();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &tmpval,
		    const_cast<char*>("val at ref pix"), &status );
    tmpval = cg_col.getAxisD();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &tmpval, 
		    const_cast<char*>("delta along axis 1"), &status );
    fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &tmpval,
		    const_cast<char*>("WCS matrix element 1 1"),&status );
    strncpy( keyval, "AXIS1MARG", FLEN_VALUE );
    fits_write_key( fp, TSTRING, const_cast<char*>("EXTNAME"), keyval,
		    const_cast<char*>("Name of extension"),
		    &status );
    fits_write_pix( fp, TDOUBLE, ifpixel, cg_col.getAxisN(), 
		    cg_col.getData(), &status );
    
    cg_col = cg.collapseAlongAxis(0);
    iaxissize[0] = static_cast<long>(cg_col.getAxisN());
    fits_create_img( fp, DOUBLE_IMG, 1, iaxissize, &status );
    strncpy( keyval, cg.getAxisLabel(1).c_str(), FLEN_VALUE - 19 );
    fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		    const_cast<char*>("Type of Data axis 1"),
		    &status);
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		    const_cast<char*>("Ref pix of axis 1"), 
		    &status );
    tmpval = cg_col.getAxisMin();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		    const_cast<char*>("val at ref pix"), &status );
    tmpval = cg_col.getAxisD();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &tmpval,
		    const_cast<char*>("delta along axis 1"), &status );
    fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &tmpval,
		    const_cast<char*>("WCS matrix element 1 1"), &status );
    strncpy( keyval, "AXIS2MARG", FLEN_VALUE );
    fits_write_key( fp, TSTRING, const_cast<char*>("EXTNAME"), keyval,
		    const_cast<char*>("Name of extension"),
		    &status );
    fits_write_pix( fp, TDOUBLE, ifpixel, cg_col.getAxisN(), 
		    cg_col.getData(), &status );
  }    

  return status;
}

int convert_to_fits3D( fitsfile* fp, const string& inputfile,
		       bool binaryfile, bool marg, int& status ) {
  if (status > 0) return status;  //No need to do anything, since FITS bad

  cosgrid3D cg;
  cg.readFile( inputfile, binaryfile );
  long axissize[3];
  axissize[0] = static_cast<long>(cg.getAxisN(0)); 
  axissize[1] = static_cast<long>(cg.getAxisN(1)); 
  axissize[2] = static_cast<long>(cg.getAxisN(2));

  fits_create_img( fp, DOUBLE_IMG, 3, axissize, &status );

  float crpix = 1, zero = 0.0, tmpval;
  char keyval[FLEN_VALUE]; //Value of key if length is in doubt.
  strncpy( keyval, cg.getAxisLabel(0).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		  const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key( fp, TFLOAT,const_cast<char*>("CRPIX1"), &crpix, 
		  const_cast<char*>("Ref pix of axis 1"), &status );
  tmpval = cg.getAxisMin(0);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		  const_cast<char*>("val at ref pix"), &status );
  tmpval = cg.getAxisD(0);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &tmpval, 
		  const_cast<char*>("delta along axis 1"), &status );
  strncpy( keyval, cg.getAxisLabel(1).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE2"), keyval,
		  const_cast<char*>("Type of Data axis 2"),&status);
  fits_write_key( fp, TFLOAT,const_cast<char*>("CRPIX2"), &crpix, 
		  const_cast<char*>("Ref pix of axis 2"), &status );
  tmpval = cg.getAxisMin(1);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL2"), &tmpval, 
		  const_cast<char*>("val at ref pix"), &status );
  tmpval = cg.getAxisD(1);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT2"), &tmpval, 
		  const_cast<char*>("delta along axis 2"), &status );
  strncpy( keyval, cg.getAxisLabel(2).c_str(), FLEN_VALUE - 19 );
  fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE3"), keyval,
		  const_cast<char*>("Type of Data axis 3"),&status);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX3"), &crpix, 
		  const_cast<char*>("Ref pix of axis 3"), &status );
  tmpval = cg.getAxisMin(2);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL3"), &tmpval,
		  const_cast<char*>("val at ref pix"), &status );
  tmpval = cg.getAxisD(2), 
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT3"), &tmpval, 
		    const_cast<char*>("delta along axis 3"), &status );
  tmpval = cg.getAxisD(0), 
    fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &tmpval, 
		    const_cast<char*>("WCS matrix element 1 1"),&status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_2"), &zero, 
		  const_cast<char*>("WCS matrix element 1 2"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_3"), &zero, 
		  const_cast<char*>("WCS matrix element 1 3"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD2_1"), &zero, 
		  const_cast<char*>("WCS matrix element 2 1"),
		  &status );
  tmpval = cg.getAxisD(1);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD2_2"), &tmpval, 
		  const_cast<char*>("WCS matrix element 2 2"), &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD2_3"), &zero, 
		  const_cast<char*>("WCS matrix element 2 3"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD3_1"), &zero, 
		  const_cast<char*>("WCS matrix element 3 1"),
		  &status );
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD3_2"), &zero, 
		  const_cast<char*>("WCS matrix element 3 2"),
		  &status );
  tmpval = cg.getAxisD(2);
  fits_write_key( fp, TFLOAT, const_cast<char*>("CD3_3"), &tmpval, 
		  const_cast<char*>("WCS matrix element 3 3"),&status );

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

  //Do marginalized 1D surfaces
  if (marg) {
    long iaxissize[1];
    long ifpixel[1] = { 1 };
    cosgrid1D cg_col;
    cg_col = cg.collapseAlongTwoAxes(1,2);
    iaxissize[0] = static_cast<long>(cg_col.getAxisN());
    fits_create_img( fp, DOUBLE_IMG, 1, iaxissize, &status );
    strncpy( keyval, cg.getAxisLabel(0).c_str(), FLEN_VALUE - 19 );
    fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		    const_cast<char*>("Type of Data axis 1"),
		    &status);
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		    const_cast<char*>("Ref pix of axis 1"), 
		    &status );
    tmpval = cg_col.getAxisMin();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		    const_cast<char*>("val at ref pix"), &status );
    tmpval = cg_col.getAxisD();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &tmpval, 
		    const_cast<char*>("delta along axis 1"), &status );
    fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &tmpval,
		    const_cast<char*>("WCS matrix element 1 1"), &status );
    strncpy( keyval, "AXIS1MARG", FLEN_VALUE );
    fits_write_key( fp, TSTRING, const_cast<char*>("EXTNAME"), keyval,
		    const_cast<char*>("Name of extension"),
		    &status );
    fits_write_pix( fp, TDOUBLE, ifpixel, cg_col.getAxisN(), 
		    cg_col.getData(), &status );
    
    cg_col = cg.collapseAlongTwoAxes(0,2);
    iaxissize[0] = static_cast<long>(cg_col.getAxisN());
    fits_create_img( fp, DOUBLE_IMG, 1, iaxissize, &status );
    strncpy( keyval, cg.getAxisLabel(1).c_str(), FLEN_VALUE - 19 );
    fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		    const_cast<char*>("Type of Data axis 1"),
		    &status);
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		    const_cast<char*>("Ref pix of axis 1"), 
		    &status );
    tmpval = cg_col.getAxisMin();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		    const_cast<char*>("val at ref pix"), &status );
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &tmpval,
		    const_cast<char*>("delta along axis 1"), &status );
    tmpval = cg_col.getAxisD();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &tmpval, 
		    const_cast<char*>("WCS matrix element 1 1"),&status );
    strncpy( keyval, "AXIS2MARG", FLEN_VALUE );
    fits_write_key( fp, TSTRING, const_cast<char*>("EXTNAME"), keyval,
		    const_cast<char*>("Name of extension"),
		    &status );
    fits_write_pix( fp, TDOUBLE, ifpixel, cg_col.getAxisN(), 
		    cg_col.getData(), &status );
    
    cg_col = cg.collapseAlongTwoAxes(0,1);
    iaxissize[0] = static_cast<long>(cg_col.getAxisN());
    fits_create_img( fp, DOUBLE_IMG, 1, iaxissize, &status );
    strncpy( keyval, cg.getAxisLabel(2).c_str(), FLEN_VALUE - 19 );
    fits_write_key( fp, TSTRING, const_cast<char*>("CTYPE1"), keyval,
		    const_cast<char*>("Type of Data axis 1"),
		    &status);
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		    const_cast<char*>("Ref pix of axis 1"), 
		    &status );
    tmpval = cg_col.getAxisMin();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRVAL1"), &tmpval, 
		    const_cast<char*>("val at ref pix"), &status );
    tmpval = cg_col.getAxisD();
    fits_write_key( fp, TFLOAT, const_cast<char*>("CRDELT1"), &tmpval, 
		    const_cast<char*>("delta along axis 1"), &status );
    fits_write_key( fp, TFLOAT, const_cast<char*>("CD1_1"), &tmpval,
		    const_cast<char*>("WCS matrix element 1 1"), &status );
    strncpy( keyval, "AXIS3MARG", FLEN_VALUE );
    fits_write_key( fp, TSTRING, const_cast<char*>("EXTNAME"), keyval,
		    const_cast<char*>("Name of extension"),
		    &status );
    fits_write_pix( fp, TDOUBLE, ifpixel, cg_col.getAxisN(), 
		    cg_col.getData(), &status );
  }

  return status;
}

int main( int argc, char **argv ) {

  int c;
  bool binaryfile, marg;
  string inputfile, outputfile;

  //Defaults
  binaryfile = false;
  marg = true;

  int option_index = 0;
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"binaryfile", no_argument, 0, 'b' },
    {"nomarg",no_argument,0,'n'},
    {0,0,0,0}
  };

  while ( ( c = getopt_long( argc, argv, "hbn",long_options,
			     &option_index ) ) != -1 )
    switch (c) {
    case 'h' :
      cerr << "NAME" << endl;
      cerr << "\tconvert_prob_to_fits -- Converts cosgrid to a FITS file" 
	   << endl;
      cerr << "SYNOPSIS" << endl;
      cerr << "\tconvert_prob_to_fits [ -b | --binaryfile ] [ -n | --nomarg ]"
	   << endl;
      cerr << "\t INFILE OUTFILE\n";
      cerr << "DESCRIPTION" << endl;
      cerr << "\tConverts a cosgrid style file INFILE to a FITS " <<
	"file OUTFILE." << endl;
      cerr << "\tOnly 1,2, and 3D surfaces currently supported." << endl;
      cerr << "\tAs usual with cfitsio, will fail if file already exists," 
	   << endl;
      cerr << "\tunless a ! is prepended to the filename.  Also supports"
	   << endl;
      cerr << "\tcompression, etc. through the normal cfitsio extended" 
	   << endl;
      cerr << "\tfilename syntax for the output file." << endl;
      cerr << "OPTIONS" << endl;
      cerr << "\t-b, --binaryfile Input file is in cosgrid binary format"
	   << endl;
      cerr << "\t-n, --nomarg Don't add marginalized distributions to output"
	   << " file" << endl;
      return 0;
      break;
    case 'b' :
      binaryfile = true;
      break;
    case 'n' :
      marg = false;
      break;
    }

  if (optind > argc - 2) {
    cerr << "ERROR in convert_prob_to_fits -- not enough arguments provided\n";
    return 1;
  }

  inputfile = string( argv[optind] );
  outputfile = string( argv[optind+1] );

  unsigned int naxes;
  try {
    naxes = NAxesInFile( inputfile, binaryfile );
  } catch (const CosFitterExcept& ex) {
    cerr << "ERROR in convert_prob_to_fits\n";
    cerr << ex << endl;
    return 2;
  }

  //Make the fits file
  int status = 0;
  fitsfile *fp;
  fits_create_file( &fp, outputfile.c_str(), &status );

  //Now read it in
  try {
    switch (naxes) {
    case 1 :
      convert_to_fits1D( fp, inputfile, binaryfile, marg, status );
      break;
    case 2 :
      convert_to_fits2D( fp, inputfile, binaryfile, marg, status );
      break;
    case 3 :
      convert_to_fits3D( fp, inputfile, binaryfile, marg, status );
      break;
    default :
      std::cerr << "Number of axes in " << inputfile << " (" 
		<< naxes << ") not"
		<< " currently supported" << endl;
      return 4;
    }
  } catch (const CosFitterExcept& ex) {
    cerr << "ERROR in convert_prob_to_fits\n";
    cerr << ex << endl;
    return 8;
  }

  fits_close_file( fp, &status );

  //Make sure everything is fine up to this point
  if (status) {
    fits_report_error(stderr, status); /* print any error message */
    return 1024;
  }

  return 0;
}
