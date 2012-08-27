//cosgrids.cc
//Implementation of cosmology grid classes

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#include <cosgrids.h>
#include <cosfitterexcept.h>
#include <utility.h>

using namespace std;

//First, a utility for reading the axisnames from
// binary files, as the writeBinaryFile stuff specifies
int readBinaryAxisName( FILE *fp, std::string& axisname ) {
  unsigned int axisnamelen;
  char *axisnamecstr;
  fread( &axisnamelen, sizeof(unsigned int), 1, fp );
  if ( axisnamelen < 1 ) return 1;
  axisnamecstr = new char[ axisnamelen + 1 ]; //Need room for \0

  fgets( axisnamecstr, axisnamelen+1, fp );//Not sure why +1 is needed, but it is
  
  axisname = string( axisnamecstr );
  
  delete[] axisnamecstr;
  return 0;
}

/*!
  \param[in] filename Name of file to check
  \param[in] binaryfile Is file in cosgrid binary format?  (def: false)
  \returns Number of axes
*/
unsigned int NAxesInFile( const std::string& filename, bool binaryfile ) {

  unsigned int naxes;
  FILE *fp;
  if (binaryfile) {
    fp = fopen( filename.c_str(), "rb" );
    if (!fp) {
      stringstream errstrng("");
      errstrng << "Couldn't open file " << filename << " for reading";
      throw CosFitterExcept("cosgrids","NAxesInFile",errstrng.str(),1);
    }
    fread( &naxes, sizeof(unsigned int), 1, fp );
    fclose(fp);
  } else {
    ifstream fl( filename.c_str() );
    if (!fl) {
      stringstream errstrng("");
      errstrng << "Couldn't open file " << filename << " for reading";
      throw CosFitterExcept("cosgrids","NAxesInFile",errstrng.str(),1);
    }
    string line;
    vector<string> words;
    getline(fl,line);
    utility::stringwords(line,words);
    naxes = atoi( words[1].c_str() ); 
    fl.close();
  }

  return naxes;

}

//////////////////////////////////////////
//              cosgrid1D               //
//////////////////////////////////////////
cosgrid1D::cosgrid1D() {
  data = 0;
  n0 = 0;
  min0 = max0 = d0 = 0.0;
}

cosgrid1D::cosgrid1D(const cosgrid1D& input) {
  n0 = input.n0;
  min0 = input.min0;
  max0 = input.max0;
  axisnames[0] = input.axisnames[0];
  spec[0] = input.spec[0];
  d0 = input.d0;
  if (input.data != 0) {
    data = new double[ n0 ];
    for (unsigned int i = 0; i < n0; ++i) data[i] = input.data[i];
  }
}

/*!
 \param[in] n0 Number of data entries along axis 0
 \param[in] min0 Minimum axis value along axis 0
 \param[in] max0 Maximum axis value along axis 0
 \param[in] d0 Step size for axis along axis 0
 \param[in] label0 Name of axis 0 parameter.  Should not have spaces in it.
 \param[in] spec0 Spec for axis 0 parameter
*/
cosgrid1D::cosgrid1D(unsigned int n0, double min0, double max0, double d0,
		     std::string label0, param_tags::paramspec spec0) :
  n0( n0 ), min0( min0 ), max0( max0 ), d0( d0 ) {
  axisnames[0] = label0;
  spec[0] = spec0;
 
  data = new double[n0];
}

/*!
  \param[in] param Information about parameter describing axis 0
*/
cosgrid1D::cosgrid1D(const param_struct& param) :
  n0( param.n ), min0( param.min ), max0( param.max ), 
  d0( param.dval ) {
  axisnames[0] = param.name;
  spec[0] = param.param_spec;
 
  data = new double[n0];
}

cosgrid1D::~cosgrid1D() {
  if (data)
    delete[] data;
}

void cosgrid1D::free() {
  n0 = 0;
  min0 = max0 = d0 = 0.0;
  if (data != 0) delete[] data;
  data = 0;
}

unsigned int cosgrid1D::getAxisIndex(const std::string& axislabel) const {
  //Returns index of axis with name axislabel
  if (axisnames[0] == axislabel) return 0;
  stringstream errstrng("");
  errstrng << "No axis found with label " << axislabel;
  throw CosFitterExcept("cosgrid1D","getAxisIndex",errstrng.str(),1);
}


/*!
  Find and return maximum value in array.
  \param[out] i0 The position of the maximum in the array
  \return The maximum value of the data array
*/

double cosgrid1D::getMaximum(unsigned int& i0) const {
  double maxval = data[0];
  i0 = 0;

  for (unsigned int i = 1; i < n0; ++i) 
    if (data[i] > maxval) {
      maxval = data[i];
      i0 = i;
    }

  return maxval;
}

/*!
  \return The sum of all of the elements in the data array
 */
double cosgrid1D::getTotal() const {
  double val = data[0];
  for (unsigned int i = 1; i < n0; ++i) 
    val += data[i];
  return val;
}  


/*!
  Normalizes the total sum of all data entries to one.
  \return The sum of all of the elements in the data array prior to 
           normalization.
 */
double cosgrid1D::normalize() {
  double val = data[0];
  for (unsigned int i = 1; i < n0; ++i) 
    val += data[i];
  if (val == 0.0) {
    string errstrng = "Total probability equal to 0 -- can't normalize";
    throw CosFitterExcept("cosgrid1D","normalize",errstrng,1);
  }
  for (unsigned int i = 0; i < n0; ++i)
    data[i] /= val;

  return val;
}  

/*!
  Serializes the contents of the data array and the axis labels to a file,
  either in binary or human readable format.
  \param[in] filename the file to write the data to
  \param[in] binary Write as binary instead of human readable
    text (def: false).
*/
int cosgrid1D::writeFile(const std::string& filename, bool binary ) const {
  if (binary) return writeBinaryFile(filename);
  return writeTextFile(filename);
}

/*!
  Reads in a data file as written by writeFile.  The user must specify
  if it is binary or text -- this code is not smart enough to figure it
  out for itself.
  \param[in] filename the file to write the data to
  \param[in] binary Read as binary instead of human readable
    text (def: false).
*/
int cosgrid1D::readFile(const std::string& filename, bool binary ) {
  if (binary) return readBinaryFile(filename);
  return readTextFile(filename);
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a human readable format.  The format of the output file is defined as:
  - The first line is axisname0, n0, min0, max0, d0
  - This is followed by n0 lines containing the data values
  \param[in] filename The file to write the data to
*/
int cosgrid1D::writeTextFile(const std::string& filename) const {
  //Outputs the probability array to a file.
  //The output format is defined as:
  //The first line is Naxes: 1
  //The second line is axisname0, n0, min0, max0, d0
  //This is followed by n0 lines containing the data
  //It is a really bad idea to let any of the axisnames have
  // spaces in them

  FILE *fp;
  fp = fopen(filename.c_str(),"w");
  if (!fp) {
    stringstream errstrng("");
    errstrng <<  "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid1D","writeTextFile",errstrng.str(),1);
  }
  
  fprintf(fp,"Naxes: 1\n");
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[0].c_str(),
	  n0,min0,max0,d0);
  
  for (unsigned int i = 0; i < n0; ++i)
    fprintf(fp,"%-15.8e\n",data[i]);

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid1D::writeTextFile
 \param[in] filename The file to read from
 */
int cosgrid1D::readTextFile(const std::string& filename) {

  if (data)
    delete[] data;

  ifstream fl( filename.c_str() );
  if (!fl) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid1D","readFile",errstrng.str(),1);
  }

  string line;
  vector<string> words;
  stringstream str("");
  
  //Read number of axes
  unsigned int naxes;
  getline(fl,line);
  utility::stringwords(line,words);
  naxes = atoi( words[1].c_str() );
  if ( naxes != 1 ) {
    fl.close();
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid1D","readFile",errstrng.str(),2);
  }
      
  //Read first axis line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[0] = words[0];
  str.str(words[1]); str.clear(); str >> n0;
  str.str(words[2]); str.clear(); str >> min0;
  str.str(words[3]); str.clear(); str >> max0;
  str.str(words[4]); str.clear(); str >> d0;

  data = new double[n0];

  //We have to clear fail flags after reads because
  // there are underflow issues with some versions of gcc --
  // that is, it can write values it can't read even if they
  // don't really underflow
  for (unsigned int i = 0; i < n0; ++i) {
    fl >> data[i];
    if (fl.fail()) fl.clear(fl.rdstate() & ~std::ios::failbit);
  }

  fl.close();

  return 0;
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a non human readable, binary, format.  The format of the output file 
  is defined as:
  - an integer specifying the number of axes.  This is always one for this
     function (duh)
  - an integer specifying the number of characters in the zero terminated
     axisname
  - the axis name
  - an integer specifying n0
  - doubles specifying min0, max0, d0
  - then n0 doubles storing the data
  Note that this will not be machine portable because of endianess and
  the sizes of various outputs.  Unlike the version which writes in text,
  here we don't care if the axisnames have spaces in them.
  \param[in] filename The file to write the data to
*/
int cosgrid1D::writeBinaryFile(const std::string& filename) const {

  FILE *fp;
  fp = fopen(filename.c_str(),"wb");
  if (!fp) {
    stringstream errstrng("");
    errstrng <<  "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid1D","writeBinaryFile",errstrng.str(),1);
  }
  
  //Write number of axes
  unsigned int naxes = 1;
  fwrite( &naxes, sizeof(unsigned int), 1, fp );

  //Write axis title
  unsigned int axisnamelength = axisnames[0].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[0].c_str(), fp );
  
  //Then extended axis information
  fwrite( &n0, sizeof(unsigned int), 1, fp );
  fwrite( &min0, sizeof(double), 1, fp );
  fwrite( &max0, sizeof(double), 1, fp );
  fwrite( &d0, sizeof(double), 1, fp );

  //And the data
  fwrite( data, sizeof(double), n0, fp );

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid1D::writeBinaryFile
 \param[in] filename The file to read from
 */
int cosgrid1D::readBinaryFile(const std::string& filename) {

  if (data)
    delete[] data;

  FILE *fp;
  fp = fopen( filename.c_str(), "rb" );
  if (!fp) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid1D","readBinaryFile",errstrng.str(),1);
  }

  //Read number of axes
  unsigned int naxes;
  fread( &naxes, sizeof(unsigned int), 1, fp );
  if ( naxes != 1 ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid1D","readBinaryFile",errstrng.str(),2);
  }
      
  //Read axis title
  int status;
  status = readBinaryAxisName( fp, axisnames[0] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid1D","readBinaryFile",errstrng.str(),4);
  }
  
  fread( &n0, sizeof(unsigned int), 1, fp );
  fread( &min0, sizeof(double), 1, fp );
  fread( &max0, sizeof(double), 1, fp );
  fread( &d0, sizeof(double), 1, fp );

  if (n0 <= 0) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "0 axis length in " << filename << endl;
    throw CosFitterExcept("cosgrid1D","readBinaryFile",errstrng.str(),8);
  }
  data = new double[n0];

  fread( data, sizeof(double), n0, fp );

  fclose(fp);

  return 0;
}

cosgrid1D& cosgrid1D::operator=(double val) {
  for (unsigned int i = 0; i < n0; ++i)
    data[i] = val;
  return *this;
}

cosgrid1D& cosgrid1D::operator=(const cosgrid1D& input) {

  //Self copy protection
  if (this == &input) return *this;

  axisnames[0] = input.axisnames[0];
  spec[0] = input.spec[0];
  if ( data == 0 ) {
    data = new double[input.n0];
  } else if ( input.n0 != n0 ) {
    delete[] data;
    data = new double[input.n0];
  }

  for (unsigned int i = 0; i < input.n0; ++i)
    data[i] = input.data[i];

  n0 = input.n0; min0 = input.min0; max0 = input.max0; d0 = input.d0;

  return *this;
}

cosgrid1D& cosgrid1D::operator+=(const cosgrid1D& input) {

  //Self addition
  if (this == &input) {
    for (unsigned int i = 0; i < n0; ++i) data[i] *= 2.0;
    return *this;
  }

  //First make sure that the array dimensions and spec codes are
  // compatible
  if (spec[0].second != input.spec[0].second) {
    stringstream errstr("");
    errstr << "Types of axis 0 not compatible in addition" << endl;
    errstr << " One has type: " << spec[0].second << " the other: "
	   << input.spec[0].second;
    errstr << " With name: " << axisnames[0] << " the other: " <<
      input.axisnames[0];
    throw CosFitterExcept("cosgrid1D","operator+=",errstr.str(),1);
  }
  
  double tol = (max0 - min0)*0.001;
  if ( (n0 != input.n0) || ( fabs(min0 - input.min0) > tol ) ||
       ( fabs(max0 - input.max0) > tol ) ) {
    throw CosFitterExcept("cosgrid1D","operator+=",
			     "Axis limits not compatible",2);
  }

  //All good
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] += input.data[i];
  }

  return *this;
}


cosgrid1D& cosgrid1D::operator*=(const cosgrid1D& input) {

  //Self multiplication
  if (this == &input) {
    for (unsigned int i = 0; i < n0; ++i) data[i] *= data[i];
    return *this;
  }

  //First make sure that the array dimensions and spec codes are
  // compatible
  if (spec[0].second != input.spec[0].second) {
    stringstream errstr("");
    errstr << "Types of axis 0 not compatible in addition" << endl;
    errstr << " One has type: " << spec[0].second << " the other: "
	   << input.spec[0].second;
    errstr << " With name: " << axisnames[0] << " the other: " <<
      input.axisnames[0];
    throw CosFitterExcept("cosgrid1D","operator*=",errstr.str(),1);
  }
  
  double tol = (max0 - min0)*0.001;
  if ( (n0 != input.n0) || ( fabs(min0 - input.min0) > tol ) ||
       ( fabs(max0 - input.max0) > tol ) ) {
    throw CosFitterExcept("cosgrid1D","operator*=",
			     "Axis limits not compatible",2);
  }

  //All good
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] *= input.data[i];
  }

  return *this;
}

cosgrid1D cosgrid1D::operator*(const cosgrid1D& input) {

  //First make sure that the array dimensions and spec codes are
  // compatible
  if (spec[0].second != input.spec[0].second) {
    stringstream errstr("");
    errstr << "Types of axis 0 not compatible in addition" << endl;
    errstr << " One has type: " << spec[0].second << " the other: "
	   << input.spec[0].second;
    errstr << " With name: " << axisnames[0] << " the other: " <<
      input.axisnames[0];
    throw CosFitterExcept("cosgrid1D","operator*",errstr.str(),1);
  }
  
  double tol = (max0 - min0)*0.001;
  if ( (n0 != input.n0) || ( fabs(min0 - input.min0) > tol ) ||
       ( fabs(max0 - input.max0) > tol ) ) {
    throw CosFitterExcept("cosgrid1D","operator*",
			     "Axis limits not compatible",2);
  }

  //All good
  cosgrid1D outcosgrid( input );
  for (unsigned int i = 0; i < input.n0; ++i) {
    outcosgrid.data[i] *= data[i];
  }

  return outcosgrid;
}


///////////////////////////////////////////////
//                cosgrid2D                  //
///////////////////////////////////////////////
cosgrid2D::cosgrid2D() {
  data = 0;
  n0 = n1 = 0;
  min0 = max0 = d0 = min1 = max1 = d1 = 0.0;
}

cosgrid2D::cosgrid2D(const cosgrid2D& input) {
  n0 = input.n0;
  min0 = input.min0;
  max0 = input.max0;
  d0 = input.d0;
  n1 = input.n1;
  min1 = input.min1;
  max1 = input.max1;
  d1 = input.d1;
  axisnames[0] = input.axisnames[0];
  axisnames[1] = input.axisnames[1];
  spec[0] = input.spec[0];
  spec[1] = input.spec[1];
  if (input.data != 0) {
    data = new double*[ n0 ];
    for (unsigned int i = 0; i < n0; ++i) {
      data[i] = new double[n1];
      for (unsigned int j = 0; j < n1; ++j) data[i][j] = input.data[i][j];
    }
  }
}

/*!
  \param[in] n0 Number of data entries along axis 0
  \param[in] min0 Minimum axis value along axis 0
  \param[in] max0 Maximum axis value along axis 0
  \param[in] d0 Step size for axis along axis 0
  \param[in] label0 Name of axis 0 parameter.  Should not have spaces in it.
  \param[in] spec0 Param spec for parameter 0
  \param[in] n1 Number of data entries along axis 1
  \param[in] min1 Minimum axis value along axis 1
  \param[in] max1 Maximum axis value along axis 1
  \param[in] d1 Step size for axis along axis 1
  \param[in] label1 Name of axis 1 parameter.  Should not have spaces in it.
  \param[in] spec1 Param spec for parameter 1
*/

cosgrid2D::cosgrid2D(unsigned int n0, double min0, double max0, double d0,
		     std::string label0, param_tags::paramspec spec0,
		     unsigned int n1, double min1, double max1, double d1, 
		     std::string label1, param_tags::paramspec spec1 ) :
  n0(n0), min0(min0), max0(max0), d0(d0), n1(n1), min1(min1), 
  max1(max1), d1(d1) {

  axisnames[0] = label0;
  axisnames[1] = label1;
  spec[0] = spec0;
  spec[1] = spec1;

  data = new double*[n0];
  for (unsigned int i = 0; i < n0; ++i) data[i] = new double[n1];
}

/*!
  \param[in] param0 Information about parameter describing axis 0
  \param[in] param1 Information about parameter describing axis 1
*/

cosgrid2D::cosgrid2D(const param_struct& param0, const param_struct& param1) :
  n0( param0.n ), min0( param0.min ), max0( param0.max ), d0( param0.dval ), 
  n1( param1.n ), min1( param1.min ), max1( param1.max ), d1( param1.dval ) {
  
  axisnames[0] = param0.name;
  axisnames[1] = param1.name;
  spec[0] = param0.param_spec;
  spec[1] = param1.param_spec;

  data = new double*[n0];
  for (unsigned int i = 0; i < n0; ++i) data[i] = new double[n1];
  
}

cosgrid2D::~cosgrid2D() {
  if (data) {
    for (unsigned int i = 0; i < n0; ++i) delete[] data[i];
    delete[] data;
  }
}

void cosgrid2D::free() {
  if (data != 0) {
    for (unsigned int i = 0; i < n0; ++i) delete[] data[i];
    delete[] data;
  }
  n0 = n1 = 0;
  min0 = max0 = d0 = min1 = max1 = d1 = 0.0;
  data = 0;
}

void cosgrid2D::transpose() {

  double **newdata = new double*[n1];
  for (unsigned int i = 0; i < n1; ++i) {
    newdata[i] = new double[n0];
    for (unsigned int j = 0; j < n0; ++j)
      newdata[i][j] = data[j][i];
  }
  for (unsigned int i = 0; i < n0; ++i) delete[] data[i];
  delete[] data;
  data = newdata;
  swap(n0,n1);
  swap(min0,min1);
  swap(max0,max1);
  swap(d0,d1);
  swap(axisnames[0],axisnames[1]);
  swap(spec[0],spec[1]);

}

/*!
  \param[in] axis Which axis to return
  \param[in] i Position along axis
*/
double cosgrid2D::getAxisVal( unsigned int axis, unsigned int i ) const {
  return getAxisVal( axis, static_cast<double>(i) );
}

/*!
  \param[in] axis Which axis to return
  \param[in] x Position along axis
*/
double cosgrid2D::getAxisVal( unsigned int axis, double x ) const {
  switch (axis) {
  case 0: return min0 + x * d0;
  case 1: return min1 + x * d1;
  default :
    stringstream errstrng("");
    errstrng << "No axis found with index " << axis;
    throw CosFitterExcept("cosgrid2D","getAxisVal",errstrng.str(),1);
  }
}

/*! Returns index of axis with the specified name, throwing a CosFitterExcpt
 if no such axis is found.
 \param[in] axislabel Name of axis to find
 \returns Index of matching axis
 */
unsigned int cosgrid2D::getAxisIndex(const std::string& axislabel) const {
  //Returns index of axis with name axislabel
  if (axisnames[0] == axislabel) return 0;
  if (axisnames[1] == axislabel) return 1;
  stringstream errstrng("");
  errstrng << "No axis found with label " << axislabel;
  throw CosFitterExcept("cosgrid2D","getAxisIndex",errstrng.str(),1);
}

/*!
  Find and return maximum value in array.
  \param[out] i0 The position of the maximum in the array along axis 0
  \param[out] i1 The position of the maximum in the array along axis 1
  \return The maximum value of the data array
*/
double cosgrid2D::getMaximum(unsigned int &i0, unsigned int& i1) const {
  //Find the index of the maximum point and return the maximum value
  double max;
  i0 = i1 = 0;
  max = data[0][0];
  
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      if (data[i][j] > max) {
	max = data[i][j];
	i0 = i; 
	i1 = j; 
      }
  return max;
}

/*!
  \return The sum of all of the elements in the data array 
*/
double cosgrid2D::getTotal() const {
  double total;
  total = 0.0;
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      total += data[i][j];
  return total;
}

/*!
  Normalizes the total sum of all data entries to one.
  \return The sum of all of the elements in the data array prior to 
           normalization.
*/
double cosgrid2D::normalize() {
  double total;

  total = 0.0;
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      total += data[i][j];

  if (total == 0.0) {
    string errstrng = "Total probability equal to 0 -- can't normalize";
    throw CosFitterExcept("cosgrid2D","normalize",errstrng,1);
  }

  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      data[i][j] /= total;

  return total;
}
 

/*!
  Serializes the contents of the data array and the axis labels to a file,
  either in binary or human readable format.
  \param[in] filename the file to write the data to
  \param[in] binary Write as binary instead of human readable
    text (def: false).
*/
int cosgrid2D::writeFile(const std::string& filename, bool binary ) const {
  if (binary) return writeBinaryFile(filename);
  return writeTextFile(filename);
}

/*!
  Reads in a data file as written by writeFile.  The user must specify
  if it is binary or text -- this code is not smart enough to figure it
  out for itself.
  \param[in] filename the file to write the data to
  \param[in] binary Read as binary instead of human readable
    text (def: false).
*/
int cosgrid2D::readFile(const std::string& filename, bool binary ) {
  if (binary) return readBinaryFile(filename);
  return readTextFile(filename);
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a human readable format.  The format of the output file is defined as:
  - The first line is Naxes: 2
  - The second line is axisname0, n0, min0, max0, d0
  - The third line is axisname1, n1, min1, max1, d1
  - This is followed by n0*n1 lines containing the data values in c++ 
     style order (i.e., the axis1 index changes the most rapidly)
  \param[in] filename The file to write the data to
*/
int cosgrid2D::writeTextFile(const std::string& filename) const {


  FILE *fp;
  fp = fopen(filename.c_str(),"w");
  if (!fp) {
    stringstream errstrng("");
    errstrng <<  "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid2D","writeTextFile",errstrng.str(),1);
  }

  fprintf(fp,"Naxes: 2\n");
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[0].c_str(),
	  n0,min0,max0,d0);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[1].c_str(),
	  n1,min1,max1,d1);
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      fprintf(fp,"%-15.8e\n",data[i][j]);
    }
  }

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid2D::writeTextFile
 \param[in] filename The file to read from
 */
int cosgrid2D::readTextFile(const std::string& filename) {
  if (data) {
    for (unsigned int i = 0; i < n0; ++i) delete[] data[i];
    delete[] data;
  }

  ifstream fl( filename.c_str() );
  if (!fl) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid2D","readFile",errstrng.str(),1);
  }

  string line;
  vector<string> words;
  stringstream str("");
  
  unsigned int naxes;
  getline(fl,line);
  utility::stringwords(line,words);
  naxes = atoi( words[1].c_str() );
  if ( naxes != 2 ) {
    fl.close();
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid2D","readFile",errstrng.str(),2);
  }

  //Read first axisname line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[0] = words[0];
  str.str(words[1]); str.clear(); str >> n0;
  str.str(words[2]); str.clear(); str >> min0;
  str.str(words[3]); str.clear(); str >> max0;
  str.str(words[4]); str.clear(); str >> d0;

  //Read second line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[1] = words[0];
  str.str(words[1]); str.clear(); str >> n1;
  str.str(words[2]); str.clear(); str >> min1;
  str.str(words[3]); str.clear(); str >> max1;
  str.str(words[4]); str.clear(); str >> d1;

  data = new double*[n0];
  for (unsigned int i = 0; i < n0; ++i) data[i] = new double[n1];

  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      fl >> data[i][j];
      //There are underflow issues with some versions of gcc
      // so we have to deliberately ignore failbit
      if (fl.fail()) fl.clear(fl.rdstate() & ~std::ios::failbit);
    }
  }

  fl.close();

  return 0;
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a non human readable, binary, format.  The format of the output file 
  is defined as:
  - an integer specifying the number of axes.
  - an integer specifying the number of characters in the zero terminated
     axisname
  - the axis name
  - an integer specifying n0
  - doubles specifying min0, max0, d0
  - a repeat of items 2-5 for the second axis
  - then n0*n1 doubles storing the data in c order
  Note that this will not be machine portable because of endianess and
  the sizes of various outputs.  Unlike the version which writes in text,
  here we don't care if the axisnames have spaces in them.
  \param[in] filename The file to write the data to
*/
int cosgrid2D::writeBinaryFile(const std::string& filename) const {

  FILE *fp;
  fp = fopen(filename.c_str(),"wb");
  if (!fp) {
    stringstream errstrng("");
    errstrng <<  "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid2D","writeBinaryFile",errstrng.str(),1);
  }
  
  //Write number of axes
  unsigned int naxes = 2;
  fwrite( &naxes, sizeof(unsigned int), 1, fp );

  //Write axis 0 title
  unsigned int axisnamelength = axisnames[0].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[0].c_str(), fp );
  
  //Then extended axis information
  fwrite( &n0, sizeof(unsigned int), 1, fp );
  fwrite( &min0, sizeof(double), 1, fp );
  fwrite( &max0, sizeof(double), 1, fp );
  fwrite( &d0, sizeof(double), 1, fp );

  //Then axis 1
  axisnamelength = axisnames[1].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[1].c_str(), fp );
  fwrite( &n1, sizeof(unsigned int), 1, fp );
  fwrite( &min1, sizeof(double), 1, fp );
  fwrite( &max1, sizeof(double), 1, fp );
  fwrite( &d1, sizeof(double), 1, fp );

  //And the data
  //It might be tempting to just try to write n0*n1 entries here,
  // but we have no actualy guarantee that different rows will
  // stored contiguously in memory.  In other words, we have to
  // loop explicitly
  for (unsigned int i = 0; i < n0; ++i ) 
    fwrite( data[i], sizeof(double), n1, fp );

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid2D::writeBinaryFile
 \param[in] filename The file to read from
 */
int cosgrid2D::readBinaryFile(const std::string& filename) {

  if (data) {
    for (unsigned int i = 0; i < n0; ++i) delete[] data[i];
    delete[] data;
  }

  FILE *fp;
  fp = fopen( filename.c_str(), "rb" );
  if (!fp) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid2D","readBinaryFile",errstrng.str(),1);
  }

  //Read number of axes
  unsigned int naxes;
  fread( &naxes, sizeof(unsigned int), 1, fp );
  if ( naxes != 2 ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid2D","readBinaryFile",errstrng.str(),2);
  }
      
  //Read axis 1 title
  unsigned int status;
  status = readBinaryAxisName( fp, axisnames[0] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid2D","readBinaryFile",errstrng.str(),4);
  }
  
  fread( &n0, sizeof(unsigned int), 1, fp );
  fread( &min0, sizeof(double), 1, fp );
  fread( &max0, sizeof(double), 1, fp );
  fread( &d0, sizeof(double), 1, fp );

  //Read axis 2 title
  status = readBinaryAxisName( fp, axisnames[1] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid2D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n1, sizeof(unsigned int), 1, fp );
  fread( &min1, sizeof(double), 1, fp );
  fread( &max1, sizeof(double), 1, fp );
  fread( &d1, sizeof(double), 1, fp );

  if ( (n0 <= 0) || ( n1 <= 0 ) ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "0 axis length in " << filename << endl;
    throw CosFitterExcept("cosgrid2D","readBinaryFile",errstrng.str(),8);
  }
  //Now read the data.  Do the allocation in a seperate step to
  // increase the chances things will be contiguous
  data = new double*[n0];
  for (unsigned int i = 0; i < n0; ++i) data[i] = new double[n1];
  for (unsigned int i = 0; i < n0; ++i ) fread( data[i], sizeof(double), n1, fp );

  fclose(fp);

  return 0;
}

cosgrid2D& cosgrid2D::operator=(double val) {
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      data[i][j] = val;
  return *this;
}

cosgrid2D& cosgrid2D::operator=(const cosgrid2D& input) {
  //Self copy protection
  if (this == &input) return *this;

  axisnames[0] = input.axisnames[0];
  axisnames[1] = input.axisnames[1];
  spec[0] = input.spec[0];
  spec[1] = input.spec[1];
  
  if (data == 0) {
    data = new double*[input.n0];
    for (unsigned int i = 0; i < input.n0; ++i) data[i] = new double[input.n1];
  } else if ( (input.n0 != n0) || (input.n1 != n1) ) {
    for (unsigned int i = 0; i < n0; ++i) delete[] data[i];
    delete[] data;
    data = new double*[input.n0];
    for (unsigned int i = 0; i < input.n0; ++i) data[i] = new double[input.n1];
  }

  for (unsigned int i = 0; i < input.n0; ++i)
    for (unsigned int j = 0; j < input.n1; ++j) 
      data[i][j] = input.data[i][j];

  n0 = input.n0; min0 = input.min0; max0 = input.max0; d0 = input.d0;
  n1 = input.n1; min1 = input.min1; max1 = input.max1; d1 = input.d1;

  return *this;
}

cosgrid2D& cosgrid2D::operator/=(double divisor) {
  for (unsigned int i = 0; i < n0; ++i) 
    for (unsigned int j = 0; j < n1; ++j) 
      data[i][j] /= divisor;
  return *this;
}

cosgrid2D& cosgrid2D::operator*=(const cosgrid2D& input) {
  //Self multiplication
  if (this == &input) {
    for (unsigned int i = 0; i < n0; ++i) 
      for (unsigned int j = 0; j < n1; ++j)
	data[i][j] *= data[i][j];
    return *this;
  }

  //First make sure that the array dimensions and spec codes are
  // compatible
  if (spec[0].second != input.spec[0].second) {
    stringstream errstr("");
    errstr << "Types of axis 0 not compatible in addition" << endl;
    errstr << " One has type: " << spec[0].second << " the other: "
	   << input.spec[0].second;
    errstr << " With name: " << axisnames[0] << " the other: " <<
      input.axisnames[0];
    throw CosFitterExcept("cosgrid2D","operator*=",errstr.str(),1);
  }
  if (spec[1].second != input.spec[1].second) {
    stringstream errstr("");
    errstr << "Types of axis 1 not compatible in addition" << endl;
    errstr << " One has type: " << spec[1].second << " the other: "
	   << input.spec[1].second;
    errstr << " With name: " << axisnames[1] << " the other: " <<
      input.axisnames[1];
    throw CosFitterExcept("cosgrid2D","operator*=",errstr.str(),1);
  }

  double tol0 = (max0 - min0)*0.001;
  if ( (n0 != input.n0) || ( fabs(min0 - input.min0) > tol0 ) ||
       ( fabs(max0 - input.max0) > tol0 ) ) {
    throw CosFitterExcept("cosgrid2D","operator*=",
			     "Axis 0 limits not compatible",2);
  }

  double tol1 = (max1 - min1)*0.001;
  if ( (n1 != input.n1) || ( fabs(min1 - input.min1) > tol1 ) ||
       ( fabs(max1 - input.max1) > tol1 ) ) {
    throw CosFitterExcept("cosgrid2D","operator*=",
			     "Axis 1 limits not compatible",2);
  }

  //All good
  for (unsigned int i = 0; i < input.n0; ++i) {
    for (unsigned int j = 0; j < input.n1; ++j) {
      data[i][j] *= input.data[i][j];
    }
  }

  return *this;
}

cosgrid2D cosgrid2D::operator*(const cosgrid2D& input) {
  //First make sure that the array dimensions and spec codes are
  // compatible
  if (spec[0].second != input.spec[0].second) {
    stringstream errstr("");
    errstr << "Types of axis 0 not compatible in addition" << endl;
    errstr << " One has type: " << spec[0].second << " the other: "
	   << input.spec[0].second;
    errstr << " With name: " << axisnames[0] << " the other: " <<
      input.axisnames[0];
    throw CosFitterExcept("cosgrid2D","operator*",errstr.str(),1);
  }
  if (spec[1].second != input.spec[1].second) {
    stringstream errstr("");
    errstr << "Types of axis 1 not compatible in addition" << endl;
    errstr << " One has type: " << spec[1].second << " the other: "
	   << input.spec[1].second;
    errstr << " With name: " << axisnames[1] << " the other: " <<
      input.axisnames[1];
    throw CosFitterExcept("cosgrid2D","operator*",errstr.str(),1);
  }

  double tol0 = (max0 - min0)*0.001;
  if ( (n0 != input.n0) || ( fabs(min0 - input.min0) > tol0 ) ||
       ( fabs(max0 - input.max0) > tol0 ) ) {
    throw CosFitterExcept("cosgrid2D","operator*",
			     "Axis 0 limits not compatible",2);
  }

  double tol1 = (max1 - min1)*0.001;
  if ( (n1 != input.n1) || ( fabs(min1 - input.min1) > tol1 ) ||
       ( fabs(max1 - input.max1) > tol1 ) ) {
    throw CosFitterExcept("cosgrid2D","operator*",
			     "Axis 1 limits not compatible",2);
  }

  //All good
  cosgrid2D outcosgrid( input );
  for (unsigned int i = 0; i < input.n0; ++i)
    for (unsigned int j = 0; j < input.n1; ++j)
      outcosgrid.data[i][j] *= data[i][j];

  return outcosgrid;
}

cosgrid2D& cosgrid2D::operator+=(const cosgrid2D& input) {
  //Self addition, reduce memory accesses
  if (this == &input) {
    for (unsigned int i = 0; i < n0; ++i) 
      for (unsigned int j = 0; j < n1; ++j)
	data[i][j] *= 2.0;
    return *this;
  }

  //First make sure that the array dimensions and spec codes are
  // compatible
  if (spec[0].second != input.spec[0].second) {
    stringstream errstr("");
    errstr << "Types of axis 0 not compatible in addition" << endl;
    errstr << " One has type: " << spec[0].second << " the other: "
	   << input.spec[0].second;
    errstr << " With name: " << axisnames[0] << " the other: " <<
      input.axisnames[0];
    throw CosFitterExcept("cosgrid2D","operator+=",errstr.str(),1);
  }
  if (spec[1].second != input.spec[1].second) {
    stringstream errstr("");
    errstr << "Types of axis 1 not compatible in addition" << endl;
    errstr << " One has type: " << spec[1].second << " the other: "
	   << input.spec[1].second;
    errstr << " With name: " << axisnames[1] << " the other: " <<
      input.axisnames[1];
    throw CosFitterExcept("cosgrid2D","operator+=",errstr.str(),1);
  }

  double tol0 = (max0 - min0)*0.001;
  if ( (n0 != input.n0) || ( fabs(min0 - input.min0) > tol0 ) ||
       ( fabs(max0 - input.max0) > tol0 ) ) {
    throw CosFitterExcept("cosgrid2D","operator+=",
			     "Axis 0 limits not compatible",2);
  }

  double tol1 = (max1 - min1)*0.001;
  if ( (n1 != input.n1) || ( fabs(min1 - input.min1) > tol1 ) ||
       ( fabs(max1 - input.max1) > tol1 ) ) {
    throw CosFitterExcept("cosgrid2D","operator+=",
			     "Axis 1 limits not compatible",2);
  }

  //All good
  for (unsigned int i = 0; i < input.n0; ++i) {
    for (unsigned int j = 0; j < input.n1; ++j) {
      data[i][j] += input.data[i][j];
    }
  }

  return *this;
}

cosgrid1D cosgrid2D::collapseAlongZero() const {
  cosgrid1D result(n1,min1,max1,d1,axisnames[1],spec[1]);
  
  double *dataptr = result.getData(); 
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    val = 0.0;
    for (unsigned int i = 0; i < n0; ++i) val += data[i][j];
    dataptr[j] = val;
  }
  
  return result;
}

cosgrid1D cosgrid2D::collapseAlongOne() const {
  cosgrid1D result(n0,min0,max0,d0,axisnames[0],spec[0]);
  
  double *dataptr = result.getData(); 
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    val = 0.0;
    for (unsigned int j = 0; j < n1; ++j) val += data[i][j];
    dataptr[i] = val;
  }
  
  return result;
}

/*!
  Collapse the grid along the specified axis by summing.
  \param[in] axisnum Index of axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid1D containing the collapsed data array
*/
cosgrid1D cosgrid2D::collapseAlongAxis(unsigned int axisnum) const {
  if (axisnum == 0) return collapseAlongZero();
  if (axisnum == 1) return collapseAlongOne();
  stringstream errstr("");
  errstr << "Invalid axis number " << axisnum;
  throw CosFitterExcept("cosgrid2D","collapseAlongAxis",errstr.str(),1);
}

/*!
  Collapse the grid along the specified axis by summing.
  \param[in] axisname Name of axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid1D containing the collapsed data array
*/
cosgrid1D cosgrid2D::collapseAlongAxis(const std::string& axisname) const {
  return collapseAlongAxis( getAxisIndex(axisname) );
}



//////////////////////////////////
//           cosgrid3D          //
//////////////////////////////////

cosgrid3D::cosgrid3D() {
  data = 0;
  n0 = n1 = n2 = 0;
  min0 = max0 = d0 = min1 = max1 = d1 = min2 = max2 = d2 = 0;
}

cosgrid3D::cosgrid3D(const cosgrid3D& input) {
  n0 = input.n0;
  min0 = input.min0;
  max0 = input.max0;
  d0 = input.d0;
  n1 = input.n1;
  min1 = input.min1;
  max1 = input.max1;
  d1 = input.d1;
  n2 = input.n2;
  min2 = input.min2;
  max2 = input.max2;
  d2 = input.d2;
  axisnames[0] = input.axisnames[0];
  axisnames[1] = input.axisnames[1];
  axisnames[2] = input.axisnames[2];
  spec[0] = input.spec[0];
  spec[1] = input.spec[1];
  spec[2] = input.spec[2];
  if (input.data != 0) {
    data = new double**[ n0 ];
    for (unsigned int i = 0; i < n0; ++i) {
      data[i] = new double*[n1];
      for (unsigned int j = 0; j < n1; ++j) {
	data[i][j] = new double[n2];
	for (unsigned int k =0; k < n2; ++k) 
	  data[i][j][k] = input.data[i][j][k];
      }
    }
  }
}


/*!
  \param[in] n0 Number of data entries along axis 0
  \param[in] min0 Minimum axis value along axis 0
  \param[in] max0 Maximum axis value along axis 0
  \param[in] d0 Step size for axis along axis 0
  \param[in] label0 Name of axis 0 parameter.  Should not have spaces in it.
  \param[in] spec0 Param spec for parameter 0
  \param[in] n1 Number of data entries along axis 1
  \param[in] min1 Minimum axis value along axis 1
  \param[in] max1 Maximum axis value along axis 1
  \param[in] d1 Step size for axis along axis 1
  \param[in] label1 Name of axis 1 parameter.  Should not have spaces in it.
  \param[in] spec1 Param spec for parameter 1
  \param[in] n2 Number of data entries along axis 2
  \param[in] min2 Minimum axis value along axis 2
  \param[in] max2 Maximum axis value along axis 2
  \param[in] d2 Step size for axis along axis 2
  \param[in] label2 Name of axis 2 parameter.  Should not have spaces in it.
  \param[in] spec2 Param spec for parameter 2
*/
cosgrid3D::cosgrid3D(unsigned int n0, double min0, double max0, double d0,
		     std::string label0, param_tags::paramspec spec0,
		     unsigned int n1, double min1, double max1, double d1, 
		     std::string label1, param_tags::paramspec spec1,
		     unsigned int n2, double min2, double max2, double d2, 
		     std::string label2, param_tags::paramspec spec2) :
  n0(n0), min0(min0), max0(max0), d0(d0), n1(n1),
  min1(min1), max1(max1), d1(d1), n2(n2), min2(min2), max2(max2), d2(d2) {

  axisnames[0] = label0;
  axisnames[1] = label1;
  axisnames[2] = label2;
  spec[0] = spec0;
  spec[1] = spec1;
  spec[2] = spec2;
  
  data = new double**[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double*[n1];
    for (unsigned int j = 0; j < n1; ++j)
      data[i][j] = new double[n2];
  }
  
}

/*!
  \param[in] param0 Information about parameter describing axis 0
  \param[in] param1 Information about parameter describing axis 1
  \param[in] param2 Information about parameter describing axis 2
*/
cosgrid3D::cosgrid3D(const param_struct& param0, const param_struct& param1,
		     const param_struct& param2) :
  n0(param0.n), min0(param0.min), max0(param0.max), d0(param0.dval), 
  n1(param1.n), min1(param1.min), max1(param1.max), d1(param1.dval), 
  n2(param2.n), min2(param2.min), max2(param2.max), d2(param2.dval) {

  axisnames[0] = param0.name;
  axisnames[1] = param1.name;
  axisnames[2] = param2.name;
  spec[0] = param0.param_spec;
  spec[1] = param1.param_spec;
  spec[2] = param2.param_spec;
  
  data = new double**[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double*[n1];
    for (unsigned int j = 0; j < n1; ++j)
      data[i][j] = new double[n2];
  }
  
}

cosgrid3D::~cosgrid3D() {
  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) delete[] data[i][j];
      delete[] data[i];
    }
    delete[] data;
  }
}

void cosgrid3D::free() {
  if (data != 0) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) delete[] data[i][j];
      delete[] data[i];
    }
    delete[] data;
  }
  n0 = n1 = n2 = 0;
  min0 = max0 = d0 = min1 = max1 = d1 = min2 = max2 = d2 = 0;
  data = 0;
}

/*!
  \param[in] axis Which axis to return
  \param[in] i Position along axis
*/
double cosgrid3D::getAxisVal( unsigned int axis, unsigned int i ) const {
  return getAxisVal( axis, static_cast<double>(i) );
}

/*!
  \param[in] axis Which axis to return
  \param[in] x Position along axis
*/
double cosgrid3D::getAxisVal( unsigned int axis, double x ) const {
  switch (axis) {
  case 0: return min0 + x * d0;
  case 1: return min1 + x * d1;
  case 2: return min2 + x * d2;
  default :
    stringstream errstrng("");
    errstrng << "No axis found with index " << axis;
    throw CosFitterExcept("cosgrid3D","getAxisVal",errstrng.str(),1);
  }
}

/*! Returns index of axis with the specified name, throwing a CosFitterExcpt
 if no such axis is found.
 \param[in] axislabel Name of axis to find
 \returns Index of matching axis
 */
unsigned int cosgrid3D::getAxisIndex(const std::string& axislabel) const {
  if (axisnames[0] == axislabel) return 0;
  if (axisnames[1] == axislabel) return 1;
  if (axisnames[2] == axislabel) return 2;
  stringstream errstrng("");
  errstrng << "No axis found with label " <<  axislabel;
  throw CosFitterExcept("cosgrid3D","getAxisIndex",errstrng.str(),1);
}

/*!
  Find and return maximum value in array.
  \param[out] i0 The position of the maximum in the array along axis 0
  \param[out] i1 The position of the maximum in the array along axis 1
  \param[out] i2 The position of the maximum in the array along axis 2
  \return The maximum value of the data array
*/
double cosgrid3D::getMaximum(unsigned int &i0, unsigned int& i1, 
			     unsigned int &i2) const {
  //Find the index of the maximum point
  double max;
  i0 = i1 = i2 = 0;
  max = data[0][0][0];
  
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	if (data[i][j][k] > max) {
	  max = data[i][j][k];
	      i0 = i; 
	      i1 = j; 
	      i2 = k; 
	}
  return max;
}

/*!
  \return The sum of all of the elements in the data array
*/
double cosgrid3D::getTotal() const {
  double total;

  total = 0.0;
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	total += data[i][j][k];
  return total;
}

/*!
  Normalizes the total sum of all data entries to one.
  \return The sum of all of the elements in the data array prior to 
           normalization.
*/
double cosgrid3D::normalize() {
  double total;

  total = 0.0;
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	total += data[i][j][k];

  if (total == 0.0) {
    std::string errstrng = "Total probability equal to 0 -- can't normalize";
    throw CosFitterExcept("cosgrid3D","normalize",errstrng,1);
  }

  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	data[i][j][k] /= total;

  return total;
}
 
/*!
  Serializes the contents of the data array and the axis labels to a file,
  either in binary or human readable format.
  \param[in] filename the file to write the data to
  \param[in] binary Write as binary instead of human readable
    text (def: false).
*/
int cosgrid3D::writeFile(const std::string& filename, bool binary ) const {
  if (binary) return writeBinaryFile(filename);
  return writeTextFile(filename);
}

/*!
  Reads in a data file as written by writeFile.  The user must specify
  if it is binary or text -- this code is not smart enough to figure it
  out for itself.
  \param[in] filename the file to write the data to
  \param[in] binary Read as binary instead of human readable
    text (def: false).
*/
int cosgrid3D::readFile(const std::string& filename, bool binary ) {
  if (binary) return readBinaryFile(filename);
  return readTextFile(filename);
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a human readable format.  The format of the output file is defined as:
  - The first line is Naxes: 3
  - The second line is axisname0, n0, min0, max0, d0
  - The third line is axisname1, n1, min1, max1, d1
  - The fourth line is axisname2, n2, min2, max2, d2
  - This is followed by n0*n1*n2 lines containing the data values in c++ 
     style order (i.e., the axis2 index changes the most rapidly)
  \param[in] filename The file to write the data to
*/
int cosgrid3D::writeTextFile(const std::string& filename) const {

  FILE *fp;
  fp = fopen(filename.c_str(),"w");
  if (!fp) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid3D","writeTextFile",errstrng.str(),1);
  }

  fprintf(fp,"Naxes: 3\n");
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[0].c_str(),
	  n0,min0,max0,d0);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[1].c_str(),
	  n1,min1,max1,d1);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[2].c_str(),
	  n2,min2,max2,d2);
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k)
	fprintf(fp,"%-15.8e\n",data[i][j][k]);
    }
  }

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid3D::writeTextFile
 \param[in] filename The file to read from
 */
int cosgrid3D::readTextFile(const std::string& filename) {
  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j)
	delete[] data[i][j];
      delete[] data[i];
    }
    delete[] data;
  }  

  ifstream fl( filename.c_str() );
  if (!fl) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid3D","readFile",errstrng.str(),1);
  }

  string line;
  vector<string> words;
  stringstream str("");
  
  unsigned int naxes;
  getline(fl,line);
  utility::stringwords(line,words);
  naxes = atoi( words[1].c_str() );
  if ( naxes != 3 ) {
    fl.close();
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid3D","readFile",errstrng.str(),2);
  }

  //Read first axisname line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[0] = words[0];
  str.str(words[1]); str.clear(); str >> n0;
  str.str(words[2]); str.clear(); str >> min0;
  str.str(words[3]); str.clear(); str >> max0;
  str.str(words[4]); str.clear(); str >> d0;

  //Read second line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[1] = words[0];
  str.str(words[1]); str.clear(); str >> n1;
  str.str(words[2]); str.clear(); str >> min1;
  str.str(words[3]); str.clear(); str >> max1;
  str.str(words[4]); str.clear(); str >> d1;

  //And third
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[2] = words[0];
  str.str(words[1]); str.clear(); str >> n2;
  str.str(words[2]); str.clear(); str >> min2;
  str.str(words[3]); str.clear(); str >> max2;
  str.str(words[4]); str.clear(); str >> d2;

  data = new double**[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double*[n1];
    for (unsigned int j = 0; j < n1; ++j) 
      data[i][j] = new double[n2];
  }

  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	fl >> data[i][j][k];
	//There are underflow issues with some versions of gcc
	// so we have to deliberately ignore failbit
	if (fl.fail()) fl.clear(fl.rdstate() & ~std::ios::failbit);
      }
    }
  }

  fl.close();

  return 0;
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a non human readable, binary, format.  The format of the output file 
  is defined as:
  - an integer specifying the number of axes.
  - an integer specifying the number of characters in the zero terminated
     axisname
  - the axis name
  - an integer specifying n0
  - doubles specifying min0, max0, d0
  - a repeat of items 2-5 for the second axis, and again for the third
  - then n0*n1*n2 doubles storing the data in c order
  Note that this will not be machine portable because of endianess and
  the sizes of various outputs.  Unlike the version which writes in text,
  here we don't care if the axisnames have spaces in them.
  \param[in] filename The file to write the data to
*/
int cosgrid3D::writeBinaryFile(const std::string& filename) const {

  FILE *fp;
  fp = fopen(filename.c_str(),"wb");
  if (!fp) {
    stringstream errstrng("");
    errstrng <<  "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid3D","writeBinaryFile",errstrng.str(),1);
  }
  
  //Write number of axes
  unsigned int naxes = 3;
  fwrite( &naxes, sizeof(unsigned int), 1, fp );

  //Write axis 0 title
  unsigned int axisnamelength = axisnames[0].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[0].c_str(), fp );
  
  //Then extended axis information
  fwrite( &n0, sizeof(unsigned int), 1, fp );
  fwrite( &min0, sizeof(double), 1, fp );
  fwrite( &max0, sizeof(double), 1, fp );
  fwrite( &d0, sizeof(double), 1, fp );

  //Then axis 1
  axisnamelength = axisnames[1].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[1].c_str(), fp );
  fwrite( &n1, sizeof(unsigned int), 1, fp );
  fwrite( &min1, sizeof(double), 1, fp );
  fwrite( &max1, sizeof(double), 1, fp );
  fwrite( &d1, sizeof(double), 1, fp );

  //Then axis 2
  axisnamelength = axisnames[2].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[2].c_str(), fp );
  fwrite( &n2, sizeof(unsigned int), 1, fp );
  fwrite( &min2, sizeof(double), 1, fp );
  fwrite( &max2, sizeof(double), 1, fp );
  fwrite( &d2, sizeof(double), 1, fp );

  //And the data
  //It might be tempting to just try to write n0*n1*n2 entries here,
  // but we have no actualy guarantee that different rows will
  // stored contiguously in memory.  In other words, we have to
  // loop explicitly
  for (unsigned int i = 0; i < n0; ++i ) 
    for ( unsigned int j = 0; j < n1; ++j )
      fwrite( data[i][j], sizeof(double), n2, fp );

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid3D::writeBinaryFile
 \param[in] filename The file to read from
 */
int cosgrid3D::readBinaryFile(const std::string& filename) {

  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j)
	delete[] data[i][j];
      delete[] data[i];
    }
    delete[] data;
  }

  FILE *fp;
  fp = fopen( filename.c_str(), "rb" );
  if (!fp) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid3D","readBinaryFile",errstrng.str(),1);
  }

  //Read number of axes
  unsigned int naxes;
  fread( &naxes, sizeof(unsigned int), 1, fp );
  if ( naxes != 3 ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid3D","readBinaryFile",errstrng.str(),2);
  }
      
  //Read axis 1 title
  int status;
  status = readBinaryAxisName( fp, axisnames[0] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid3D","readBinaryFile",errstrng.str(),4);
  }
  
  fread( &n0, sizeof(unsigned int), 1, fp );
  fread( &min0, sizeof(double), 1, fp );
  fread( &max0, sizeof(double), 1, fp );
  fread( &d0, sizeof(double), 1, fp );

  //Read axis 2
  status = readBinaryAxisName( fp, axisnames[1] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid3D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n1, sizeof(unsigned int), 1, fp );
  fread( &min1, sizeof(double), 1, fp );
  fread( &max1, sizeof(double), 1, fp );
  fread( &d1, sizeof(double), 1, fp );

  //Read axis 3
  status = readBinaryAxisName( fp, axisnames[2] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid3D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n2, sizeof(unsigned int), 1, fp );
  fread( &min2, sizeof(double), 1, fp );
  fread( &max2, sizeof(double), 1, fp );
  fread( &d2, sizeof(double), 1, fp );

  if ( (n0 <= 0) || ( n1 <= 0 ) || ( n2 <= 0 ) ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "0 axis length in " << filename << endl;
    throw CosFitterExcept("cosgrid3D","readBinaryFile",errstrng.str(),8);
  }
  //Now read the data.  Do the allocation in a seperate step to
  // increase the chances things will be contiguous
  data = new double**[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double*[n1];
    for ( unsigned int j = 0; j < n1; ++j )
      data[i][j] = new double[n2];
  }
  for (unsigned int i = 0; i < n0; ++i ) 
    for (unsigned int j = 0; j < n1; ++j)
      fread( data[i][j], sizeof(double), n2, fp );

  fclose(fp);

  return 0;
}

cosgrid3D& cosgrid3D::operator=(double val) {
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	data[i][j][k] = val;
  return *this;
}

cosgrid3D& cosgrid3D::operator=(const cosgrid3D& input) {
  //Self copy protection
  if (this == &input) return *this;

  axisnames[0] = input.axisnames[0];
  axisnames[1] = input.axisnames[1];
  axisnames[2] = input.axisnames[2];
  spec[0] = input.spec[0];
  spec[1] = input.spec[1];
  spec[2] = input.spec[2];

  if (data == 0) {
    data = new double**[input.n0];
    for (unsigned int i = 0; i < input.n0; ++i) {
      data[i] = new double*[input.n1];
      for (unsigned int j = 0; j < input.n1; ++j) 
	data[i][j] = new double[input.n2];
    }
  } else if ( (input.n0 != n0) || (input.n1 != n1) || (input.n2 != n2) ) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j)
	delete[] data[i][j];
      delete[] data[i];
    }
    delete[] data;
    data = new double**[input.n0];
    for (unsigned int i = 0; i < input.n0; ++i) {
      data[i] = new double*[input.n1];
      for (unsigned int j = 0; j < input.n1; ++j) 
	data[i][j] = new double[input.n2];
    }
  }

  for (unsigned int i = 0; i < input.n0; ++i) 
    for (unsigned int j = 0; j < input.n1; ++j)
      for (unsigned int k = 0; k < input.n2; ++k)
	data[i][j][k] = input.data[i][j][k];
  
  n0 = input.n0; min0 = input.min0; max0 = input.max0; d0 = input.d0;
  n1 = input.n1; min1 = input.min1; max1 = input.max1; d1 = input.d1;
  n2 = input.n2; min2 = input.min2; max2 = input.max2; d2 = input.d2;

  return *this;
}

cosgrid3D& cosgrid3D::operator/=(double divisor) {
  //Divides all entries in cosgrid by specified value
  for (unsigned int i = 0; i < n0; ++i) 
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k) 
	data[i][j][k] /= divisor;
  return *this;
}

cosgrid3D& cosgrid3D::operator*=(const cosgrid3D& input) {
  double **ptr2, *ptr1;

  //Self multiplication
  if (this == &input) {
    for (unsigned int i = 0; i < n0; ++i) {
      ptr2 = data[i];
      for (unsigned int j = 0; j < n1; ++j) {
	ptr1 = ptr2[j];
	for (unsigned int k = 0; k < n2; ++k)
	  ptr1[k] *= ptr1[k];
      }
    }
    return *this;
  }

  //First make sure that the array dimensions and spec codes are
  // compatible
  if (spec[0].second != input.spec[0].second) {
    stringstream errstr("");
    errstr << "Types of axis 0 not compatible in multiplication" << endl;
    errstr << " One has type: " << spec[0].second << " the other: "
	   << input.spec[0].second;
    errstr << " With name: " << axisnames[0] << " the other: " <<
      input.axisnames[0];
    throw CosFitterExcept("cosgrid3D","operator*=",errstr.str(),1);
  }
  if (spec[1].second != input.spec[1].second) {
    stringstream errstr("");
    errstr << "Types of axis 1 not compatible in multiplication" << endl;
    errstr << " One has type: " << spec[1].second << " the other: "
	   << input.spec[1].second;
    errstr << " With name: " << axisnames[1] << " the other: " <<
      input.axisnames[1];
    throw CosFitterExcept("cosgrid3D","operator*=",errstr.str(),1);
  }
  if (spec[2].second != input.spec[2].second) {
    stringstream errstr("");
    errstr << "Types of axis 2 not compatible in multiplication" << endl;
    errstr << " One has type: " << spec[2].second << " the other: "
	   << input.spec[2].second;
    errstr << " With name: " << axisnames[2] << " the other: " <<
      input.axisnames[2];
    throw CosFitterExcept("cosgrid3D","operator*=",errstr.str(),1);
  }

  double tol = (max0 - min0)*0.001;
  if ( (n0 != input.n0) || ( fabs(min0 - input.min0) > tol ) ||
       ( fabs(max0 - input.max0) > tol ) ) {
    throw CosFitterExcept("cosgrid3D","operator*=",
			     "Axis 0 limits not compatible",2);
  }
  tol = (max1 - min1)*0.001;
  if ( (n1 != input.n1) || ( fabs(min1 - input.min1) > tol ) ||
       ( fabs(max1 - input.max1) > tol ) ) {
    throw CosFitterExcept("cosgrid3D","operator*=",
			     "Axis 1 limits not compatible",2);
  }
  tol = (max2 - min2)*0.001;
  if ( (n2 != input.n2) || ( fabs(min2 - input.min2) > tol ) ||
       ( fabs(max2 - input.max2) > tol ) ) {
    throw CosFitterExcept("cosgrid3D","operator*=",
			     "Axis 2 limits not compatible",2);
  }

  //All good
  double **i_ptr2, *i_ptr1;
  for (unsigned int i = 0; i < input.n0; ++i) {
    ptr2 = data[i]; i_ptr2 = input.data[i];
    for (unsigned int j = 0; j < input.n1; ++j) {
      ptr1 = ptr2[j]; i_ptr1 = i_ptr2[j];
      for (unsigned int k = 0; k < input.n2; ++k)
	ptr1[k] *= i_ptr1[k];
    }
  }

  return *this;
}

cosgrid3D& cosgrid3D::operator+=(const cosgrid3D& input) {
  //Self addition, reduce memory accesses
  if (this == &input) {
    for (unsigned int i = 0; i < n0; ++i) 
      for (unsigned int j = 0; j < n1; ++j)
	for (unsigned int k = 0; k < n2; ++k)
	  data[i][j][k] *= 2.0;
    return *this;
  }

  //First make sure that the array dimensions and param codes are
  // compatible
  if (spec[0].second != input.spec[0].second) {
    stringstream errstr("");
    errstr << "Types of axis 0 not compatible in addition" << endl;
    errstr << " One has type: " << spec[0].second << " the other: "
	   << input.spec[0].second;
    errstr << " With name: " << axisnames[0] << " the other: " <<
      input.axisnames[0];
    throw CosFitterExcept("cosgrid3D","operator+=",errstr.str(),1);
  }
  if (spec[1].second != input.spec[1].second) {
    stringstream errstr("");
    errstr << "Types of axis 1 not compatible in addition" << endl;
    errstr << " One has type: " << spec[1].second << " the other: "
	   << input.spec[1].second;
    errstr << " With name: " << axisnames[1] << " the other: " <<
      input.axisnames[1];
    throw CosFitterExcept("cosgrid3D","operator+=",errstr.str(),1);
  }
  if (spec[2].second != input.spec[2].second) {
    stringstream errstr("");
    errstr << "Types of axis 2 not compatible in addition" << endl;
    errstr << " One has type: " << spec[2].second << " the other: "
	   << input.spec[2].second;
    errstr << " With name: " << axisnames[2] << " the other: " <<
      input.axisnames[2];
    throw CosFitterExcept("cosgrid3D","operator+=",errstr.str(),1);
  }
  
  double tol0 = (max0 - min0)*0.001;
  if ( (n0 != input.n0) || ( fabs(min0 - input.min0) > tol0 ) ||
       ( fabs(max0 - input.max0) > tol0 ) ) {
    throw CosFitterExcept("cosgrid3D","operator+=",
			     "Axis 0 limits not compatible",2);
  }

  double tol1 = (max1 - min1)*0.001;
  if ( (n1 != input.n1) || ( fabs(min1 - input.min1) > tol1 ) ||
       ( fabs(max1 - input.max1) > tol1 ) ) {
    throw CosFitterExcept("cosgrid3D","operator+=",
			     "Axis 1 limits not compatible",2);
  }

  double tol2 = (max2 - min2)*0.001;
  if ( (n2 != input.n2) || ( fabs(min2 - input.min2) > tol2 ) ||
       ( fabs(max2 - input.max2) > tol2 ) ) {
    throw CosFitterExcept("cosgrid3D","operator+=",
			     "Axis 2 limits not compatible",2);
  }

  //All good
  for (unsigned int i = 0; i < input.n0; ++i) 
    for (unsigned int j = 0; j < input.n1; ++j) 
      for (unsigned int k = 0; k < input.n2; ++k)
	data[i][j][k] += input.data[i][j][k];

  return *this;
}


cosgrid2D cosgrid3D::collapseAlongZero() const {
  //Collapse 3D grid along the zeroth axis
  
  cosgrid2D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2]);
  
  double **dataptr = result.getData(); 
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int k = 0; k < n2; ++k) {
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i) val += data[i][j][k];
      dataptr[j][k] = val;
    }
  }
  
  return result;
}

cosgrid2D cosgrid3D::collapseAlongOne() const {
  //Collapse 3D grid along the 1st axis
  
  cosgrid2D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n2,min2,max2,d2,axisnames[2],spec[1]);
  
  double **dataptr = result.getData();
  double val;
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int k = 0; k < n2; ++k) {
      val = 0.0;
      for (unsigned int j = 0; j < n1; ++j) val += data[i][j][k];
      dataptr[i][k] = val;
    }
  }
  
  return result;
}

cosgrid2D cosgrid3D::collapseAlongTwo() const {
  //Collapse 3D grid along the 2nd axis
  
  cosgrid2D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1]);
 
  double **dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      val = 0.0;
      for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k];
      dataptr[i][j] = val;
    }
  }

  return result;
}

cosgrid1D cosgrid3D::collapseAlongZeroOne() const {
  //Collapse 3D grid along the 0th and 1st axis
  
  cosgrid1D result(n2,min2,max2,d2,axisnames[2],spec[2]);
 
  double *dataptr = result.getData();
  double val;
  for (unsigned int k = 0; k < n2; ++k) {
    val = 0.0;
    for (unsigned int i = 0; i < n0; ++i)
      for (unsigned int j = 0; j < n1; ++j) val += data[i][j][k];
    dataptr[k] = val;
  }
  return result;
}

cosgrid1D cosgrid3D::collapseAlongZeroTwo() const {
  //Collapse 3D grid along the 0th and 2nd axis
  
  cosgrid1D result(n1,min1,max1,d1,axisnames[1],spec[1]);
 
  double *dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    val = 0.0;
    for (unsigned int i = 0; i < n0; ++i)
      for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k];
    dataptr[j] = val;
  }

  return result;
}

cosgrid1D cosgrid3D::collapseAlongOneTwo() const {
  //Collapse 3D grid along the 1st and 2nd axis
  
  cosgrid1D result(n0,min0,max0,d0,axisnames[0],spec[0]);
 
  double *dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    val = 0.0;
    for (unsigned int j = 0; j < n1; ++j)
      for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k];
    dataptr[i] = val;
  }

  return result;
}

/*!
  Collapse the grid along the specified axis by summing.
  \param[in] axisnum Index of axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid2D containing the collapsed data array
*/
cosgrid2D cosgrid3D::collapseAlongAxis(unsigned int axisnum) const {
  if (axisnum == 0) return collapseAlongZero();
  if (axisnum == 1) return collapseAlongOne();
  if (axisnum == 2) return collapseAlongTwo();
  stringstream errstr("");
  errstr << "Invalid axis number " << axisnum;
  throw CosFitterExcept("cosgrid3D","collapseAlongAxis",errstr.str(),1);
}

/*!
  Collapse the grid along the specified axis by summing.
  \param[in] axisname Nameof axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid2D containing the collapsed data array
*/
cosgrid2D cosgrid3D::collapseAlongAxis(const std::string& axisname) const {
  return collapseAlongAxis( getAxisIndex(axisname) );
}

/*!
  Collapse the grid along the specified axes by summing.
  \param[in] axisnum1 Index of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisnum2 Index of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid1D containing the collapsed data array
*/
cosgrid1D cosgrid3D::collapseAlongTwoAxes(unsigned int axisnum1, 
					  unsigned int axisnum2) const {
  if (axisnum1 == 0) {
    if (axisnum2 == 1) return collapseAlongZeroOne();
    if (axisnum2 == 2) return collapseAlongZeroTwo();
  } else if (axisnum1 == 1) {
    if (axisnum2 == 0) return collapseAlongZeroOne();
    if (axisnum2 == 2) return collapseAlongOneTwo();
  } else if (axisnum1 == 2) {
    if (axisnum2 == 0) return collapseAlongZeroTwo();
    if (axisnum2 == 1) return collapseAlongOneTwo();
  }
  stringstream errstr("");
  errstr << "Can't collapse along axes: " << axisnum1 << " " << axisnum2;
  throw CosFitterExcept("cosgrid3D","collapseAlongTwoAxes",
			   errstr.str(),1);

}

/*!
  Collapse the grid along the specified axes by summing.
  \param[in] axisname1 Name of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisname2 Name of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid1D containing the collapsed data array
*/
cosgrid1D cosgrid3D::collapseAlongTwoAxes(const std::string& axisname1,
					  const std::string& axisname2) const {
  //Collapse grid along 2 axes
  return collapseAlongTwoAxes( getAxisIndex(axisname1),
			       getAxisIndex(axisname2) );

}



//////////////////////////////////
//           cosgrid4D          //
//////////////////////////////////

cosgrid4D::cosgrid4D() {
  data = 0;
  n0 = n1 = n2 = n3 = 0;
  min0 = max0 = d0 = min1 = max1 = d1 = min2 = max2 = d2 = 0.0;
  min3 = max3 = d3 = 0.0;
}


cosgrid4D::cosgrid4D(const cosgrid4D& input) {
  n0 = input.n0; min0 = input.min0; max0 = input.max0; d0 = input.d0;
  n1 = input.n1; min1 = input.min1; max1 = input.max1; d1 = input.d1;
  n2 = input.n2; min2 = input.min2; max2 = input.max2; d2 = input.d2;
  n3 = input.n3; min3 = input.min3; max3 = input.max3; d3 = input.d3;
  axisnames[0] = input.axisnames[0];
  axisnames[1] = input.axisnames[1];
  axisnames[2] = input.axisnames[2];
  axisnames[3] = input.axisnames[3];
  spec[0] = input.spec[0];
  spec[1] = input.spec[1];
  spec[2] = input.spec[2];
  spec[3] = input.spec[3];

  if (input.data != 0) {
    data = new double***[ n0 ];
    for (unsigned int i = 0; i < n0; ++i) {
      data[i] = new double**[n1];
      for (unsigned int j = 0; j < n1; ++j) {
	data[i][j] = new double*[n2];
	for (unsigned int k = 0; k < n2; ++k) { 
	  data[i][j][k] = new double[n3];
	  for (unsigned int m = 0; m < n3; ++m) 
	    data[i][j][k][m] = input.data[i][j][k][m];
	}
      }
    }
  }
}

/*!
  \param[in] n0 Number of data entries along axis 0
  \param[in] min0 Minimum axis value along axis 0
  \param[in] max0 Maximum axis value along axis 0
  \param[in] d0 Step size for axis along axis 0
  \param[in] label0 Name of axis 0 parameter.  Should not have spaces in it.
  \param[in] spec0 Param spec for parameter 0
  \param[in] n1 Number of data entries along axis 1
  \param[in] min1 Minimum axis value along axis 1
  \param[in] max1 Maximum axis value along axis 1
  \param[in] d1 Step size for axis along axis 1
  \param[in] label1 Name of axis 1 parameter.  Should not have spaces in it.
  \param[in] spec1 Param spec for parameter 1
  \param[in] n2 Number of data entries along axis 2
  \param[in] min2 Minimum axis value along axis 2
  \param[in] max2 Maximum axis value along axis 2
  \param[in] d2 Step size for axis along axis 2
  \param[in] label2 Name of axis 2 parameter.  Should not have spaces in it.
  \param[in] spec2 Param spec for parameter 2
  \param[in] n3 Number of data entries along axis 3
  \param[in] min3 Minimum axis value along axis 3
  \param[in] max3 Maximum axis value along axis 3
  \param[in] d3 Step size for axis along axis 3
  \param[in] label3 Name of axis 3 parameter.  Should not have spaces in it.
  \param[in] spec3 Param spec for parameter 3
*/
cosgrid4D::cosgrid4D(unsigned int n0, double min0, double max0, double d0,
		     std::string label0, param_tags::paramspec spec0,
		     unsigned int n1, double min1, double max1,  double d1, 
		     std::string label1, param_tags::paramspec spec1,
		     unsigned int n2, double min2, double max2, double d2, 
		     std::string label2, param_tags::paramspec spec2, 
		     unsigned int n3, double min3, double max3, double d3, 
		     std::string label3, param_tags::paramspec spec3) :
  n0(n0), min0(min0), max0(max0), d0(d0), n1(n1),
  min1(min1), max1(max1), d1(d1), n2(n2), min2(min2), max2(max2), d2(d2),
  n3(n3), min3(min3), max3(max3), d3(d3) {

  axisnames[0] = label0;
  axisnames[1] = label1;
  axisnames[2] = label2;
  axisnames[3] = label3;
  spec[0] = spec0;
  spec[1] = spec1;
  spec[2] = spec2;
  spec[3] = spec3;

  data = new double***[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double**[n1];
    for (unsigned int j = 0; j < n1; ++j) {
      data[i][j] = new double*[n2];
      for (unsigned int k = 0; k < n2; ++k) {
	data[i][j][k] = new double[n3];
      }
    }
  }
}

/*!
  \param[in] param0 Information about parameter describing axis 0
  \param[in] param1 Information about parameter describing axis 1
  \param[in] param2 Information about parameter describing axis 2
  \param[in] param3 Information about parameter describing axis 3
*/
cosgrid4D::cosgrid4D(const param_struct& param0, const param_struct& param1,
		     const param_struct& param2, const param_struct& param3) :
  n0(param0.n), min0(param0.min), max0(param0.max), d0(param0.dval), 
  n1(param1.n), min1(param1.min), max1(param1.max), d1(param1.dval), 
  n2(param2.n), min2(param2.min), max2(param2.max), d2(param2.dval),
  n3(param3.n), min3(param3.min), max3(param3.max), d3(param3.dval) {

  axisnames[0] = param0.name;
  axisnames[1] = param1.name;
  axisnames[2] = param2.name;
  axisnames[3] = param3.name;
  spec[0] = param0.param_spec;
  spec[1] = param1.param_spec;
  spec[2] = param2.param_spec;
  spec[3] = param3.param_spec;

  data = new double***[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double**[n1];
    for (unsigned int j = 0; j < n1; ++j) {
      data[i][j] = new double*[n2];
      for (unsigned int k = 0; k < n2; ++k) {
	data[i][j][k] = new double[n3];
      }
    }
  }
}
   

cosgrid4D::~cosgrid4D() {
  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k)
	  delete[] data[i][j][k];
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }
}

void cosgrid4D::free() {
  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k)
	  delete[] data[i][j][k];
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }
  n0 = n1 = n2 = n3 = 0;
  min0 = max0 = d0 = min1 = max1 = d1 = min2 = max2 = d2 = 0.0;
  min3 = max3 = d3 = 0.0;
  data = 0;
}

/*!
  \param[in] axis Which axis to return
  \param[in] i Position along axis
*/
double cosgrid4D::getAxisVal( unsigned int axis, unsigned int i ) const {
  return getAxisVal( axis, static_cast<double>(i) );
}

/*!
  \param[in] axis Which axis to return
  \param[in] x Position along axis
*/
double cosgrid4D::getAxisVal( unsigned int axis, double x ) const {
  switch (axis) {
  case 0: return min0 + x * d0;
  case 1: return min1 + x * d1;
  case 2: return min2 + x * d2;
  case 3: return min3 + x * d3;
  default :
    stringstream errstrng("");
    errstrng << "No axis found with index " << axis;
    throw CosFitterExcept("cosgrid4D","getAxisVal",errstrng.str(),1);
  }
}

/*! Returns index of axis with the specified name, throwing a CosFitterExcpt
 if no such axis is found.
 \param[in] axislabel Name of axis to find
 \returns Index of matching axis
*/
unsigned int cosgrid4D::getAxisIndex(const std::string& axislabel) const {
  //Returns index of axis with name axislabel
  if (axisnames[0] == axislabel) return 0;
  if (axisnames[1] == axislabel) return 1;
  if (axisnames[2] == axislabel) return 2;
  if (axisnames[3] == axislabel) return 3;
  stringstream errstrng("");
  errstrng << "No axis found with label " << axislabel;
  throw CosFitterExcept("cosgrid4D","getAxisIndex",errstrng.str(),1);
}

/*!
  \param[out] i0 Index of maximum along axis 0
  \param[out] i1 Index of maximum along axis 1
  \param[out] i2 Index of maximum along axis 2
  \param[out] i3 Index of maximum along axis 3
  \returns The maximum value
*/
double cosgrid4D::getMaximum(unsigned int &i0, unsigned int& i1, 
			     unsigned int &i2, unsigned int &i3) const {
  double max;
  i0 = i1 = i2 = i3 = 0;
  max = data[0][0][0][0];
  
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h)
	  if (data[i][j][k][h] > max) {
	    max = data[i][j][k][h];
	    i0 = i; i1 = j; i2 = k; i3 = h; 
	  }
  return max;
}

/*!
  \returns The total of all of the elements in the grid
*/
double cosgrid4D::getTotal() const {
  double total;
  total = 0.0;
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h)
	  total += data[i][j][k][h];
  return total;
}

/*!
  \returns The previous total of all of the elements in the grid prior
     to normalization
*/
double cosgrid4D::normalize() {
  double total;

  total = 0.0;
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h)
	  total += data[i][j][k][h];

  if (total == 0.0) {
    string errstrng = "Total value equal to 0 -- can't normalize";
    throw CosFitterExcept("cosgrid4D","normalize",errstrng,1);
  }

  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h)
	  data[i][j][k][h] /= total;

  return total;
}

/*!
  Serializes the contents of the data array and the axis labels to a file,
  either in binary or human readable format.
  \param[in] filename the file to write the data to
  \param[in] binary Write as binary instead of human readable
    text (def: false).
*/
int cosgrid4D::writeFile(const std::string& filename, bool binary ) const {
  if (binary) return writeBinaryFile(filename);
  return writeTextFile(filename);
}

/*!
  Reads in a data file as written by writeFile.  The user must specify
  if it is binary or text -- this code is not smart enough to figure it
  out for itself.
  \param[in] filename the file to write the data to
  \param[in] binary Read as binary instead of human readable
    text (def: false).
*/
int cosgrid4D::readFile(const std::string& filename, bool binary ) {
  if (binary) return readBinaryFile(filename);
  return readTextFile(filename);
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a human readable format.  The format of the output file is defined as:
  - The first line is Naxes: 4
  - The second line is axisname0, n0, min0, max0, d0
  - The third line is axisname1, n1, min1, max1, d1
  - The fourth line is axisname2, n2, min2, max2, d2
  - The fifth line is axisname3, n3, min3, max3, d3
  - This is followed by n0*n1*n2*n3 lines containing the data values in c++ 
     style order (i.e., the axis3 index changes the most rapidly)
  \param[in] filename The file to write the data to
*/
int cosgrid4D::writeTextFile(const std::string& filename) const {
  FILE *fp;
  fp = fopen(filename.c_str(),"w");
  if (!fp) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid4D","writeTextFile",errstrng.str(),1);
  }

  fprintf(fp,"Naxes: 4\n");
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[0].c_str(),
	  n0,min0,max0,d0);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[1].c_str(),
	  n1,min1,max1,d1);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[2].c_str(),
	  n2,min2,max2,d2);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[3].c_str(),
	  n3,min3,max3,d3);
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	for (unsigned int h = 0; h < n3; ++h) {
	  fprintf(fp,"%-15.8e\n",data[i][j][k][h]);
	}
      }
    }
  }

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid4D::writeTextFile
 \param[in] filename The file to read from
 */
int cosgrid4D::readTextFile(const std::string& filename) {

  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k)
	  delete[] data[i][j][k];
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }
  
  ifstream fl( filename.c_str() );
  if (!fl) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid4D","readFile",errstrng.str(),1);
  }

  string line;
  vector<string> words;
  stringstream str("");
  
  unsigned int naxes;
  getline(fl,line);
  utility::stringwords(line,words);
  naxes = atoi( words[1].c_str() );
  if ( naxes != 4 ) {
    fl.close();
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid4D","readFile",errstrng.str(),2);
  }

  //Read first line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[0] = words[0];
  str.str(words[1]); str.clear(); str >> n0;
  str.str(words[2]); str.clear(); str >> min0;
  str.str(words[3]); str.clear(); str >> max0;
  str.str(words[4]); str.clear(); str >> d0;

  //Read second line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[1] = words[0];
  str.str(words[1]); str.clear(); str >> n1;
  str.str(words[2]); str.clear(); str >> min1;
  str.str(words[3]); str.clear(); str >> max1;
  str.str(words[4]); str.clear(); str >> d1;

  //And third
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[2] = words[0];
  str.str(words[1]); str.clear(); str >> n2;
  str.str(words[2]); str.clear(); str >> min2;
  str.str(words[3]); str.clear(); str >> max2;
  str.str(words[4]); str.clear(); str >> d2;

  //And four
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[3] = words[0];
  str.str(words[1]); str.clear(); str >> n3;
  str.str(words[2]); str.clear(); str >> min3;
  str.str(words[3]); str.clear(); str >> max3;
  str.str(words[4]); str.clear(); str >> d3;

  data = new double***[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double**[n1];
    for (unsigned int j = 0; j < n1; ++j) {
      data[i][j] = new double*[n2];
      for (unsigned int k = 0; k < n2; ++k) {
	data[i][j][k] = new double[n3];
      }
    }
  }
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	for (unsigned int h = 0; h < n3; ++h) {
	  fl >> data[i][j][k][h];
	  //There are underflow issues with some versions of gcc
	  // so we have to deliberately ignore failbit
	  if (fl.fail()) fl.clear(fl.rdstate() & ~std::ios::failbit);
	}
      }
    }
  }
  
  fl.close();
  
  return 0;
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a non human readable, binary, format.  The format of the output file 
  is defined as:
  - an integer specifying the number of axes.
  - an integer specifying the number of characters in the zero terminated
     axisname
  - the axis name
  - an integer specifying n0
  - doubles specifying min0, max0, d0
  - a repeat of items 2-5 for the second, third, and fourth
  - then n0*n1*n2*n3 doubles storing the data in c order
  Note that this will not be machine portable because of endianess and
  the sizes of various outputs.  Unlike the version which writes in text,
  here we don't care if the axisnames have spaces in them.
  \param[in] filename The file to write the data to
*/
int cosgrid4D::writeBinaryFile(const std::string& filename) const {

  FILE *fp;
  fp = fopen(filename.c_str(),"wb");
  if (!fp) {
    stringstream errstrng("");
    errstrng <<  "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid4D","writeBinaryFile",errstrng.str(),1);
  }
  
  //Write number of axes
  unsigned int naxes = 4;
  fwrite( &naxes, sizeof(unsigned int), 1, fp );

  //Write axis 0 title
  unsigned int axisnamelength = axisnames[0].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[0].c_str(), fp );
  
  //Then extended axis information
  fwrite( &n0, sizeof(unsigned int), 1, fp );
  fwrite( &min0, sizeof(double), 1, fp );
  fwrite( &max0, sizeof(double), 1, fp );
  fwrite( &d0, sizeof(double), 1, fp );

  //Then axis 1
  axisnamelength = axisnames[1].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[1].c_str(), fp );
  fwrite( &n1, sizeof(unsigned int), 1, fp );
  fwrite( &min1, sizeof(double), 1, fp );
  fwrite( &max1, sizeof(double), 1, fp );
  fwrite( &d1, sizeof(double), 1, fp );

  //Then axis 2
  axisnamelength = axisnames[2].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[2].c_str(), fp );
  fwrite( &n2, sizeof(unsigned int), 1, fp );
  fwrite( &min2, sizeof(double), 1, fp );
  fwrite( &max2, sizeof(double), 1, fp );
  fwrite( &d2, sizeof(double), 1, fp );

  //Then axis 3
  axisnamelength = axisnames[3].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[3].c_str(), fp );
  fwrite( &n3, sizeof(unsigned int), 1, fp );
  fwrite( &min3, sizeof(double), 1, fp );
  fwrite( &max3, sizeof(double), 1, fp );
  fwrite( &d3, sizeof(double), 1, fp );

  //And the data
  //It might be tempting to just try to write n0*n1*n2*n3 entries here,
  // but we have no actualy guarantee that different rows will
  // stored contiguously in memory.  In other words, we have to
  // loop explicitly
  for (unsigned int i = 0; i < n0; ++i ) 
    for ( unsigned int j = 0; j < n1; ++j )
      for ( unsigned int k = 0; k < n2; ++k )
	fwrite( data[i][j][k], sizeof(double), n3, fp );

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid4D::writeBinaryFile
 \param[in] filename The file to read from
 */
int cosgrid4D::readBinaryFile(const std::string& filename) {

  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k)
	  delete[] data[i][j][k];
	delete [] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }

  FILE *fp;
  fp = fopen( filename.c_str(), "rb" );
  if (!fp) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid4D","readBinaryFile",errstrng.str(),1);
  }

  //Read number of axes
  unsigned int naxes;
  fread( &naxes, sizeof(unsigned int), 1, fp );
  if ( naxes != 4 ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid4D","readBinaryFile",errstrng.str(),2);
  }
      
  //Read axis 1 title
  unsigned int status;
  status = readBinaryAxisName( fp, axisnames[0] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid4D","readBinaryFile",errstrng.str(),4);
  }
  
  fread( &n0, sizeof(unsigned int), 1, fp );
  fread( &min0, sizeof(double), 1, fp );
  fread( &max0, sizeof(double), 1, fp );
  fread( &d0, sizeof(double), 1, fp );

  //Read axis 2
  status = readBinaryAxisName( fp, axisnames[1] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid4D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n1, sizeof(unsigned int), 1, fp );
  fread( &min1, sizeof(double), 1, fp );
  fread( &max1, sizeof(double), 1, fp );
  fread( &d1, sizeof(double), 1, fp );

  //Read axis 3
  status = readBinaryAxisName( fp, axisnames[2] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid4D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n2, sizeof(unsigned int), 1, fp );
  fread( &min2, sizeof(double), 1, fp );
  fread( &max2, sizeof(double), 1, fp );
  fread( &d2, sizeof(double), 1, fp );

  //Read axis 4
  status = readBinaryAxisName( fp, axisnames[3] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid4D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n3, sizeof(unsigned int), 1, fp );
  fread( &min3, sizeof(double), 1, fp );
  fread( &max3, sizeof(double), 1, fp );
  fread( &d3, sizeof(double), 1, fp );

  if ( (n0 <= 0) || ( n1 <= 0 ) || ( n2 <= 0 ) || ( n3 <= 0 ) ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "0 axis length in " << filename << endl;
    throw CosFitterExcept("cosgrid4D","readBinaryFile",errstrng.str(),8);
  }
  //Now read the data.  Do the allocation in a seperate step to
  // increase the chances things will be contiguous
  data = new double***[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double**[n1];
    for ( unsigned int j = 0; j < n1; ++j ) {
      data[i][j] = new double*[n2];
      for (unsigned int k = 0; k < n2; ++k ) 
	data[i][j][k] = new double[n3];
    }
  }
  for (unsigned int i = 0; i < n0; ++i ) 
    for (unsigned int j = 0; j < n1; ++j)
      for (unsigned int k = 0; k < n2; ++k)
	fread( data[i][j][k], sizeof(double), n3, fp );

  fclose(fp);

  return 0;
}

cosgrid4D& cosgrid4D::operator=(double val) {
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int m = 0; m < n3; ++m)
	  data[i][j][k][m] = val;
  return *this;
}

cosgrid4D& cosgrid4D::operator=(const cosgrid4D& input) {
  //Self copy protection
  if (this == &input) return *this;

  axisnames[0] = input.axisnames[0];
  axisnames[1] = input.axisnames[1];
  axisnames[2] = input.axisnames[2];
  axisnames[3] = input.axisnames[3];
  spec[0] = input.spec[0];
  spec[1] = input.spec[1];
  spec[2] = input.spec[2];
  spec[3] = input.spec[3];

  if (data == 0) {
    data = new double***[input.n0];
    for (unsigned int i = 0; i < input.n0; ++i) {
      data[i] = new double**[input.n1];
      for (unsigned int j = 0; j < input.n1; ++j) {
	data[i][j] = new double*[input.n2];
	for (unsigned int k = 0; k < input.n2; ++k) 
	  data[i][j][k] = new double[input.n3];
      }
    }
  } else if ( (input.n0 != n0) || (input.n1 != n1) ||
	      (input.n2 != n2) || (input.n3 != n3) ) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k)
	  delete[] data[i][j][k];
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
    data = new double***[input.n0];
    for (unsigned int i = 0; i < input.n0; ++i) {
      data[i] = new double**[input.n1];
      for (unsigned int j = 0; j < input.n1; ++j) {
	data[i][j] = new double*[input.n2];
	for (unsigned int k = 0; k < input.n2; ++k) 
	  data[i][j][k] = new double[input.n3];
      }
    }
  }

  for (unsigned int i = 0; i < input.n0; ++i) 
    for (unsigned int j = 0; j < input.n1; ++j) 
      for (unsigned int k = 0; k < input.n2; ++k) 
	for (unsigned int h = 0; h < input.n3; ++h) 
	  data[i][j][k][h] = input.data[i][j][k][h];
  
  n0 = input.n0; min0 = input.min0; max0 = input.max0; d0 = input.d0;
  n1 = input.n1; min1 = input.min1; max1 = input.max1; d1 = input.d1;
  n2 = input.n2; min2 = input.min2; max2 = input.max2; d2 = input.d2;
  n3 = input.n3; min3 = input.min3; max3 = input.max3; d3 = input.d3;

  return *this;
}

cosgrid4D& cosgrid4D::operator/=(double divisor) {
  for (unsigned int i = 0; i < n0; ++i) 
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k) 
	for (unsigned int m = 0; m < n3; ++m) 
	    data[i][j][k][m] /= divisor;
  return *this;
}

cosgrid3D cosgrid4D::collapseAlongZero() const {

  cosgrid3D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3]);
  
  double ***dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int k = 0; k < n2; ++k) {
      for (unsigned int h = 0; h < n3; ++h) {
	val = 0.0;
	for (unsigned int i = 0; i < n0; ++i) val += data[i][j][k][h];
	dataptr[j][k][h] = val;
      }
    }
  }
  
  return result;
}

cosgrid3D cosgrid4D::collapseAlongOne() const {
  
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double ***dataptr = result.getData();
  double val;
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int k = 0; k < n2; ++k) {
      for (unsigned int h = 0; h < n3; ++h) {
	val = 0.0;
	for (unsigned int j = 0; j < n1; ++j) val += data[i][j][k][h];
	dataptr[i][k][h] = val;
      }
    }
  }
  
  return result;
}

cosgrid3D cosgrid4D::collapseAlongTwo() const {
  
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double ***dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int h = 0; h < n3; ++h) {
	val = 0.0;
	for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k][h];
	dataptr[i][j][h] = val;
      }
    }
  }

  return result;
}
 
cosgrid3D cosgrid4D::collapseAlongThree() const {
  
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2]);
  
  double ***dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	val = 0.0;
	for (unsigned int h = 0; h < n3; ++h) val += data[i][j][k][h];
	dataptr[i][j][k] = val;
      }
    }
  }
  
  return result;
}
 
cosgrid2D cosgrid4D::collapseAlongZeroOne() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int k = 0; k < n2; ++k) {
    for (unsigned int h = 0; h < n3; ++h) {
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i)
	for (unsigned int j = 0; j < n1; ++j) val += data[i][j][k][h];
      dataptr[k][h] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid4D::collapseAlongZeroTwo() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n3,min3,max3,d3,axisnames[3],spec[2]);
  
  double **dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int h = 0; h < n3; ++h) {
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i)
	for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k][h];
      dataptr[j][h] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid4D::collapseAlongZeroThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int k = 0; k < n2; ++k) {
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i)
	for (unsigned int h = 0; h < n3; ++h) val += data[i][j][k][h];
      dataptr[j][k] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid4D::collapseAlongOneTwo() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int h = 0; h < n3; ++h) {
      val = 0.0;
      for (unsigned int j = 0; j < n1; ++j)
	for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k][h];
      dataptr[i][h] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid4D::collapseAlongOneThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n2,min2,max2,d2,axisnames[2],spec[1]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int k = 0; k < n2; ++k) {
      val = 0.0;
      for (unsigned int j = 0; j < n1; ++j)
	for (unsigned int h = 0; h < n3; ++h) val += data[i][j][k][h];
      dataptr[i][k] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid4D::collapseAlongTwoThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      val = 0.0;
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h) val += data[i][j][k][h];
      dataptr[i][j] = val;
    }
  }

  return result;
}


cosgrid1D cosgrid4D::collapseAlongZeroOneTwo() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid1D result(n3,min3,max3,d3,axisnames[3],spec[3]);

  double *dataptr = result.getData();
  double val;
  for (unsigned int h = 0; h < n3; ++h) {
    val = 0.0;
    for (unsigned int i = 0; i < n0; ++i)
      for (unsigned int j = 0; j < n1; ++j) 
	for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k][h];
    dataptr[h] = val;
  }
  return result;
}

cosgrid1D cosgrid4D::collapseAlongZeroOneThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid1D result(n2,min2,max2,d2,axisnames[2],spec[2]);

  double *dataptr = result.getData();
  double val;
  for (unsigned int k = 0; k < n2; ++k) {
    val = 0.0;
    for (unsigned int i = 0; i < n0; ++i)
      for (unsigned int j = 0; j < n1; ++j) 
	for (unsigned int h = 0; h < n3; ++h) val += data[i][j][k][h];
    dataptr[k] = val;
  }
  return result;
}

cosgrid1D cosgrid4D::collapseAlongZeroTwoThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid1D result(n1,min1,max1,d1,axisnames[1],spec[1]);

  double *dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    val = 0.0;
    for (unsigned int i = 0; i < n0; ++i)
      for (unsigned int k = 0; k < n2; ++k) 
	for (unsigned int h = 0; h < n3; ++h) val += data[i][j][k][h];
    dataptr[j] = val;
  }
  return result;
}

cosgrid1D cosgrid4D::collapseAlongOneTwoThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid1D result(n0,min0,max0,d0,axisnames[0],spec[0]);

  double *dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    val = 0.0;
    for (unsigned int j = 0; j < n1; ++j)
      for (unsigned int k = 0; k < n2; ++k) 
	for (unsigned int h = 0; h < n3; ++h) val += data[i][j][k][h];
    dataptr[i] = val;
  }
  return result;
}

/*!
  Collapse the grid along the specified axis by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisnum Index of axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid3D containing the collapsed data array
*/
cosgrid3D cosgrid4D::collapseAlongAxis(unsigned int axisnum) const {
  if (axisnum == 0) return collapseAlongZero();
  if (axisnum == 1) return collapseAlongOne();
  if (axisnum == 2) return collapseAlongTwo();
  if (axisnum == 3) return collapseAlongThree();
  stringstream errstrng("");
  errstrng << "Invalid axis number " << axisnum;
  throw CosFitterExcept("cosgrid4D","collapseAlongAxis",errstrng.str(),1);
}
 
/*!
  Collapse the grid along the specified axis by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisname Nameof axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid3D containing the collapsed data array
*/
cosgrid3D cosgrid4D::collapseAlongAxis(const std::string& axisname) const {
  return collapseAlongAxis( getAxisIndex(axisname) );
}

/*!
  Collapse the grid along the specified axes by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisnum1 Index of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisnum2 Index of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid2D containing the collapsed data array
*/
cosgrid2D cosgrid4D::collapseAlongTwoAxes(unsigned int axisnum1,
					  unsigned int axisnum2) const {
  if (axisnum1 == 0) {
    if (axisnum2 == 1) return collapseAlongZeroOne();
    if (axisnum2 == 2) return collapseAlongZeroTwo();
    if (axisnum2 == 3) return collapseAlongZeroThree();
  } else if (axisnum1 == 1) {
    if (axisnum2 == 0) return collapseAlongZeroOne();
    if (axisnum2 == 2) return collapseAlongOneTwo();
    if (axisnum2 == 3) return collapseAlongOneThree();
  } else if (axisnum1 == 2) {
    if (axisnum2 == 0) return collapseAlongZeroTwo();
    if (axisnum2 == 1) return collapseAlongOneTwo();
    if (axisnum2 == 3) return collapseAlongTwoThree();
  } else if (axisnum1 == 3) {
    if (axisnum2 == 0) return collapseAlongZeroThree();
    if (axisnum2 == 1) return collapseAlongOneThree();
    if (axisnum2 == 2) return collapseAlongTwoThree();
  }
  stringstream errstr("");
  errstr << "Can't collapse along axes: " << axisnum1 << " " << axisnum2;
  throw CosFitterExcept("cosgrid4D","collapseAlongTwoAxes",
			   errstr.str(),1);

}

/*!
  Collapse the grid along the specified axes by summing.
  \param[in] axisname1 Name of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisname2 Name of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid2D containing the collapsed data array
*/
cosgrid2D cosgrid4D::collapseAlongTwoAxes(const std::string& axisname1,
					  const std::string& axisname2) const {
  return collapseAlongTwoAxes( getAxisIndex(axisname1),
			       getAxisIndex(axisname2) );

}

/*!
  Collapse the grid along the specified axes by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisnum1 Index of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisnum2 Index of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisnum3 Index of 3rd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid1D containing the collapsed data array
*/
cosgrid1D cosgrid4D::collapseAlongThreeAxes(unsigned int axisnum1,
					    unsigned int axisnum2,
					    unsigned int axisnum3) const {
  if (axisnum1 == 0) {
    if (axisnum2 == 1) {
      if (axisnum3 == 2) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongZeroOneThree();
    }
    if (axisnum2 == 2) {
      if (axisnum3 == 1) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongZeroTwoThree();
    }
    if (axisnum2 == 3) {
      if (axisnum3 == 1) return collapseAlongZeroOneThree();
      if (axisnum3 == 2) return collapseAlongZeroTwoThree();
    }
  } else if (axisnum1 == 1) {
    if (axisnum2 == 0) {
      if (axisnum3 == 2) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongZeroOneThree();
    }
    if (axisnum2 == 2) {
      if (axisnum3 == 0) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongOneTwoThree();
    }
    if (axisnum2 == 3) {
      if (axisnum3 == 0) return collapseAlongZeroOneThree();
      if (axisnum3 == 2) return collapseAlongOneTwoThree();
    }
  } else if (axisnum1 == 2) {
    if (axisnum2 == 0) {
      if (axisnum3 == 1) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongZeroTwoThree();
    }
    if (axisnum2 == 1) {
      if (axisnum3 == 0) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongOneTwoThree();
    }
    if (axisnum2 == 3) {
      if (axisnum3 == 0) return collapseAlongZeroTwoThree();
      if (axisnum3 == 1) return collapseAlongOneTwoThree();
    }
  } else if (axisnum1 == 3) {
    if (axisnum2 == 0) {
      if (axisnum3 == 1) return collapseAlongZeroOneThree();
      if (axisnum3 == 2) return collapseAlongZeroTwoThree();
    }
    if (axisnum2 == 1) {
      if (axisnum3 == 0) return collapseAlongZeroOneThree();
      if (axisnum3 == 2) return collapseAlongOneTwoThree();
    }
    if (axisnum2 == 2) {
      if (axisnum3 == 0) return collapseAlongZeroTwoThree();
      if (axisnum3 == 1) return collapseAlongOneTwoThree();
    }
  }
  stringstream errstr("");
  errstr << "Can't collapse along axes: " << axisnum1 << " " << axisnum2 <<
    " " << axisnum3;
  throw CosFitterExcept("cosgrid4D","collapseAlongThreeAxes",
			   errstr.str(),1);

}

/*!
  Collapse the grid along the specified axes by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisname1 Name of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisname2 Name of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisname3 Name of 3rd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid1D containing the collapsed data array
*/
cosgrid1D cosgrid4D::collapseAlongThreeAxes(const std::string& axisname1,
					    const std::string& axisname2,
					    const std::string& axisname3) const {
  return collapseAlongThreeAxes( getAxisIndex(axisname1),
				 getAxisIndex(axisname2),
				 getAxisIndex(axisname3));

}

//////////////////////////////////
//           cosgrid5D          //
//////////////////////////////////

cosgrid5D::cosgrid5D() {
  data = 0;
  n0 = n1 = n2 = n3 = n4 = 0;
  min0 = max0 = d0 = min1 = max1 = d1 = min2 = max2 = d2 = 0.0;
  min3 = max3 = d3 = min4 = max4 = d4 = 0.0;
}

cosgrid5D::cosgrid5D(const cosgrid5D& input) {
  n0 = input.n0; min0 = input.min0; max0 = input.max0; d0 = input.d0;
  n1 = input.n1; min1 = input.min1; max1 = input.max1; d1 = input.d1;
  n2 = input.n2; min2 = input.min2; max2 = input.max2; d2 = input.d2;
  n3 = input.n3; min3 = input.min3; max3 = input.max3; d3 = input.d3;
  n4 = input.n4; min4 = input.min4; max4 = input.max4; d4 = input.d4;
  axisnames[0] = input.axisnames[0];
  axisnames[1] = input.axisnames[1];
  axisnames[2] = input.axisnames[2];
  axisnames[3] = input.axisnames[3];
  axisnames[4] = input.axisnames[4];
  spec[0] = input.spec[0];
  spec[1] = input.spec[1];
  spec[2] = input.spec[2];
  spec[3] = input.spec[3];
  spec[4] = input.spec[4];
  if (input.data != 0) {
    data = new double****[ n0 ];
    for (unsigned int i = 0; i < n0; ++i) {
      data[i] = new double***[n1];
      for (unsigned int j = 0; j < n1; ++j) {
	data[i][j] = new double**[n2];
	for (unsigned int k = 0; k < n2; ++k) { 
	  data[i][j][k] = new double*[n3];
	  for (unsigned int m = 0; m < n3; ++m) {
	    data[i][j][k][m] = new double[n4];
	    for (unsigned int n = 0; n < n4; ++n) 
	      data[i][j][k][m][n] = input.data[i][j][k][m][n];
	  }
	}
      }
    }
  }
}

/*!
  \param[in] n0 Number of data entries along axis 0
  \param[in] min0 Minimum axis value along axis 0
  \param[in] max0 Maximum axis value along axis 0
  \param[in] d0 Step size for axis along axis 0
  \param[in] label0 Name of axis 0 parameter.  Should not have spaces in it.
  \param[in] spec0 Param spec for parameter 0
  \param[in] n1 Number of data entries along axis 1
  \param[in] min1 Minimum axis value along axis 1
  \param[in] max1 Maximum axis value along axis 1
  \param[in] d1 Step size for axis along axis 1
  \param[in] label1 Name of axis 1 parameter.  Should not have spaces in it.
  \param[in] spec1 Param spec for parameter 1
  \param[in] n2 Number of data entries along axis 2
  \param[in] min2 Minimum axis value along axis 2
  \param[in] max2 Maximum axis value along axis 2
  \param[in] d2 Step size for axis along axis 2
  \param[in] label2 Name of axis 2 parameter.  Should not have spaces in it.
  \param[in] spec2 Param spec for parameter 2
  \param[in] n3 Number of data entries along axis 3
  \param[in] min3 Minimum axis value along axis 3
  \param[in] max3 Maximum axis value along axis 3
  \param[in] d3 Step size for axis along axis 3
  \param[in] label3 Name of axis 3 parameter.  Should not have spaces in it.
  \param[in] spec3 Param spec for parameter 3
  \param[in] n4 Number of data entries along axis 4
  \param[in] min4 Minimum axis value along axis 4
  \param[in] max4 Maximum axis value along axis 4
  \param[in] d4 Step size for axis along axis 4
  \param[in] label4 Name of axis 4 parameter.  Should not have spaces in it.
  \param[in] spec4 Param spec for parameter 4
*/
cosgrid5D::cosgrid5D(unsigned int n0, double min0, double max0, double d0,
		     std::string label0, param_tags::paramspec spec0,
		     unsigned int n1, double min1, double max1, double d1, 
		     std::string label1, param_tags::paramspec spec1,
		     unsigned int n2, double min2, double max2, double d2, 
		     std::string label2, param_tags::paramspec spec2, 
		     unsigned int n3, double min3, double max3, double d3, 
		     std::string label3, param_tags::paramspec spec3,
		     unsigned int n4, double min4, double max4, double d4, 
		     std::string label4, param_tags::paramspec spec4) :
  n0(n0), min0(min0), max0(max0), d0(d0), n1(n1), 
  min1(min1), max1(max1), d1(d1), n2(n2), min2(min2), max2(max2), d2(d2),
  n3(n3), min3(min3), max3(max3), d3(d3), n4(n4), min4(min4), max4(max4), 
  d4(d4) {

  axisnames[0] = label0;
  axisnames[1] = label1;
  axisnames[2] = label2;
  axisnames[3] = label3;
  axisnames[4] = label4;
  spec[0] = spec0;
  spec[1] = spec1;
  spec[2] = spec2;
  spec[3] = spec3;
  spec[4] = spec4;

  data = new double****[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double***[n1];
    for (unsigned int j = 0; j < n1; ++j) {
      data[i][j] = new double**[n2];
      for (unsigned int k = 0; k < n2; ++k) {
	data[i][j][k] = new double*[n3];
	for (unsigned int m = 0; m < n3; ++m)
	  data[i][j][k][m] = new double[n4];
      }
    }
  }
}

/*!
  \param[in] param0 Information about parameter describing axis 0
  \param[in] param1 Information about parameter describing axis 1
  \param[in] param2 Information about parameter describing axis 2
  \param[in] param3 Information about parameter describing axis 3
  \param[in] param4 Information about parameter describing axis 4
*/
cosgrid5D::cosgrid5D(const param_struct& param0, const param_struct& param1,
		     const param_struct& param2, const param_struct& param3,
		     const param_struct& param4) :
  n0(param0.n), min0(param0.min), max0(param0.max), d0(param0.dval), 
  n1(param1.n), min1(param1.min), max1(param1.max), d1(param1.dval), 
  n2(param2.n), min2(param2.min), max2(param2.max), d2(param2.dval),
  n3(param3.n), min3(param3.min), max3(param3.max), d3(param3.dval),
  n4(param4.n), min4(param4.min), max4(param4.max), d4(param4.dval) {

  axisnames[0] = param0.name;
  axisnames[1] = param1.name;
  axisnames[2] = param2.name;
  axisnames[3] = param3.name;
  axisnames[4] = param4.name;
  spec[0] = param0.param_spec;
  spec[1] = param1.param_spec;
  spec[2] = param2.param_spec;
  spec[3] = param3.param_spec;
  spec[4] = param4.param_spec;

  data = new double****[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double***[n1];
    for (unsigned int j = 0; j < n1; ++j) {
      data[i][j] = new double**[n2];
      for (unsigned int k = 0; k < n2; ++k) {
	data[i][j][k] = new double*[n3];
	for (unsigned int m = 0; m < n3; ++m)
	  data[i][j][k][m] = new double[n4];
      }
    }
  }
}

cosgrid5D::~cosgrid5D() {
  if (data != 0) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k) {
	  for (unsigned int m = 0; m < n3; ++m)
	    delete[] data[i][j][k][m];
	  delete[] data[i][j][k];
	}
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }
}

void cosgrid5D::free() {
  if (data != 0) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k) {
	  for (unsigned int m = 0; m < n3; ++m)
	    delete[] data[i][j][k][m];
	  delete[] data[i][j][k];
	}
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }
  n0 = n1 = n2 = n3 = n4 = 0;
  min0 = max0 = d0 = min1 = max1 = d1 = min2 = max2 = d2 = 0.0;
  min3 = max3 = d3 = min4 = max4 = d4 = 0.0;
  data = 0;
}

/*!
  \param[in] axis Which axis to return
  \param[in] i Position along axis
*/
double cosgrid5D::getAxisVal( unsigned int axis, unsigned int i ) const {
  return getAxisVal( axis, static_cast<double>(i) );
}

/*!
  \param[in] axis Which axis to return
  \param[in] x Position along axis
*/
double cosgrid5D::getAxisVal( unsigned int axis, double x ) const {
  switch (axis) {
  case 0: return min0 + x * d0;
  case 1: return min1 + x * d1;
  case 2: return min2 + x * d2;
  case 3: return min3 + x * d3;
  case 4: return min4 + x * d4;
  default :
    stringstream errstrng("");
    errstrng << "No axis found with index " << axis;
    throw CosFitterExcept("cosgrid5D","getAxisVal",errstrng.str(),1);
  }
}


/*!
  \param[in] axislabel Name of axis
  \returns Index (0,1,2,3,4) of matching axis
*/
unsigned int cosgrid5D::getAxisIndex(const std::string& axislabel) const {
  //Returns index of axis with name axislabel
  if (axisnames[0] == axislabel) return 0;
  if (axisnames[1] == axislabel) return 1;
  if (axisnames[2] == axislabel) return 2;
  if (axisnames[3] == axislabel) return 3;
  if (axisnames[4] == axislabel) return 4;
  stringstream errstrng("");
  errstrng << "No axis found with label " << axislabel;
  throw CosFitterExcept("cosgrid5D","getAxisIndex",errstrng.str(),1);
}

/*!
  \param[out] i0 Index of maximum along axis 0
  \param[out] i1 Index of maximum along axis 1
  \param[out] i2 Index of maximum along axis 2
  \param[out] i3 Index of maximum along axis 3
  \param[out] i4 Index of maximum along axis 4
  \returns The maximum value
*/
double cosgrid5D::getMaximum(unsigned int &i0, unsigned int& i1, 
			     unsigned int &i2, unsigned int &i3, 
			     unsigned int &i4) const {
  double max;
  i0 = i1 = i2 = i3 = i4 = 0;
  max = data[0][0][0][0][0];
  
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h)
	  for (unsigned int m = 0; m < n4; ++m)
	    if (data[i][j][k][h][m] > max) {
	      max = data[i][j][k][h][m];
	      i0 = i; 
	      i1 = j; 
	      i2 = k; 
	      i3 = h; 
	      i4 = m;
	    }
  return max;
}

/*!
  \returns The total of all of the elements in the grid
*/
double cosgrid5D::getTotal() const {
  double total;
  total = 0.0;
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h)
	  for (unsigned int m = 0; m < n4; ++m)
	    total += data[i][j][k][h][m];
  return total;
}

/*!
  \returns The previous total of all of the elements in the grid
            prior to normalization
*/
double cosgrid5D::normalize() {
  double total;

  total = 0.0;
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h)
	  for (unsigned int m = 0; m < n4; ++m)
	    total += data[i][j][k][h][m];

  if (total == 0.0) {
    std::string errstrng = "Total value equal to 0 -- can't normalize";
    throw CosFitterExcept("cosgrid5D","normalize",errstrng,1);
  }

  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int h = 0; h < n3; ++h)
	  for (unsigned int m = 0; m < n4; ++m)
	    data[i][j][k][h][m] /= total;

  return total;
}
 
/*!
  Serializes the contents of the data array and the axis labels to a file,
  either in binary or human readable format.
  \param[in] filename the file to write the data to
  \param[in] binary Write as binary instead of human readable
    text (def: false).
*/
int cosgrid5D::writeFile(const std::string& filename, bool binary ) const {
  if (binary) return writeBinaryFile(filename);
  return writeTextFile(filename);
}

/*!
  Reads in a data file as written by writeFile.  The user must specify
  if it is binary or text -- this code is not smart enough to figure it
  out for itself.
  \param[in] filename the file to write the data to
  \param[in] binary Read as binary instead of human readable
    text (def: false).
*/
int cosgrid5D::readFile(const std::string& filename, bool binary ) {
  if (binary) return readBinaryFile(filename);
  return readTextFile(filename);
}	      

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a human readable format.  The format of the output file is defined as:
  - The first line is Naxes: 5
  - The second line is axisname0, n0, min0, max0, d0
  - The third line is axisname1, n1, min1, max1, d1
  - The fourth line is axisname2, n2, min2, max2, d2
  - The fifth line is axisname3, n3, min3, max3, d3
  - The sixth line is axisname4, n4, min4, max4, d4
  - This is followed by n0*n1*n2*n3*n4 lines containing the data values in c++ 
     style order (i.e., the axis4 index changes the most rapidly)
  \param[in] filename The file to write the data to
*/
int cosgrid5D::writeTextFile(const std::string& filename) const {
  FILE *fp;
  fp = fopen(filename.c_str(),"w");
  if (!fp) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid5D","writeTextFile",errstrng.str(),1);
  }

  fprintf(fp,"Naxes: 5\n");
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[0].c_str(),
	  n0,min0,max0,d0);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[1].c_str(),
	  n1,min1,max1,d1);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[2].c_str(),
	  n2,min2,max2,d2);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[3].c_str(),
	  n3,min3,max3,d3);
  fprintf(fp,"%s %d %11.8f %11.8f %11.8f\n",axisnames[4].c_str(),
	  n4,min4,max4,d4);
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	for (unsigned int m = 0; m < n3; ++m) {
	  for (unsigned int n = 0; n < n4; ++n)
	    fprintf(fp,"%-15.8e\n",data[i][j][k][m][n]);
	}
      }
    }
  }

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid5D::writeTextFile
 \param[in] filename The file to read from
 */

int cosgrid5D::readTextFile(const std::string& filename) {

  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k) {
	  for (unsigned int m = 0; m < n3; ++m)
	    delete[] data[i][j][k][m];
	  delete[] data[i][j][k];
	}
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }

  ifstream fl( filename.c_str() );
  if (!fl) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid5D","readFile",errstrng.str(),1);
  }

  string line;
  vector<string> words;
  stringstream str("");

  unsigned int naxes;
  getline(fl,line);
  utility::stringwords(line,words);
  naxes = atoi( words[1].c_str() );
  if ( naxes != 5 ) {
    fl.close();
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid1D","readFile",errstrng.str(),2);
  }
  
  //Read first axis line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[0] = words[0];
  str.str(words[1]); str.clear(); str >> n0;
  str.str(words[2]); str.clear(); str >> min0;
  str.str(words[3]); str.clear(); str >> max0;
  str.str(words[4]); str.clear(); str >> d0;

  //Read second line
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[1] = words[0];
  str.str(words[1]); str.clear(); str >> n1;
  str.str(words[2]); str.clear(); str >> min1;
  str.str(words[3]); str.clear(); str >> max1;
  str.str(words[4]); str.clear(); str >> d1;

  //And third
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[2] = words[0];
  str.str(words[1]); str.clear(); str >> n2;
  str.str(words[2]); str.clear(); str >> min2;
  str.str(words[3]); str.clear(); str >> max2;
  str.str(words[4]); str.clear(); str >> d2;

  //And four
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[3] = words[0];
  str.str(words[1]); str.clear(); str >> n3;
  str.str(words[2]); str.clear(); str >> min3;
  str.str(words[3]); str.clear(); str >> max3;
  str.str(words[4]); str.clear(); str >> d3;

  //five
  getline(fl,line);
  utility::stringwords(line,words);
  axisnames[4] = words[0];
  str.str(words[1]); str.clear(); str >> n4;
  str.str(words[2]); str.clear(); str >> min4;
  str.str(words[3]); str.clear(); str >> max4;
  str.str(words[4]); str.clear(); str >> d4;

  data = new double****[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double***[n1];
    for (unsigned int j = 0; j < n1; ++j) {
      data[i][j] = new double**[n2];
      for (unsigned int k = 0; k < n2; ++k) {
	data[i][j][k] = new double*[n3];
	for (unsigned int m = 0; m < n3; ++m)
	  data[i][j][k][m] = new double[n4];
      }
    }
  }
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	for (unsigned int m = 0; m < n3; ++m) {
	  for (unsigned int n = 0; n < n4; ++n) {
	    fl >> data[i][j][k][m][n];
	    //There are underflow issues with some versions of gcc
	    // so we have to deliberately ignore failbit
	    if (fl.fail()) fl.clear(fl.rdstate() & ~std::ios::failbit);
	  }
	}
      }
    }
  }
  
  fl.close();
  
  return 0;
}

/*!
  Serializes the contents of the data array and the axis labels to a file
  in a non human readable, binary, format.  The format of the output file 
  is defined as:
  - an integer specifying the number of axes.
  - an integer specifying the number of characters in the zero terminated
     axisname
  - the axis name
  - an integer specifying n0
  - doubles specifying min0, max0, d0
  - a repeat of items 2-5 for the second-fifth axes
  - then n0*n1*n2*n3*n4 doubles storing the data in c order
  Note that this will not be machine portable because of endianess and
  the sizes of various outputs.  Unlike the version which writes in text,
  here we don't care if the axisnames have spaces in them.
  \param[in] filename The file to write the data to
*/
int cosgrid5D::writeBinaryFile(const std::string& filename) const {

  FILE *fp;
  fp = fopen(filename.c_str(),"wb");
  if (!fp) {
    stringstream errstrng("");
    errstrng <<  "Couldn't open file " << filename << " for writing";
    throw CosFitterExcept("cosgrid5D","writeBinaryFile",errstrng.str(),1);
  }
  
  //Write number of axes
  unsigned int naxes = 5;
  fwrite( &naxes, sizeof(unsigned int), 1, fp );

  //Write axis 0 title
  unsigned int axisnamelength = axisnames[0].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[0].c_str(), fp );
  
  //Then extended axis information
  fwrite( &n0, sizeof(unsigned int), 1, fp );
  fwrite( &min0, sizeof(double), 1, fp );
  fwrite( &max0, sizeof(double), 1, fp );
  fwrite( &d0, sizeof(double), 1, fp );

  //Then axis 1
  axisnamelength = axisnames[1].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[1].c_str(), fp );
  fwrite( &n1, sizeof(unsigned int), 1, fp );
  fwrite( &min1, sizeof(double), 1, fp );
  fwrite( &max1, sizeof(double), 1, fp );
  fwrite( &d1, sizeof(double), 1, fp );

  //Then axis 2
  axisnamelength = axisnames[2].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[2].c_str(), fp );
  fwrite( &n2, sizeof(unsigned int), 1, fp );
  fwrite( &min2, sizeof(double), 1, fp );
  fwrite( &max2, sizeof(double), 1, fp );
  fwrite( &d2, sizeof(double), 1, fp );

  //Then axis 3
  axisnamelength = axisnames[3].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[3].c_str(), fp );
  fwrite( &n3, sizeof(unsigned int), 1, fp );
  fwrite( &min3, sizeof(double), 1, fp );
  fwrite( &max3, sizeof(double), 1, fp );
  fwrite( &d3, sizeof(double), 1, fp );

  //Then axis 4
  axisnamelength = axisnames[4].size();
  fwrite( &axisnamelength, sizeof(unsigned int), 1, fp );
  fputs( axisnames[4].c_str(), fp );
  fwrite( &n4, sizeof(unsigned int), 1, fp );
  fwrite( &min4, sizeof(double), 1, fp );
  fwrite( &max4, sizeof(double), 1, fp );
  fwrite( &d4, sizeof(double), 1, fp );

  //And the data
  //It might be tempting to just try to write n0*n1*n2*n3 entries here,
  // but we have no actualy guarantee that different rows will
  // stored contiguously in memory.  In other words, we have to
  // loop explicitly
  for (unsigned int i = 0; i < n0; ++i ) 
    for ( unsigned int j = 0; j < n1; ++j )
      for ( unsigned int k = 0; k < n2; ++k )
	for (unsigned int m = 0; m < n3; ++m )
	  fwrite( data[i][j][k][m], sizeof(double), n4, fp );

  fclose(fp);

  return 0;
}

/*!
  Reads in the data from a file in the format specified by 
 cosgrid5D::writeBinaryFile
 \param[in] filename The file to read from
 */
int cosgrid5D::readBinaryFile(const std::string& filename) {

  if (data) {
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k) {
	  for (unsigned int m = 0; m < n3; ++m)
	    delete[] data[i][j][k][m];
	  delete[] data[i][j][k];
	}
	delete [] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
  }

  FILE *fp;
  fp = fopen( filename.c_str(), "rb" );
  if (!fp) {
    stringstream errstrng("");
    errstrng << "Couldn't open file " << filename << " for reading";
    throw CosFitterExcept("cosgrid5D","readBinaryFile",errstrng.str(),1);
  }

  //Read number of axes
  unsigned int naxes;
  fread( &naxes, sizeof(unsigned int), 1, fp );
  if ( naxes != 5 ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Wrong number of axes in file " << filename << " : " <<
      naxes;
    throw CosFitterExcept("cosgrid5D","readBinaryFile",errstrng.str(),2);
  }
      
  //Read axis 1 title
  unsigned int status;
  status = readBinaryAxisName( fp, axisnames[0] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid3D","readBinaryFile",errstrng.str(),4);
  }
  
  fread( &n0, sizeof(unsigned int), 1, fp );
  fread( &min0, sizeof(double), 1, fp );
  fread( &max0, sizeof(double), 1, fp );
  fread( &d0, sizeof(double), 1, fp );

  //Read axis 2
  status = readBinaryAxisName( fp, axisnames[1] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid5D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n1, sizeof(unsigned int), 1, fp );
  fread( &min1, sizeof(double), 1, fp );
  fread( &max1, sizeof(double), 1, fp );
  fread( &d1, sizeof(double), 1, fp );

  //Read axis 3
  status = readBinaryAxisName( fp, axisnames[2] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid5D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n2, sizeof(unsigned int), 1, fp );
  fread( &min2, sizeof(double), 1, fp );
  fread( &max2, sizeof(double), 1, fp );
  fread( &d2, sizeof(double), 1, fp );

  //Read axis 4
  status = readBinaryAxisName( fp, axisnames[3] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid5D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n3, sizeof(unsigned int), 1, fp );
  fread( &min3, sizeof(double), 1, fp );
  fread( &max3, sizeof(double), 1, fp );
  fread( &d3, sizeof(double), 1, fp );

  //Read axis 5
  status = readBinaryAxisName( fp, axisnames[4] );
  if (status) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "Difficulty reading axis title string in file " <<
      filename << endl;
    throw CosFitterExcept("cosgrid5D","readBinaryFile",errstrng.str(),4);
  }
  fread( &n4, sizeof(unsigned int), 1, fp );
  fread( &min4, sizeof(double), 1, fp );
  fread( &max4, sizeof(double), 1, fp );
  fread( &d4, sizeof(double), 1, fp );

  if ( (n0 <= 0) || ( n1 <= 0 ) || ( n2 <= 0 ) || ( n3 <= 0 ) ||
       (n4 <= 0) ) {
    fclose(fp);
    stringstream errstrng("");
    errstrng << "0 axis length in " << filename << endl;
    throw CosFitterExcept("cosgrid5D","readBinaryFile",errstrng.str(),8);
  }
  //Now read the data.  Do the allocation in a seperate step to
  // increase the chances things will be contiguous
  data = new double****[n0];
  for (unsigned int i = 0; i < n0; ++i) {
    data[i] = new double***[n1];
    for ( unsigned int j = 0; j < n1; ++j ) {
      data[i][j] = new double**[n2];
      for (unsigned int k = 0; k < n2; ++k ) {
	data[i][j][k] = new double*[n3];
	for (unsigned int m = 0; m < n3; ++m )
	  data[i][j][k][m] = new double[n4];
      }
    }
  }
  for (unsigned int i = 0; i < n0; ++i ) 
    for (unsigned int j = 0; j < n1; ++j)
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int m = 0; m < n3; ++m)
	  fread( data[i][j][k][m], sizeof(double), n4, fp );

  fclose(fp);

  return 0;
}

cosgrid5D& cosgrid5D::operator=(double val) {
  for (unsigned int i = 0; i < n0; ++i)
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k)
	for (unsigned int m = 0; m < n3; ++m)
	  for (unsigned int n = 0; n < n4; ++n)
	    data[i][j][k][m][n] = val;
  return *this;
}

cosgrid5D& cosgrid5D::operator=(const cosgrid5D& input) {
  //Self copy protection
  if (this == &input) return *this;

  axisnames[0] = input.axisnames[0];
  axisnames[1] = input.axisnames[1];
  axisnames[2] = input.axisnames[2];
  axisnames[3] = input.axisnames[3];
  axisnames[4] = input.axisnames[4];
  spec[0] = input.spec[0];
  spec[1] = input.spec[1];
  spec[2] = input.spec[2];
  spec[3] = input.spec[3];
  spec[4] = input.spec[4];

  if (data == 0) {
    data = new double****[input.n0];
    for (unsigned int i = 0; i < input.n0; ++i) {
      data[i] = new double***[input.n1];
      for (unsigned int j = 0; j < input.n1; ++j) {
	data[i][j] = new double**[input.n2];
	for (unsigned int k = 0; k < input.n2; ++k) {
	  data[i][j][k] = new double*[input.n3];
	  for (unsigned int m = 0; m < input.n3; ++m) {
	    data[i][j][k][m] = new double[input.n4];
	  }
	}
      }
    }
  } else if ( (input.n0 != n0) || (input.n1 != n1) ||
	      (input.n2 != n2) || (input.n3 != n3) || 
	      (input.n4 != n4) ) {      
    for (unsigned int i = 0; i < n0; ++i) {
      for (unsigned int j = 0; j < n1; ++j) {
	for (unsigned int k = 0; k < n2; ++k) {
	  for (unsigned int m = 0; m < n3; ++m) 
	    delete[] data[i][j][k][m];
	  delete[] data[i][j][k];
	}
	delete[] data[i][j];
      }
      delete[] data[i];
    }
    delete[] data;
    
    data = new double****[input.n0];
    for (unsigned int i = 0; i < input.n0; ++i) {
      data[i] = new double***[input.n1];
      for (unsigned int j = 0; j < input.n1; ++j) {
	data[i][j] = new double**[input.n2];
	for (unsigned int k = 0; k < input.n2; ++k) {
	  data[i][j][k] = new double*[input.n3];
	  for (unsigned int m = 0; m < input.n3; ++m) {
	    data[i][j][k][m] = new double[input.n4];
	  }
	}
      }
    }
    
  }
  
  for (unsigned int i = 0; i < input.n0; ++i) 
    for (unsigned int j = 0; j < input.n1; ++j) 
      for (unsigned int k = 0; k < input.n2; ++k) 
	for (unsigned int m = 0; m < input.n3; ++m) 
	  for (unsigned int n = 0; n < input.n4; ++n) 
	    data[i][j][k][m][n] = input.data[i][j][k][m][n];

  n0 = input.n0; min0 = input.min0; max0 = input.max0; d0 = input.d0;
  n1 = input.n1; min1 = input.min1; max1 = input.max1; d1 = input.d1;
  n2 = input.n2; min2 = input.min2; max2 = input.max2; d2 = input.d2;
  n3 = input.n3; min3 = input.min3; max3 = input.max3; d3 = input.d3;
  n4 = input.n4; min4 = input.min4; max4 = input.max4; d4 = input.d4;

  
  return *this;
}

/*!
  \param[in] divisor Number to divide by
  \returns Reference to this object
 */
cosgrid5D& cosgrid5D::operator/=(double divisor) {
  for (unsigned int i = 0; i < n0; ++i) 
    for (unsigned int j = 0; j < n1; ++j) 
      for (unsigned int k = 0; k < n2; ++k) 
	for (unsigned int m = 0; m < n3; ++m) 
	  for (unsigned int n = 0; n < n4; ++n) 
	    data[i][j][k][m][n] /= divisor;
  return *this;
}
	    

cosgrid4D cosgrid5D::collapseAlongZero() const {

  cosgrid4D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3],
		   n4,min4,max4,d4,axisnames[4],spec[4]);
  
  double ****dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int k = 0; k < n2; ++k) {
      for (unsigned int m = 0; m < n3; ++m) {
	for (unsigned int n = 0; n < n4; ++n) {
	  val = 0.0;
	  for (unsigned int i = 0; i < n0; ++i) val += data[i][j][k][m][n];
	  dataptr[j][k][m][n] = val;
	}
      }
    }
  }
  
  return result;
}

cosgrid4D cosgrid5D::collapseAlongOne() const {

  cosgrid4D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3],
		   n4,min4,max4,d4,axisnames[4],spec[4]);
  
  double ****dataptr = result.getData();
  double val;
  
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int k = 0; k < n2; ++k) {
      for (unsigned int m = 0; m < n3; ++m) {
	for (unsigned int n = 0; n < n4; ++n) {
	  val = 0.0;
	  for (unsigned int j = 0; j < n1; ++j) val += data[i][j][k][m][n];
	  dataptr[i][k][m][n] = val;
	}
      }
    }
  }
  
  return result;
}

cosgrid4D cosgrid5D::collapseAlongTwo() const {
  
  cosgrid4D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1],
		   n3,min3,max3,d3,axisnames[3],spec[3],
		   n4,min4,max4,d4,axisnames[4],spec[4]);
  
  double ****dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int m = 0; m < n3; ++m) {
	for (unsigned int n = 0; n < n4; ++n) {
	  val = 0.0;
	  for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k][m][n];
	  dataptr[i][j][m][n] = val;
	}
      }
    }
  }

  return result;
}
 
cosgrid4D cosgrid5D::collapseAlongThree() const {
  
  cosgrid4D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n4,min4,max4,d4,axisnames[4],spec[4]);
  
  double ****dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	for (unsigned int n = 0; n < n4; ++n) {
	  val = 0.0;
	  for (unsigned int m = 0; m < n3; ++m) val += data[i][j][k][m][n];
	  dataptr[i][j][k][n] = val;
	}
      }
    }
  }
  
  return result;
}

cosgrid4D cosgrid5D::collapseAlongFour() const {
  
  cosgrid4D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3]);
  
  double ****dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	for (unsigned int m = 0; m < n3; ++m) {
	  val = 0.0;
	  for (unsigned int n = 0; n < n4; ++n) val += data[i][j][k][m][n];
	  dataptr[i][j][k][m] = val;
	}
      }
    }
  }
  
  return result;
}
 
cosgrid3D cosgrid5D::collapseAlongZeroOne() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3],
		   n4,min4,max4,d4,axisnames[4],spec[4]);

  double ***dataptr = result.getData();
  double val;
  for (unsigned int k = 0; k < n2; ++k) {
    for (unsigned int m = 0; m < n3; ++m) {
      for (unsigned int n = 0; n < n4; ++n) {
	val = 0.0;
	for (unsigned int i = 0; i < n0; ++i)
	  for (unsigned int j = 0; j < n1; ++j) val += data[i][j][k][m][n];
	dataptr[k][m][n] = val;
      }
    }
  }

  return result;
}

cosgrid3D cosgrid5D::collapseAlongZeroTwo() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n3,min3,max3,d3,axisnames[3],spec[3],
		   n4,min4,max4,d4,axisnames[4],spec[4]);
  
  double ***dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int m = 0; m < n3; ++m) {
      for (unsigned int n = 0; n < n4; ++n) {
	val = 0.0;
	for (unsigned int i = 0; i < n0; ++i)
	  for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k][m][n];
	dataptr[j][m][n] = val;
      }
    }
  }
  
  return result;
}

cosgrid3D cosgrid5D::collapseAlongZeroThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n4,min4,max4,d4,axisnames[4],spec[4]);
  
  double ***dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int k = 0; k < n2; ++k) {
      for (unsigned int n = 0; n < n4; ++n) {
	val = 0.0;
	for (unsigned int i = 0; i < n0; ++i)
	  for (unsigned int m = 0; m < n3; ++m) val += data[i][j][k][m][n];
	dataptr[j][k][n] = val;
      }
    }
  }

  return result;
}

cosgrid3D cosgrid5D::collapseAlongZeroFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double ***dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int k = 0; k < n2; ++k) {
      for (unsigned int m = 0; m < n3; ++m) {
	val = 0.0;
	for (unsigned int i = 0; i < n0; ++i)
	  for (unsigned int n = 0; n < n4; ++n) val += data[i][j][k][m][n];
	dataptr[j][k][m] = val;
      }
    }
  }


  return result;
}

cosgrid3D cosgrid5D::collapseAlongOneTwo() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n3,min3,max3,d3,axisnames[3],spec[3],
		   n4,min4,max4,d4,axisnames[4],spec[4]);

  double ***dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int m = 0; m < n3; ++m) {
      for (unsigned int n = 0; n < n4; ++n) {
	val = 0.0;
	for (unsigned int j = 0; j < n1; ++j)
	  for (unsigned int k = 0; k < n2; ++k) val += data[i][j][k][m][n];
	dataptr[i][m][n] = val;
      }
    }
  }

  return result;
}

cosgrid3D cosgrid5D::collapseAlongOneThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n4,min4,max4,d4,axisnames[4],spec[4]);

  double ***dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int k = 0; k < n2; ++k) {
      for (unsigned int n = 0; n < n4; ++n) {
	val = 0.0;
	for (unsigned int j = 0; j < n1; ++j)
	  for (unsigned int m = 0; m < n3; ++m) val += data[i][j][k][m][n];
	dataptr[i][k][n] = val;
      }
    }
  }

  return result;
}

cosgrid3D cosgrid5D::collapseAlongOneFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double ***dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int k = 0; k < n2; ++k) {
      for (unsigned int m = 0; m < n3; ++m) {
	val = 0.0;
	for (unsigned int j = 0; j < n1; ++j)
	  for (unsigned int n = 0; n < n4; ++n) val += data[i][j][k][m][n];
	dataptr[i][k][m] = val;
      }
    }
  }

  return result;
}

cosgrid3D cosgrid5D::collapseAlongTwoThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1],
		   n4,min4,max4,d4,axisnames[4],spec[4]);
  
  double ***dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int n = 0; n < n4; ++n) {
	val = 0.0;
	for (unsigned int k = 0; k < n2; ++k)
	  for (unsigned int m = 0; m < n3; ++m) val += data[i][j][k][m][n];
	dataptr[i][j][n] = val;
      }
    }
  }

  return result;
}
 
cosgrid3D cosgrid5D::collapseAlongTwoFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double ***dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int m = 0; m < n3; ++m) {
	val = 0.0;
	for (unsigned int k = 0; k < n2; ++k)
	  for (unsigned int n = 0; n < n4; ++n) val += data[i][j][k][m][n];
	dataptr[i][j][m] = val;
      }
    }
  }

  return result;
}

cosgrid3D cosgrid5D::collapseAlongThreeFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid3D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2]);

  double ***dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) {
      for (unsigned int k = 0; k < n2; ++k) {
	val = 0.0;
	for (unsigned int m = 0; m < n3; ++m)
	  for (unsigned int n = 0; n < n4; ++n) val += data[i][j][k][m][n];
	dataptr[i][j][k] = val;
      }
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongZeroOneTwo() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n3,min3,max3,d3,axisnames[3],spec[3],
		   n4,min4,max4,d4,axisnames[4],spec[4]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int m = 0; m < n3; ++m) {
    for (unsigned int n = 0; n < n4; ++n) { 
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i) 
	for (unsigned int j = 0; j < n1; ++j) 
	  for (unsigned int k = 0; k < n2; ++k) 
	    val += data[i][j][k][m][n];
      dataptr[m][n] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongZeroOneThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n2,min2,max2,d2,axisnames[2],spec[2],
		   n4,min4,max4,d4,axisnames[4],spec[4]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int k = 0; k < n2; ++k) {
    for (unsigned int n = 0; n < n4; ++n) { 
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i) 
	for (unsigned int j = 0; j < n1; ++j) 
	  for (unsigned int m = 0; m < n3; ++m) 
	    val += data[i][j][k][m][n];
      dataptr[k][n] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongZeroOneFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n2,min2,max2,d2,axisnames[2],spec[2],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int k = 0; k < n2; ++k) {
    for (unsigned int m = 0; m < n3; ++m) { 
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i) 
	for (unsigned int j = 0; j < n1; ++j) 
	  for (unsigned int n = 0; n < n4; ++n) 
	    val += data[i][j][k][m][n];
      dataptr[k][m] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongZeroTwoThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n4,min4,max4,d4,axisnames[4],spec[4]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int n = 0; n < n4; ++n) { 
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i) 
	for (unsigned int k = 0; k < n2; ++k) 
	  for (unsigned int m = 0; m < n3; ++m) 
	    val += data[i][j][k][m][n];
      dataptr[j][n] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongZeroTwoFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int m = 0; m < n3; ++m) { 
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i) 
	for (unsigned int k = 0; k < n2; ++k) 
	  for (unsigned int n = 0; n < n4; ++n) 
	    val += data[i][j][k][m][n];
      dataptr[j][m] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongZeroThreeFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int j = 0; j < n1; ++j) {
    for (unsigned int k = 0; k < n2; ++k) { 
      val = 0.0;
      for (unsigned int i = 0; i < n0; ++i) 
	for (unsigned int m = 0; m < n3; ++m) 
	  for (unsigned int n = 0; n < n4; ++n) 
	    val += data[i][j][k][m][n];
      dataptr[j][k] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongOneTwoThree() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n4,min4,max4,d4,axisnames[4],spec[4]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int n = 0; n < n4; ++n) { 
      val = 0.0;
      for (unsigned int j = 0; j < n1; ++j) 
	for (unsigned int k = 0; k < n2; ++k) 
	  for (unsigned int m = 0; m < n3; ++m) 
	    val += data[i][j][k][m][n];
      dataptr[i][n] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongOneTwoFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n3,min3,max3,d3,axisnames[3],spec[3]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int m = 0; m < n3; ++m) { 
      val = 0.0;
      for (unsigned int j = 0; j < n1; ++j) 
	for (unsigned int k = 0; k < n2; ++k) 
	  for (unsigned int n = 0; n < n4; ++n) 
	    val += data[i][j][k][m][n];
      dataptr[i][m] = val;
    }
  }

  return result;
}

cosgrid2D cosgrid5D::collapseAlongOneThreeFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n0,min0,max0,d0,axisnames[0],spec[0],
		   n2,min2,max2,d2,axisnames[2],spec[2]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int k = 0; k < n2; ++k) { 
      val = 0.0;
      for (unsigned int j = 0; j < n1; ++j) 
	for (unsigned int m = 0; m < n3; ++m)
	  for (unsigned int n = 0; n < n4; ++n) 
	    val += data[i][j][k][m][n];
      dataptr[i][k] = val;
    }
  }
  
  return result;
}

cosgrid2D cosgrid5D::collapseAlongTwoThreeFour() const {
  //This could be constructed from multiple calls to
  // other functions, but for efficiency purposes it
  // is preferrable to do this explicitly
  cosgrid2D result(n1,min1,max1,d1,axisnames[1],spec[1],
		   n2,min2,max2,d2,axisnames[2],spec[2]);

  double **dataptr = result.getData();
  double val;
  for (unsigned int i = 0; i < n0; ++i) {
    for (unsigned int j = 0; j < n1; ++j) { 
      val = 0.0;
      for (unsigned int k = 0; k < n2; ++k) 
	for (unsigned int m = 0; m < n3; ++m)
	  for (unsigned int n = 0; n < n4; ++n) 
	    val += data[i][j][k][m][n];
      dataptr[i][j] = val;
    }
  }

  return result;
}

/*!
  Collapse the grid along the specified axis by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisnum Index of axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid4D containing the collapsed data array
*/
cosgrid4D cosgrid5D::collapseAlongAxis(unsigned int axisnum) const {
  if (axisnum == 0) return collapseAlongZero();
  if (axisnum == 1) return collapseAlongOne();
  if (axisnum == 2) return collapseAlongTwo();
  if (axisnum == 3) return collapseAlongThree();
  if (axisnum == 4) return collapseAlongFour();
  stringstream errstrng("");
  errstrng << "Invalid axis number " << axisnum;
  throw CosFitterExcept("cosgrid5D","collapseAlongAxis",errstrng.str(),1);
}
 
/*!
  Collapse the grid along the specified axis by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisname Nameof axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid4D containing the collapsed data array
*/
cosgrid4D cosgrid5D::collapseAlongAxis(const std::string& axisname) const {
  return collapseAlongAxis( getAxisIndex(axisname) );
}

/*!
  Collapse the grid along the specified axes by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisnum1 Index of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisnum2 Index of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid3D containing the collapsed data array
*/

cosgrid3D cosgrid5D::collapseAlongTwoAxes(unsigned int axisnum1,
					  unsigned int axisnum2) const {
  if (axisnum1 == 0) {
    if (axisnum2 == 1) return collapseAlongZeroOne();
    if (axisnum2 == 2) return collapseAlongZeroTwo();
    if (axisnum2 == 3) return collapseAlongZeroThree();
    if (axisnum2 == 4) return collapseAlongZeroFour();
  } else if (axisnum1 == 1) {
    if (axisnum2 == 0) return collapseAlongZeroOne();
    if (axisnum2 == 2) return collapseAlongOneTwo();
    if (axisnum2 == 3) return collapseAlongOneThree();
    if (axisnum2 == 4) return collapseAlongOneFour();
  } else if (axisnum1 == 2) {
    if (axisnum2 == 0) return collapseAlongZeroTwo();
    if (axisnum2 == 1) return collapseAlongOneTwo();
    if (axisnum2 == 3) return collapseAlongTwoThree();
    if (axisnum2 == 4) return collapseAlongTwoFour();
  } else if (axisnum1 == 3) {
    if (axisnum2 == 0) return collapseAlongZeroThree();
    if (axisnum2 == 1) return collapseAlongOneThree();
    if (axisnum2 == 2) return collapseAlongTwoThree();
    if (axisnum2 == 4) return collapseAlongThreeFour();
  } else if (axisnum1 == 4) {
    if (axisnum2 == 0) return collapseAlongZeroFour();
    if (axisnum2 == 1) return collapseAlongOneFour();
    if (axisnum2 == 2) return collapseAlongTwoFour();
    if (axisnum2 == 3) return collapseAlongThreeFour();
  }
  stringstream errstr("");
  errstr << "Can't collapse along axes: " << axisnum1 << " " << axisnum2;
  throw CosFitterExcept("cosgrid5D","collapseAlongTwoAxes",
			   errstr.str(),1);

}

/*!
  Collapse the grid along the specified axes by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisname1 Name of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisname2 Name of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid3D containing the collapsed data array
*/
cosgrid3D cosgrid5D::collapseAlongTwoAxes(const std::string& axisname1,
					  const std::string& axisname2) const {
  //Collapse grid along 2 axes
  return collapseAlongTwoAxes( getAxisIndex(axisname1),
			       getAxisIndex(axisname2) );

}

/*!
  Collapse the grid along the specified axes by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisnum1 Index of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisnum2 Index of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisnum3 Index of 3rd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid2D containing the collapsed data array
*/
cosgrid2D cosgrid5D::collapseAlongThreeAxes(unsigned int axisnum1,
					    unsigned int axisnum2,
					    unsigned int axisnum3) const {
  if (axisnum1 == 0) {
    if (axisnum2 == 1) {
      if (axisnum3 == 2) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongZeroOneThree();
      if (axisnum3 == 4) return collapseAlongZeroTwoFour();
    }
    if (axisnum2 == 2) {
      if (axisnum3 == 1) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongZeroTwoThree();
      if (axisnum3 == 4) return collapseAlongZeroTwoFour();
    }
    if (axisnum2 == 3) {
      if (axisnum3 == 1) return collapseAlongZeroOneThree();
      if (axisnum3 == 2) return collapseAlongZeroTwoThree();
      if (axisnum3 == 4) return collapseAlongZeroThreeFour();
    }
    if (axisnum2 == 4) {
      if (axisnum3 == 1) return collapseAlongZeroOneFour();
      if (axisnum3 == 2) return collapseAlongZeroTwoFour();
      if (axisnum3 == 3) return collapseAlongZeroThreeFour();
    }
  } else if (axisnum1 == 1) {
    if (axisnum2 == 0) {
      if (axisnum3 == 2) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongZeroOneThree();
      if (axisnum3 == 4) return collapseAlongZeroOneFour();
    }
    if (axisnum2 == 2) {
      if (axisnum3 == 0) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongOneTwoThree();
      if (axisnum3 == 4) return collapseAlongOneTwoFour();
    }
    if (axisnum2 == 3) {
      if (axisnum3 == 0) return collapseAlongZeroOneThree();
      if (axisnum3 == 2) return collapseAlongOneTwoThree();
      if (axisnum3 == 4) return collapseAlongOneThreeFour();
    }
    if (axisnum2 == 4) {
      if (axisnum3 == 0) return collapseAlongZeroOneFour();
      if (axisnum3 == 2) return collapseAlongOneTwoFour();
      if (axisnum3 == 3) return collapseAlongOneThreeFour();
    }
  } else if (axisnum1 == 2) {
    if (axisnum2 == 0) {
      if (axisnum3 == 1) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongZeroTwoThree();
      if (axisnum3 == 4) return collapseAlongZeroTwoFour();
    }
    if (axisnum2 == 1) {
      if (axisnum3 == 0) return collapseAlongZeroOneTwo();
      if (axisnum3 == 3) return collapseAlongOneTwoThree();
      if (axisnum3 == 4) return collapseAlongOneTwoFour();
    }
    if (axisnum2 == 3) {
      if (axisnum3 == 0) return collapseAlongZeroTwoThree();
      if (axisnum3 == 1) return collapseAlongOneTwoThree();
      if (axisnum3 == 4) return collapseAlongTwoThreeFour();
    }
    if (axisnum2 == 4) {
      if (axisnum3 == 0) return collapseAlongZeroTwoFour();
      if (axisnum3 == 1) return collapseAlongOneTwoFour();
      if (axisnum3 == 3) return collapseAlongTwoThreeFour();
    }
  } else if (axisnum1 == 3) {
    if (axisnum2 == 0) {
      if (axisnum3 == 1) return collapseAlongZeroOneThree();
      if (axisnum3 == 2) return collapseAlongZeroTwoThree();
      if (axisnum3 == 4) return collapseAlongZeroThreeFour();
    }
    if (axisnum2 == 1) {
      if (axisnum3 == 0) return collapseAlongZeroOneThree();
      if (axisnum3 == 2) return collapseAlongOneTwoThree();
      if (axisnum3 == 4) return collapseAlongOneThreeFour();
    }
    if (axisnum2 == 2) {
      if (axisnum3 == 0) return collapseAlongZeroTwoThree();
      if (axisnum3 == 1) return collapseAlongOneTwoThree();
      if (axisnum3 == 4) return collapseAlongTwoThreeFour();
    }
    if (axisnum2 == 4) {
      if (axisnum3 == 0) return collapseAlongZeroThreeFour();
      if (axisnum3 == 1) return collapseAlongOneThreeFour();
      if (axisnum3 == 2) return collapseAlongTwoThreeFour();
    } 
  } else if (axisnum1 == 4) {
    if (axisnum2 == 0) {
      if (axisnum3 == 1) return collapseAlongZeroOneFour();
      if (axisnum3 == 2) return collapseAlongZeroTwoFour();
      if (axisnum3 == 3) return collapseAlongZeroThreeFour();
    }
    if (axisnum2 == 1) {
      if (axisnum3 == 0) return collapseAlongZeroOneFour();
      if (axisnum3 == 2) return collapseAlongOneTwoFour();
      if (axisnum3 == 4) return collapseAlongOneThreeFour();
    }
    if (axisnum2 == 2) {
      if (axisnum3 == 0) return collapseAlongZeroTwoFour();
      if (axisnum3 == 1) return collapseAlongOneTwoFour();
      if (axisnum3 == 4) return collapseAlongTwoThreeFour();
    }
    if (axisnum2 == 3) {
      if (axisnum3 == 0) return collapseAlongZeroThreeFour();
      if (axisnum3 == 1) return collapseAlongOneThreeFour();
      if (axisnum3 == 2) return collapseAlongTwoThreeFour();
    }
  }
  stringstream errstr("");
  errstr << "Can't collapse along axes: " << axisnum1 << " " << axisnum2 <<
    " " << axisnum3;
  throw CosFitterExcept("cosgrid5D","collapseAlongThreeAxes",
			   errstr.str(),1);

}

/*!
  Collapse the grid along the specified axes by summing.
  Note that the efficiency of this operation varies considerably
  depending on which axies are summed over.  Generally, summing over
  later indicides is far faster than earlier ones.

  \param[in] axisname1 Name of 1st axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisname2 Name of 2nd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \param[in] axisname3 Name of 3rd axis.  A CosFitterExcpt is thrown if this is
              out of range
  \returns A cosgrid2D containing the collapsed data array
*/
cosgrid2D cosgrid5D::collapseAlongThreeAxes(const std::string& axisname1,
					    const std::string& axisname2,
					    const std::string& axisname3) const {
  //Collapse grid along 3 axes
  return collapseAlongThreeAxes( getAxisIndex(axisname1),
				 getAxisIndex(axisname2),
				 getAxisIndex(axisname3));

}
