//cosgrids.h
//Cosmology Grid classes

#ifndef __cosgrids__
#define __cosgrids__

#include <string>

#include "paramfile.h"

/*! \defgroup cosgrids Cosmological Grids

\brief Tabulated grids of values as functions of parameters

Set of classes meant to represent some tabulated quantity as a function
of various cosmological parameters.  These are essentially multi-dimensional
arrays with built in support for collapsing the array by summing across
any of the axes.  The axes have names, which can be used as arguments
to the various summing commands.  Thus, one can execute commands like
\code
new_cosgrid2D = cosgrid4D_example.sumAlongTwoAxes("ParameterA","ParameterB");
\endcode
which results in the two dimensions corresponding to ParameterA and ParameterB
being summed over and the resulting (2D) grid being returned.

It would be tempting to do this as a base class and derived classes,
but experiments showed that the performance penalty was unacceptable.
*/

unsigned int NAxesInFile(const std::string& filename, 
			 bool binaryfile = false); //!< Returns the number of axes present in a given serialized cosgrid file

/*! 
  \brief Value grid for one cosmological parameter

  \ingroup cosgrids
  
  Grid class for 1D case.
  
*/

class cosgrid1D {
 private:
  std::string axisnames[1];  //!< Name of cosmological parameter axis
  param_tags::paramspec spec[1]; //!< Spec of parameter string

  int writeTextFile(const std::string&) const; //!< Write this data structure to the specified file as text
  int readTextFile(const std::string&); //!< Read in a serialized version of this object from the specified file as text
  int writeBinaryFile(const std::string&) const; //!< Write this data structure to the specified file in binary format
  int readBinaryFile(const std::string&); //!< Read in serialized binary version of data as written by writeBinaryFile

  double *data; //!< Holds values
  unsigned int n0; //!< Number of entries in data
  double min0; //!< Minimum axis value
  double max0; //!< Maximum axis value
  double d0; //!< Step size along axis

 public:
  cosgrid1D();  //!< Default constructor.  Not very useful
  cosgrid1D(const param_struct&); //!< Constructor
  cosgrid1D(unsigned int, double, double, double, std::string, 
	    param_tags::paramspec); //!< Constructor 
  cosgrid1D(const cosgrid1D&); //!< Copy constructor
  ~cosgrid1D(); //!< Destructor.
  
  void free(); //!< Frees memory

  double* getData() { return data; } //!< Bad form, but needed for the sake of efficiency
  double const* getData() const { return data; } //!< Bad form as above

  double getVal(unsigned int i) const { return data[i]; } //!< Non modifying data access
  double getAxisVal(unsigned int i) const { return min0 + static_cast<double>(i)*d0; } //!< Get axis value at index
  double getAxisVal(double x) const { return min0 + x*d0; } //!< Get axis value at index

  double getMaximum(unsigned int& i0) const; //!< Find maximum value
  double getTotal() const; //!< Total value in grid
  double normalize(); //!< Normalizes the values in the grid to unity.

  unsigned int getAxisIndex(const std::string&) const; //!< Get index of axis corresponding to label, if any
  std::string getAxisLabel() const { return axisnames[0];} //!< Get name of axis
  param_tags::paramspec getAxisSpec() const { return spec[0]; } //!< Get parameter spec

  unsigned int getAxisN() const { return n0; } //!< Return axis extent
  double getAxisMin() const { return min0; } //!< Return min of axis
  double getAxisMax() const { return max0; } //!< Return max of axis
  double getAxisD() const { return d0; } //!< Return delta of axis

  int writeFile(const std::string& filename, bool binary=false) const; //!< Write data to file
  int readFile(const std::string& filename, bool binary=false); //!< Read serialized version of object

  cosgrid1D& operator=(double); //!< Sets all elements to the value
  cosgrid1D& operator+=(const cosgrid1D&); //!< Adds the contents
  cosgrid1D& operator=(const cosgrid1D&); //!< Copies over all data and values
  cosgrid1D& operator*=(const cosgrid1D&); //!< Multiply in place
  cosgrid1D operator*(const cosgrid1D&); //!< Multiplication

  double& operator[](unsigned int i) { return data[i]; } //!< Element access
  const double& operator[](unsigned int i) const { return data[i]; } //!< Constant element access
};
 

/*! 
  \brief Value grid for two cosmological parameters
  
  \ingroup cosgrids

  Grid class for 2D case.

*/

class cosgrid2D {
 private:
  std::string axisnames[2];  //!< Names of axes
  param_tags::paramspec spec[2]; //!< Spec of parameter

  cosgrid1D collapseAlongZero() const; //!< Collapses the grid along the zero axis
  cosgrid1D collapseAlongOne() const; //!< Collapses the grid along the one axis

  int writeTextFile(const std::string&) const; //!< Write this data structure to the specified file as text
  int readTextFile(const std::string&); //!< Read in a serialized version of this object from the specified file as text
  int writeBinaryFile(const std::string&) const; //!< Write this data structure to the specified file in binary format
  int readBinaryFile(const std::string&); //!< Read in serialized binary version of data as written by writeBinaryFile

  double **data; //!< Value array
  unsigned int n0; //!< Number of entries along axis 0
  double min0; //!< Minimum value of axis 0
  double max0; //!< Maximum value of axis 0
  double d0; //!< Step size along axis 
  unsigned int n1; //!< Number of entriea along axis 1
  double min1; //!< Minimum value along axis 1
  double max1; //!< Maximum value along axis 1
  double d1; //!< Step size along axis 1

 public:
  cosgrid2D();  //!< Default constructor
  cosgrid2D(const param_struct&, const param_struct&); //!< Constructor
  cosgrid2D(unsigned int,double,double,double,std::string, 
	    param_tags::paramspec, unsigned int,double,double,double,
	    std::string, param_tags::paramspec);//!< Constructor
  cosgrid2D(const cosgrid2D&); //!< Copy constructor
  ~cosgrid2D(); //!< Destructor
    
  void free(); //!< Frees memory

  void transpose(); //!< Transposes the array

  double** getData() { return data; } //!< Bad form, but needed for efficiency
  double const* const* getData() const { return data; } //!< As above

  double getVal(unsigned int i, unsigned int j) const { return data[i][j]; } //!< Non-modifying data access
  double getAxisVal(unsigned int axis, unsigned int i) const; //!< Get axis value at index
  double getAxisVal(unsigned int axis, double i) const; //!< Get axis value at index

  double getMaximum(unsigned int& i0, unsigned int& i1) const; //!< Returns maximum value in array as well as the position of this maximum
  double getTotal() const; //!< Total value in grid
  double normalize(); //!< Normalizes the grid to unity.

  unsigned int getAxisIndex( const std::string& ) const; //!< Returns index of axis with specified name
  std::string getAxisLabel(unsigned int index) const { return axisnames[index]; } //!< Returns name of axis with specified index
  param_tags::paramspec getAxisSpec(unsigned int index) const { return spec[index]; } //!< Get parameter spec

  unsigned int getAxisN(unsigned int index) const { 
    switch (index) {
    case 0: return n0;
    case 1: return n1;
    default: return 0;
    }
  } //!< Return axis extent

  double getAxisMin(unsigned int index) const { 
    switch (index) {
    case 0: return min0; 
    case 1: return min1;
    default: return 0.0;
    }
  } //!< Return min of axis

  double getAxisMax(unsigned int index) const { 
    switch (index) {
    case 0: return max0;
    case 1: return max1;
    default: return 0.0;
    }
  } //!< Return max of axis
  double getAxisD(unsigned int index) const { 
    switch (index) {
    case 0 : return d0;
    case 1 : return d1;
    default: return 0.0;
    }
  } //!< Return delta of axis


  int writeFile(const std::string& filename, bool binary=false) const; //!< Write data to file
  int readFile(const std::string& filename, bool binary=false); //!< Read serialized version of object

  cosgrid1D collapseAlongAxis(unsigned int axisnum) const;  //!< Collapses the grid along the axis with the specified index
  cosgrid1D collapseAlongAxis(const std::string& axisname) const; //!< Collapses the grid along the axis with the specified label

  cosgrid2D& operator=(double); //!< Sets all elements to value
  cosgrid2D& operator=(const cosgrid2D&); //!< Copies over all data and values
  cosgrid2D& operator/=(double);  //!< Divides all data entries by value

  cosgrid2D& operator*=(const cosgrid2D&); //!< Multiply in place
  cosgrid2D operator*(const cosgrid2D&); //!< Multiplication

  cosgrid2D& operator+=(const cosgrid2D&); //!< Add two grids

  double* operator[](unsigned int i) { return data[i]; } //!< Element access
  double const* operator[](unsigned int i) const { return data[i]; } //!< Constant element access
  
};

/*! 
  \brief Value grid for three cosmological parameters
  
  \ingroup cosgrids
  Grid class for 3D case.
*/
class cosgrid3D {
 private:
  std::string axisnames[3]; //!< Names of axes
  param_tags::paramspec spec[3]; //!< Spec of parameter string

  cosgrid2D collapseAlongZero() const; //!< Collapse the grid along axis 0
  cosgrid2D collapseAlongOne() const; //!< Collapse the grid along axis 1
  cosgrid2D collapseAlongTwo() const; //!< Collapse the grid along axis 2

  cosgrid1D collapseAlongZeroOne() const; //!< Collapse the grid along axes 0 and 1
  cosgrid1D collapseAlongZeroTwo() const; //!< Collapse the grid along axes 0 and 2
  cosgrid1D collapseAlongOneTwo() const; //!< Collapse the grid along axes 1 and 2
  int writeTextFile(const std::string&) const; //!< Write this data structure to the specified file as text
  int readTextFile(const std::string&); //!< Read in a serialized version of this object from the specified file as text
  int writeBinaryFile(const std::string&) const; //!< Write this data structure to the specified file in binary format
  int readBinaryFile(const std::string&); //!< Read in serialized binary version of data as written by writeBinaryFile

  double ***data; //!< Value array

  unsigned int n0; //!< Number of entries along axis 0
  double min0; //!< Minimum value of axis 0
  double max0; //!< Maximum value of axis 0
  double d0; //!< Step size along axis 0
  unsigned int n1; //!< Number of entries along axis 1
  double min1; //!< Minimum value of axis 1
  double max1; //!< Maximum value of axis 1
  double d1; //!< Step size along axis 1
  unsigned int n2; //!< Number of entries along axis 2
  double min2; //!< Minimum value of axis 2
  double max2; //!< Maximum value of axis 2
  double d2; //!< Step size along axis 2

 public:
  cosgrid3D();  //!< Default constructor
  cosgrid3D(const param_struct&, const param_struct&,
	    const param_struct&); //!< Constructor
  cosgrid3D(unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,
	    param_tags::paramspec); //!< Constructor
  cosgrid3D(const cosgrid3D&); //!< Copy constructor
  ~cosgrid3D(); //!< Destructor
  
  void free(); //!< Frees memory

  double*** getData() { return data; } //!< Bad form, but needed for efficiency
  double const* const* const* getData() const { return data; } //!< Bad form as above

  double getVal(unsigned int i, unsigned int j, unsigned int k) const { return data[i][j][k]; } //!< Non-modifying data access
  double getAxisVal(unsigned int axis, unsigned int i) const; //!< Get axis value at index
  double getAxisVal(unsigned int axis, double i) const; //!< Get axis value at index

  unsigned int getAxisIndex(const std::string&) const; //!< Returns index of axis with specified label
  std::string getAxisLabel(unsigned int index) const { return axisnames[index]; }  //!< Returns name of axis with specified index
  param_tags::paramspec getAxisSpec(unsigned int index) const { return spec[index]; } //!< Get parameter spec

  double getMaximum(unsigned int& i0, unsigned int& i1, unsigned int &i2) const; //!< Returns maximum value in array as well as the position of this maximum
  double getTotal() const; //!< Total value in grid
  double normalize(); //!< Normalizes the grid to unity.

  unsigned int getAxisN(unsigned int index) const { 
    switch (index) {
    case 0: return n0;
    case 1: return n1;
    case 2: return n2;
    default: return 0;
    }
  } //!< Return axis extent

  double getAxisMin(unsigned int index) const { 
    switch (index) {
    case 0: return min0; 
    case 1: return min1;
    case 2: return min2;
    default: return 0.0;
    }
  } //!< Return min of axis

  double getAxisMax(unsigned int index) const { 
    switch (index) {
    case 0: return max0;
    case 1: return max1;
    case 2: return max2;
    default: return 0.0;
    }
  } //!< Return max of axis
  double getAxisD(unsigned int index) const { 
    switch (index) {
    case 0: return d0;
    case 1: return d1;
    case 2: return d2;
    default: return 0.0;
    }
  } //!< Return delta of axis


  cosgrid2D collapseAlongAxis(unsigned int) const; //!< Collapses the grid along a single axis with the specified index
  cosgrid2D collapseAlongAxis(const std::string& axisname) const; //!< Collapses the grid along a single axis with the specified label

  cosgrid1D collapseAlongTwoAxes(unsigned int,unsigned int) const; //!< Collapses the grid along two axes with the specified indicies
  cosgrid1D collapseAlongTwoAxes(const std::string&, const std::string&) const; //!< Collapses the grid along two axes with the specified labels

  int writeFile(const std::string& filename, bool binary=false) const; //!< Write data to file
  int readFile(const std::string& filename, bool binary=false); //!< Read serialized version of object
 
  cosgrid3D& operator=(double); //!< Sets all elements to value
  cosgrid3D& operator=(const cosgrid3D&);  //!< Copies over all data and values
  cosgrid3D& operator/=(double);  //!< Divides all data entries by value

  cosgrid3D& operator*=(const cosgrid3D&); //!< Multiply in place
  cosgrid3D& operator+=(const cosgrid3D&); //!< Add two grids in place

  double** operator[](unsigned int i) { return data[i]; } //!< Element access
  double const* const* operator[](unsigned int i) const { return data[i]; } //!< Constant element access

};

/*! 
  \brief Value grid for four cosmological parameters
  
  \ingroup cosgrids
  Grid class for 4D case 
 */
class cosgrid4D {
 private:
  std::string axisnames[4]; //!< Names of axes
  param_tags::paramspec spec[4]; //!< Spec of parameter string

  double ****data;  //!< Value array

  unsigned int n0; //!< Number of entries along axis 0
  double min0; //!< Minimum value of axis 0
  double max0; //!< Maximum value of axis 0
  double d0; //!< Step size along axis 0
  unsigned int n1; //!< Number of entries along axis 1
  double min1; //!< Minimum value of axis 1
  double max1; //!< Maximum value of axis 1
  double d1; //!< Step size along axis 1
  unsigned int n2; //!< Number of entries along axis 2
  double min2; //!< Minimum value of axis 2
  double max2; //!< Maximum value of axis 2
  double d2; //!< Step size along axis 2
  unsigned int n3; //!< Number of entries along axis 3
  double min3; //!< Minimum value of axis 3
  double max3; //!< Maximum value of axis 3
  double d3; //!< Step size along axis 3

  cosgrid3D collapseAlongZero() const; //!< Collapses the grid along the zero axis
  cosgrid3D collapseAlongOne() const; //!< Collapses the grid along the one axis
  cosgrid3D collapseAlongTwo() const; //!< Collapses the grid along the two axis
  cosgrid3D collapseAlongThree() const;  //!< Collapses the grid along the three axis
  cosgrid2D collapseAlongZeroOne() const; //!< Collapses the grid along the zero and one axes
  cosgrid2D collapseAlongZeroTwo() const; //!< Collapses the grid along the zero and two axes
  cosgrid2D collapseAlongZeroThree() const; //!< Collapses the grid along the zero and three axes
  cosgrid2D collapseAlongOneTwo() const; //!< Collapses the grid along the one and two axes
  cosgrid2D collapseAlongOneThree() const; //!< Collapses the grid along the zero and three axes
  cosgrid2D collapseAlongTwoThree() const; //!< Collapses the grid along the two and three axes

  cosgrid1D collapseAlongZeroOneTwo() const; //!< Collapse the grid along zero, one, two axes
  cosgrid1D collapseAlongZeroOneThree() const; //!< Collapse the grid along zero, one, three axes
  cosgrid1D collapseAlongZeroTwoThree() const; //!< Collapse the grid along zero, two, three axes
  cosgrid1D collapseAlongOneTwoThree() const; //!< Collapse the grid along zero, two, three axes

  int writeTextFile(const std::string&) const; //!< Write this data structure to the specified file as text
  int readTextFile(const std::string&); //!< Read in a serialized version of this object from the specified file as text
  int writeBinaryFile(const std::string&) const; //!< Write this data structure to the specified file in binary format
  int readBinaryFile(const std::string&); //!< Read in serialized binary version of data as written by writeBinaryFile  

 public:
  cosgrid4D(); //!< Default constructor
  cosgrid4D(unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,
	    param_tags::paramspec); //!< Constructor
  cosgrid4D(const param_struct&, const param_struct&,
	    const param_struct&, const param_struct&); //!< Constructor
  cosgrid4D(const cosgrid4D&); //!< Copy constructor
  ~cosgrid4D(); //!< Destructor
  

  void free(); //!< Frees memory

  double**** getData() { return data; } //!< Bad form, but necessary for efficiency
  double const* const* const* const* getData() const { return data; } //!< Bad form as above

  double getVal(unsigned int i, unsigned int j, unsigned int k, unsigned int h) const { return data[i][j][k][h]; }  //!< Non-modifying data access
  double getAxisVal(unsigned int axis, unsigned int i) const; //!< Get axis value at index  
  double getAxisVal(unsigned int axis, double i) const; //!< Get axis value at index

  unsigned int getAxisIndex(const std::string&) const; //!< Returns index of axis with specified label
  std::string getAxisLabel(unsigned int index) const { return axisnames[index]; }  //!< Returns name of axis with specified index 
  param_tags::paramspec getAxisSpect(unsigned int index) const { return spec[index]; } //!< Get parameter spec

  unsigned int getAxisN(unsigned int index) const { 
    switch (index) {
    case 0: return n0;
    case 1: return n1;
    case 2: return n2;
    case 3: return n3;
    default: return 0;
    }
  } //!< Return axis extent

  double getAxisMin(unsigned int index) const { 
    switch (index) {
    case 0: return min0; 
    case 1: return min1;
    case 2: return min2;
    case 3: return min3;
    default: return 0.0;
    }
  } //!< Return min of axis

  double getAxisMax(unsigned int index) const { 
    switch (index) {
    case 0: return max0;
    case 1: return max1;
    case 2: return max2;
    case 3: return max3;
    default: return 0.0;
    }
  } //!< Return max of axis
  double getAxisD(unsigned int index) const { 
    switch (index) {
    case 0: return d0;
    case 1: return d1;
    case 2: return d2;
    case 3: return d3;
    default: return 0.0;
    }
  } //!< Return delta of axis

  double getMaximum(unsigned int& i0, unsigned int& i1, unsigned int &i2, unsigned int &i3) const; //!< Returns maximum value in array as well as the position of this maximum
  double getTotal() const; //!< Total value in grid
  double normalize(); //!< Normalizes the grid to unity.

  cosgrid3D collapseAlongAxis(unsigned int) const; //!< Collapses the grid along a single axis with the specified index
  cosgrid3D collapseAlongAxis(const std::string&) const; //!< Collapses the grid along two axes with the specified labels

  cosgrid2D collapseAlongTwoAxes(unsigned int,unsigned int) const; //!< Collapses the grid along two axes with the specified indicies
  cosgrid2D collapseAlongTwoAxes(const std::string&, const std::string&) const; //!< Collapses the grid along two axes with the specified labels

  cosgrid1D collapseAlongThreeAxes(unsigned int,unsigned int,unsigned int) const; //!< Collapses the grid along three axes with the specified indicies
  cosgrid1D collapseAlongThreeAxes(const std::string&, const std::string&, const std::string&) const; //!< Collapses the grid along three axes with the specified labels

  int writeFile(const std::string& filename, bool binary=false) const; //!< Write data to file
  int readFile(const std::string& filename, bool binary=false); //!< Read serialized version of object

  cosgrid4D& operator=(double); //!< Sets all elements to value
  cosgrid4D& operator=(const cosgrid4D&); //!< Copies over all data and values
  cosgrid4D& operator/=(double); //!< Divides all data entries by value

  double*** operator[](unsigned int i) { return data[i]; } //!< Element access
  double const* const* const* operator[](unsigned int i) const { return data[i]; } //!< Constant element access
};

/*! 
  \brief Value grid for five cosmological parameters
  
  \ingroup cosgrids
  Grid class for 5D case 
 */
class cosgrid5D {
 private:
  std::string axisnames[5]; //!< Names of axes
  param_tags::paramspec spec[5]; //!< Spec of parameter string

  double *****data; //!< Value array
  unsigned int n0; //!< Number of entries along axis 0
  double min0; //!< Minimum value of axis 0
  double max0; //!< Maximum value of axis 0
  double d0; //!< Step size along axis 0
  unsigned int n1; //!< Number of entries along axis 1
  double min1; //!< Minimum value of axis 1
  double max1; //!< Maximum value of axis 1
  double d1; //!< Step size along axis 1
  unsigned int n2; //!< Number of entries along axis 2
  double min2; //!< Minimum value of axis 2
  double max2; //!< Maximum value of axis 2
  double d2; //!< Step size along axis 2
  unsigned int n3; //!< Number of entries along axis 3
  double min3; //!< Minimum value of axis 3
  double max3; //!< Maximum value of axis 3
  double d3; //!< Step size along axis 3
  unsigned int n4; //!< Number of entries along axis 4
  double min4; //!< Minimum value of axis 4
  double max4; //!< Maximum value of axis 4
  double d4; //!< Step size along axis 4

  cosgrid4D collapseAlongZero() const; //!< Collapses the grid along the zero axis
  cosgrid4D collapseAlongOne() const; //!< Collapses the grid along the one axis
  cosgrid4D collapseAlongTwo() const; //!< Collapses the grid along the two axis
  cosgrid4D collapseAlongThree() const; //!< Collapses the grid along the three axis
  cosgrid4D collapseAlongFour() const; //!< Collapses the grid along the four axis

  cosgrid3D collapseAlongZeroOne() const; //!< Collapses the grid along the zero and one axes
  cosgrid3D collapseAlongZeroTwo() const; //!< Collapses the grid along the zero and two axes
  cosgrid3D collapseAlongZeroThree() const; //!< Collapses the grid along the zero and three axes
  cosgrid3D collapseAlongZeroFour() const; //!< Collapses the grid along the zero and four axes
  cosgrid3D collapseAlongOneTwo() const; //!< Collapses the grid along the one and two axes
  cosgrid3D collapseAlongOneThree() const; //!< Collapses the grid along the one and three axes
  cosgrid3D collapseAlongOneFour() const; //!< Collapses the grid along the one and four axes
  cosgrid3D collapseAlongTwoThree() const; //!< Collapses the grid along the two and three
  cosgrid3D collapseAlongTwoFour() const; //!< Collapses the grid along the two and four axes
  cosgrid3D collapseAlongThreeFour() const; //!< Collapses the grid along the three and four axes

  cosgrid2D collapseAlongZeroOneTwo() const; //!< Collapse the grid along zero, one, two axes
  cosgrid2D collapseAlongZeroOneThree() const; //!< Collapse the grid along zero, one, three axes
  cosgrid2D collapseAlongZeroOneFour() const; //!< Collapse the grid along zero, one, four axes
  cosgrid2D collapseAlongZeroTwoThree() const; //!< Collapse the grid along zero, two, three axes
  cosgrid2D collapseAlongZeroTwoFour() const; //!< Collapse the grid along zero, two, four axes
  cosgrid2D collapseAlongZeroThreeFour() const; //!< Collapse the grid along zero, three, four axes  
  cosgrid2D collapseAlongOneTwoThree() const; //!< Collapse the grid along one, two, three axes
  cosgrid2D collapseAlongOneTwoFour() const; //!< Collapse the grid along one, two, three axes
  cosgrid2D collapseAlongOneThreeFour() const; //!< Collapse the grid along one, three, four axes
  cosgrid2D collapseAlongTwoThreeFour() const; //!< Collapse the grid along two, three, four axes

  int writeTextFile(const std::string&) const; //!< Write this data structure to the specified file as text
  int readTextFile(const std::string&); //!< Read in a serialized version of this object from the specified file as text
  int writeBinaryFile(const std::string&) const; //!< Write this data structure to the specified file in binary format
  int readBinaryFile(const std::string&); //!< Read in serialized binary version of data as written by writeBinaryFile

 public:
  cosgrid5D(); //!< Default constructor
  cosgrid5D(unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,param_tags::paramspec,
	    unsigned int,double,double,double,std::string,
	    param_tags::paramspec); //!< Constructor
  cosgrid5D(const param_struct&, const param_struct&,
	    const param_struct&, const param_struct&,
	    const param_struct&); //!< Constructor
  cosgrid5D(const cosgrid5D&); //!< Copy constructor
  ~cosgrid5D(); //!< Destructor
  
  void free(); //!< Frees memory
 
  double***** getData() { return data; } //!< Bad form, but necessary for efficiency
  double const* const* const* const* const* getData() const { return data; } //!< Bad form as above

  double getVal(unsigned int i, unsigned int j, unsigned int k, unsigned int m, unsigned int n) const 
    { return data[i][j][k][m][n]; } //!< Non-modifying data access
  double getAxisVal(unsigned int axis, unsigned int i) const; //!< Get axis value at index  
  double getAxisVal(unsigned int axis, double i) const; //!< Get axis value at index

  unsigned int getAxisIndex(const std::string&) const; //!< Returns index of axis with specified label
  std::string getAxisLabel(unsigned int index) const { return axisnames[index]; }  //!< Returns name of axis with specified index 
  param_tags::paramspec getAxisSpec(unsigned int index) const { return spec[index]; } //!< Get parameter spec

  unsigned int getAxisN(unsigned int index) const { 
    switch (index) {
    case 0: return n0;
    case 1: return n1;
    case 2: return n2;
    case 3: return n3;
    case 4: return n4;
    default: return 0;
    }
  } //!< Return axis extent

  double getAxisMin(unsigned int index) const { 
    switch (index) {
    case 0: return min0; 
    case 1: return min1;
    case 2: return min2;
    case 3: return min3;
    case 4: return min4;
    default: return 0.0;
    }
  } //!< Return min of axis

  double getAxisMax(unsigned int index) const { 
    switch (index) {
    case 0: return max0;
    case 1: return max1;
    case 2: return max2;
    case 3: return max3;
    case 4: return max4;
    default: return 0.0;
    }
  } //!< Return max of axis
  double getAxisD(unsigned int index) const { 
    switch (index) {
    case 0: return d0;
    case 1: return d1;
    case 2: return d2;
    case 3: return d3;
    case 4: return d4;
    default: return 0.0;
    }
  } //!< Return delta of axis


  double getMaximum(unsigned int& i0, unsigned int& i1, unsigned int &i2, unsigned int &i3, unsigned int &i4) const; //!< Returns maximum value in array as well as the position of this maximum
  double getTotal() const; //!< Total value in grid
  double normalize(); //!< Normalizes the grid to unity.

  cosgrid4D collapseAlongAxis(unsigned int) const; //!< Collapses the grid along the axis with the specified index
  cosgrid4D collapseAlongAxis(const std::string&) const; //!< Collapses the grid along the axis with the specified label

  cosgrid3D collapseAlongTwoAxes(unsigned int,unsigned int) const; //!< Collapses the grid along two axes with the specified indicies
  cosgrid3D collapseAlongTwoAxes(const std::string&, const std::string&) const; //!< Collapses the grid along two axes with the specified labels

  cosgrid2D collapseAlongThreeAxes(unsigned int,unsigned int,unsigned int) const; //!< Collapses the grid along three axes with the specified indicies
  cosgrid2D collapseAlongThreeAxes(const std::string&, const std::string&, const std::string&) const; //!< Collapses the grid along three axes with the specified labels

  int writeFile(const std::string& filename, bool binary=false) const; //!< Write data to file
  int readFile(const std::string& filename, bool binary=false); //!< Read serialized version of object

  cosgrid5D& operator=(double); //!< Sets all elements to value
  cosgrid5D& operator=(const cosgrid5D&); //!< Copies over all data and values
  cosgrid5D& operator/=(double); //!< Divides all data entries by value

  double**** operator[](unsigned int i) { return data[i]; } //!< Element access
  double const* const* const* const* operator[](unsigned int i) const { return data[i]; } //!< Constant element access

};

#endif
