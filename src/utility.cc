#include <algorithm>
#include <iterator>
#include <sstream>

#include "utility.h"
#include "cosfitterexcept.h"
#include "paramfile.h"

using namespace std;

/*!
  \param[in] ins The string which is converted
  \return The uppercase version
*/
//This is disturbingly difficult to do in c++
//The locale stuff doesn't work on strings properly, so
// we have to go through C
std::string utility::uppercase(const std::string &ins) {
  unsigned int n;
  n = ins.size();
  char *chr = new char[ n+1 ];
  unsigned int ncopy = ins.copy( chr, n );
  for (unsigned int i=0; i < ncopy; ++i) chr[i]=toupper(chr[i]);
  chr[ncopy] = '\0';  //Or bad things happen
  std::string retstr(chr);
  delete[] chr;
  return retstr;
}

/*! 
  Breaks an input string up into a vector of string, which correspond
  to the input string split on spaces.  Cheerfully stolen from Rob
  Knop.
*/
void utility::stringwords(const std::string &ins,
			  std::vector<std::string> &words) {
  string s,tmp;
  unsigned int i,p;
  int first,last;

  s = ins;

  // Trim spaces from beginning and end

  first=s.find_first_not_of(" ");
  if (first==-1) {
    s="";
  } else {
    last=s.find_last_not_of(" ");
    s=s.substr(first,last-first+1);
  }

  words.clear();

  p=s.find_first_not_of(" \t\r\n");
  if (p>=s.size()) return;
  while ((i=s.find_first_of(" \t\r\n",p))) {
    tmp=s.substr(p,i-p);
    words.push_back(tmp);
    p=s.find_first_not_of(" \t\r\n",i);
    if (p>=s.size()) return;
  }
  tmp=s.substr(p);
  words.push_back(tmp);
}

/*!
  Given three data points, fits a parabola and determines the x value
  of the maximum for that parabola.  Assumes that x1,x2,x3 are equally
  spaced in ascending order
*/
double utility::parabolamax(double x1, double y1, double x2, 
			   double y2, double x3, double y3) {
  //double a;
  double b,c;
  
  //cout << "x1: " << x1 << " y1: " << y1 << endl;
  //cout << "x2: " << x2 << " y2: " << y2 << endl;
  //cout << "x3: " << x3 << " y3: " << y3 << endl;  

  //a = ( (x1*x2*y3 + x2*x3*y1) / 2.0 - x1*x3*y2 );
  b = ( x1*y2 + x3*y2 - ( x2*y1 + x3*y1 + x1*y3 + x2*y3)/2.0 );
  c = ( ( y1 + y3 ) / 2.0 - y2 );

  if (c == 0) {
    string errstrng = "Singular value for c";
    throw CosFitterExcept("cosfitter","parabolamax",errstrng,40);
  }

  //cout << "a: " << a << " b: " << b << " c: " << c << endl;

  return (- b/ (2.0*c));
}

/*!
  Uses parabolamax to estimate the abcissa of the maximum value of the
  input cosgrid1D, given the maximum index
*/
int utility::parabestimate(const cosgrid1D& c1, unsigned int maxindx, 
			   double& abcissa) {

  if ( (maxindx < 1) || (maxindx > (c1.getAxisN() - 2)) ) {
    std::stringstream errstr("");
    errstr << "For parameter: " << c1.getAxisLabel() << std::endl;
    errstr << " Can't estimate parabolic max -- maximum is too close " << 
      "to edge" << std::endl;
    errstr << " maxindx: " << maxindx << " Npoints: " << c1.getAxisN() 
	   << std::endl;
    errstr << " Using maximum value" << std::endl;
    std::cerr << errstr.str();
    abcissa = static_cast<double>(maxindx);
    return 1;
  }

  double const* data = c1.getData();
  double x1,x2,x3,y1,y2,y3;

  x1 = static_cast<double>(maxindx - 1);
  y1 = data[ maxindx - 1 ];
  x2 = static_cast<double>(maxindx);
  y2 = data[ maxindx ];
  x3 = static_cast<double>(maxindx + 1);
  y3 = data[ maxindx + 1 ];

  //cout << "x1: " << x1 << " y1: " << y1 << endl;
  //cout << "x2: " << x2 << " y2: " << y2 << endl;
  //cout << "x3: " << x3 << " y3: " << y3 << endl;

  abcissa = parabolamax(x1,y1,x2,y2,x3,y3);
  return 0;
}

/*!
  \param[in] axis Abscissa values for interpolating.  Must be ordered!
  \param[in] values Data values
  \param[in] x Point to interpolate to
  \param[out] status Status variable, 0 on success.  1 if the value is
    off the beginning, and 2 if it's off the end
  \returns The interpolant

*/
double utility::linterp( std::vector<double>& axis,
			 std::vector<double>& values,
			 double x, int& status ) {
  //We can't make axis and values constant because distance
  // doesn't support const iterators.

  status = 0;
  
  //We will use equal_range to find where the value is
  //After this call, range.first is guaranteed to point at an
  // element that is greater than or equal to x, while range.second
  // points at one which is greater than x.  Or, if we ran off the end,
  // they point at .end() or .begin().
  pair < vector<double>::iterator, vector<double>::iterator > range;
  range = equal_range( axis.begin(), axis.end(), x );

  if ( range.first == axis.begin() ) {
    //First element in array is larger than x
    status = 1;
    return values[0];
  }

  if ( range.first == axis.end() ) {
    //Last element is smaller than x
    status = 2;
    return values[ values.size() - 1 ];
  }

  if ( *range.first == x ) {
    //We actually found the exact value!
    //We have to see if there is a range of exact values
    if ( distance( range.first, range.second ) > 1 ) {
      //There is a range.  It's unclear what to do here, so we
      // shall take the approach of returning the average value
      int i1 = static_cast<int>( distance( axis.begin(), range.first ) );
      int i2 = static_cast<int>( distance( axis.begin(), range.second ) );
      double retval = values[i1];
      for (int i = i1+1; i < i2; ++i) retval += values[i];
      status = 0;
      return retval / static_cast<double>( i2 - i1 );
    } else {
      //This isn't a range, just a single value
      status = 0;
      int i1 = static_cast<int>( distance( axis.begin(), range.first ) );
      return values[i1];
    }
  }

  //This should be the most common case
  //Get index location.  We want i1 to point at something less than x,
  // and i2 to point at something larger.  We've already made sure
  // that i1 will be well defined
  int i1 = static_cast<int>( distance( axis.begin(), range.first ) ) - 1;
  int i2 = static_cast<int>( distance( axis.begin(), range.second ) );

  //Linearly interpolate to get the location of the middle
  double x1, x2, val1, val2, retval;
  x1 = axis[i1];
  x2 = axis[i2];
  val1 = values[i1];
  val2 = values[i2];
  retval = (val2 - val1)/(x2-x1) * ( x - x1 ) + val1;

  status = 0;
  return retval;

}


/*!
  Little utility function to update a vector with values
  of the cosmological parameters in the order that getLumDist
  understands.

  \param[in] code The parameter code from param_tags
  \param[in] fparam The fitparam holding info about the fit
  \param[out] vec The vector to update
  \param[in] val The value to insert

  vec is in the order \f$\Omega_m, \Omega_{DE}, w_0, w_a\f$.
*/
void utility::updateStandardCosVec( param_tags::paramcodes code, 
				    const fitparam& fparam,
				    std::vector<double>& vec,
				    double val ) {
  if ( vec.size() != 4 ) vec.resize(4);
  switch ( code ) {
  case param_tags::omegam :
      vec[0] = val;
      if (fparam.fixcurv) vec[1] = 1.0 - fparam.ocurv - val;
      break;
  case param_tags::omegade :
    vec[1] = val;
    break;
  case param_tags::w0 :
    vec[2] = val;
    break;
  case param_tags::wa :
    vec[3] = val;
    break;
  default :
    //Do nothing -- don't know parameter
    break;
  }
}
