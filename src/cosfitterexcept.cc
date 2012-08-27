#include <cmath>
#include <iostream>
#include <cosfitterexcept.h>

using namespace std;

//Cheerfully stolen from deepexcpt, written by Rob Knop

/*!
  Does the work of the constructors, except for setting the flags
  \param[in] inclass Class generating exception
  \param[in] inmethod Method generating exception
  \param[in] inerrstr Error string 
  \param[in] inerr Error number
*/
void CosFitterExcept::init(const std::string& inclass,
			   const std::string& inmethod,
			   const std::string& inerrstr,int inerr)
{
  errclass=inclass;
  errmethod=inmethod;
  errstr=inerrstr;
  errnum=inerr;
}


CosFitterExcept::CosFitterExcept() {
  init("","","",0);
  classset = methodset = strset = errset = false;
}

/*!
  Most basic error, specifying only the error message
 */
CosFitterExcept::CosFitterExcept(const std::string errstr) {
  init("","",errstr,0);
  strset = true;
}

/*!
  Error with error string and error number
*/
CosFitterExcept::CosFitterExcept(const std::string errstr,int err) {
  init("","",errstr,err);
  strset = errset = true;
}

/*!
  Error with error string, class and method generating exception
*/
CosFitterExcept::CosFitterExcept(const std::string errclass,
				 const std::string errmethod,
				 const std::string errstr) {
  init(errclass,errmethod,errstr,0);
  classset = methodset = strset = true;
}

/*!
  Full error specification: error std::string, number, class, and method.
 */
CosFitterExcept::CosFitterExcept(const std::string errclass,
				 const std::string errmethod,
				 const std::string errstr,int err) {
  init(errclass,errmethod,errstr,err);
  classset = methodset = strset = errset = true;
}

/*
  Provides output capabilities, and is smart enough not to output
  information that hasn't been set.
 */
ostream& operator<<(std::ostream& os, const CosFitterExcept& err) {
  os << "Error encountered" << endl;
  if (err.classset) os << "Error Class/Namespace: " << err.errclass << endl;
  if (err.methodset) os << "Method: " << err.errmethod << endl;
  if (err.strset) os << "Error Message: " << err.errstr << endl;
  if (err.errset) os << "Error Code: " << err.errnum << endl;
  return os;
}
