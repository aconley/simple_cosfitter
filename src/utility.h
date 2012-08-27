//utility.h

#ifndef __cosfitter_utility__
#define __cosfitter_utility__

#include <string>
#include <vector>

#include <cosgrids.h>
#include <param_tags.h>

/*!
  \brief Utility functions
*/
namespace utility {

  std::string uppercase(const std::string &ins);  //!< Converts string to uppercase

  void stringwords(const std::string& ins, std::vector<std::string>& words); //!< Breaks up an input string into a vector

  double parabolamax(double, double, double, double, double, double); //!< Fits a parabola and find the max
  int parabestimate(const cosgrid1D&, unsigned int, double&); //!< Make a parabolic estimate for the maximum of a cosgrid1D

  double linterp( std::vector<double>&, std::vector<double>&,
		  double, int& ); //!< Linearly interpolates on ordered vector

  //Utility vector updater
  void updateStandardCosVec( param_tags::paramcodes, const fitparam&,
			     std::vector<double>&, double ); //!< Maps codes to positions in vec

}

#endif
