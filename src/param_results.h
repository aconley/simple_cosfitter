//param_results.h

#ifndef __param_results__
#define __param_results__

#include <string>
#include <utility>
#include <vector>
#include <map>

#include "param_tags.h"
#include "cosgrids.h"

/*!
  \brief Structure for holding estimates for a single paramter from a fit.

  \ingroup fitter
*/

struct param_results {
  std::string name; //!< Name of parameter, from param_tags
  param_tags::paramspec spec; //!< Specification of parameter
  param_struct::fittypes fit; //!< What type of fit was performed on this variable
  double value; //!< Single value to report
  double error; //!< Single error to report
  double mostlikelyval; //!< Best value (highest likelihood, lowest \f$\chi^2\f$)
  double margpeak1D; //!< Peak of 1D marginalized value
  double expecval; //!< Expectation value
  double errorlow; //!< Error on low side of expecval
  double errorhigh; //!< Error on high side of expecval

  param_results(); //!< Constructor

  void getMarginalized1DVals(const cosgrid1D&); //!< Calculate the marginalized value from the 1D probability distribution

  double getVal() const; //!< Returns a representative value


};

void GetParamVec( const cosgrid1D& c, std::vector<param_results>& retvec); //!< Fill vector with results of fit
void GetParamVec( const cosgrid2D& c, std::vector<param_results>& retvec); //!< Fill vector with results of fit
void GetParamVec( const cosgrid3D& c, std::vector<param_results>& retvec); //!< Fill vector with results of fit

void GetParamMap( const cosgrid1D& c, std::map<param_tags::paramcodes,
		    param_results>& retvec); //!< Fill vector with results of fit
void GetParamMap( const cosgrid2D& c, std::map<param_tags::paramcodes,
		    param_results>& retvec); //!< Fill vector with results of fit
void GetParamMap( const cosgrid3D& c, std::map<param_tags::paramcodes,
		    param_results>& retvec); //!< Fill vector with results of fit

void LoadFixedParams( const std::map< param_tags::paramcodes,
		      param_struct>& params,
		      std::map<param_tags::paramcodes, 
		      param_results>& retmap ); //!< Loads the fixed parameter specs into the results_map


std::ostream& operator<<(std::ostream &s, const param_results& p); //!< Ouput operator


#endif
