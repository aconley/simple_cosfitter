#include <iostream>

#include <cosfitterexcept.h>
#include <param_results.h>
#include <paramfile.h>
#include <utility.h>

using namespace std;

/////////////////////////////////////////////////
//            param_results                    //
/////////////////////////////////////////////////
param_results::param_results() {
  name = string("unknown");
  spec = std::pair< param_tags::paramcodes, param_tags::paramtype>
    ( param_tags::unknown, param_tags::unset );
  fit = param_struct::not_set;
  value = 0.0;
  error = 0.0;
  mostlikelyval = 0.0;
  expecval = 0.0;
  errorlow = 0.0;
  errorhigh = 0.0;
}

/*!
  Which value is returned depends on what type of fit was performed.
*/
double param_results::getVal() const {
  switch (fit) {
  case param_struct::not_set :
    throw CosFitterExcept("param_results","getVal",
			  "Unknown fit type for parameter",1);
    break;
  case param_struct::fixed :
    return value;
    break;
  case param_struct::analytic :
    return value;
  case param_struct::loop :
    return expecval;
    break;
  case param_struct::dependent :
    return value;
    break;
  }
  return 0.0; //Should never get here
}

void param_results::getMarginalized1DVals(const cosgrid1D& c) {

  //Minimum number of points to find error
  const unsigned int minpoints_for_error=5; 

  const double conflevel = 0.682690; //!< 1 sigma

  name = c.getAxisLabel();
  spec = c.getAxisSpec();
  fit = param_struct::loop;

  double meanval, cdat, axisval, totprob;
  vector< double > cumprob;

  if (c.getAxisN() == 0) 
    throw CosFitterExcept("param_results","getMarginalized1DVals",
			  "No data in input cosgrid",1);

  unsigned int maxindex;
  c.getMaximum(maxindex);
  mostlikelyval = c.getAxisVal( maxindex );

  meanval = c.getAxisMin() * c[0];

  if (c.getAxisN() <= minpoints_for_error) {
    //We can't calculate errors or anything
    if (c[0] != 0.0) expecval = meanval / c[0];
    return;
  }

  unsigned int n0 = c.getAxisN();
  cumprob.resize( n0 );
  cumprob[0] = c[0];
  totprob = cumprob[0];
  axisval = c.getAxisMin();
  double dval = c.getAxisD();
  for (unsigned int i = 1; i < n0; ++i) {
    axisval += dval;
    cdat = c[i];
    totprob += cdat;
    cumprob[i] = totprob;

    meanval += cdat * axisval;
  }
  meanval /= totprob;
  
  expecval = meanval;
    
  //Use inverse linear interpolation to find the location of this value
  vector<double> axis( n0 );
  for (unsigned int i = 0; i < n0; ++i) 
    axis[i] = c.getAxisVal( i );
  int status;
  double expecloc;
  //Get the axis location of the expectation value by inverse 
  // interpolating (i.e., with cumprob as the axis)
  expecloc = utility::linterp( cumprob, axis, 0.5*totprob, status );
  //Make sure interpolation worked
  if (status != 0) return;

  //These are the values we will use to find the errors
  double highexpecval = totprob * (0.5 + conflevel / 2.0);
  double lowexpecval = totprob * (0.5 - conflevel / 2.0);
  
  //Make sure we have room to do the interpolation
  if ( cumprob[ n0 - 1 ] < highexpecval ) return;
  if ( cumprob[ 0 ] > lowexpecval ) return;

  //Now we hunt 
  double lowexpecloc, highexpecloc;
  lowexpecloc = utility::linterp( cumprob, axis, lowexpecval, status );
  if (status != 0) return;
  highexpecloc = utility::linterp( cumprob, axis, highexpecval, status );
  if (status != 0) return;

  errorlow = expecloc - lowexpecloc;
  errorhigh = highexpecloc - expecloc;

}

/*!
  \ingroup fitter

  \param[in] c Cosgrid to get results from
  \param[out] retvec Vector of results
*/
void GetParamVec( const cosgrid1D& c, std::vector<param_results>& retvec) {
  int st;
  unsigned int indexmax;
  double abcissa;
  retvec.reserve(1);
  retvec[0].getMarginalized1DVals( c );
  c.getMaximum( indexmax );
  st = utility::parabestimate( c, indexmax, abcissa );
  if (st != 0) {
    retvec[0].margpeak1D = c.getAxisVal( indexmax );
  } else {
    retvec[0].margpeak1D = c.getAxisVal( abcissa );
  }

  unsigned int i0;
  c.getMaximum(i0);
  retvec[0].mostlikelyval = c.getAxisVal(i0);

}

/*!
  \ingroup fitter

  \param[in] c Cosgrid to get results from
  \param[out] retvec Vector of results
*/
void GetParamVec( const cosgrid2D& c, std::vector<param_results>& retvec) {
  int st;
  unsigned int indexmax;
  double abcissa;
  retvec.reserve(2);

  for (unsigned int i = 0; i < 2; ++i) {
    cosgrid1D marg = c.collapseAlongAxis( (i + 1) % 2 );
    retvec[i].getMarginalized1DVals( marg );
    marg.getMaximum( indexmax );
    st = utility::parabestimate( marg, indexmax, abcissa );
    if (st != 0) {
      retvec[i].margpeak1D = marg.getAxisVal( indexmax );
    } else {
      retvec[i].margpeak1D = marg.getAxisVal( abcissa );
    }
  }

  unsigned int i0, i1;
  c.getMaximum(i0,i1);
  retvec[0].mostlikelyval = c.getAxisVal(0,i0);
  retvec[1].mostlikelyval = c.getAxisVal(1,i1);
}

/*!
  \ingroup fitter

  \param[in] c Cosgrid to get results from
  \param[out] retvec Vector of results
*/
void GetParamVec( const cosgrid3D& c, std::vector<param_results>& retvec) {
  int st;
  unsigned int indexmax;
  double abcissa;
  retvec.reserve(3);

  cosgrid2D marg2 = c.collapseAlongAxis(0);

  cosgrid1D marg = marg2.collapseAlongAxis( 1 );
  retvec[1].getMarginalized1DVals( marg );
  marg.getMaximum( indexmax );
  st = utility::parabestimate( marg, indexmax, abcissa );
  if (st != 0) {
    retvec[1].margpeak1D = marg.getAxisVal( indexmax );
  } else {
    retvec[1].margpeak1D = marg.getAxisVal( abcissa );
  }

  marg = marg2.collapseAlongAxis( 0 );
  retvec[2].getMarginalized1DVals( marg );
  marg.getMaximum( indexmax );
  st = utility::parabestimate( marg, indexmax, abcissa );
  if (st != 0) {
    retvec[2].margpeak1D = marg.getAxisVal( indexmax );
  } else {
    retvec[2].margpeak1D = marg.getAxisVal( abcissa );
  }

  marg = c.collapseAlongTwoAxes( 1, 2 );
  retvec[0].getMarginalized1DVals( marg );
  marg.getMaximum( indexmax );
  st = utility::parabestimate( marg, indexmax, abcissa );
  if (st != 0) {
    retvec[0].margpeak1D = marg.getAxisVal( indexmax );
  } else {
    retvec[0].margpeak1D = marg.getAxisVal( abcissa );
  }

  unsigned int i0, i1, i2;
  c.getMaximum(i0,i1,i2);
  retvec[0].mostlikelyval = c.getAxisVal(0,i0);
  retvec[1].mostlikelyval = c.getAxisVal(1,i1);
  retvec[2].mostlikelyval = c.getAxisVal(2,i2);

}

/*!
  \ingroup fitter

  \param[in] c Cosgrid to get results from
  \param[out] retmap Map of results
*/
void GetParamMap( const cosgrid1D& c, 
		    std::map<param_tags::paramcodes, param_results>& retmap) {
  int st;
  unsigned indexmax;
  double abcissa;

  param_results temp;

  temp.getMarginalized1DVals( c );
  c.getMaximum( indexmax );
  st = utility::parabestimate( c, indexmax, abcissa );
  if (st != 0) {
    temp.margpeak1D = c.getAxisVal( indexmax );
  } else {
    temp.margpeak1D = c.getAxisVal( abcissa );
  }

  unsigned int i0;
  c.getMaximum(i0);
  temp.mostlikelyval = c.getAxisVal(i0);

  retmap[ temp.spec.first ] = temp;

}

/*!
  \ingroup fitter

  \param[in] c Cosgrid to get results from
  \param[out] retmap Map of results
*/
void GetParamMap( const cosgrid2D& c, 
		    std::map<param_tags::paramcodes, param_results>& retmap) {
  int st;
  unsigned indexmax;
  double abcissa;
  param_results temp;

  for (unsigned int i = 0; i < 2; ++i) {
    cosgrid1D marg = c.collapseAlongAxis( (i + 1) % 2 );
    temp.getMarginalized1DVals( marg );
    marg.getMaximum( indexmax );
    st = utility::parabestimate( marg, indexmax, abcissa );
    if (st != 0) {
      temp.margpeak1D = marg.getAxisVal( indexmax );
    } else {
      temp.margpeak1D = marg.getAxisVal( abcissa );
    }
    retmap[ temp.spec.first ] = temp;
  }

  unsigned int i0, i1;
  c.getMaximum(i0,i1);
  retmap[ c.getAxisSpec(0).first ].mostlikelyval = c.getAxisVal(0,i0);
  retmap[ c.getAxisSpec(1).first ].mostlikelyval = c.getAxisVal(1,i1);
}

/*!
  \ingroup fitter

  \param[in] c Cosgrid to get results from
  \param[out] retmap Map of results
*/
void GetParamMap( const cosgrid3D& c, 
		    std::map< param_tags::paramcodes, param_results >& 
		    retmap ) {
  int st;
  unsigned int indexmax;
  double abcissa;
  param_results temp;

  cosgrid2D marg2 = c.collapseAlongAxis(0);

  cosgrid1D marg = marg2.collapseAlongAxis( 1 );
  temp.getMarginalized1DVals( marg );
  marg.getMaximum( indexmax );
  st = utility::parabestimate( marg, indexmax, abcissa );
  if (st != 0) {
    temp.margpeak1D = marg.getAxisVal( indexmax );
  } else {
    temp.margpeak1D = marg.getAxisVal( abcissa );
  }
  retmap[ temp.spec.first ] = temp;

  marg = marg2.collapseAlongAxis( 0 );
  temp.getMarginalized1DVals( marg );
  marg.getMaximum( indexmax );
  st = utility::parabestimate( marg, indexmax, abcissa );
  if (st != 0) {
    temp.margpeak1D = marg.getAxisVal( indexmax );
  } else {
    temp.margpeak1D = marg.getAxisVal( abcissa );
  }
  retmap[ temp.spec.first ] = temp;

  marg = c.collapseAlongTwoAxes( 1, 2 );
  temp.getMarginalized1DVals( marg );
  marg.getMaximum( indexmax );
  st = utility::parabestimate( marg, indexmax, abcissa );
  if (st != 0) {
    temp.margpeak1D = marg.getAxisVal( indexmax );
  } else {
    temp.margpeak1D = marg.getAxisVal( abcissa );
  }
  retmap[ temp.spec.first ] = temp; 

  unsigned int i0, i1, i2;
  c.getMaximum(i0,i1,i2);
  retmap[ c.getAxisSpec(0).first ].mostlikelyval = c.getAxisVal(0,i0);
  retmap[ c.getAxisSpec(1).first ].mostlikelyval = c.getAxisVal(1,i1);
  retmap[ c.getAxisSpec(2).first ].mostlikelyval = c.getAxisVal(2,i2);

}

std::ostream& operator<<(std::ostream &s, const param_results& p) {
  s << "Parameter: ";
  s << p.name << endl;
  s << " Fit type: ";
  switch ( p.fit ) {
  case param_struct::not_set :
    s << "Unknown" << endl;
    s << " Value: " << p.value;
    if (p.error != 0.0) s << " +- " << p.error;
    break;
  case param_struct::fixed:
    s << "Fixed" << endl;
    s << " Value: " << p.value;
    break;
  case param_struct::analytic :
    s << "Analytic best fit estimate" << endl;
    s << " Value: " << p.value;
    if (p.error != 0.0) s << " +- " << p.error;
    break;
  case param_struct::loop :
    s << "Loop" << endl;
    s << " Most likely value: " << p.mostlikelyval << endl;
    s << " 1D marginalized peak: " << p.margpeak1D << endl;
    s << " 1D marginalized expectation value: " << p.expecval;
    if (p.errorhigh != 0.0) s << " +" << p.errorhigh << " -" <<
      p.errorlow;
    break;
  case param_struct::dependent :
    s << "Dependent variable" << endl;
    s << " Value: " << p.value;
    break;
  }
  return s;
}

/*!
  \param[in] params Map of parameter specifications to be processed
  \param[out] retmap Map that the fixed params are added to

  This is how we propogate the parameters we held fixed in the fit
  into the results.  retmap is appended to, so it can already
  have entries in it.  
 */
void LoadFixedParams( const std::map< param_tags::paramcodes,
		      param_struct>& params,std::map<param_tags::paramcodes, 
		      param_results>& retmap ) {
  param_results presult;
  std::map< param_tags::paramcodes, param_struct >::const_iterator it;
  for ( it = params.begin(); it != params.end(); ++it)
    if ( it->second.fit == param_struct::fixed ) {
      presult.name = it->second.name;
      presult.spec = it->second.param_spec;
      presult.fit = it->second.fit;
      presult.value = it->second.fixval;
      retmap[ it->first ] = presult;
    }
}
