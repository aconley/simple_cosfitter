//lumdist_array.h

#ifndef __lumdist_array__
#define __lumdist_array__

#include <vector>

#include "paramfile.h"
#include "cosgrids.h"
#include "snedata.h"
#include "lumdist.h"
#include "cosgrids.h"

/*!
  \brief Builds an array of luminosity distances, 0D case
  \ingroup fitter
*/
class lumdist_array0D {
private :
  unsigned int nz;   //!<Number of redshifts

  std::vector<double> z_hel;  //!< \f$z_{\mbox{hel}}\f$
  std::vector<double> z_cmb;  //!< \f$z_{\mbox{cmb}}\f$
  double* data; //!< Holds actual luminosity distances
  bool good;       //!< Tells you if the lumdist calc was good 
public :
  lumdist_array0D(); //!< Default constructor
  lumdist_array0D( const SNeData&, const fitparam& ); //!< Constructor
  lumdist_array0D( double, double, double, double, 
		   const SNeData& ); //!< Constructor
  ~lumdist_array0D(); //!< Destructor

  void free(); //!< Frees memory


  int loadData( const SNeData&, const fitparam&, bool verbose=false ); //!< Loads data
  int loadData( double, double, double, double, const SNeData&, 
		bool verbose=false ); //!< Loads data

  /*! \brief Check if lumdist is good*/
  bool isGood() const { return good; }

  /*! \brief Return array of lum dists*/
  double const* getLumdists() const {
    return data;
  }

};

/*!
  \brief Builds an array of luminosity distances, 1D case
  \ingroup fitter
*/
class lumdist_array1D {
private :
  unsigned int n0;   //!<Extent along first cosmological parameter
  unsigned int nz;   //!<Number of redshifts

  std::vector<double> z_hel;  //!< \f$z_{\mbox{hel}}\f$
  std::vector<double> z_cmb;  //!< \f$z_{\mbox{cmb}}\f$
  double** data; //!< Holds actual luminosity distances
  bool *good;       //!< Tells you if the lumdist calc was good 
public :
  lumdist_array1D(); //!< Default constructor
  lumdist_array1D( const cosgrid1D&, const SNeData&,
		   const fitparam& ); //!< Constructor
  ~lumdist_array1D(); //!< Destructor

  void free(); //!< Frees memory

  //!< Get the dimensions
  unsigned int getNAxis(unsigned int i) const;

  int loadData( const cosgrid1D&, const SNeData&,
		const fitparam&, bool verbose=false ); //!< Loads data

  /*! \brief Check if lumdist is good for a particular set of params */
  bool isGood(unsigned int i) const { return good[i]; }

  /*! \brief Return array of lum dists*/
  double const* getLumdists(unsigned int i) const {
    return data[i];
  }

};

/*!
  \brief Builds an array of luminosity distances, 2D case
  \ingroup fitter
*/
class lumdist_array2D {
private :
  unsigned int n0;   //!<Extent along first cosmological parameter
  unsigned int n1;   //!<Extend along second cosmological parameter
  unsigned int nz;   //!<Number of redshifts

  std::vector<double> z_hel;  //!< \f$z_{\mbox{hel}}\f$
  std::vector<double> z_cmb;  //!< \f$z_{\mbox{cmb}}\f$
  double*** data; //!< Holds actual luminosity distances
  bool **good;       //!< Tells you if the lumdist calc was good 
public :
  lumdist_array2D(); //!< Default constructor
  lumdist_array2D( const cosgrid2D&, const SNeData&,
		   const fitparam& ); //!< Constructor
  ~lumdist_array2D(); //!< Destructor

  void free(); //!< Frees memory

  int loadData( const cosgrid2D&, const SNeData&,
		const fitparam&, bool verbose=false ); //!< Loads data
  
  /*! \brief Check if lumdist is good for a particular set of params */
  bool isGood(unsigned int i, unsigned int j) const { return good[i][j]; }

  /*! \brief Return array of lum dists*/
  double const* getLumdists(unsigned int i, unsigned int j) const {
    return data[i][j];
  }

};

/*!
  \brief Builds an array of luminosity distances, 3D case
  \ingroup fitter
*/
class lumdist_array3D {
private :
  unsigned int n0;   //!<Extent along first cosmological parameter
  unsigned int n1;   //!<Extend along second cosmological parameter
  unsigned int n2;   //!<Extend along third cosmological parameter
  unsigned int nz;   //!<Number of redshifts

  std::vector<double> z_hel;  //!< \f$z_{\mbox{hel}}\f$
  std::vector<double> z_cmb;  //!< \f$z_{\mbox{cmb}}\f$
  double**** data; //!< Holds actual luminosity distances
  bool ***good;       //!< Tells you if the lumdist calc was good 
public :
  lumdist_array3D(); //!< Default constructor
  lumdist_array3D( const cosgrid3D&, const SNeData&,
		   const fitparam& ); //!< Constructor
  ~lumdist_array3D(); //!< Destructor

  void free(); //!< Frees memory

  int loadData( const cosgrid3D&, const SNeData&,
		const fitparam&, bool verbose=false ); //!< Loads data
  
  /*! \brief Check if lumdist is good for a particular set of params */
  bool isGood(unsigned int i, unsigned int j, unsigned int k) const 
  { return good[i][j][k]; }

  /*! \brief Return array of lum dists*/
  double const* getLumdists(unsigned int i, unsigned int j,
			    unsigned int k) const 
    { return data[i][j][k]; }

};

#endif
