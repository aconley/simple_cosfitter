#ifndef __snedata__
#define __snedata__

#include <vector>
#include <string>
#include <iterator>
#include <set>

#include <covmatrix.h>
#include <covmatrix_woodbury.h>

/*!
  \defgroup snedata Supernova Data
  \brief Collection of information about supernovae

  These classes provide a mechanism for holding lists of SN parameters --
  that is, the data the cosmological parameters are fit to.
*/

/*!
  \brief Data structure specialized for holding supernova data.

  \ingroup snedata
*/
struct SNeDataEntry {
  std::string name; //!< Name of supernova
  double zcmb; //!< Redshift of supernova in CMB frame
  double zhel; //!< Heliocentric redshift of supernova
  double var_z; //!< Redshift variance
  double widthpar; //!< Width-luminosity parameter (stretch, delta m_15, etc.)
  double var_widthpar; //!< Variance of widthpar
  double mag; //!< Magnitude
  double var_mag; //!< Variance of magnitude
  double colourpar; //!< Colour related parameter ( E(B-V), A_V, etc. )
  double var_colourpar; //!< Variance of colourpar
  double cov_mag_widthpar; //!< Covariance between magnitude and widthpar
  double cov_mag_colourpar; //!< Covariance between magnitude and colourpar
  double cov_widthpar_colourpar; //!< Covariance between widthpar and colourpar
  double thirdpar;  //!< Third parameter
  double var_thirdpar; //!< Variance of third par
  unsigned int scriptmset; //!< Which scriptm this is in
  unsigned int dataset; //!< Subset identifier

  /*! \brief Constructor */
  SNeDataEntry(std::string SNNAME="",double ZCMB=0,double ZHEL=0,double VARZ=0,
	       double WIDTHPAR=0, double VARWIDTHPAR=0,double MAG=0,
	       double VARMAG=0, double COLOURPAR=0, double VARCOLOURPAR=0, 
	       double COV_MAG_WIDTHPAR=0, double COV_MAG_COLOURPAR=0, 
	       double COV_WIDTHPAR_COLOURPAR=0, double THIRDPAR=0,
	       double VARTHIRDPAR=0, unsigned int SCRIPTMSET=0, 
	       unsigned int DATASET=0) :
    name(SNNAME),zcmb(ZCMB),zhel(ZHEL),var_z(VARZ),widthpar(WIDTHPAR),
    var_widthpar(VARWIDTHPAR),mag(MAG),var_mag(VARMAG),colourpar(COLOURPAR),
    var_colourpar(VARCOLOURPAR),cov_mag_widthpar(COV_MAG_WIDTHPAR),
    cov_mag_colourpar(COV_MAG_COLOURPAR),
    cov_widthpar_colourpar(COV_WIDTHPAR_COLOURPAR),
    thirdpar(THIRDPAR), var_thirdpar(VARTHIRDPAR),
    scriptmset(SCRIPTMSET), dataset(DATASET)
  { };

  /*! \brief Copy Constructor */
  SNeDataEntry(const SNeDataEntry& inval) {
    name = inval.name; zcmb = inval.zcmb; zhel = inval.zhel; 
    var_z = inval.var_z; widthpar = inval.widthpar; 
    var_widthpar = inval.var_widthpar; mag = inval.mag;
    var_mag = inval.var_mag; colourpar = inval.colourpar; 
    var_colourpar = inval.var_colourpar; 
    cov_mag_widthpar = inval.cov_mag_widthpar;
    cov_mag_colourpar = inval.cov_mag_colourpar;
    cov_widthpar_colourpar = inval.cov_widthpar_colourpar;
    thirdpar = inval.thirdpar; var_thirdpar = inval.var_thirdpar;
    scriptmset = inval.scriptmset; dataset = inval.dataset;
  }

  /*! \brief Returns corrected magnitude */
  double getCorrectedMag(double, double) const;

  /*! \brief Returns error on corrected magnitude */
  double getCorrectedMagError(double, double) const;

  /*! \brief Returns variance of corrected magnitude */
  double getCorrectedMagVar(double, double) const;

  unsigned int getDataSet() const { return dataset; } //!< Returns dataset number
  void setScriptmSet(double);
  unsigned int getScriptmSet() const { return scriptmset; }

  bool isCovValid() const; //!< Check if covariance info is valid

  /*! \brief Copy operator */
  SNeDataEntry& operator=(const SNeDataEntry& inval) {
    if (this == &inval) return *this; //Self copy protection
    name = inval.name; zcmb = inval.zcmb; zhel = inval.zhel; 
    var_z = inval.var_z; widthpar = inval.widthpar; 
    var_widthpar = inval.var_widthpar; mag = inval.mag;
    var_mag = inval.var_mag; colourpar = inval.colourpar; 
    var_colourpar = inval.var_colourpar; 
    cov_mag_widthpar = inval.cov_mag_widthpar;
    cov_mag_colourpar = inval.cov_mag_colourpar;
    cov_widthpar_colourpar = inval.cov_widthpar_colourpar;
    thirdpar = inval.thirdpar; var_thirdpar = inval.var_thirdpar;
    scriptmset = inval.scriptmset; dataset = inval.dataset;
    return *this;
  }

};


/*!
  \brief Sorting class for SNeDataEntry based on zcmb.
  \ingroup snedata
*/
class SNeSortByZcmb {
 public:
  /*! \brief Comparison operator on zcmb */
  int operator()(const SNeDataEntry& sn1, const SNeDataEntry& sn2) {
    return sn1.zcmb < sn2.zcmb;
  }
};


/*!
  \brief Collection of SNeDataEntry objects

  \ingroup snedata

  All data points are stored in the vector<SNeDataEntry> points, which you can
  access through the [] operator.  This class also provides serialization
  capability (i.e., it can write the data to a file and read it back).
*/
class SNeData {
 private:
  std::string Name; //!< Just a handy name, potentially useful for debugging
  std::vector<SNeDataEntry> points; //!< The data as SNeDataEntry's

  //Covariance matrix stuff
  bool isdiagonal; //!< True if the errors are diagonal
  bool havefullcov; //!< Do we have a covariance matrix in the full form?
  bool havewoodburycov; //!< Do we have a covariance matrix in the Woodbury form?
  covMatrixSNe covmatrix; //!< Stores the covariance matrix info, non-Woodbury form
  covMatrixWoodburySNe woodbury_covmatrix; //!< Stores a covariance matrix in the woodbury form

 public:
  
  typedef std::vector<SNeDataEntry>::iterator iterator;  //!< Forward iterator
  typedef std::vector<SNeDataEntry>::const_iterator const_iterator; //!< Forward iterator on constant
  typedef std::vector<SNeDataEntry>::reverse_iterator reverse_iterator; //!< Reverse iterator
  typedef std::vector<SNeDataEntry>::const_reverse_iterator const_reverse_iterator; //!< Reverse iterator on constant

  SNeData(); //!< Default constructor
  explicit SNeData(const std::string& name); //!< No data read
  SNeData(const std::string& name, const std::string& fileName); //!< Create SNeData object and read in data from file 
  void readData(const std::string&, bool verbose=false); //!< Read in the data from file 

  void readCovData(const std::vector< std::string >&); //!< Read in full style covariance matricies
  void readWoodburyCovData(const std::string&, const std::string&,
			   const std::string&); //!< Read in a woodbury style covariance matrix

  void writeData(const std::string&) const; //!< Write the data to a file
  void setName(const std::string& name) { Name = name; } //!< Set name of object
  std::string getName() const { return Name; } //!< Return the name of the object

  bool areErrorsDiagonal() const { return isdiagonal; }  //!< Returns true if there are no covariance terms
  bool haveFullCovMatrix() const { return havefullcov; } //!< Do we have a full style covariance matrix
  bool haveWoodburyCovMatrix() const { return havewoodburycov; } //!< Do we have a woodbury style covariance matrix

  covMatrixSNe getCovMatrixSNe() const; //!< Return a copy of the covariance matrix
  const covMatrixSNe& getCovMatrixSNeRef() const { return covmatrix; } //!< Return a constant reference to the Woodbury cov matrix
  covMatrixWoodburySNe getWoodburyCovMatrixSNe() const; //!< Return a copy of the covariance matrix
  const covMatrixWoodburySNe& getWoodburyCovMatrixSNeRef() const { return woodbury_covmatrix; } //!< Return a constant reference to the Woodbury cov matrix

  //Iterators
  iterator begin() { return points.begin(); } //!< Returns read/write iterator pointing to first element
  const_iterator begin() const { return points.begin(); } //!< Returns read only iterator pointing to first element
  iterator end() { return points.end(); } //!< Returns read/write iterator pointing to one past the last element
  const_iterator end() const { return points.end(); } //!< Returns read only iterator pointing to one past the last element

  //Reverse iterators
  reverse_iterator rbegin() { return points.rbegin(); } //!< Returns read/write reverse iterator pointing to last element
  const_reverse_iterator rbegin() const { return points.rbegin(); } //!< Returns read only reverse iterator pointing to last element
  reverse_iterator rend() { return points.rend(); } //!< Returns read/write reverse iterator pointing to one before first element
  const_reverse_iterator rend() const { return points.rend(); } //!< Returns read only reverse iterator pointing to one before first element

  //Size/capacity stuff
  unsigned int capacity() const { return points.capacity(); } //!< Space allocated (but not necessarily filled) in points
  unsigned int size() const { return points.size(); } //!< Number of elements in points
  bool empty() const { return points.empty(); } //!< Is points empty?
  void resize(unsigned int newsize) { 
    points.resize(newsize); 
    if (havefullcov) covmatrix.resize( newsize );
    if (havewoodburycov) woodbury_covmatrix.resize( newsize, 
						    woodbury_covmatrix.getNu());
  } //!< Resize points
  void reserve(unsigned int newsize) { points.reserve(newsize); } //!< Allocate but don't initialize more space in points
  
  //sorting
  void zcmbsort(); //!< Sort supernovae by zcmb

  //concatenation
  SNeData& operator+=(const SNeData&); //!< Concatenate list to end of current one

  //indexing
  SNeData operator[](const std::vector<unsigned int>&) const; //!< Return new list indexed from old
  SNeData operator[](const std::vector<std::string>&) const; //!< Return new list indexed by supernova name
  SNeDataEntry& operator[](const int i) { return points[i]; } //!< Unchecked subscript
  const SNeDataEntry& operator[](const int i) const { return points[i]; } //!< Unchecked subscript
  SNeDataEntry& at(const int i) { return points.at(i); } //!<Checked subscript
  const SNeDataEntry& at(const int i) const { return points.at(i); } //!< Checked subscript

  //Removal
  void remove(const std::vector<int>&); //!< Removes specified elements in plac
  void remove(const std::vector<std::string>&, bool strict=true); //!< Removes specified elements in place using names
  SNeData copy_remove(const std::vector<int>&) const; //!< Returns a new list with indexed entries removed
  SNeData copy_remove(const std::vector<std::string>&, bool strict=true) const; //!< Returns a new list with entries with names in argument removed

  void setScriptmSet(double);
  std::set<unsigned int> getDataSetList() const; //!< Returns the datasets present

  bool isCovValid() const; //!< Check if covariance info is valid for all SN
  bool isCovValid(unsigned int&) const; //!< Check if covariance info is valid for all SN

};

/*!
  \brief Sorting class for SNeDataEntry based on zcmb.

  \ingroup snedata

  This differs from SNeSortByZcmb in that it actually sorts
  an auxillary array based on the contents of a SN list
  instead of actually sorting the entries themselves.  This is
  useful when dealing with the covariance matrix.
*/
class SNeSortByZcmbIndex {
private:
  SNeData& data;  //!< SNeData to sort by
public :
  SNeSortByZcmbIndex( SNeData& d ) : data( d ) { }  //!< Constructor
  /*!
    \brief Comparison operator
   */
  int operator()( unsigned int i, unsigned int j ) {
    return data[i].zcmb < data[j].zcmb;
  }
};

#endif
