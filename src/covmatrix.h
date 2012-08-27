//covmatrix.h

#include <vector>
#include <fstream>

#ifndef __covmatrix__
#define __covmatrix__

/*!
  \defgroup linalg Linear Algebra
   
  \brief Code assosciated with linear algebra, primarily
    as it relates to covariance matricies.
*/

/*!
  \brief Code to handle covariance matricies
  \ingroup linalg

  More generally, this is specialized for handling
  double precision symmetric positive definite matricies.
  Based closely on the TNT3 beta code by Roldan Pozo.
  The code comes with a default inversion implementation using
  Cholesky decomposition (since covariance matricies are symmetric
  positive definite).  However, there are hooks during compilation
  for linking in LAPACK routines (ATLAS and MKL) which are highly
  recommended if you are thinking of using this code.

  Run-time error checking is mostly avoided for performance
  reasons, so be careful.  Some routines take a status variable,
  similar to the one used in cfitsio.  If a routine is passed
  a status variable that is non-zero, it will return immediately.
  If the routine fails somehow, it will set status to something
  that isn't 0.

*/
class covMatrix {
 private :
  unsigned int n_; //!< Number of rows and columns
  unsigned int size_; //!< Number of elements

  //The internal represenatation is as a single block rather than
  // the c-style because this makes everything allocate contiguously,
  // and because it is necessary for interfacing with external LAPACK
  // libraries
  //Note that the matrix is symmetric so row order versus column order
  // doesn't matter
  double *v_; //!< Holds elements
  double** row_;  //!< row pointers for efficient access
  
  /*! \brief Internal helper function for creating the array of row ptrs */
  void initialize(unsigned int M);

  /*! \brief Internal copying operator */
  void copy(const double* v);

  /*! \brief Sets all matrix elements to specified value */
  void set(const double&);

  /*! \brief Free memory of array */
  void destroy();

  //Routines having to do with the included inversion routine
  /*! \brief Factors matrix into Cholesky form in place */
  void factorCholesky(int& status);

  /*! \brief Inverts the matrix if it's in lower-triangular form*/
  void invertLowerTriangular();

  /*! \brief Multiply a lower triangular matrix by it's transpose */
  void lowerTriangularTransposeMult();

 public :

  /*! \brief Default constructor */
  covMatrix();
  
  /* ! \brief Copy constructor */
  covMatrix(const covMatrix&);

  /* ! \brief Constructor for given size (and value) */
  explicit covMatrix(unsigned int, double val = 0.0);

  /*! \brief Make from C array */
  explicit covMatrix(unsigned int, double**);

  /*! \brief Destructor */
  ~covMatrix() { destroy(); }

  /*! \brief Return the total number of elements */
  unsigned int getSize() const { return size_; }
  /*! \brief Return dimensions */
  unsigned int getNRows() const { return n_; }
  /*! \brief Return dimensions */
  unsigned int getNCols() const { return n_; }

  /*! \brief Change size of matrix, reallocating memory if necessary. */
  covMatrix& resize(unsigned int);
  /*! \brief Resize and set to value */
  covMatrix& resize(unsigned int,double);

  /*! \brief Re-order the elements */
  covMatrix& reorder(const std::vector<unsigned int>& );
  /*! \brief Select elements */
  void select( covMatrix&, const std::vector<unsigned int>&) const;
  /*! \brief Append another covmatrix to the current one */
  covMatrix& append(const covMatrix&);
  /*! \brief Append zeros to curren covmatrix */
  covMatrix& appendzero(unsigned int);

  /*! \brief Clear out elements */
  void clear() { resize(0); }

#ifdef COVFITS
  covMatrix& readFits(const std::string&); //!< Read from fits file
#endif  

  /*! \brief Single subscripting operator */
  inline double* operator[](unsigned int i) {
    return row_[i];
  }
  /*! \brief Single subscripting operator */
  inline double* operator[](int i) {
    return row_[i];
  }

  /*! \brief Single subscripting operator */
  inline double const* operator[](unsigned int i) const {
    return row_[i];
  }
  /*! \brief Single subscripting operator */
  inline double const* operator[](int i) const {
    return row_[i];
  }

  /*! \brief Get axis to raw data internals. Dangerous, but good for performance*/
  inline double* getData() { return v_; }
  /*! \brief Get axis to raw data internals. Dangerous, but good for performance*/
  inline double const* getData() const { return v_; }

  double** operator()() { return row_; } //!< Access to row pointers
  double const* const* operator()() const { return row_; } //!< Const access to row pointers

  /*! \brief Get the total of all the elements */
  double getTotal() const;
  /*! \brief Get the sum of the rows, also returning the total */
  double getRowTotals( std::vector<double>& ) const;

  covMatrix& operator=(const covMatrix& B); //!< Copy one matrix to another
  covMatrix& operator=(double val); //!< Set all elements to value

  /*! \brief Get the diagonal elements */
  std::vector<double> getDiag() const;
  /*! \brief Set the diagonal elements */
  void setDiag( const std::vector<double>& );

  /*! \brief Matrix/Matrix multiplication: C = A*B */
  covMatrix& mult(const covMatrix&A, const covMatrix&B, int& status);

  /*! \brief Matrix/Vector multiplication: \f$A \cdot \vec{b}\f$ */
  std::vector<double> mult( const std::vector<double> b, int& status ) const;

  /*! \brief Matrix/Vector multiplication: \f$A \cdot \vec{b}\f$ */
  std::vector<float> mult( const std::vector<float> b, int& status ) const;

  /*! \brief Add a scalar times another matrix to this one */
  covMatrix& scalarMultAndAdd( double scal, const covMatrix& A, int& status );

  /*! \brief Matrix addition */
  covMatrix& add(const covMatrix&A, int& status );

  /*! \brief Matrix scaling */
  covMatrix& operator*=(double);

  /*! \brief Matrix addition */
  covMatrix& operator+=(const covMatrix&);

  /*! \brief Invert in place*/
  covMatrix& invert(int& status);
  /*! \brief Invert in place, also returning log of determinant*/
  covMatrix& invert(int& status, double& logdet);
  /*! \brief Invert in place, also providing Cholesky decomposed matrix*/
  covMatrix& invert(int& status, covMatrix& triang);

};

/*! \brief Output operator */
std::ostream& operator<<(std::ostream &s, const covMatrix& A);
/*! \brief Input operator */
std::istream& operator>>(std::istream &s, covMatrix& A);

//Linear algebra
/*! \brief Matrix-Matrix multiplication: C = A*B */
covMatrix operator*( const covMatrix &A, const covMatrix& B);
/*! \brief Matrix-vector multiplication */
std::vector<double> operator*(const covMatrix& A, 
			      const std::vector<double>& b);


/*!
  \brief Interface between covariance matricies and SN data
  \ingroup linalg
*/

class covMatrixSNe {
 private:
  unsigned int nsn_; //!< Number of SN
  bool haveMagCov; //!< Do we have a mag covariance matrix
  bool haveWidthCov; //!< Do we have a width covariance matrix
  bool haveColourCov; //!< Do we have a colour covariance matrix
  bool haveMagWidthCov; //!< Do we have a mag-width covariance matrix
  bool haveMagColourCov; //!< Do we have a mag-colour covariance matrix
  bool haveWidthColourCov; //!< Do we have a width-colour covariance matrix

  covMatrix mag_covmatrix; //!< Stores the magnitude covariance matrix
  covMatrix width_covmatrix; //!< Stores the width covariance matrix
  covMatrix colour_covmatrix; //!< Stores the colour covariance matrix
  covMatrix magwidth_covmatrix; //!< Stores the mag-width covariance matrix
  covMatrix magcolour_covmatrix; //!< Stores the mag-colour covariance matrix
  covMatrix widthcolour_covmatrix; //!< Stores the width-colour covariance matrix

  void readCovData(const std::string& FileName,
		   covMatrix& covmat); //!< Reading function for cov matricies
  void writeCovData(const std::string& FileName,
		    const covMatrix& covmat) const; //!< Writing function for cov matricies

 public:
  /*! \brief Default constructor */
  covMatrixSNe();
  
  /*! \brief Full filenames constructor */
  covMatrixSNe( const std::string&, const std::string&, const std::string&, 
		const std::string&, const std::string&, const std::string& );

  /*! \brief Makes the current covmat have the same covmats set
    as the argument, initializing to the given size as needed */
  covMatrixSNe& matchCovmats( const covMatrixSNe&, unsigned int );

  bool haveAnyCovMatrix() const { return haveMagCov || haveWidthCov || 
      haveColourCov || haveMagWidthCov || haveMagColourCov || 
      haveWidthColourCov; } //!< Do we have a full style covariance matrix
  bool haveMagCovMatrix() const { return haveMagCov; } //!< Do we have a mag covariance matrix?
  bool haveWidthCovMatrix() const { return haveWidthCov; } //!< Do we have a width covariance matrix
  bool haveColourCovMatrix() const { return haveColourCov; } //!< Do we have a colour covariance matrix
  bool haveMagWidthCovMatrix() const { return haveMagWidthCov; } //!< Do we have a mag-width covariance matrix
  bool haveMagColourCovMatrix() const { return haveMagColourCov; } //!< Do we have a mag-colour covariance matrix
  bool haveWidthColourCovMatrix() const { return haveWidthColourCov; } //!< Do we have a width-colour covariance matrix

  unsigned int getNsn() const { return nsn_; } //!< Returns number of SN

  covMatrixSNe& resize(unsigned int); //!< Resize the covariance matricies
  /*! \brief Re-order the elements */
  covMatrixSNe& reorder(const std::vector<unsigned int>& );
  /*! \brief Select elements */
  void select( covMatrixSNe&, const std::vector<unsigned int>&) const;
  /*! \brief Append a covmatrix to the current one(s) */
  covMatrixSNe& append(const covMatrixSNe&);
  /*! \brief Append a covmatrix with all zero entries to the current one(s) */
  covMatrixSNe& appendzero(unsigned int);

  //Reads
  void readCovData(const std::vector<std::string>&); //!< Read in all (or a subset) of covariance matricies
  void readMagCovData(const std::string&); //!< Read in the mag covariance matrix from a file 
  void readWidthCovData(const std::string&); //!< Read in the width covariance matrix from a file 
  void readColourCovData(const std::string&); //!< Read in the colour covariance matrix from a file 
  void readMagWidthCovData(const std::string&); //!< Read in the mag-width covariance matrix from a file 
  void readMagColourCovData(const std::string&); //!< Read in the mag-colour covariance matrix from a file 
  void readWidthColourCovData(const std::string&); //!< Read in the width-colour covariance matrix from a file 

  //Writes
  void writeCovData( const std::vector<std::string>&) const; //!< Write all covmatricies
  void writeMagCovData(const std::string&) const; //!< Write the mag covariance matrix to a file
  void writeWidthCovData(const std::string&) const; //!< Write the width covariance matrix to a file
  void writeColourCovData(const std::string&) const; //!< Write the colour covariance matrix to a file
  void writeMagWidthCovData(const std::string&) const; //!< Write the mag-width covariance matrix to a file
  void writeMagColourCovData(const std::string&) const; //!< Write the mag-colour covariance matrix to a file
  void writeWidthColourCovData(const std::string&) const; //!< Write the width-colour covariance matrix to a file

  //Getting covaraince matrix elements and the like
  covMatrix getMagCovMatrix() const; //!< Returns a copy of the mag covariance matrix
  const covMatrix& getMagCovMatrixRef() const { return mag_covmatrix; } //!< Returns a reference to the mag covariance matrix
  double getMagCovElement(unsigned int,unsigned int) const; //!< Returns the corresponding element of the mag covariance matrix
  unsigned int getMagCovMatrixNumRows() const { return mag_covmatrix.getNRows(); } //!< Number of rows in the mag covariance matrix

  //width cov matrix
  covMatrix getWidthCovMatrix() const; //!< Return a copy of the width covariance matrix
  const covMatrix& getWidthCovMatrixRef() const { return width_covmatrix; } //!< Return a reference to the width covariance matrix
  double getWidthCovElement(unsigned int,unsigned int) const; //!< Returns the corresponding element of the width covariance matrix
  unsigned int getWidthCovMatrixNumRows() const { return width_covmatrix.getNRows(); } //!< Number of rows in the width covariance matrix

  //colour cov matrix
  covMatrix getColourCovMatrix() const; //!< Return a copy of the colour covariance matrix
  const covMatrix& getColourCovMatrixRef() const { return colour_covmatrix; } //!< Return a reference to the colour covariance matrix
  double getColourCovElement(unsigned int,unsigned int) const; //!< Returns the corresponding element of the colour covariance matrix
  unsigned int getColourCovMatrixNumRows() const { return colour_covmatrix.getNRows(); } //!< Number of rows in the colour covariance matrix

  //mag-width
  covMatrix getMagWidthCovMatrix() const; //!< Returns a copy of the mag-width covariance matrix
  const covMatrix& getMagWidthCovMatrixRef() const { return magwidth_covmatrix; } //!< Returns a reference to the mag-width covariance matri
  double getMagWidthCovElement(unsigned int,unsigned int) const; //!< Returns the corresponding element of the mag-width covariance matrix
  unsigned int getMagWidthCovMatrixNumRows() const { return magwidth_covmatrix.getNRows(); } //!< Number of rows in the mag-width covariance matrix

  //mag-colour
  covMatrix getMagColourCovMatrix() const; //!< Returns a copy of the mag-colour covariance matrix
  double getMagColourCovElement(unsigned int,unsigned int) const; //!< Returns the corresponding element of the mag-colour covariance matrix
  const covMatrix& getMagColourCovMatrixRef() const { return magcolour_covmatrix; } //!< Returns a reference to the mag-colour covariance matrix
  unsigned int getMagColourCovMatrixNumRows() const { return magcolour_covmatrix.getNRows(); } //!< Number of rows in the mag-colour covariance matrix

  //width-colour
  covMatrix getWidthColourCovMatrix() const; //!< Returns a copy of the width-colour covariance matrix
  const covMatrix& getWidthColourCovMatrixRef() const { return widthcolour_covmatrix; } //!< Returns a reference to the width-colour covariance matrix
  double getWidthColourCovElement(unsigned int,unsigned int) const; //!< Returns the corresponding element of the width-colour covariance matrix
  unsigned int getWidthColourCovMatrixNumRows() const { return widthcolour_covmatrix.getNRows(); } //!< Number of rows in the width-colour covariance matrix

  //Combined covariance matrix, given alpha, beta
  void getCombinedCovMatrix(covMatrix&,double,double) const; //!<Fill in the combined \f$\alpha, \beta\f$ and \f$m_B\f$ cov matrix
  double getCombinedCovElement(double,double,unsigned int,unsigned int) const; //!< Get element of combined covariance matrix
};

#endif
