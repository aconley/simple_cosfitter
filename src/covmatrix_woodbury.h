#include <vector>
#include <fstream>

#ifndef __covmatrix_woodbury__
#define __covmatrix_woodbury__

#if USEACCELERATE

#include <Accelerate/Accelerate.h>

#elif ! (USEMKL || USEATLAS)
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#endif

#include "covmatrix.h"

/*!
  \brief Covariance matricies that are of the form that
         allows them to be inverted using Sherman-Morrison-Woodbury
  \ingroup linalg

  Specifically, we must have a \f$N \times N\f$ covariance matrix that is 
  comprised of diagonal terms \f$\left( \mathbf{D} \right)\f$
  plus the products of a limited number of rank N updates
  \f[
     \mathbf{V} = \mathbf{D} + \sum^M_i \mathbf{u}_i \otimes \mathbf{u}_i
  \f]
  for the \f$M\f$ vectors \f$\mathbf{u}_i\f$, each of which has \f$N\f$
  elements.  Instead of doing these one at a time, we use the block
  matrix form (Sherman-Morrison-Woodbury) and use the following relation
  to do the inverse:
  \f[
    \mathbf{V}^{-1} = \mathbf{D}^{-1} - \left[ \mathbf{D}^{-1}
      \cdot \mathbf{U} \cdot \left( \mathbf{1} +
      \mathbf{U}^T \cdot \mathbf{D}^{-1} \cdot \mathbf{U} \right)^{-1}
      \cdot \mathbf{U}^T \cdot \mathbf{D}^{-1} \right]
  \f]
  The required inverse is now only of a \f$M \times M\f$ matrix,
  so this is only of interest if \f$M < N\f$.  The \f$M\f$ columns of 
  \f$\mathbf{U}\f$ are just the \f$\mathbf{u}\f$.

*/
class covMatrixWoodbury {
private :

  unsigned int n_; //!< Number of diagonal elements
  unsigned int m_; //!< Number of columns of \f$\mathbf{U}\f$ (U is n_ by m_)
 
  std::vector<double> diag_; //!< Diagonal elements \f$\mathbf{D}\f$
  
  //Internal matricies for computations
#if USEMKL || USEATLAS 
  //The internal representation of U is as a single block
  // rather than as a C-style matrix.  This is to ensure that
  // the storage is allocated contiguously, and to allow interfacing
  // with external LAPACK libraries.
  double *U; //<! \f$\mathbf{U}\f$, size n_ by m_
  double *working; //!< A working array of size m_ by n_
  double *working2; //!< A working array of size n_ by m_
  double *workingsq; //!< A working array of size m_ by m_
  int *ipiv; //!< Pivot array for LU decomp
#elif USEACCELERATE
  //Like the above, but ipiv is a different type
  double *U; //<! \f$\mathbf{U}\f$, size n_ by m_
  double *working; //!< A working array of size m_ by n_
  double *working2; //!< A working array of size n_ by m_
  double *workingsq; //!< A working array of size m_ by m_
  __CLPK_integer *ipiv; //!< Pivot array for LU decomp
#else
  //Here the representation is as a gsl_matrix
  gsl_matrix *U; //<! \f$\mathbf{U}\f$, size n_ by m_
  gsl_matrix *working; //!< A working array of size m_ by n_
  gsl_matrix *working2; //!< A working array of size n_ by m_
  gsl_matrix *workingsq; //!< A working array of size m_ by m_
  gsl_matrix *inverse; //!< The gsl matrix version of the inverse
  gsl_permutation *ipiv; //!< Pivot array for LU decomp
#endif

  /*! \brief Initialize memory */
  void initialize(unsigned int, unsigned int);

  /*! \brief Free memory of array */
  void destroy();

public :

  /*! \brief Default constructor */
  covMatrixWoodbury();

  /* ! \brief Copy constructor */
  covMatrixWoodbury(const covMatrixWoodbury&);

  /* ! \brief Constructor for given size (and diagonal value) */
  explicit covMatrixWoodbury(unsigned int n, unsigned int m_,
			     double diagval = 0.0);

  /* ! \brief Constructor given diagonal and U) */
  explicit covMatrixWoodbury(const std::vector<double>& diag,
			     unsigned int m, const double *const u);

  /*! \brief Destructor */
  ~covMatrixWoodbury() { destroy(); }

  /*! \brief Return the total number of elements */
  unsigned int getSize() const { return n_ * n_; }
  /*! \brief Return dimensions */
  unsigned int getNRows() const { return n_; }
  /*! \brief Return dimensions */
  unsigned int getNCols() const { return n_; }
  /*! \brief Return number of columns of \f$\mathbf{U}\f$ */
  unsigned int getUNCols() const { return m_; }

  /*! \brief Change size of matrix, reallocating memory if necessary. */
  covMatrixWoodbury& resize(unsigned int, unsigned int);
  /*! \brief Resize and set diag to value.  Clears out U */
  covMatrixWoodbury& resize(unsigned int,unsigned int, double);

  /*! \brief Clear out elements */
  void clear() { destroy(); }

  /*! \brief Access Diagonal elements */
  double getDiagElement( unsigned int i ) const { return diag_[i]; }
  
  /*! \brief Set all diagonal elements to value */
  void setDiag( double value ) { diag_.assign( n_, value ); }
  /*! \brief Set diagonal elements to a copy of vector */
  void setDiag( const std::vector<double> d );

  /*! \brief Fill column of U with vector*/
  void setUColumn( unsigned int, const std::vector<double>&) ;

  /*! \brief Get diagonal vector copy */
  void getDiag( std::vector<double>& d ) const { d = diag_; }
  /*! \brief Get column of U*/
  void getUColumn( unsigned int, std::vector<double>&) const;

  /*! \brief Re-order the elements */
  void reorder( const std::vector<unsigned int>& ); 

  /*! \brief Select elements */
  void select( covMatrixWoodbury&, const std::vector<unsigned int>&) const;

  /*! \brief Build and return the actual covariance matrix */
  void getCovMatrix( covMatrix& cov ) const;

  /*! \brief Build and return the inverse covaraince matrix */
  void getInvCovMatrix( covMatrix& invcov, int& status ) const;

  /*! \brief Read from file */
  void readFile( const std::string& filename );

  /*! \brief Write to file */
  void writeFile( const std::string& filename ) const;

  /*! \brief Copy one matrix to another */
  covMatrixWoodbury& operator=(const covMatrixWoodbury& B); 

  /*! \brief Add two covmatricies */
  covMatrixWoodbury& add(const covMatrixWoodbury&,int&);

  /*! \brief Add two covmatricies */
  covMatrixWoodbury& operator+=(const covMatrixWoodbury&);

#if USEMKL || USEATLAS || USEACCELERATE 
  double* getU() { return U; } //!< Dangerous -- used for efficiency
#else
  gsl_matrix* getU() { return U; } //!< Dangerous -- used for efficiency
#endif

};

/*!
  \brief An interface between SN data and a Woodbury form covariance matrix

  \ingroup linalg
*/
class covMatrixWoodburySNe {
 private:
  unsigned int nsn_; //!< Number of SN
  unsigned int nu_; //!< Number of u vectors
  
  //Note that we don't store the diagonal!
#if USEMKL || USEATLAS || USEACCELERATE
  //We have to store 3 nsn by nu matricies holding the m_B, s, c
  // components of V
  double *U0; //!< \f$\mathbf{U}\f$ for \f$m_B\f$, size nsn by nu
  double *Ua; //!< \f$\mathbf{U}\f$ for \f$s\f$, size nsn by nu
  double *Ub; //!< \f$\mathbf{U}\f$ for \f$\mathcal C\f$, size nsn by nu
#else
  //Here the representations are as gsl_matrix-es
  gsl_matrix *U0; //!< \f$\mathbf{U}\f$ for \f$m_B\f$, size nsn by nu
  gsl_matrix *Ua; //!< \f$\mathbf{U}\f$ for \f$s\f$, size nsn by nu
  gsl_matrix *Ub; //!< \f$\mathbf{U}\f$ for \f$\mathcal C\f$, size nsn by nu
  gsl_matrix *tmp; //!< Helper matrix, size nsn by nu
#endif

  mutable covMatrixWoodbury covWood; //!< Storage space for current matrix in woodbury form

  void initialize(unsigned int, unsigned int); //!< Initialize memory
  void destroy(); //!< Frees the memory

 public:
  /*! \brief Default constructor */
  covMatrixWoodburySNe();

  /*! \brief Constructor with size*/
  covMatrixWoodburySNe(unsigned int, unsigned int);

  /*! \brief Read from files constructor*/
  covMatrixWoodburySNe( const std::string&, const std::string&, 
			const std::string& );

  /*! \brief Destructor */
  ~covMatrixWoodburySNe() { destroy(); }

  unsigned int getNsn() const { return nsn_; } //!< Get number of SN
  unsigned int getNu() const { return nu_; } //!< Get number of u vectors
  
  /*! \brief Change the size, reallocating memory if necessary */
  covMatrixWoodburySNe& resize(unsigned int, unsigned int);

  /*! \brief Clear out elements */
  void clear() { destroy(); }

  /*! \brief Builds the woodbury style covariance matrix */
  covMatrixWoodbury& buildWoodburyCovMatrix(const std::vector<double>&, 
					    double, double) const;

  /*! \brief Returns the normal form covariance matrix */
  void getCovMatrix( covMatrix&, const std::vector<double>&, 
		     double, double) const;

  /*! \brief Returns the inverse covariance matrix */
  void getInvCovMatrix( covMatrix&, const std::vector<double>&, 
			double, double, int&) const;

  /*! \brief Read from file */
  void readFile( const std::string&, const std::string&, const std::string& );

  /*! \brief Write to file */
  void writeFile( const std::string&, const std::string&, 
		  const std::string& ) const;

  /*! \brief Re-order the elements */
  void reorder( const std::vector<unsigned int>& ); 

  /*! \brief Select elements */
  void select( covMatrixWoodburySNe&, const std::vector<unsigned int>&) const;

  /*! \brief Add two covmatricies */
  covMatrixWoodburySNe& add(const covMatrixWoodburySNe&,int&);

  /*! \brief Add two covmatricies */
  covMatrixWoodburySNe& operator+=(const covMatrixWoodburySNe&);

  /*! \brief Append a covmatrix to the current one */
  covMatrixWoodburySNe& append(const covMatrixWoodburySNe&,int&);

  /*! \brief Copy one matrix to another */
  covMatrixWoodburySNe& operator=(const covMatrixWoodburySNe&); 

};

#endif
