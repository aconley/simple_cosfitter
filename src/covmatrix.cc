//covmatrix.cc

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include "covmatrix.h"
#include "cosfitterexcept.h"

#if USEMKL
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#elif USEATLAS
extern "C" {
  #include <clapack.h>
  #include <cblas.h>
} 
#elif USEACCELERATE
#include <Accelerate/Accelerate.h>
#endif

#ifdef COVFITS
#include <fitsio.h>
#endif

//We avoid exceptions and most run-time error checking for 
//performance reasons, so be careful.

/*!
  \param[in] N Number of rows and columns
 */
void covMatrix::initialize(unsigned int N) {
  size_ = N*N;
  n_ = N;
  
  if ( N == 0 ) {
    v_ = NULL;
    row_ = NULL;
  } else {
    v_ = new double[size_]; 
    row_ = new double*[N];
    double* p = v_;              
    for (unsigned int i=0; i<N; i++) {
      row_[i] = p;
      p += N ;
    }      
  }
}

/*!
  \param[in] v Data to copy into internal array
*/
void covMatrix::copy(const double* v) {
  for (unsigned int i = 0; i < size_; ++i)
    v_[i] = v[i];
}

/*! 
  \param[in] val Sets all elements to this value
*/
void covMatrix::set(const double& val) {
  for (unsigned int i = 0; i < size_; ++i)
    v_[i] = val;
}

void covMatrix::destroy() {     
  /* do nothing, if no memory has been previously allocated */
  if (v_ == NULL) return ;
  
  /* if we are here, then matrix was previously allocated */
  if (v_ != NULL) delete [] (v_);     
  if (row_ != NULL) delete [] (row_);
  v_ = NULL;
  row_ = NULL;
  n_ = size_ = 0;
}

/*!
  This is internal because the user shouldn't be able to interface
  with the matrix when it's in lower triangular form.  Also, if you
  pass it something that isn't lower triangular, bad things will happen.
*/
void covMatrix::invertLowerTriangular() {
  double sum, *rowptr;
  for (unsigned int i = 0; i < n_; ++i) {
    rowptr = row_[i];
    //Do diagonal
    rowptr[i] = 1.0 / rowptr[i];
    for (unsigned int j = 0; j < i; ++j)
      row_[j][i] = 0.0;
    for (unsigned int j = i+1; j < n_; ++j) {
      rowptr = row_[j];
      sum = 0.0;
      for (unsigned int k = i; k < j; ++k)
	sum -= rowptr[k] * row_[k][i];
      rowptr[i] = sum / rowptr[j];
    }
  }
    
}

/*
  Factors the matrix into it's lower diagonal Cholesky form
  such that \f$A = L \cdot L^T\f$.  If the process fails,
  possibly because the matrix is not symmetric positive definite,
  then status will be non-zero on return.

  \param status 0 on success, -1 if non-spd
 */
void covMatrix::factorCholesky(int& status) {
  if (status) return;
  double d,s, *ptri, *ptrj;
  bool isspd;
  isspd = true;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    d = ptri[i];
    for (unsigned int j = 0; j < i; ++j)
      d -= ptri[j]*ptri[j];
    isspd = isspd && ( d > 0 );
    d = sqrt( d > 0.0 ? d : 0.0 );
    ptri[i] = d;
    for (unsigned int j = i+1; j < n_; ++j) {
      ptrj = row_[j];
      s = ptri[j];
      for (unsigned int k = 0; k < i; ++k)
	s -= ptri[k]*ptrj[k];
      ptrj[i] = ( d > 0 ? s/d : 0.0 );
    }
  }
  if ( isspd ) status = 0; else status = -1;
}
    
/*!
  Form the product of a lower triangular matrix and it's transpose,
  overwriting the current matrix.

  Note that the result is always symmetric.  This is internal
  because the user shouldn't interface with the lower triangular
  forms.
*/
void covMatrix::lowerTriangularTransposeMult() {

  double aii, sum, val;

  for (unsigned int i = 0; i < n_; ++i) {
    aii = v_[n_*i + i];
    if ( i < n_-1 ) {
      //Diagonal term
      sum = aii * aii;
      for ( unsigned int j = i+1; j < n_; ++j ) {
	val = v_[n_ * j + i];
	sum += val * val;
      }
      v_[ n_ * i + i ] = sum;

      //Do the lower part only, reflect later
      double *ptrk;
      for ( unsigned int j = 0; j < i; ++j) {
	//We have to handle the k=i case seperately,
	// since we already adjusted it's value
	ptrk = row_[i];
	sum = aii * ptrk[j];
	for ( unsigned int k = i+1; k < n_; ++k) {
	  ptrk = row_[k];
	  sum += ptrk[i] * ptrk[j];
	}
	v_[ n_ * i + j ] = sum;
      }
    } else {					
      //Bottom row, simple form
      double *rowptr;
      rowptr = row_[n_-1];
      for (unsigned int k = 0; k < n_; ++k)
	rowptr[k] *= aii;
    }
  }

  //Now reflect to upper half
  for (unsigned int i = 0; i < n_; ++i)
    for (unsigned int j = i+1; j < n_; ++j)
      v_[n_ * i + j] = v_[n_ * j + i];      
}

covMatrix::covMatrix() : n_(0), size_(0), v_(NULL), row_(NULL) { };

covMatrix::covMatrix(const covMatrix &A) {
  initialize(A.n_);
  copy(A.v_);
}

covMatrix::covMatrix(unsigned int M, double value) {
  initialize(M);
  set(value);
}

covMatrix::covMatrix(unsigned int N, double** arr) {
  initialize(N);
  unsigned int cntr;
  double *rowptr;
  for (unsigned int i = 0; i < N; ++i) {
    rowptr = arr[i];
    cntr = n_ * i;
    for (unsigned int j = 0; j < N; ++j) 
      v_[ cntr + j ] = rowptr[j];
  }
}

/*
  This operations occurs in place, i.e. when resizing to 
  a new matrix, original matrix elements
  are <b>NOT</b> retained.  Instead, one must explicit create
  a new matrix of this size and manually copy the elements, e.g.
  <pre>
  
  covMatrix B(N);
  
  unsigned int min_N = N < A.getNRows() ? N : A.getNRows();
  for (unsigned int i=0; i<=min_N; i++)
  for (unsigned int j=0; j<=min_N; j++)
  B[i][j] = A[i][j];
  
  A.destroy();
  </pre>
  
  \param[in] N new size
*/
covMatrix& covMatrix::resize(unsigned int N) {
  if ( n_ == N ) {
    return *this;
  }
  destroy();
  initialize(N);
  return *this;
}

covMatrix& covMatrix::resize(unsigned int N, double value) {
  if ( n_ == N ) {
    set(value);
    return *this;
  }
  destroy();
  initialize(N);
  set(value);
  return *this;
}

/*!
  There is no requirement that the indexing vector have unique
  elements, just that the size doesn't change.  So this has to
  be done with a full copy.
*/
covMatrix& covMatrix::reorder(const std::vector<unsigned int>& idx ) {
    if (idx.size() != n_)
      throw CosFitterExcept("covMatrix","reorder",
			    "index array not same size as data array",1);
    
    double *newv;
    unsigned int locold, locnew;
    newv = new double[ size_ ];
    //Recall matrix is symmetric, so ignore row/col ordering issues
    for (unsigned int i = 0; i < n_; ++i) {
      locnew = i*n_;
      locold = idx[i]*n_;
      for (unsigned int j = 0; j < n_; ++j)
	newv[locnew + j] = v_[locold+idx[j]];
    }
    delete[] v_;
    v_ = newv;
    //Fix the row pointers.  Note we don't have to re-alloc
    double* p = v_; 
    for (unsigned int i=0; i < n_ ; i++) {
      row_[i] = p;
      p += n_ ;
    }      
    return *this;
}

/*!
  \param[out] newmat Matrix with output information
  \param[in] idx  Index array to select with.  The boundaries are
                  not checked
*/
void covMatrix::select( covMatrix& newmat,
			const std::vector<unsigned int>& idx ) const {
  unsigned int n = idx.size();
  if (newmat.n_ != n ) newmat.resize( n );
  //Recall matrix is symmetric, so ignore row/col ordering issues
  unsigned int locnew, locold;
  for (unsigned int i = 0; i < n; ++i) {
    locnew = i*n;
    locold = idx[i] * n_;
    for (unsigned int j = 0; j < n; ++j)
      newmat.v_[locnew + j] = v_[locold+idx[j]];
  }
}

covMatrix& covMatrix::append(const covMatrix& B) {
  //Allocate space for new one
  unsigned int newsize;
  double *newv;
  newsize = n_ + B.n_;
  newv = new double[ newsize*newsize ];

  unsigned int locold, locnew;
  //Stick in current matrix
  for (unsigned int i = 0; i < n_; ++i) {
    locnew = i*newsize;
    locold = i*n_;
    for (unsigned int j = 0; j < n_; ++j)
      newv[locnew + j] = v_[locold + j];
  }
  //Add B matrix
  for (unsigned int i = 0; i < B.n_; ++i) {
    locnew = (i+n_)*newsize;
    locold = i*n_;
    for (unsigned int j = 0; j < B.n_; ++j)
      newv[locnew + j + n_] = B.v_[locold + j];
  }

  //Set diagonal blocks to zero
  for (unsigned int i = 0; i < n_; ++i)
    for (unsigned int j = n_; j < newsize; ++j)
      newv[i*newsize + j] = 0.0;
  for (unsigned int i = n_; i < newsize; ++i)
    for (unsigned int j = 0; j < n_; ++j)
      newv[i*newsize + j] = 0.0;

  delete[] v_;
  delete[] row_;
  v_ = newv;
  n_ = newsize;
  size_ = newsize * newsize;
  row_ = new double*[size_];
  double* p = v_;              
  for (unsigned int i=0; i<n_; i++) {
    row_[i] = p;
    p += n_;
  }    
  return *this;
}

covMatrix& covMatrix::appendzero(unsigned int n) {
  //Allocate space for new one
  unsigned int newsize;
  double *newv;
  newsize = n_ + n;
  newv = new double[ newsize*newsize ];

  unsigned int locold, locnew;
  //Stick in current matrix
  for (unsigned int i = 0; i < n_; ++i) {
    locnew = i*newsize;
    locold = i*n_;
    for (unsigned int j = 0; j < n_; ++j)
      newv[locnew + j] = v_[locold + j];
  }
  //Add zeros everywhere else
  for (unsigned int i = 0; i < newsize; ++i)
    for (unsigned int j = n_; j < newsize; ++j)
      newv[i*n_+j] = 0.0;
  for (unsigned int i = n_; i < newsize; ++i)
    for (unsigned int j = 0; j < n_; ++j)
      newv[i*n_+j] = 0.0;
  delete[] v_;
  delete[] row_;
  v_ = newv;
  n_ = newsize;
  size_ = newsize * newsize;
  row_ = new double*[size_];
  double* p = v_;              
  for (unsigned int i=0; i<n_; i++) {
    row_[i] = p;
    p += n_;
  }    
  return *this;
}

double covMatrix::getTotal() const {
  if (size_ == 0) return 0.0;
  double sum = v_[0];
  for (unsigned int i = 1; i < size_; ++i)
    sum += v_[i];
  return sum;
}

double covMatrix::getRowTotals( std::vector<double>& vec) const {
  if (size_ == 0) {
    vec.clear();
    return 0.0;
  }
  if (vec.size() != n_) vec.resize(n_);
  double sum, rowsum, *rowptr;
  sum = 0.0;
  for (unsigned int i = 0; i < n_; ++i) {
    rowptr = row_[i];
    rowsum = rowptr[0];
    for (unsigned int j = 1; j < n_; ++j)
      rowsum += rowptr[j];
    sum += rowsum;
    vec[i] = rowsum;
  }
  return sum;
}

covMatrix& covMatrix::operator=(const covMatrix& B) {
  //Self assignment check
  if (v_ == B.v_) return *this;
  if ( n_ == B.n_ ) //Don't need to realloc memory
    copy(B.v_);
  else {
    resize(B.n_);
    copy(B.v_);
  }
  return *this;
}

covMatrix& covMatrix::operator=(double val) {
  set(val);
  return *this;
}

std::vector<double> covMatrix::getDiag() const {
  std::vector<double> diag(n_);
  for (unsigned int i=0; i < n_; ++i)
    diag[i] = row_[i][i];
  return diag;
}

/*!
  Fails if vector is wrong size, but fails silently, so be careful
*/
void covMatrix::setDiag(const std::vector<double>& diag) {
  if ( diag.size() != n_ ) return;
  for (unsigned int i=0; i < n_; ++i)
    row_[i][i] = diag[i];
}

/*!
  Replaces the contents of the current matrix with
  \f$ A \cdot B \f$.

  \param[in] A One of the matricies
  \param[in] B The other matricix
  \param     status zero on successful calculation
*/
covMatrix& covMatrix::mult(const covMatrix&A, const covMatrix&B, int& status) {
  if ( status ) return *this;
  if ( v_ == A.v_ || v_ == B.v_ ) {
    //We can't handle this
    status = 2;
    return *this;
  }

  unsigned int n_A;
  n_A = A.getNRows();
  if ( n_A != B.getNRows() ) {
    status = 1;
    return *this;
  }
  if ( n_ != n_A ) resize( n_A );

#if USEMKL || USEATLAS || USEACCELERATE
  int in = static_cast<int>(n_);
  CBLAS_ORDER order = CblasRowMajor;
  CBLAS_TRANSPOSE trans = CblasNoTrans;
  cblas_dgemm(order,trans,trans,in,in,in,1.0,
	      A.v_,in,B.v_,in,0.0,v_,in );
#else
  double sum;
  const double* row_i;
  const double* col_k;
  for (unsigned int i = 0; i < n_; ++i) {
    row_i = A.row_[i];
    for (unsigned int k = 0; k < n_; ++k) {
      col_k = &(B[0][k]);
      sum = 0;
      for (unsigned int j=0; j< n_; ++j) {
	sum  += *row_i * *col_k;
	row_i++;
	col_k += n_;
      }
      row_[i][k] = sum; 
    }
  }
#endif

  return *this;
}

/*!
  Returns \f$ A \cdot \vec{b}\f$ where A is the current matrix

  \param[in] b The input vector
  \param     status zero on successful calculation

  The arithmetic is done in double precision internally.
*/
std::vector<double> covMatrix::mult( const std::vector<double> b,
				    int& status ) const {
  if ( status ) return std::vector<double>();
  unsigned int n_b;
  n_b = b.size();
  if ( n_b != n_ ) {
    status = 1;
    return std::vector<double>();
  }
  std::vector<double> tmp(n_);
  double sum, *rowptr;
  for (unsigned int i=0; i<n_; ++i) {
    sum = 0;
    rowptr = row_[i];
    for (unsigned int j=0; j<n_; ++j)
      sum +=  rowptr[j] * b[j];
    tmp[i] = sum;
  }
  status = 0;
  return tmp;
}


/*!
  Returns \f$ A \cdot \vec{b}\f$, where A is the current matrix

  \param[in] b      The input vector
  \param     status zero on successful calculation, should be zero
                     on input.
*/
std::vector<float> covMatrix::mult( const std::vector<float> b,
				    int& status ) const {
  if ( status ) return std::vector<float>();
  unsigned int n_b;
  n_b = b.size();
  if ( n_b != n_ ) {
    status = 1;
    return std::vector<float>();
  }
  std::vector<float> tmp(n_);
  double sum, *rowptr;
  for (unsigned int i=0; i<n_; ++i) {
    sum = 0;
    rowptr = row_[i];
    for (unsigned int j=0; j<n_; ++j)
      sum +=  rowptr[j] * b[j];
    tmp[i] = static_cast<float>(sum); 
  }
  status = 0;
  return tmp;
}

covMatrix& covMatrix::scalarMultAndAdd( double scal, const covMatrix&A,
					int& status ) {
  if (status) return *this;
  if ( n_ != A.n_ ) {
    status = 1;
    return *this;
  }
  double *rowptr, *Arowptr;
  for (unsigned int i = 0; i < n_; ++i) {
    rowptr = row_[i];
    Arowptr = A.row_[i];
    for (unsigned int j = 0; j < n_; ++j)
      rowptr[j] += scal * Arowptr[j];
  }
  status = 0;
  return *this;
}

covMatrix& covMatrix::add(const covMatrix& A, int& status) {
  if (status) return *this;
  if ( v_ == A.v_ ) {
    //Self addition
    for (unsigned int i = 0; i < size_; ++i)
      v_[i] += v_[i];
    status = 0;
    return *this;
  }
  if ( A.getNRows() != n_ ) {
    status = 1;
    return *this;
  }
  for (unsigned int i = 0; i < size_; ++i)
    v_[i] += A.v_[i];
  status = 0;
  return *this;
}

covMatrix& covMatrix::operator*=(double val) {
  for (unsigned int i = 0; i < size_; ++i)
    v_[i] *= val;
  return *this;
}

/*!
  Invert using Cholesky decomposition
  \param[out] status Returns 0 on success 
*/
covMatrix& covMatrix::invert(int& status) {
  if (status) return *this;
  
#if USEMKL
  char uplo = 'U';
  int in = static_cast<int>(n_);
  dpotrf( &uplo, &in, v_, &in, &status );
  //std::cerr << "status from dpotrf: " << status << std::endl;
  if (status) return *this;
  dpotri( &uplo, &in, v_, &in, &status );
  if (status) return *this;
  //std::cerr << "status from dpotri: " << status << std::endl;
  //After this, only the lower half is correct.  Normally this
  // would be the upper half, since uplo='U', but CLAPACK
  // thinks we are in Column Major order
  //Anyways, reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }
#elif USEATLAS
  ATLAS_ORDER Order = CblasRowMajor;
  ATLAS_UPLO Uplo = CblasLower;
  status = clapack_dpotrf( Order, Uplo, n_, v_, n_ );
  if (status) return *this;
  status = clapack_dpotri( Order, Uplo, n_, v_, n_ );
  if (status) return *this;
  //Only the lower half is correct -- reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }
#elif USEACCELERATE
  char Uplo = 'U';
  __CLPK_integer N = static_cast<__CLPK_integer>(n_);
  __CLPK_integer STATUS = static_cast<__CLPK_integer>(status);
  dpotrf_(&Uplo, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;
  dpotri_(&Uplo, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }  
#else
  status = 0;
  factorCholesky(status);
  if (status) return *this;
  invertLowerTriangular();
  lowerTriangularTransposeMult();
#endif
  return *this;
}


/*!
  Invert using Cholesky decomposition
  \param[out] status Returns 0 on success 
  \param[out] logdet The log of the determinant
*/
covMatrix& covMatrix::invert(int& status, double& logdet) {
  if (status) return *this;
  
  logdet = 0.0;

#if USEMKL
  char uplo = 'U';
  int in = static_cast<int>(n_);
  dpotrf( &uplo, &in, v_, &in, &status );
  //std::cerr << "status from dpotrf: " << status << std::endl;
  if (status) return *this;

  logdet = log(row_[0][0]);
  for (unsigned int i = 0; i < n_; ++i)
    logdet += log(row_[i][i]);
  logdet *= 2.0;

  dpotri( &uplo, &in, v_, &in, &status );
  if (status) return *this;
  //std::cerr << "status from dpotri: " << status << std::endl;
  //After this, only the lower half is correct.  Normally this
  // would be the upper half, since uplo='U', but CLAPACK
  // thinks we are in Column Major order
  //Anyways, reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }
#elif USEATLAS
  ATLAS_ORDER Order = CblasRowMajor;
  ATLAS_UPLO Uplo = CblasLower;
  status = clapack_dpotrf( Order, Uplo, n_, v_, n_ );
  if (status) return *this;

  logdet = log(row_[0][0]);
  for (unsigned int i = 0; i < n_; ++i)
    logdet += log(row_[i][i]);
  logdet *= 2.0;

  status = clapack_dpotri( Order, Uplo, n_, v_, n_ );
  if (status) return *this;
  //Only the lower half is correct -- reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }
#elif USEACCELERATE
  char Uplo = 'U';
  __CLPK_integer N = static_cast<__CLPK_integer>(n_);
  __CLPK_integer STATUS = static_cast<__CLPK_integer>(status);
  dpotrf_(&Uplo, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;

  logdet = log(row_[0][0]);
  for (unsigned int i = 0; i < n_; ++i)
    logdet += log(row_[i][i]);
  logdet *= 2.0;

  dpotri_(&Uplo, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }  
#else
  status = 0;
  factorCholesky(status);
  if (status) return *this;

  logdet = log(row_[0][0]);
  for (unsigned int i = 0; i < n_; ++i)
    logdet += log(row_[i][i]);
  logdet *= 2.0;

  invertLowerTriangular();
  lowerTriangularTransposeMult();
#endif
  return *this;
}


/*!
  Invert using Cholesky decomposition
  \param[out] status Returns 0 on success
  \param[out] triang Returns the inverse of the lower triangular decomposition
*/
covMatrix& covMatrix::invert(int& status, covMatrix& triang) {
  if (status) return *this;
  
#if USEMKL
  char uplo = 'U';
  char diag = 'N';
  int in = static_cast<int>(n_);
  dpotrf( &uplo, &in, v_, &in, &status );
  if (status) return *this;
  dtrtri( &uplo, &diag, &in, v_, &in, &status );
  if (status) return *this;
  triang = *this;
  dlauum(&uplo, &in, v_, &in, &status );
  if (status) return *this;
  //After this, only the lower half is correct.  Normally this
  // would be the upper half, since uplo='U', but CLAPACK
  // thinks we are in Column Major order
  //Anyways, reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }
#elif USEATLAS
  ATLAS_ORDER Order = CblasRowMajor;
  ATLAS_UPLO Uplo = CblasLower;
  ATLAS_DIAG Diag = CblasNonUnit;
  status = clapack_dpotrf( Order, Uplo, n_, v_, n_ );
  if (status) return *this;
  status = clapack_dtrtri( Order, Uplo, Diag, n_, v_, n_);
  if (status) return *this;
  triang = *this;
  status = clapack_dlauum( Order, Uplo, n_, v_, n_ );
  if (status) return *this;
  //Only the lower half is correct -- reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }
#elif USEACCELERATE
  char Uplo = 'U';
  char diag = 'N';
  __CLPK_integer N = static_cast<__CLPK_integer>(n_);
  __CLPK_integer STATUS = static_cast<__CLPK_integer>(status);

  //Factor
  dpotrf_(&Uplo, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;

  //Invert triangular
  dtrtri_(&Uplo, &diag, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;

  //Get copy of inverted triangular
  triang = *this;

  //Multiply to get full inverse
  dlauum_(&Uplo, &N, v_, &N, &STATUS);
  status = static_cast<int>(STATUS);
  if (status) return *this;

  //Reflect
  double *ptri;
  for (unsigned int i = 0; i < n_; ++i) {
    ptri = row_[i];
    for (unsigned int j = i+1; j < n_; ++j) 
      ptri[j] = v_[ n_*j+i ];
  }  
#else
  status = 0;
  factorCholesky(status);
  triang = *this;
  if (status) return *this;
  invertLowerTriangular();
  lowerTriangularTransposeMult();
#endif
  return *this;
}

covMatrix& covMatrix::operator+=(const covMatrix& A) {
  int status;
  status = 0;
  return add( A, status );
}

std::ostream& operator<<(std::ostream &s, const covMatrix& A) {
  unsigned int N = A.getNRows();
  double const* rowptr;
  s << N << std::endl;
  for (unsigned int i=0; i<N; ++i) {
    rowptr = A[i];
    for (unsigned int j=0; j<N; ++j)
      s << rowptr[j] << " ";
    s << std::endl;
  }
  return s;
}

std::istream& operator>>(std::istream &s, covMatrix &A) {
  unsigned int N;
  s >> N;

  if ( !( N == A.getNRows() ) )
    A.resize(N);

  double *rowptr;
  for (unsigned int i=0; i<N; ++i) {
    rowptr = A[i];
    for (unsigned int j=0; j<N; ++j)
      s >> rowptr[j];
  }
  return s;
}

#ifdef COVFITS
covMatrix& covMatrix::readFits(const std::string& FileName) {
  int status=0;
  fitsfile *fptr;

  fits_open_file(&fptr, FileName.c_str(), READONLY, &status);
  if (status) {
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    clear();
    return *this;
  }

  int naxis;
  long naxes[2];
  fits_get_img_dim(fptr, &naxis, &status);
  fits_get_img_size(fptr, 2, naxes, &status);
  if (status || naxis != 2) 
    throw CosFitterExcept("covMatrix","readFits",
			  "Input fits file is not 2D",1);
  if (naxes[0] != naxes[1])
    throw CosFitterExcept("covMatrix","readFits",
			  "Input fits file is not square",2);
  resize(naxes[0]);

  //Now read the actual data
  long fpixel[2];
  fpixel[0] = 1;
  for (unsigned int i = 1; i <= static_cast<unsigned int>(naxes[1]); ++i) {
    fpixel[1] = i;
    if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], 0, 
		      v_ +  naxes[0]*(i-1), 0, &status))
      break;   /* jump out of loop on error */
  }
  fits_close_file(fptr,&status);
  if (status) {
    fits_report_error(stderr,status);
    clear();
    return *this;
  } 
  return *this;
}
#endif

covMatrix operator*( const covMatrix&A, const covMatrix& B) {
  int status;
  status = 0;
  covMatrix tmp( A.getNRows() );
  return tmp.mult(A,B,status);
}

std::vector<double> operator*(const covMatrix&A, const std::vector<double>&b) {
  int status;
  status = 0;
  return A.mult( b, status );
}

/////////////////////////////////////////////////////////////////

void covMatrixSNe::readCovData(const std::string& FileName,
			       covMatrix& covmat) {

#ifdef COVFITS
  bool have_fits;
  have_fits = false;
  //See if we have a fits file
  std::string::size_type dotpos=FileName.rfind(".");
  if (dotpos == std::string::npos) have_fits=false; else
    if ( FileName.substr(dotpos) == ".fits" ) have_fits=true; else
      have_fits=false;
  
  if (have_fits) {
    covmat.readFits(FileName);
    return;
  }
#endif
  std::ifstream fl(FileName.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Data file " << FileName << " not found";
    throw CosFitterExcept("covMatrixSNe","readCovData",errstrng.str(),2);
  }
  fl >> covmat;
  fl.close();

}


/*!
  Writes a covaraince matrix.  Right now this just uses the overloaded
  operator<< for covMatrix
  \param[in] FileName File to write to
  \param[in] covmat  Covariance matrix to write
*/
void covMatrixSNe::writeCovData(const std::string& FileName, 
			   const covMatrix& covmat) const {
  std::ofstream fl(FileName.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Can't open " << FileName << " to write";
    throw CosFitterExcept("covMatrixSNe","writeCovData",errstrng.str(),2);
  }
  fl << covmat;
  fl.close();
}

covMatrixSNe::covMatrixSNe() : nsn_(0), haveMagCov(false), haveWidthCov(false),
			       haveColourCov(false), haveMagWidthCov(false),
			       haveMagColourCov(false), 
			       haveWidthColourCov(false){ }

covMatrixSNe::covMatrixSNe(const std::string& file0, 
			   const std::string& filea, 
			   const std::string& fileb, 
			   const std::string& file0a, 
			   const std::string& file0b, 
			   const std::string& fileab) {
  readMagCovData(file0);
  readWidthCovData(filea);
  readColourCovData(fileb);
  readMagWidthCovData(file0a);
  readMagColourCovData(file0b);
  readWidthColourCovData(fileab);
}

/*!
  New covmatric are set to zero.
 */
covMatrixSNe& 
covMatrixSNe::matchCovmats( const covMatrixSNe& B, unsigned int n ) {
  if ( haveAnyCovMatrix() && n != nsn_ )
    throw CosFitterExcept("covMatrixSNe","matchCovmats",
			  "Covmat has entries and requested size doesn't match",1);
  nsn_ = n;
  if ( B.haveMagCov ) {
    haveMagCov = true;
    mag_covmatrix.resize(n,0.0);
  }
  if ( B.haveWidthCov ) {
    haveWidthCov = true;
    width_covmatrix.resize(n,0.0);
  }
  if ( B.haveColourCov ) {
    haveColourCov = true;
    colour_covmatrix.resize(n,0.0);
  }
  if ( B.haveMagWidthCov ) {
    haveMagWidthCov = true;
    magwidth_covmatrix.resize(n,0.0);
  }
  if ( B.haveMagColourCov ) {
    haveMagColourCov = true;
    magcolour_covmatrix.resize(n,0.0);
  }
  if ( B.haveWidthColourCov ) {
    haveWidthColourCov = true;
    widthcolour_covmatrix.resize(n,0.0);
  }
  return *this;
}

covMatrixSNe& covMatrixSNe::resize( unsigned int n ) {
  if (haveMagCov) mag_covmatrix.resize(n);
  if (haveWidthCov) width_covmatrix.resize(n);
  if (haveColourCov) colour_covmatrix.resize(n);
  if (haveMagWidthCov) magwidth_covmatrix.resize(n);
  if (haveMagColourCov) magcolour_covmatrix.resize(n);
  if (haveWidthColourCov) widthcolour_covmatrix.resize(n);
  nsn_ = n;
  return *this;
}

covMatrixSNe& covMatrixSNe::reorder(const std::vector<unsigned int>& idx ) {
  if (haveMagCov) mag_covmatrix.reorder(idx);
  if (haveWidthCov) width_covmatrix.reorder(idx);
  if (haveColourCov) colour_covmatrix.reorder(idx);
  if (haveMagWidthCov) magwidth_covmatrix.reorder(idx);
  if (haveMagColourCov) magcolour_covmatrix.reorder(idx);
  if (haveWidthColourCov) widthcolour_covmatrix.reorder(idx);
  return *this;
}

void covMatrixSNe::select( covMatrixSNe& in, 
			   const std::vector<unsigned int>& idx ) const {
  in.resize( idx.size() );
  in.matchCovmats( *this, idx.size() );
  if (haveMagCov) mag_covmatrix.select(in.mag_covmatrix,idx);
  if (haveWidthCov) width_covmatrix.select(in.width_covmatrix,idx);
  if (haveColourCov) colour_covmatrix.select(in.colour_covmatrix,idx);
  if (haveMagWidthCov) magwidth_covmatrix.select(in.magwidth_covmatrix,idx);
  if (haveMagColourCov) magcolour_covmatrix.select(in.magcolour_covmatrix,idx);
  if (haveWidthColourCov) 
    widthcolour_covmatrix.select(in.widthcolour_covmatrix,idx);
}

covMatrixSNe& covMatrixSNe::append(const covMatrixSNe& ap) {
  //This is slightly tricker -- we have to account for the possibility
  // that either this one or the one we are appending doesn't
  // have a given cov matrix
  unsigned int n = ap.nsn_;

  if (haveMagCov) {
    if (!ap.haveMagCov) mag_covmatrix.appendzero(n); else 
      mag_covmatrix.append(ap.mag_covmatrix);
  } else if (ap.haveMagCov) {
    mag_covmatrix.resize( nsn_, 0.0 );
    haveMagCov = true;
    mag_covmatrix.append(ap.mag_covmatrix);
  }

  if (haveWidthCov) {
    if (!ap.haveWidthCov) width_covmatrix.appendzero(n); else 
      width_covmatrix.append(ap.width_covmatrix);
  } else if (ap.haveWidthCov) {
    width_covmatrix.resize( nsn_, 0.0 );
    haveWidthCov = true;
    width_covmatrix.append(ap.width_covmatrix);
  }

  if (haveColourCov) {
    if (!ap.haveColourCov) colour_covmatrix.appendzero(n); else 
      colour_covmatrix.append(ap.colour_covmatrix);
  } else if (ap.haveColourCov) {
    colour_covmatrix.resize( nsn_, 0.0 );
    haveColourCov = true;
    colour_covmatrix.append(ap.colour_covmatrix);
  }

  if (haveMagWidthCov) {
    if (!ap.haveMagWidthCov) magwidth_covmatrix.appendzero(n); else 
      magwidth_covmatrix.append(ap.magwidth_covmatrix);
  } else if (ap.haveMagWidthCov) {
    magwidth_covmatrix.resize( nsn_, 0.0 );
    haveMagWidthCov = true;
    magwidth_covmatrix.append(ap.magwidth_covmatrix);
  }

  if (haveMagColourCov) {
    if (!ap.haveMagColourCov) magcolour_covmatrix.appendzero(n); else 
      magcolour_covmatrix.append(ap.magcolour_covmatrix);
  } else if (ap.haveMagColourCov) {
    magcolour_covmatrix.resize( nsn_, 0.0 );
    haveMagColourCov = true;
    magcolour_covmatrix.append(ap.magcolour_covmatrix);
  }

  if (haveWidthColourCov) {
    if (!ap.haveWidthColourCov) widthcolour_covmatrix.appendzero(n); else 
      widthcolour_covmatrix.append(ap.widthcolour_covmatrix);
  } else if (ap.haveWidthColourCov) {
    widthcolour_covmatrix.resize( nsn_, 0.0 );
    haveWidthColourCov = true;
    widthcolour_covmatrix.append(ap.widthcolour_covmatrix);
  }

  return *this;
}


covMatrixSNe& covMatrixSNe::appendzero(unsigned int n) {
  if (haveMagCov) mag_covmatrix.appendzero(n);
  if (haveWidthCov) width_covmatrix.appendzero(n);
  if (haveColourCov) colour_covmatrix.appendzero(n);
  if (haveMagWidthCov) magwidth_covmatrix.appendzero(n);
  if (haveMagColourCov) magcolour_covmatrix.appendzero(n);
  if (haveWidthColourCov) widthcolour_covmatrix.appendzero(n);
  return *this;
}

/*!
  \param[in] FileNames Vector of string filenames, matched by position
    to the six cov matricies.  An empty string will result in a skipped read.
*/
void covMatrixSNe::readCovData(const std::vector<std::string>& FileNames) {
  unsigned int nentries;
  nentries = FileNames.size();
  if ( nentries > 0 && FileNames[0] != "" )
    readMagCovData( FileNames[0] );
  if ( nentries > 1 && FileNames[1] != "" )
    readWidthCovData( FileNames[1] );
  if ( nentries > 2 && FileNames[2] != "" )
    readColourCovData( FileNames[2] );
  if ( nentries > 3 && FileNames[3] != "" )
    readMagWidthCovData( FileNames[3] );
  if ( nentries > 4 && FileNames[4] != "" )
    readMagColourCovData( FileNames[4] );
  if ( nentries > 5 && FileNames[5] != "" )
    readWidthColourCovData( FileNames[5] );
}


/*!
  Reads in the magnitude covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void covMatrixSNe::readMagCovData(const std::string& FileName) { 
  readCovData( FileName, mag_covmatrix );
  unsigned int n;
  n = mag_covmatrix.getNRows();
  if ( n != 0 ) {
    if (nsn_ != 0 && n != nsn_ )
      throw CosFitterExcept("covMatrixSNe","readMagCovData",
			    "new mag cov is not same size as pre-existing data",1);
    nsn_ = n;
    haveMagCov = true;
  }
}

/*!
  Reads in the width covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void covMatrixSNe::readWidthCovData(const std::string& FileName) { 
  readCovData( FileName, width_covmatrix );
  unsigned int n;
  n = width_covmatrix.getNRows();
  if ( n != 0 ) {
    if (nsn_ != 0 && n != nsn_ )
      throw CosFitterExcept("covMatrixSNe","readWidthCovData",
			    "new width cov is not same size as pre-existing data",1);
    nsn_ = n;
    haveWidthCov = true;
  }
}

/*!
  Reads in the colour covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void covMatrixSNe::readColourCovData(const std::string& FileName) { 
  readCovData( FileName, colour_covmatrix );
  unsigned int n;
  n = colour_covmatrix.getNRows();
  if ( n != 0 ) {
    if (nsn_ != 0 && n != nsn_ )
      throw CosFitterExcept("covMatrixSNe","readColourCovData",
			    "new colour cov is not same size as pre-existing data",1);
    nsn_ = n;
    haveColourCov = true;
  }
}

/*!
  Reads in the mag-width covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void covMatrixSNe::readMagWidthCovData(const std::string& FileName) { 
  readCovData( FileName, magwidth_covmatrix );
  unsigned int n;
  n = magwidth_covmatrix.getNRows();
  if ( n != 0 ) {
    if (nsn_ != 0 && n != nsn_ )
      throw CosFitterExcept("covMatrixSNe","readMagWidthCovData",
			    "new mag-width cov is not same size as pre-existing data",1);
    nsn_ = n;
    haveMagWidthCov = true;
  }
}


/*!
  Reads in the mag-colour covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void covMatrixSNe::readMagColourCovData(const std::string& FileName) { 
  readCovData( FileName, magcolour_covmatrix );
  unsigned int n;
  n = magcolour_covmatrix.getNRows();
  if ( n != 0 ) {
    if (nsn_ != 0 && n != nsn_ )
      throw CosFitterExcept("covMatrixSNe","readMagColourCovData",
			    "new mag-colour cov is not same size as pre-existing data",1);
    nsn_ = n;
    haveMagColourCov = true;
  }
}

/*!
  Reads in the width-colour covariance matrix.  Right now this just uses the 
  overloaded operator>> for covMatrix.  Replaces any existing one.
  \param[in] FileName The file to read the matrix from
*/
void covMatrixSNe::readWidthColourCovData(const std::string& FileName) { 
  readCovData( FileName, widthcolour_covmatrix );
  unsigned int n;
  n = widthcolour_covmatrix.getNRows();
  if ( n != 0 ) {
    if (nsn_ != 0 && n != nsn_ )
      throw CosFitterExcept("covMatrixSNe","readWidthColourCovData",
			    "new width-colour cov is not same size as pre-existing data",1);
    nsn_ = n;
    haveWidthColourCov = true;
  }
}

/*!
  \param[in] FileNames Vector of string filenames, matched by position
    to the six cov matricies.  An empty string will result in a skipped write.
*/
void covMatrixSNe::writeCovData(const std::vector<std::string>& FileNames) 
  const { 
  unsigned int nentries;
  nentries = FileNames.size();
  if ( haveMagCov && nentries > 0 && FileNames[0] != "" )
    writeCovData( FileNames[0], mag_covmatrix );
  if ( haveMagCov && nentries > 1 && FileNames[1] != "" )
    writeCovData( FileNames[1], width_covmatrix );
  if ( haveMagCov && nentries > 2 && FileNames[2] != "" )
    writeCovData( FileNames[2], colour_covmatrix );
  if ( haveMagCov && nentries > 3 && FileNames[3] != "" )
    writeCovData( FileNames[3], magwidth_covmatrix );
  if ( haveMagCov && nentries > 4 && FileNames[4] != "" )
    writeCovData( FileNames[4], magcolour_covmatrix );
  if ( haveMagCov && nentries > 5 && FileNames[5] != "" )
    writeCovData( FileNames[5], widthcolour_covmatrix );
}


/*!
  Writes the mag covariance matrix.  
*/
void covMatrixSNe::writeMagCovData(const std::string& FileName) const { 
  if (!haveMagCov) return;
  writeCovData(FileName,mag_covmatrix);
}

/*!
  Writes the width covariance matrix.  
*/
void covMatrixSNe::writeWidthCovData(const std::string& FileName) const { 
  if (!haveWidthCov) return;
  writeCovData(FileName,width_covmatrix);
}

/*!
  Writes the colour covariance matrix.  
*/
void covMatrixSNe::writeColourCovData(const std::string& FileName) const { 
  if (!haveColourCov) return;
  writeCovData(FileName,colour_covmatrix);
}

/*!
  Writes the mag-width covariance matrix.  
*/
void covMatrixSNe::writeMagWidthCovData(const std::string& FileName) const { 
  if (!haveMagWidthCov) return;
  writeCovData(FileName,magwidth_covmatrix);
}

/*!
  Writes the mag-colour covariance matrix.  
*/
void covMatrixSNe::writeMagColourCovData(const std::string& FileName) const { 
  if (!haveMagColourCov) return;
  writeCovData(FileName,magcolour_covmatrix);
}

/*!
  Writes the width-colour covariance matrix.  
*/
void covMatrixSNe::writeWidthColourCovData(const std::string& FileName) const { 
  if (!haveWidthColourCov) return;
  writeCovData(FileName,widthcolour_covmatrix);
}

/*!
  Use getMagCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of mag_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix covMatrixSNe::getMagCovMatrix() const {
  if (!haveMagCov) return covMatrix();
  return mag_covmatrix;
}

double covMatrixSNe::getMagCovElement(unsigned int i, unsigned int j) const {
  if (!haveMagCov) return 0.0;
  return mag_covmatrix[i][j];
}

/*!
  Use getWidthCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of width_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix covMatrixSNe::getWidthCovMatrix() const {
  if (!haveWidthCov) return covMatrix();
  return width_covmatrix;
}

double covMatrixSNe::getWidthCovElement(unsigned int i, unsigned int j) const {
  if (!haveWidthCov) return 0.0;
  return width_covmatrix[i][j];
}

/*!
  Use getColourCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of colour_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix covMatrixSNe::getColourCovMatrix() const {
  if (!haveColourCov) return covMatrix();
  return colour_covmatrix;
}

double 
covMatrixSNe::getColourCovElement(unsigned int i, unsigned int j) const {
  if (!haveColourCov) return 0.0;
  return colour_covmatrix[i][j];
}

/*!
  Use getMagWidthCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of magwidth_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix covMatrixSNe::getMagWidthCovMatrix() const {
  if (!haveMagWidthCov) return covMatrix();
  return magwidth_covmatrix;
}

double 
covMatrixSNe::getMagWidthCovElement(unsigned int i, unsigned int j) const {
  if (!haveMagWidthCov) return 0.0;
  return magwidth_covmatrix[i][j];
}

/*!
  Use getMagColourCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of magcolour_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix covMatrixSNe::getMagColourCovMatrix() const {
  if (!haveMagColourCov) return covMatrix();
  return magcolour_covmatrix;
}

double 
covMatrixSNe::getMagColourCovElement(unsigned int i, unsigned int j) const {
  if (!haveMagColourCov) return 0.0;
  return magcolour_covmatrix[i][j];
}

/*!
  Use getWidthColourCovmatrixRef to get a reference if you want to avoid
  copying.  This is a lot safer, but much slower.

  \return A copy of widthcolour_covmatrix.  If none is present, an empty
    matrix is returned.
*/
covMatrix covMatrixSNe::getWidthColourCovMatrix() const {
  if (!haveWidthColourCov) return covMatrix();
  return widthcolour_covmatrix;
}

double 
covMatrixSNe::getWidthColourCovElement(unsigned int i, unsigned int j) const {
  if (!haveWidthColourCov) return 0.0;
  return widthcolour_covmatrix[i][j];
}

/*!
  \param[out] cov  The combined covariance matrix
  \param[in] alpha  \f$\alpha\f$
  \param[in] beta   \f$\beta\f$

*/
void covMatrixSNe::getCombinedCovMatrix(covMatrix& cov,
					double alpha,double beta) const {
  int status = 0;
  if (cov.getNRows() != nsn_) cov.resize(nsn_);

  if (!haveAnyCovMatrix()) {
    //Return a zero cov matrix
    cov = 0;
    return;
  }
  if (haveMagCov) cov = mag_covmatrix; else 
    cov = covMatrix(nsn_,0.0);
  if (haveWidthCov && alpha != 0.0) 
    cov.scalarMultAndAdd( alpha * alpha, width_covmatrix, status );
  if (haveColourCov && beta != 0.0 ) 
    cov.scalarMultAndAdd( beta * beta, colour_covmatrix, status );
  if (haveMagWidthCov && alpha != 0.0) 
    cov.scalarMultAndAdd( 2.0 * alpha, magwidth_covmatrix, status );
  if (haveMagColourCov && beta != 0.0) 
    cov.scalarMultAndAdd( - 2.0 * beta, magcolour_covmatrix, status );
  if (haveWidthColourCov && alpha != 0.0 && beta != 0.0) 
    cov.scalarMultAndAdd( - 2.0 * alpha * beta, widthcolour_covmatrix, 
			     status );
  if (status) throw CosFitterExcept("covMatrixSNe","getCombinedCovMatrix",
				    "Error combining cov matricies",1);
}

/*!
  \param[in] alpha \f$\alpha\f$
  \param[in] beta  \f$\beta\f$
  \param[in] i     Element index
  \param[in] j     Element index
  \return The corresponding element of the combined covariance matrix

  //This doesn't work if you have a woodbury cov matrix
 */
double 
covMatrixSNe::getCombinedCovElement(double alpha, double beta,
				    unsigned int i,unsigned int j) const 
{
  if (! haveAnyCovMatrix() ) return 0.0;

  double val;
  val = 0.0;
  if (haveMagCov) val += mag_covmatrix[i][j];
  if (haveWidthCov) val += alpha * width_covmatrix[i][j];
  if (haveColourCov)  val += beta * colour_covmatrix[i][j];
  if (haveMagWidthCov) val += 2.0 * alpha * magwidth_covmatrix[i][j];
  if (haveMagColourCov) val -= 2.0 * beta * magcolour_covmatrix[i][j];
  if (haveWidthColourCov) 
    val -= 2.0 * alpha * beta * widthcolour_covmatrix[i][j];
  return val;
}
