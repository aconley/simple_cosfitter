//covmatrix.cc

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

//Linear algebra stuff
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

#else

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#endif

#include <covmatrix_woodbury.h>
#include <cosfitterexcept.h>


//We avoid exceptions and most run-time error checking for 
//performance reasons, so be careful.
/*!
  Sets all the elements of the vectors to zero.
*/
void covMatrixWoodbury::initialize(unsigned int n, unsigned int m) {
  if (n == 0 || m == 0) {
    destroy();
    return;
  }

  //Deal with the diagonal
  if ( n < n_ ) {
    diag_.clear();
  } 
  if ( n != n_ ) {
    diag_.resize( n );
  }
  //Clear values
  diag_.assign( n, 0.0 );

  //Now deal with U and the working matrices
  if ( n != n_  || m != m_ ) {
#if USEMKL || USEATLAS || USEACCELERATE
    if (U != NULL) delete[] U;
    if (working != NULL) delete[] working;
    if (working2 != NULL) delete[] working2;

    U = new double[ n * m ];
    for (unsigned int i = 0; i < n*m; ++i) U[i]=0.0;
    working = new double[ m * n ];
    working2 = new double[ n * m ];
#else
    if (U != NULL) gsl_matrix_free( U );
    if (working != NULL) gsl_matrix_free( working );
    if (working2 != NULL) gsl_matrix_free( working2 );
    if (n != n_ && inverse != NULL) gsl_matrix_free( inverse );

    U = gsl_matrix_alloc( n, m );
    gsl_matrix_set_zero(U);
    working = gsl_matrix_alloc( m, n );
    working2 = gsl_matrix_alloc( n, m );
    if (n != n_) inverse = gsl_matrix_alloc( n, n );
#endif
  }

  if ( m != m_ ) {
#if USEMKL || USEATLAS 
    if ( workingsq != NULL ) delete[] workingsq;
    if ( ipiv != NULL ) delete[] ipiv;
    workingsq = new double[ m * m ];
    ipiv = new int[ m ];
#elif USEACCELERATE
    if ( workingsq != NULL ) delete[] workingsq;
    if ( ipiv != NULL ) delete[] ipiv;
    workingsq = new double[ m * m ];
    ipiv = new __CLPK_integer[ m ];
#else
    if ( workingsq != NULL ) gsl_matrix_free( workingsq );
    if ( ipiv != NULL ) gsl_permutation_free( ipiv );
    workingsq = gsl_matrix_alloc( m, m );
    gsl_permutation_alloc( m );
#endif
  }

  n_ = n;
  m_ = m;
}

void covMatrixWoodbury::destroy() {
  diag_.clear();

#if USEMKL || USEATLAS || USEACCELERATE 
  if (U != NULL) delete[] U;
  if (working != NULL) delete[] working;
  if (working2 != NULL) delete[] working2;
  if (workingsq != NULL) delete[] workingsq;
  if (ipiv != NULL) delete[] ipiv;
#else
  if (U != NULL) gsl_matrix_free( U );
  if (working != NULL) gsl_matrix_free( working );
  if (working2 != NULL) gsl_matrix_free( working2 );
  if (workingsq != NULL) gsl_matrix_free( workingsq );
  if (inverse != NULL) gsl_matrix_free( inverse );
  if (ipiv != NULL) gsl_permutation_free( ipiv );
  inverse = NULL;
#endif

  U = working = working2 = workingsq = NULL;
  ipiv = NULL;

  n_ = 0;
  m_ = 0;
}


covMatrixWoodbury::covMatrixWoodbury() {
  n_ = 0;
  m_ = 0;
  U = working = working2 = workingsq = NULL;  
#if ! (USEMKL || USEATLAS || USEACCELERATE)
  inverse = NULL;
#endif
  ipiv = NULL;

}

/*!
  \param[in] A  Matrix to copy
*/
covMatrixWoodbury::covMatrixWoodbury(const covMatrixWoodbury& A) {
  U = working = working2 = workingsq = NULL;
#if ! (USEMKL || USEATLAS || USEACCELERATE)
  inverse = NULL;
#endif
  ipiv = NULL;
  initialize(A.n_,A.m_);
  std::copy( A.diag_.begin(), A.diag_.end(), diag_.begin() );
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < n_ * m_; ++i) U[i] = A.U[i];
#else
  gsl_matrix_memcpy (U, A.U);
#endif
}

covMatrixWoodbury::covMatrixWoodbury(unsigned int n, unsigned int m,
				     double diagval) {
  U = working = working2 = workingsq = NULL;
#if ! (USEMKL || USEATLAS || USEACCELERATE)
  inverse = NULL;
#endif
  ipiv = NULL;
  initialize(n,m);
  setDiag( diagval );
}

/*!
  We don't check to make sure the vectors are the same size,
  so be careful.  The input u should be in column major order
  (i.e., Fortran convention)
*/
covMatrixWoodbury::covMatrixWoodbury(const std::vector<double>& diag,
				     unsigned int m, 
				     const double *const u) {
  U = working = working2 = workingsq = NULL;
  ipiv = NULL;
  initialize( diag.size(), m );
  std::copy( diag.begin(), diag.end(), diag_.begin() );
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < n_ * m_; ++i) U[i] = u[i];
#else
  inverse = NULL;
  for (unsigned int j = 0; j < m_; ++j)
    for (unsigned int i = 0; i < n_; ++i)
      gsl_matrix_set( U, i, j, u[i + j*n_] );  //Column major
#endif
}

/*!
  \param[in] d Vector to set diagonal to

  This is not checked for performance reasons, so make sure
  that the input vector is the right size...
*/
void covMatrixWoodbury::setDiag( const std::vector<double> d ) {
  std::copy( d.begin(), d.end(), diag_.begin() );
}

/*!
  This doesn't propogate boundary checks for performance reasons, so be careful
*/
void covMatrixWoodbury::setUColumn(unsigned int m, 
				   const std::vector<double>& invec) {
  if ( n_ == 0 || m_ == 0 ) return;
  if ( m < 0 || m > m_ ) return;
  if ( invec.size() != n_ ) return;
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < n_; ++i)
    U[ i + m*n_ ] = invec[i];
#else
  for (unsigned int i = 0; i < n_; ++i)
    gsl_matrix_set(U,i,m,invec[i]);
#endif
}

/*!
  This doesn't propogate boundary checks for performance reasons, so be careful
*/
void covMatrixWoodbury::getUColumn(unsigned int m, 
				   std::vector<double>& invec) const {
  if ( n_ == 0 || m_ == 0 ) return;
  if ( m < 0 || m > m_ ) return;
  invec.resize(n_);
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < n_; ++i)
    invec[i] = U[ i + m*n_ ];
#else
  for (unsigned int i = 0; i < n_; ++i)
    invec[i] = gsl_matrix_get(U,i,m);
#endif
}


/*!
  This operations occurs in place, i.e. when resizing to 
  a new matrix, original matrix elements
  are <b>NOT</b> retained.  Instead, one must explicit create
  a new matrix of this size and manually copy the elements.
  However, if you tell it to resize to the same size, the
  old elements are left alone.
*/
covMatrixWoodbury& covMatrixWoodbury::resize(unsigned int N, unsigned int M) {
  if ( n_ == N && m_ == M ) {
    return *this;
  }
  initialize(N,M);
  return *this;
}

covMatrixWoodbury& covMatrixWoodbury::resize(unsigned int N, unsigned int M,
					     double diagval) {
  if ( n_ == N && m_ == M ) {
    setDiag( diagval );
#if USEMKL || USEATLAS || USEACCELERATE
    for (unsigned int i = 0; i < n_ * m_; ++i) U[i] = 0.0;
#else
    gsl_matrix_set_zero(U);
#endif
    return *this;
  }
  initialize(N,M);
  setDiag( diagval );
  return *this;
}

/*!
  \param[in] idx  Index array to reorder with.  The boundaries are
                  not checked, nor are repeats.  Must be of length n_.
*/
void covMatrixWoodbury::reorder( const std::vector<unsigned int>& idx ) {
  if (idx.size() != n_)
    throw CosFitterExcept("covMatrixWoodbury","reorder",
			  "index array not same size as data arrays",1);

  std::vector<double> working;
  working.resize(n_);

  //Do the diagonal
  for (unsigned int i = 0; i < n_; ++i)
    working[i] = diag_[ idx[i] ];
  diag_ = working;

  //Now do U, which is less fun.  We are reordering the first index
#if USEMKL || USEATLAS || USEACCELERATE
  double *Uworking;
  Uworking = new double[ n_ * m_ ];
  for (unsigned int j = 0; j < m_; ++j)
    for (unsigned int i = 0; i < n_; ++i)
      Uworking[ i + j*n_ ] = U[ idx[i] + j*n_ ];
  delete[] U;
  U = Uworking;
#else
  //We can't use swaps because uniqueness is not required for the reorder
  gsl_vector *vworking;
  gsl_matrix *Uworking;
  vworking = gsl_vector_alloc( m_ );
  Uworking = gsl_matrix_alloc( n_, m_ );
  for (unsigned int i = 0; i < n_; ++i) {
    gsl_matrix_get_row(vworking, U, idx[i]);
    gsl_matrix_set_row(Uworking,i,vworking);
  }
  gsl_matrix_free(U);
  U = Uworking;
  gsl_vector_free( vworking );
#endif 
}

/*!
  \param[out] newmat Matrix with output information
  \param[in] idx  Index array to select with.  The boundaries are
                  not checked
*/
void covMatrixWoodbury::select( covMatrixWoodbury& newmat,
				const std::vector<unsigned int>& idx ) const {
  unsigned int n = idx.size();
  if (newmat.n_ != n || newmat.m_ != m_ )
    newmat.resize( n, m_ );

  //Do the diagonal
  for (unsigned int i = 0; i < n; ++i)
    newmat.diag_[i] = diag_[ idx[i] ];
  
    //Now do U, which is less fun.  We are selecting on the first index
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int j = 0; j < m_; ++j)
    for (unsigned int i = 0; i < n; ++i)
      newmat.U[ i + j*n ] = U[ idx[i] + j*n_ ];
#else
  //We can't use swaps because uniqueness is not required for the reorder
  gsl_vector *vworking;
  vworking = gsl_vector_alloc( m_ );
  for (unsigned int i = 0; i < n; ++i) {
    gsl_matrix_get_row(vworking, U, idx[i]);
    gsl_matrix_set_row(newmat.U,i,vworking);
  }
  gsl_vector_free( vworking );
#endif 

}

void covMatrixWoodbury::getCovMatrix( covMatrix& cov ) const {
  if ( cov.getNRows() != n_ ) cov.resize( n_ );
  cov = 0.0;

#if USEMKL || USEATLAS || USEACCELERATE
  int in = static_cast<int>(n_);
  int im = static_cast<int>(m_);

  cov.setDiag( diag_ );
  cblas_dgemm( CblasColMajor, CblasNoTrans, CblasTrans, in, in, im, 1.0,
	       U, in, U, in, 1.0, cov.getData(), in );
#else
  gsl_matrix *covmat;
  covmat = gsl_matrix_alloc( n_, n_ );
  for (unsigned int i=0; i < n_; ++i)
    gsl_matrix_set( covmat, i, i, diag_[i] );  //Column major

  gsl_blas_dgemm(CblasNoTrans, CblasTrans,
		 1.0, U, U, 1.0, covmat );
  //Copy into covmatrix
  double *covdata = cov.getData();
  for (unsigned int i=0; i < n_; ++i)
    for (unsigned int j=0; j < n_; ++j)
      covdata[ i + j*n_ ] = gsl_matrix_get( covmat, i, j );
  gsl_matrix_free(covmat);
  for (unsigned int i = 0; i < n_; ++i)
    cov[i][i] += diag_[i];
#endif



}

void covMatrixWoodbury::getInvCovMatrix( covMatrix& invcov, 
					 int& status) const {
  if (status) return;

  if ( invcov.getNRows() != n_ ) invcov.resize( n_ );

#if USEMKL || USEATLAS || USEACCELERATE
  //We need int versions for calls, since they don't accept
  // unsigned ints
  int in, im;
  in = static_cast<int>(n_);
  im = static_cast<int>(m_);

  //Convenience pointer
  double *dptr, *dptr2;
#else
  double *elemptr;
#endif

  //Invert the diagonal
  std::vector<double> invdiag( n_ );
  for (unsigned int i = 0; i < n_; ++i)
    invdiag[i] = 1.0 / diag_[i];

  //Now form the product U^T * D^-1 and store in the matrix
  // working.  This is relatively easy becauce D is, well, diagonal
  // D has n by n elements, U is n by m, so U^T is m by n
  // and the product is therefore m by n
  //std::cerr << "Forming U^T D^-1" << std::endl;
#if USEMKL || USEATLAS || USEACCELERATE
  //for (unsigned int i = 0; i < m_; ++i)
  //  for (unsigned int j = 0; j < n_; ++j)
  //    working[ i + j*m_ ] = invdiag[j] * U[ j + i*n_ ];
  for (unsigned int i = 0; i < m_; ++i) {
    dptr = U + i * n_;
    for (unsigned int j = 0; j < n_; ++j) {
      //working[i,j] = invdiag[j] * U[j,i]
      working[ i + j * m_ ] = invdiag[j] * dptr[j];
    }
  }
#else
  for (unsigned int j = 0; j < m_; ++j)
    for (unsigned int k = 0; k < n_; ++k) 
      gsl_matrix_set( working, j, k, invdiag[k] * 
		      gsl_matrix_get( U, k, j ) );
#endif
  /*
  std::cerr << "U^T D^-1" << std::endl;
  for (unsigned int i = 0; i < m_; ++i) {
    for (unsigned int j = 0; j < n_-1; ++j)
      std::cout << working[i+j*m_] << " ";
    std::cout << working[i+(n_-1)*m_] << std::endl;
  }
  */


  //Next form U^T * D^-1 * U and store in workingsq (m_ x m_)
  // U^T * D^-1 is m_ x n_ (stored in working)
  // U is n_ x m_
  // We already calculated U^T D^-1, so just right multiply by U
  //Thus, for gemm, m = m_, n = m_, k = n_
  //std::cerr << "Forming U^T D^-1 U" << std::endl;
#if USEMKL || USEATLAS || USEACCELERATE
  //CBLAS version
  cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
	       im, im, in, 1.0, working, im, U, in, 0.0,
	       workingsq, im );
#else
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, working, U,
		  0.0, workingsq );
#endif

  //Add the identity matrix to form 1 + U^T * D^-1 * U
  //in workingsq
  //std::cerr << "Adding identity" << std::endl;
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < m_; ++i)
    workingsq[i*m_+ i] += 1.0;
#else
  for (unsigned int i = 0; i < m_; ++i) {
    elemptr = gsl_matrix_ptr( workingsq, i, i );
    *elemptr += 1.0;
  }
#endif
  /*
  std::cerr << "1 + U^T D^-1 U" << std::endl;
  for (unsigned int i = 0; i < m_; ++i) {
    for (unsigned int j = 0; j < m_-1; ++j)
      std::cout << workingsq[i+j*m_] << " ";
    std::cout << workingsq[i+(m_-1)*m_] << std::endl;
  }
  */

  //Form (1 + U^T * D^-1 * U)^-1 * (U^T D^-1)
  //Instead of actually inverting, we can 
  // solve the system (1 + U^T * D^-1 * U) * C = (U^T D^-1)
  // for C. That is, WORKINGSQ * C = WORKING
  // Store result (C) in WORKING (m_ x n_)		      
  //We use LU decomposition to do this
  //In the GSL case we actually invert because it doesn't
  // have a matrix solver for LU decomposition
  //std::cerr << "Forming (1 + U^T * D^-1 * U)^-1 * (U^T D^-1)" << std::endl;
  if ( m_ == 1 ) {
    //This special case causes problems for the linear algebra routines
    //But is pretty darn easy to handle.  Obviously, we don't bother
    // to factor here.
#if USEMKL || USEATLAS || USEACCELERATE
    double invval = 1.0 / workingsq[0];
    for (unsigned int i = 0; i < n_; ++i)
      working[i] *= invval;
#else
    double invval = 1.0 / gsl_matrix_get(workingsq,0,0);
    gsl_matrix_scale(working,invval);
#endif

  } else {
#if USEMKL
    //Factor
    dgetrf( &im, &im, workingsq, &im, ipiv, &status);
    if (status) return;
    //Then solve
    char trans = 'N';
    dgetrs( &trans, &im, &in, workingsq, &im, ipiv, working,
	    &im, &status );
    if (status) return;
#elif USEATLAS
    //Factor
    status = clapack_dgetrf( CblasColMajor, im, im, workingsq, im, ipiv );
    if (status) return;
    //Solve
    status = clapack_dgetrs( CblasColMajor, CblasNoTrans, im, in, workingsq,
			     im, ipiv, working, im);
    if (status) return;
#elif USEACCELERATE
    __CLPK_integer N = static_cast<__CLPK_integer>(n_);
    __CLPK_integer M = static_cast<__CLPK_integer>(m_);
    __CLPK_integer STATUS = static_cast<__CLPK_integer>(status);
    dgetrf_( &M, &N, workingsq, &M, ipiv, &STATUS );
    if (STATUS) return;
    char trans = 'N';
    dgetrs_( &trans, &M, &N, workingsq, &M, ipiv, working, &M, &STATUS );
    status = static_cast<int>(STATUS);
#else
    //Factor
    int sgnum;
    status = gsl_linalg_LU_decomp( workingsq, ipiv, &sgnum );
    if (status) return;
    
    //Invert
    gsl_matrix *gsl_inverse = gsl_matrix_alloc( m_, m_ );
    status = gsl_linalg_LU_invert( workingsq, ipiv, gsl_inverse );
    if (status) return;
    
    //Now form the product and copy it back into working
    gsl_matrix* gsl_prod = gsl_matrix_alloc( m_, n_ );
    status = gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0,
			     gsl_inverse, working, 0.0, gsl_prod );
    gsl_matrix_memcpy( working, gsl_prod );
    gsl_matrix_free( gsl_inverse );
    gsl_matrix_free( gsl_prod );
#endif
  }
  /*
  std::cerr << "(1 + U^T D^-1 U)^-1 (U^T D^-1)" << std::endl;
  for (unsigned int i = 0; i < m_; ++i) {
    for (unsigned int j = 0; j < n_-1; ++j)
      std::cout << working[i+j*m_] << " ";
    std::cout << working[i+(n_-1)*m_] << std::endl;
  }
  */


  //Now calculate D^-1 U and store in working2
  //std::cerr << "Forming D^-1 U" << std::endl;
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int j = 0; j < m_; ++j) {
    dptr = working2 + j * n_;
    dptr2 = U + j * n_;
    for (unsigned int i = 0; i < n_; ++i)
      // working2[i][j] = invdiag[i] * U[i][j];
      dptr[i] = invdiag[i] * dptr2[i];
  }
#else
  for (unsigned int j = 0; j < m_; ++j)
    for (unsigned int i = 0; i < n_; ++i)
      gsl_matrix_set( working2, i, j, invdiag[i] *
		      gsl_matrix_get( U, i, j ) );
#endif

  //Form - D^-1 U (in working2) times what is now in working, and store
  // in the inverse covariance matrix
  // D^-1 * U is n_ by m_
  // working is m_ by n_
#if USEMKL || USEATLAS || USEACCELERATE
  //CBLAS version
  // so m = n_, n = n_, k = m_
  cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
	       in, in, im, -1.0, working2, in, working,
	       im, 0.0, invcov.getData(), in );
#else
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0,
		  working2, working, 0.0, inverse );
  //Now copy into invcov
  double *rowptr;
  for (unsigned int i = 0; i < n_; ++i) {
    rowptr = invcov[i];
    for (unsigned int j = 0; j < n_; ++j)
      rowptr[j] = gsl_matrix_get( inverse, i, j );
  }
#endif

  //Now, finally, update the diagonal
  for (unsigned int i = 0; i < n_; ++i)
    invcov[i][i] += invdiag[i];  
}

void covMatrixWoodbury::readFile( const std::string& FileName ) {
  unsigned int n, m;
  std::ifstream fl(FileName.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Data file " << FileName << " not found";
    throw CosFitterExcept("covmatrix_woodbury","readFile",errstrng.str(),2);
  }
  fl >> n >> m;
  initialize(n,m);
  //Note that diag isn't in the file, since it's provided by
  // the main data file
  for (unsigned int i = 0; i < n_; ++i)
    diag_[i] = 0;

#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < n_; ++i)
    for (unsigned int j = 0; j < m_; ++j)
      fl >> U[ i + j*n_ ];
#else
  double tmp;
  for (unsigned int i = 0; i < n_; ++i)
    for (unsigned int j = 0; j < m_; ++j) {
      fl >> tmp;
      gsl_matrix_set(U,i,j,tmp);
    }
#endif

  fl.close();
}

void covMatrixWoodbury::writeFile( const std::string& FileName ) const {
  std::ofstream fl(FileName.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Couldn't open " << FileName;
    throw CosFitterExcept("covmatrix_woodbury","writeFile",errstrng.str(),2);
  }
  fl << n_ << " " << m_ << std::endl;
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < n_; ++i) {
    for (unsigned int j = 0; j < m_-1; ++j)
      fl << U[ i + j * n_ ] << " ";
    fl << U[ i + (m_-1) * n_ ] << std::endl;
  }
#else
  for (unsigned int i = 0; i < n_; ++i) {
    for (unsigned int j = 0; j < m_-1; ++j)
      fl << gsl_matrix_get(U,i,j) << " ";
    fl << gsl_matrix_get(U,i,m_-1) << std::endl;
  }
#endif
  fl.close();
}

covMatrixWoodbury& covMatrixWoodbury::operator=(const covMatrixWoodbury& B) {
  if ( this == &B ) return *this; //Self assignment
  if ( n_ != B.n_ || m_ != B.m_ )
    resize( B.n_, B.m_ );
  if ( n_ != B.n_ )
    diag_.resize( B.n_ );
  
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < m_ * n_; ++i)
    U[i] = B.U[i];
#else
  gsl_matrix_memcpy(U, B.U);
#endif
  
  return *this;
}

covMatrixWoodbury& 
covMatrixWoodbury::add(const covMatrixWoodbury& A, int& status) {
  if (status) return *this;

  if (&A == this) {
    //Self addition
    for (unsigned int i = 0; i < n_; ++i)
      diag_[i] *= 2.0;
#if USEMKL || USEATLAS || USEACCELERATE
    if ( U != NULL ) for (unsigned int i = 0; i < n_*m_; ++i)
      U[i] *= 2.0;
#else
    gsl_matrix_scale(U,2.0);
#endif
    return *this;
  }

  //Size check
  if (A.n_ != n_ || A.m_ != m_) {
    status = 1;
    return *this;
  }

  for (unsigned int i = 0; i < n_; ++i)
    diag_[i] += A.diag_[i];

#if USEMKL || USEATLAS || USEACCELERATE
  if ( U != NULL ) for (unsigned int i = 0; i < n_*m_; ++i)
    U[i] += A.U[i];
#else
  if ( U != NULL ) gsl_matrix_add(U,A.U);
#endif
  return *this;
}

covMatrixWoodbury& covMatrixWoodbury::operator+=(const covMatrixWoodbury& A) {
  int status;
  status = 0;
  return add(A,status);
}

/////////////////////////////////////////////////////////////////////
/*!
  Sets all the elements of the vectors to zero.
*/
void covMatrixWoodburySNe::initialize(unsigned int n, unsigned int m) {
  if (n == 0 || m == 0) {
    destroy();
    return;
  }
  
  if ( nsn_ != n || nu_ != m ) {
#if USEMKL || USEATLAS || USEACCELERATE
    if (U0 != NULL) delete[] U0;
    if (Ua != NULL) delete[] Ua;
    if (Ub != NULL) delete[] Ub;

    U0 = new double[ n * m ];
    for (unsigned int i = 0; i < n*m; ++i) U0[i]=0.0;
    Ua = new double[ n * m ];
    for (unsigned int i = 0; i < n*m; ++i) Ua[i]=0.0;
    Ub = new double[ n * m ];
    for (unsigned int i = 0; i < n*m; ++i) Ub[i]=0.0;
#else
    if (U0 != NULL) gsl_matrix_free( U0 );
    if (Ua != NULL) gsl_matrix_free( Ua );
    if (Ub != NULL) gsl_matrix_free( Ub );
    if (tmp != NULL) gsl_matrix_free( tmp );
    U0 = gsl_matrix_alloc( n, m );
    Ua = gsl_matrix_alloc( n, m );
    Ub = gsl_matrix_alloc( n, m );
    tmp = gsl_matrix_alloc( n, m );
    gsl_matrix_set_zero(U0);
    gsl_matrix_set_zero(Ua);
    gsl_matrix_set_zero(Ub);
#endif
  } 

  nsn_ = n;
  nu_ = m;
}

void covMatrixWoodburySNe::destroy() {

#if USEMKL || USEATLAS || USEACCELERATE 
  if (U0 != NULL) delete[] U0;
  if (Ua != NULL) delete[] Ua;
  if (Ub != NULL) delete[] Ub;
#else
  if (U0 != NULL) gsl_matrix_free( U0 );
  if (Ua != NULL) gsl_matrix_free( Ua );
  if (Ub != NULL) gsl_matrix_free( Ub );
  if (tmp != NULL) gsl_matrix_free( tmp );
  tmp = NULL;
#endif
  U0 = Ua = Ub = NULL;
  nsn_ = 0;
  nu_ = 0;
}

covMatrixWoodburySNe::covMatrixWoodburySNe() {
  nsn_ = 0;
  nu_ = 0;
  U0 = Ua = Ub = NULL;
#if !(USEMKL || USEATLAS || USEACCELERATE)
  tmp = NULL;
#endif
}

covMatrixWoodburySNe::covMatrixWoodburySNe(unsigned int n, unsigned int m) {
  U0 = Ua = Ub = NULL;
#if !(USEMKL || USEATLAS || USEACCELERATE)
  tmp = NULL;
#endif
  initialize(n,m);
}

covMatrixWoodburySNe::covMatrixWoodburySNe( const std::string& infile0,
					    const std::string& infilea,
					    const std::string& infileb) {
  U0 = Ua = Ub = NULL;
#if !(USEMKL || USEATLAS || USEACCELERATE)
  tmp = NULL;
#endif
  readFile( infile0, infilea, infileb );
}

covMatrixWoodburySNe& covMatrixWoodburySNe::resize(unsigned int n, 
						   unsigned int m) {
  if ( nsn_ == n && nu_ == m ) {
    return *this;
  }
  initialize(n,m);
  return *this;
}

/*!
  \param[in] diag The diagonal elements of the cov matrix
  \param[in] alpha The stretch correction parameter -- set to 0 to ignore
  \param[in] beta The colour correction parameter -- set to 0 to ignore
  \returns A reference to the Woodbury Cov matrix
*/
covMatrixWoodbury& 
covMatrixWoodburySNe::buildWoodburyCovMatrix(const std::vector<double>& diag, 
					     double alpha, double beta) const {

  //Size checks
  if (diag.size() != nsn_)
     throw CosFitterExcept("covMatrixWoodburySNe","buildWoodburyCovMatrix",
			   "diag array not matching length",1); 
  if (nsn_ != covWood.getNRows() || nu_ != covWood.getUNCols() )
    covWood.resize( nsn_, nu_ );

  //Set diagonal
  covWood.setDiag( diag );

  //Now make combined U
#if USEMKL || USEATLAS || USEACCELERATE 
  double *covU;
  covU = covWood.getU();  //Kinda evil
  for (unsigned int i = 0; i < nsn_ * nu_; ++i)
    covU[i] = U0[i];
  if (alpha != 0.0) for (unsigned int i = 0; i < nsn_ * nu_; ++i)
    covU[i] += alpha*Ua[i];
  if (beta != 0.0) for (unsigned int i = 0; i < nsn_ * nu_; ++i)
    covU[i] -= beta*Ub[i];
#else
  //Perhaps tmp should be pre-allocated
  gsl_matrix *covU;
  covU = covWood.getU();  //Kinda evil
  gsl_matrix_memcpy(covU, U0);
  if (alpha != 0) {
    gsl_matrix_memcpy(tmp, Ua);
    gsl_matrix_scale(tmp, alpha);
    gsl_matrix_add(covU,tmp);
  }
  if (beta != 0) {
    gsl_matrix_memcpy(tmp, Ub);
    gsl_matrix_scale(tmp, -1.0*beta);
    gsl_matrix_add(covU,tmp);
  }
#endif
  return covWood;
}

/*!
  \param[out] cov The covariance matrix in standard form
  \param[in] diag The diagonal elements of the cov matrix
  \param[in] alpha The stretch correction parameter -- set to 0 to ignore
  \param[in] beta The colour correction parameter -- set to 0 to ignore
  \returns A reference to the Woodbury Cov matrix
*/
void covMatrixWoodburySNe::getCovMatrix(covMatrix& cov,
					const std::vector<double>& diag, 
					double alpha, double beta) const {  
  buildWoodburyCovMatrix( diag, alpha, beta );
  cov.resize( nsn_, nsn_ );
  covWood.getCovMatrix( cov );
}

/*!
  \param[out] invcov The inverse covariance matrix in normal form
  \param[in] diag The diagonal elements of the cov matrix
  \param[in] alpha The stretch correction parameter -- set to 0 to ignore
  \param[in] beta The colour correction parameter -- set to 0 to ignore
  \param[out] status Status variable -- 0 if everything is copacetic
  \returns A reference to the Woodbury Cov matrix
*/
void covMatrixWoodburySNe::getInvCovMatrix(covMatrix& invcov,
					   const std::vector<double>& diag, 
					   double alpha, double beta,
					   int& status) const {  
  if (status) return;
  buildWoodburyCovMatrix( diag, alpha, beta );
  invcov.resize( nsn_, nsn_ );
  covWood.getInvCovMatrix( invcov, status );
}

/*!
  \param[in] filename0 File to read \f$m_B\f$ matrix from
  \param[in] filenamea File to read \f$s\f$ matrix from
  \param[in] filenameb File to read \f$\mathcal C\f$ matrix from
*/
void covMatrixWoodburySNe::readFile( const std::string& filename0,
				     const std::string& filenamea,
				     const std::string& filenameb) {
  unsigned int nsn, nu;

  std::ifstream fl(filename0.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Data file 0: " << filename0 << " not found";
    throw CosFitterExcept("covMatrixWoodburySNe","readFile",errstrng.str(),2);
  }
  fl >> nsn >> nu;

  initialize(nsn,nu);

#if USEMKL || USEATLAS || USEACCELERATE 
  for (unsigned int i = 0; i < nsn_; ++i)
    for (unsigned int j = 0; j < nu_; ++j)
      fl >> U0[ i + j*nsn_ ];
#else
  double tmp;
  for (unsigned int i = 0; i < nsn_; ++i)
    for (unsigned int j = 0; j < nu_; ++j) {
      fl >> tmp;
      gsl_matrix_set(U0,i,j,tmp);
    }
#endif
  fl.close();

  //Next Ua
  unsigned int ntest, mtest;
  fl.open(filenamea.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Data file a: " << filenamea << " not found";
    throw CosFitterExcept("covMatrixWoodburySNe","readFile",errstrng.str(),2);
  }
  fl >> ntest >> mtest;
  if (ntest != nsn_ || mtest != nu_) {
    std::stringstream errstrng("");
    errstrng << "Matrix in " << filenamea << " different size than in "
	     << filenameb;
    throw CosFitterExcept("covMatrixWoodburySNe","readFile",
			  errstrng.str(),4);
  }

#if USEMKL || USEATLAS || USEACCELERATE 
  for (unsigned int i = 0; i < nsn_; ++i)
    for (unsigned int j = 0; j < nu_; ++j)
      fl >> Ua[ i + j*nsn_ ];
#else
  for (unsigned int i = 0; i < nsn_; ++i)
    for (unsigned int j = 0; j < nu_; ++j) {
      fl >> tmp;
      gsl_matrix_set(Ua,i,j,tmp);
    }
#endif
  fl.close();

  //And Ub
  fl.open(filenameb.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Data file b: " << filenameb << " not found";
    throw CosFitterExcept("covMatrixWoodburySNe","readFile",errstrng.str(),2);
  }
  fl >> ntest >> mtest;
  if (ntest != nsn_ || mtest != nu_) {
    std::stringstream errstrng("");
    errstrng << "Matrix in " << filenameb << " different size than in "
	     << filenameb;
    throw CosFitterExcept("covMatrixWoodburySNe","readFile",
			  errstrng.str(),4);
  }
#if USEMKL || USEATLAS || USEACCELERATE 
  for (unsigned int i = 0; i < nsn_; ++i)
    for (unsigned int j = 0; j < nu_; ++j)
      fl >> Ub[ i + j*nsn_ ];
#else
  for (unsigned int i = 0; i < nsn_; ++i)
    for (unsigned int j = 0; j < nu_; ++j) {
      fl >> tmp;
      gsl_matrix_set(Ub,i,j,tmp);
    }
#endif
  fl.close();
}

/*!
  \param[in] filename0 File to write \f$m_B\f$ matrix to
  \param[in] filenamea File to write \f$s\f$ matrix to
  \param[in] filenameb File to write \f$\mathcal C\f$ matrix to
*/
void covMatrixWoodburySNe::writeFile( const std::string& filename0,
				      const std::string& filenamea,
				      const std::string& filenameb ) const {
  std::ofstream fl(filename0.c_str());
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Couldn't open " << filename0;
    throw CosFitterExcept("covMatrixWoodburySNe","writeFile",errstrng.str(),2);
  }
  fl << nsn_ << " " << nu_ << std::endl;
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < nsn_; ++i) {
    for (unsigned int j = 0; j < nu_-1; ++j)
      fl << U0[ i + j * nsn_ ] << " ";
    fl << U0[ i + (nu_-1) * nsn_ ] << std::endl;
  }
#else
  for (unsigned int i = 0; i < nsn_; ++i) {
    for (unsigned int j = 0; j < nu_-1; ++j)
      fl << gsl_matrix_get(U0,i,j) << " ";
    fl << gsl_matrix_get(U0,i,nu_) << std::endl;
  }
#endif
  fl.close();

  fl.open( filenamea.c_str() );
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Couldn't open " << filenamea;
    throw CosFitterExcept("covMatrixWoodburySNe","writeFile",errstrng.str(),2);
  }
  fl << nsn_ << " " << nu_ << std::endl;
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < nsn_; ++i) {
    for (unsigned int j = 0; j < nu_-1; ++j)
      fl << Ua[ i + j * nsn_ ] << " ";
    fl << Ua[ i + (nu_-1) * nsn_ ] << std::endl;
  }
#else
  for (unsigned int i = 0; i < nsn_; ++i) {
    for (unsigned int j = 0; j < nu_-1; ++j)
      fl << gsl_matrix_get(Ua,i,j) << " ";
    fl << gsl_matrix_get(Ua,i,nu_) << std::endl;
  }
#endif
  fl.close();

  fl.open( filenameb.c_str() );
  if (!fl) {
    std::stringstream errstrng("");
    errstrng << "Couldn't open " << filenameb;
    throw CosFitterExcept("covMatrixWoodburySNe","writeFile",errstrng.str(),2);
  }
  fl << nsn_ << " " << nu_ << std::endl;
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < nsn_; ++i) {
    for (unsigned int j = 0; j < nu_-1; ++j)
      fl << Ub[ i + j * nsn_ ] << " ";
    fl << Ub[ i + (nu_-1) * nsn_ ] << std::endl;
  }
#else
  for (unsigned int i = 0; i < nsn_; ++i) {
    for (unsigned int j = 0; j < nu_-1; ++j)
      fl << gsl_matrix_get(Ub,i,j) << " ";
    fl << gsl_matrix_get(Ub,i,nu_) << std::endl;
  }
#endif
  fl.close();
}

/*!
  \param[in] idx  Index array to reorder with.  The boundaries are
                  not checked, nor are repeats.  Must be of length n_.
*/
void covMatrixWoodburySNe::reorder( const std::vector<unsigned int>& idx ) {
  if (idx.size() != nsn_)
    throw CosFitterExcept("covMatrixWoodburySNe","reorder",
			  "index array not same size as data arrays",1);

#if USEMKL || USEATLAS || USEACCELERATE
  double *Uworking;
  Uworking = new double[ nsn_ * nu_ ];
  for (unsigned int j = 0; j < nu_; ++j)
    for (unsigned int i = 0; i < nsn_; ++i)
      Uworking[ i + j*nsn_ ] = U0[ idx[i] + j*nsn_ ];
  delete[] U0;
  U0 = Uworking;
  Uworking = new double[ nsn_ * nu_ ];
  for (unsigned int j = 0; j < nu_; ++j)
    for (unsigned int i = 0; i < nsn_; ++i)
      Uworking[ i + j*nsn_ ] = Ua[ idx[i] + j*nsn_ ];
  delete[] Ua;
  Ua = Uworking;
  Uworking = new double[ nsn_ * nu_ ];
  for (unsigned int j = 0; j < nu_; ++j)
    for (unsigned int i = 0; i < nsn_; ++i)
      Uworking[ i + j*nsn_ ] = Ub[ idx[i] + j*nsn_ ];
  delete[] Ub;
  Ub = Uworking;
#else
  gsl_vector *vworking;
  gsl_matrix *Uworking;
  vworking = gsl_vector_alloc( nu_ );
  Uworking = gsl_matrix_alloc( nsn_, nu_ );
  for (unsigned int i = 0; i < nsn_; ++i) {
    gsl_matrix_get_row(vworking, U0, idx[i]);
    gsl_matrix_set_row(Uworking,i,vworking);
  }
  gsl_matrix_free(U0);
  U0 = Uworking;
  Uworking = gsl_matrix_alloc( nsn_, nu_ );
  for (unsigned int i = 0; i < nsn_; ++i) {
    gsl_matrix_get_row(vworking, Ua, idx[i]);
    gsl_matrix_set_row(Uworking,i,vworking);
  }
  gsl_matrix_free(Ua);
  Ua = Uworking;
  Uworking = gsl_matrix_alloc( nsn_, nu_ );
  for (unsigned int i = 0; i < nsn_; ++i) {
    gsl_matrix_get_row(vworking, Ub, idx[i]);
    gsl_matrix_set_row(Uworking,i,vworking);
  }
  gsl_matrix_free(Ub);
  Ub = Uworking;
  gsl_vector_free( vworking );
#endif 
}

/*!
  \param[out] newmat Matrix with output information
  \param[in] idx  Index array to select with.  The boundaries are
                  not checked
*/
void covMatrixWoodburySNe::select( covMatrixWoodburySNe& newmat,
				   const std::vector<unsigned int>& idx ) 
  const {

  unsigned int n = idx.size();
  if (newmat.nsn_ != nsn_ || newmat.nu_ != nu_ )
    newmat.resize( n, nu_ );

#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int j = 0; j < nu_; ++j)
    for (unsigned int i = 0; i < n; ++i)
      newmat.U0[ i + j*n ] = U0[ idx[i] + j*nsn_ ];
  for (unsigned int j = 0; j < nu_; ++j)
    for (unsigned int i = 0; i < n; ++i)
      newmat.Ua[ i + j*n ] = Ua[ idx[i] + j*nsn_ ];
  for (unsigned int j = 0; j < nu_; ++j)
    for (unsigned int i = 0; i < n; ++i)
      newmat.Ub[ i + j*n ] = Ub[ idx[i] + j*nsn_ ];
#else
  gsl_vector *vworking;
  vworking = gsl_vector_alloc( nu_ );
  for (unsigned int i = 0; i < n; ++i) {
    gsl_matrix_get_row(vworking, U0, idx[i]);
    gsl_matrix_set_row(newmat.U0,i,vworking);
  }
  for (unsigned int i = 0; i < n; ++i) {
    gsl_matrix_get_row(vworking, Ua, idx[i]);
    gsl_matrix_set_row(newmat.Ua,i,vworking);
  }
  for (unsigned int i = 0; i < n; ++i) {
    gsl_matrix_get_row(vworking, Ub, idx[i]);
    gsl_matrix_set_row(newmat.Ub,i,vworking);
  }
  gsl_vector_free( vworking );
#endif 

}

covMatrixWoodburySNe& 
covMatrixWoodburySNe::add(const covMatrixWoodburySNe& A, int& status) {
  if (status) return *this;
  if (&A == this) {
    //Self addition
#if USEMKL || USEATLAS || USEACCELERATE
    if ( U0 != NULL ) for (unsigned int i = 0; i < nsn_*nu_; ++i)
      U0[i] *= 2.0;
    if ( Ua != NULL ) for (unsigned int i = 0; i < nsn_*nu_; ++i)
      Ua[i] *= 2.0;
    if ( Ub != NULL ) for (unsigned int i = 0; i < nsn_*nu_; ++i)
      Ub[i] *= 2.0;
#else
    if (U0 != NULL) gsl_matrix_scale(U0,2.0);
    if (Ua != NULL) gsl_matrix_scale(Ua,2.0);
    if (Ub != NULL) gsl_matrix_scale(Ub,2.0);
#endif
    return *this;
  }

  //Size check
  if (A.nsn_ != nsn_ || A.nu_ != nu_) {
    status = 1;
    return *this;
  }

#if USEMKL || USEATLAS || USEACCELERATE
  if ( U0 != NULL ) for (unsigned int i = 0; i < nsn_*nu_; ++i)
    U0[i] += A.U0[i];
  if ( Ua != NULL ) for (unsigned int i = 0; i < nsn_*nu_; ++i)
    Ua[i] += A.Ua[i];
  if ( Ub != NULL ) for (unsigned int i = 0; i < nsn_*nu_; ++i)
    Ub[i] += A.Ub[i];
#else
  if ( U0 != NULL ) gsl_matrix_add(U0,A.U0);
  if ( Ua != NULL ) gsl_matrix_add(Ua,A.Ua);
  if ( Ub != NULL ) gsl_matrix_add(Ub,A.Ub);
#endif
  return *this;
}

covMatrixWoodburySNe& 
covMatrixWoodburySNe::operator+=(const covMatrixWoodburySNe& A) {
  int status;
  status = 0;
  return add(A,status);
}

covMatrixWoodburySNe& 
covMatrixWoodburySNe::append(const covMatrixWoodburySNe& cov, int& status) {
  if (status) return *this;
  
  if ( nu_ != cov.nu_ ) {
    status = 1;
    return *this;
  }

#if USEMKL || USEATLAS || USEACCELERATE
  unsigned int newnsn_;
  double* Utmp;
  newnsn_ = nsn_ + cov.nsn_;
  Utmp = new double[ newnsn_ * nu_ ];
  for (unsigned int j = 0; j < nu_; ++j) {
    for (unsigned int i = 0; i < nsn_; ++i)
      Utmp[i + j*newnsn_] = U0[i + j*nsn_];
    for (unsigned int i = 0; i < cov.nsn_; ++i)
      Utmp[ i + nsn_ + j*newnsn_] = cov.U0[i + j*cov.nsn_];
  }
  delete[] U0;
  U0 = Utmp;
  Utmp = new double[ newnsn_ * nu_ ];
  for (unsigned int j = 0; j < nu_; ++j) {
    for (unsigned int i = 0; i < nsn_; ++i)
      Utmp[i + j*newnsn_] = Ua[i + j*nsn_];
    for (unsigned int i = 0; i < cov.nsn_; ++i)
      Utmp[ i + nsn_ + j*newnsn_] = cov.Ua[i + j*cov.nsn_];
  }
  delete[] Ua;
  Ua = Utmp;
  Utmp = new double[ newnsn_ * nu_ ];
  for (unsigned int j = 0; j < nu_; ++j) {
    for (unsigned int i = 0; i < nsn_; ++i)
      Utmp[i + j*newnsn_] = Ub[i + j*nsn_];
    for (unsigned int i = 0; i < cov.nsn_; ++i)
      Utmp[ i + nsn_ + j*newnsn_] = cov.Ub[i + j*cov.nsn_];
  }
  delete[] Ub;
  Ub = Utmp;
#else
  unsigned int newnsn_;
  gsl_matrix* Utmp;
  newnsn_ = nsn_ + cov.nsn_;
  Utmp = gsl_matrix_alloc( newnsn_, nu_ );
  for (unsigned int j = 0; j < nu_; ++j) {
    for (unsigned int i = 0; i < nsn_; ++i)
      gsl_matrix_set( Utmp, i, j, gsl_matrix_get( U0, i, j ) );
    for (unsigned int i = 0; i < cov.nsn_; ++i)
      gsl_matrix_set( Utmp, i+nsn_, j, gsl_matrix_get( cov.U0, i, j ) );
  }
  gsl_matrix_free(U0);
  U0 = Utmp;
  Utmp = gsl_matrix_alloc( newnsn_, nu_ );
  for (unsigned int j = 0; j < nu_; ++j) {
    for (unsigned int i = 0; i < nsn_; ++i)
      gsl_matrix_set( Utmp, i, j, gsl_matrix_get( Ua, i, j ) );
    for (unsigned int i = 0; i < cov.nsn_; ++i)
      gsl_matrix_set( Utmp, i+nsn_, j, gsl_matrix_get( cov.Ua, i, j ) );
  }
  gsl_matrix_free(Ua);
  Ua = Utmp;
  Utmp = gsl_matrix_alloc( newnsn_, nu_ );
  for (unsigned int j = 0; j < nu_; ++j) {
    for (unsigned int i = 0; i < nsn_; ++i)
      gsl_matrix_set( Utmp, i, j, gsl_matrix_get( Ub, i, j ) );
    for (unsigned int i = 0; i < cov.nsn_; ++i)
      gsl_matrix_set( Utmp, i+nsn_, j, gsl_matrix_get( cov.Ub, i, j ) );
  }
  gsl_matrix_free(Ub);
  Ub = Utmp;
#endif
  nsn_ = newnsn_;
  return *this;
}

covMatrixWoodburySNe& 
covMatrixWoodburySNe::operator=(const covMatrixWoodburySNe& B) {
  if ( this == &B ) return *this; //Self assignment
  if ( nsn_ != B.nsn_ || nu_ != B.nu_ )
    resize( B.nsn_, B.nu_ );
  
#if USEMKL || USEATLAS || USEACCELERATE
  for (unsigned int i = 0; i < nsn_ * nu_; ++i)
    U0[i] = B.U0[i];
  for (unsigned int i = 0; i < nsn_ * nu_; ++i)
    Ua[i] = B.Ua[i];
  for (unsigned int i = 0; i < nsn_ * nu_; ++i)
    Ub[i] = B.Ub[i];
#else
  gsl_matrix_memcpy(U0, B.U0);
  gsl_matrix_memcpy(Ua, B.Ua);
  gsl_matrix_memcpy(Ub, B.Ub);
#endif
  
  return *this;
}
