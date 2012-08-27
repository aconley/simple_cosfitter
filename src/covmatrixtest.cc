#include <iostream>
#include <cmath>
#include <ctime>

#include <covmatrix.h>
#include <covmatrix_woodbury.h>

int main( int argc, char **argv ) {

  //The first order is just to test whether or not the inversions
  // are working.
  const unsigned int ntest = 3;
  int status = 0;
  double test_matrix[ntest][ntest] =
    { { 100.0, 15.0, 0.01 }, { 15.0, 2.3, 0.01 }, { 0.01, 0.01, 1.0 } };
  double test_inverse[ntest][ntest] = 
    { { 0.46064652, -3.0043271, 0.025436805},
      { -3.0043271, 20.028961, -0.17024633 },
      { 0.025436805, -0.17024633, 1.0014481 } };
  std::cout << "Initial matrix to invert: " << std::endl;
  for (unsigned int i = 0; i < ntest; ++i) {
    for (unsigned int j = 0; j < ntest; ++j)
      printf("%10.5f ",test_matrix[i][j]);
    printf("\n");
  }
  covMatrix covmatrix_test( ntest );
  for (unsigned int i = 0; i < ntest; ++i)
    for (unsigned int j = 0; j < ntest; ++j)
      covmatrix_test[i][j] = test_matrix[i][j];
  covmatrix_test.invert( status );
  if (status) {
    std::cerr << "Error: matrix inversion failed with status: "
	      << status << std::endl;
    return status;
  }
  std::cout << "covMatrix" << std::endl;
  for (unsigned int i = 0; i < ntest; ++i)
    for (unsigned int j = 0; j < ntest; ++j)
      printf("%d %d: %10.5f should be: %10.5f |diff| %10.5f\n", i,j,
	     covmatrix_test[i][j],test_inverse[i][j],
	     fabs(covmatrix_test[i][j]-test_inverse[i][j]));
 
  std::cout << "Beginning timing test" << std::endl;
  const unsigned int nsn = 300;
  const unsigned int nrep = 500;
  clock_t starttime, endtime;

  covMatrix covmatrix_c( nsn ), covmatrix_init( nsn );
  for (unsigned int i = 0; i < nsn; ++i) {
    covmatrix_init[i][i] = 1.0;
    for (unsigned int j = 0; j < i; ++j )
      covmatrix_init[i][j] = 0.01;
    for (unsigned int j = i+1; j < nsn; ++j )
      covmatrix_init[i][j] = 0.01;
  }
  starttime = clock();
  for (unsigned int i = 0; i < nrep; ++i) {
    covmatrix_c = covmatrix_init; 
    covmatrix_c.invert(status);
    if ( status ) {
      printf("Problem with inversion. status: %d\n",status);
      exit(1);
    }
  }
  endtime = clock();
  std::cout << "covMatrix inversion: " << nsn << " nrep: " << nrep
	    << " time: "
	    << (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) 
	    << " sec" << std::endl;
  
  //Now we want to test the woodbury thing
  //We test this on a 3x3 matrix
  std::cout << "Testing Woodbury" << std::endl;

  const unsigned int wood_size = 3;
  std::vector<double> diag(wood_size,0.0);
  double offdiag = 0.01, vecelem = sqrt(offdiag);
  diag[0] = 1.0 - offdiag; diag[1] = 2.0-offdiag; diag[2] = 3.0-offdiag;
  double uvec[wood_size];
  for (unsigned int i = 0; i < wood_size; ++i) uvec[i] = vecelem;
  covMatrixWoodbury woodbury_test( diag, 1, uvec );

  //This should have the same values as the above
  covMatrix woodbury_cov( wood_size );
  woodbury_cov[0][0] = 1; woodbury_cov[0][1] = 0.01; woodbury_cov[0][2] = 0.01;
  woodbury_cov[1][0] = 0.01; woodbury_cov[1][1] = 2; woodbury_cov[1][2] = 0.01;
  woodbury_cov[2][0] = 0.01; woodbury_cov[2][1] = 0.01; woodbury_cov[2][2] = 3;
  
  //See if we get back the covariance matrix we expect
  status = 0;
  std::cout << "Checking if cov matrix of expected form" << std::endl;
  covMatrix woodbury_cov_construct( wood_size );
  woodbury_test.getCovMatrix( woodbury_cov_construct );
  const double tol = 1e-5;
  for (unsigned int i = 0; i < wood_size; ++i)
    for (unsigned int j = 0; j < wood_size; ++j)
      if (fabs(woodbury_cov_construct[i][j] - woodbury_cov[i][j]) > tol) {
	status = 1;
	std::cerr << "Error in woodbury cov construction" << std::endl;
	std::cerr << " For index: " << i << " " << j << " expected " 
		  << woodbury_cov[i][j] << " got " << woodbury_cov_construct[i][j]
		  << " diff " << fabs(woodbury_cov_construct[i][j] - 
				      woodbury_cov[i][j]) << std::endl;
      }
  if (status) {
    std::cerr << "Aborting further tests" << std::endl;
    exit(1);
  }
  std::cout << "success" << std::endl;

  //Now test the inversion
  covMatrix woodbury_inv( wood_size );
  woodbury_inv[0][0] = 1.00008; woodbury_inv[0][1] = -0.00498383; 
  woodbury_inv[0][2] = -0.00331700;
  woodbury_inv[1][0] = -0.00498383; woodbury_inv[1][1] = 0.500033; 
  woodbury_inv[1][2] = -0.00165016;
  woodbury_inv[2][0] = -0.00331700; woodbury_inv[2][1] = -0.00165016; 
  woodbury_inv[2][2] = 0.333350;

  covMatrix woodbury_inv_test( wood_size );
  status = 0;
  std::cout << "Testing woodbury inversion" << std::endl;
  woodbury_test.getInvCovMatrix( woodbury_inv_test, status );
  if (status) std::cerr << "Woodbury inversion failed" << std::endl; else {
    for (unsigned int i = 0; i < wood_size; ++i)
      for (unsigned int j = 0; j < wood_size; ++j)
	if (fabs(woodbury_inv_test[i][j] - woodbury_inv[i][j]) > tol) {
	  status = 1;
	  std::cerr << "Error in woodbury inverse" << std::endl;
	  std::cerr << " For index: " << i << " " << j << " expected " 
		    << woodbury_inv[i][j] << " got " << woodbury_inv_test[i][j]
		    << " diff " << fabs(woodbury_inv_test[i][j] - 
					woodbury_inv[i][j]) << std::endl;
	}
  }
  if (status) {
    std::cout << "Aborting further tests" << std::endl;
    exit(2);
  }
  std::cout << "success" << std::endl;

  //Do that again but set m = 5, since the single vector is a special
  // case
  const unsigned int bigm = 5;
  std::vector<double> uvec2( wood_size, vecelem / 
			     sqrt( static_cast<double>(bigm) ) );
  woodbury_test.resize(wood_size,bigm);
  for (unsigned int i = 0; i < bigm; ++i)
    woodbury_test.setUColumn( i, uvec2 );
  woodbury_test.setDiag( diag );
  woodbury_test.getCovMatrix( woodbury_cov_construct );
  for (unsigned int i = 0; i < wood_size; ++i)
    for (unsigned int j = 0; j < wood_size; ++j)
      if (fabs(woodbury_cov_construct[i][j] - woodbury_cov[i][j]) > tol) {
	status = 1;
	std::cerr << "Error in second woodbury cov construction" << std::endl;
	std::cerr << " For index: " << i << " " << j << " expected " 
		  << woodbury_cov[i][j] << " got " << woodbury_cov_construct[i][j]
		  << " diff " << fabs(woodbury_cov_construct[i][j] - 
				      woodbury_cov[i][j]) << std::endl;
      }
  if (status) {
    std::cerr << "Aborting further tests" << std::endl;
    exit(1);
  }
  std::cout << "success" << std::endl;
  std::cout << "Testing second woodbury inversion" << std::endl;
  woodbury_test.getInvCovMatrix( woodbury_inv_test, status );
  if (status) std::cerr << "Woodbury inversion failed" << std::endl; else {
    for (unsigned int i = 0; i < wood_size; ++i)
      for (unsigned int j = 0; j < wood_size; ++j)
	if (fabs(woodbury_inv_test[i][j] - woodbury_inv[i][j]) > tol) {
	  status = 1;
	  std::cerr << "Error in second woodbury inverse" << std::endl;
	  std::cerr << " For index: " << i << " " << j << " expected " 
		    << woodbury_inv[i][j] << " got " << woodbury_inv_test[i][j]
		    << " diff " << fabs(woodbury_inv_test[i][j] - 
					woodbury_inv[i][j]) << std::endl;
	}
  }
  if (status) {
    std::cout << "Aborting further tests" << std::endl;
    exit(2);
  }
  std::cout << "success" << std::endl;

  //Ok -- timing tests.  Bump the cross size up to 50 systematics,
  // since that is roughly what we expect
  const unsigned int bigm2 = 50;
  woodbury_test.resize( nsn, bigm2 );
  woodbury_inv_test.resize( nsn ); //To recieve the inverse cov matrix
  diag.resize( nsn );
  for (unsigned int i = 0; i < nsn; ++i)
    diag[i] = 1 + static_cast<double>(i) * 0.01;
  uvec2.resize(nsn);
  uvec2.assign( nsn, vecelem / sqrt( static_cast<double>(bigm2) ) );
  woodbury_test.setDiag( diag );
  for (unsigned int i = 0; i < bigm2; ++i)
    woodbury_test.setUColumn( i, uvec2 );

  std::cout << "Beginning Woodbury timing test" << std::endl;
  starttime = clock();
  status = 0;
  for (unsigned int i = 0; i < nrep; ++i) {
    woodbury_test.getInvCovMatrix( woodbury_inv_test, status);
    if ( status ) {
      printf("Problem with inversion. status: %d\n",status);
      exit(1);
    }
  }
  endtime = clock();
  std::cout << "covMatrixWoodbury inversion: " << nsn << " by " << bigm2
	    << " nrep: " << nrep << " time: "
	    << (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) 
	    << " sec" << std::endl;
  
  //Don't have to do this either, but it's a test to make sure they work
  woodbury_test.clear();
  woodbury_inv_test.clear();
  return 0;
}
