//cosgridtest.cc

//Tests filling a cosgrid4 object, then writing it to a file
// cosgridtest.txt

#include <iostream>
#include <cmath>
#include <assert.h>

#include <cosgrids.h>
#include <param_results.h>

using namespace std;

int main(int argc, char** argv) {

  int n0,n1,n2,n3;
  double min0,max0,d0,min1,max1,d1,min2,max2,d2,min3,max3,d3;
  param_tags::paramspec psp;
  psp = make_pair( param_tags::unknown, param_tags::unset );

  n0 = 5; min0 = 0.0; d0 = 0.1; max0 = min0+(n0-1)*d0;
  n1 = 60; min1 = 1.0; d1 = 0.3; max1 = min1+(n1-1)*d1;
  n2 = 3; min2 = 0.0; d2 = 0.01; max2 = min2+(n2-1)*d2;
  n3 = 4; min3 = 0.0; d3 = 0.01; max3 = min3+(n3-1)*d3;

  cosgrid2D cg2(n2,min2,max2,d2,"axis0",psp,
		n3,min3,max3,d3,"axis1",psp);
  for (int i = 0; i < n2; ++i)
    for (int j = 0; j < n3; ++j)
      cg2[i][j] = j;

  cosgrid1D cg1 = cg2.collapseAlongAxis("axis0");
  assert( cg1[0] == 0 );
  assert( cg1[2] == 6 );

  //Test binary write
  cg1.writeFile("cosgridtest.dat",true);

  //Try reading again
  cg1.readFile("cosgridtest.dat",true);
  assert( cg1[0] == 0 );
  assert( cg1[2] == 6 );

  cosgrid3D cg3(n0,min0,max0,d0,"axis0",psp,
		n2,min2,max2,d2,"axis2",psp,
		n3,min3,max3,d3,"axis3",psp);
  for (int i = 0; i < n0; ++i)
    for (int j = 0; j < n2; ++j)
      for (int k = 0; k < n3; ++k)
	cg3[i][j][k] = k;

  cosgrid2D cg2a = cg3.collapseAlongAxis("axis0");
  assert( cg2a[0][0] == 0 );
  assert( cg2a[2][3] == 15 );


  cosgrid4D cg4(n0,min0,max0,d0,"axis0",psp,
		n1,min1,max1,d1,"axis1",psp,
		n2,min2,max2,d2,"axis2",psp,
		n3,min3,max3,d3,"axis3",psp);

  for (int i = 0; i < n0; ++i)
    for (int j = 0; j < n1; ++j)
      for (int k = 0; k < n2; ++k)
	for (int m = 0; m < n3; ++m)
	  cg4[i][j][k][m] = 1000*i + 10*j + 0.1*k + 0.001*m;
	
  cg4.writeFile("cosgridtest.txt");

  cg4.writeFile("cosgridtest.dat",true);
  cg4.readFile("cosgridtest.dat",true);
  assert( cg4[1][0][0][0] == 1000 );
  assert( cg4[1][2][0][0] == 1020 );


  //Now check our ability to marginalize
  cout << "Marginalization test" << endl;
  cosgrid1D gaussian( 501, -2.5, 2.5, 0.01, "axis0", psp );
  double xval;
  for (unsigned int i = 0; i < gaussian.getAxisN(); ++i) {
    xval = gaussian.getAxisVal(i);
    gaussian[i] = gaussian.getAxisN()/sqrt(2 * 3.14159 * 0.4 * 0.4 ) *
      exp( -0.5 * ( xval - 0.1 )*( xval - 0.1 ) / ( 0.4 * 0.4 ) );
  }

  param_results res;
  res.getMarginalized1DVals( gaussian );
  res.fit = param_struct::loop;
  cout << res << endl;
  
  return 0;
}
