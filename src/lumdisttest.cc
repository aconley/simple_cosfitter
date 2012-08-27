//lumdisttest

//Tests the luminosity distance function

#include <iostream>
#include <cstdio>
#include <ctime>

#include "lumdist.h"

using namespace std;

int main() {
  const int nsn = 9;
  double zhelio[]={ 0.016, 0.043, 0.05, 0.1, 0.13, 0.23, 0.56, 0.58, 0.92 };
  //int nsn = 2;
  //double zhel[]={0.016,0.043};
  //double realval2[]={ -8.97080004958959343, -6.80992261417432854};
  const double dz = 0.0001;
  double om,ol,w,wa;
  vector<double> dl(nsn); 
  //om = 0.2, ol=0.8
  double realval1[] = {-8.95017976263079618, -6.75550527762053932, 
  		       -6.41584459945, -4.82700939383, 
		       -4.20958120237, -2.82386529515, 
		       -0.512923928526, -0.4179260305256,
		       0.8512265247497 };
  //om = 1.0, ol = 0.0
  double realval2[] = {-8.97080004958959343, -6.80992261417432942, 
		       -6.47882488631776443, -4.94887551402269121,
		       -4.36494863920477627, -3.08188620085880149,
		       -1.03104528142059704, -0.948706379384076204,
		       0.144317058846084562 };
  //om = 0.2, ow = 0.8, w=-0.6
  double realval3[] = { -8.9584197292238592, -6.7772117417371760,
			-6.4409544195873858, -4.8754129485227136,
			-4.2711368187489045, -2.9251290036090802,
			-0.7082128477783773, -0.6174321426342105,
			0.5985051669605238 };

  //om = 0.5, ow = 0.5, w = -0.6
  double realval4[] = { -8.9630871934567420, -6.7896520774076965,
			-6.4553888464563478, -4.9038407585255914,
			-4.3077524148281485, -2.9879389330821375,
			-0.8464353662220475, -0.7597300378829626,
			0.3942911932673751 };

  SNeData sne("Test");
  sne.resize( nsn );
  SNeDataEntry sndata;
  for (int i = 0; i < nsn; ++i) {
    sndata.zhel = zhelio[i];
    sndata.zcmb = zhelio[i];
    sne[i] = sndata;
  }

  lumdist lm;
  lm.setDz( dz );

  om = 0.2; ol = 0.8;
  lm.getLumDist(sne,dl,om,ol,-1,0.0,lumdist::TRAPEZOIDAL);
  //lm.lumDistSum(sne,dl,om,ol,dz);
  cout << "Results for Omega_m=" << om << " Omega_l=" << ol << endl;
  printf("%-8s %-8s %-8s\n","z","dl","diff");
  for (int i = 0 ; i < nsn; ++i) {
    printf("%8.5f %8.5f %8.5f\n",zhelio[i],dl[i],dl[i]-realval1[i]);
  }

  om = 1.0; ol = 0.0;
  lm.getLumDist(sne,dl,om,ol,-1,0.0,lumdist::TRAPEZOIDAL);
  //lm.lumDistSum(sne,dl,om,ol, dz);
  cout << "Results for Omega_m=" << om << " Omega_l=" << ol << endl;
  printf("%-8s %-8s %-8s\n","z","dl","diff");
  for (int i = 0 ; i < nsn; ++i) {
    printf("%8.5f %8.5f %8.5f\n",zhelio[i],dl[i],dl[i]-realval2[i]);
  }

  om = 0.2; ol = 0.8; w = -1.0;
  lm.getLumDist(sne,dl,om,ol,w,0.0,lumdist::TRAPEZOIDAL);
  //lm.lumDistSumW(sne,dl,w,om,ol,dz);
  cout << "Results for Omega_m=" << om << " Omega_l=" << ol << 
    " w=" << w << endl;
  printf("%-8s %-8s %-8s\n","z","dl","diff");
  for (int i = 0 ; i < nsn; ++i) {
    printf("%8.5f %8.5f %8.5f\n",zhelio[i],dl[i],dl[i]-realval1[i]);
  }

  om = 0.2; ol = 0.8; w = -0.6;
  lm.getLumDist(sne,dl,om,ol,w,0.0,lumdist::TRAPEZOIDAL);
  //lm.lumDistSumW(sne,dl,w,om,ol,dz);
  cout << "Results for Omega_m=" << om << " Omega_l=" << ol << 
    " w=" << w << endl;
  printf("%-8s %-8s %-8s\n","z","dl","diff");
  for (int i = 0 ; i < nsn; ++i) {
    printf("%8.5f %8.5f %8.5f\n",zhelio[i],dl[i],dl[i]-realval3[i]);
  }
  
  om = 0.5; ol = 0.5; w = -0.6;
  lm.getLumDist(sne,dl,om,ol,w,0.0,lumdist::TRAPEZOIDAL);
  //lm.lumDistSumW(sne,dl,w,om,ol,dz);
  cout << "Results for Omega_m=" << om << " Omega_l=" << ol << 
    " w=" << w << endl;
  printf("%-8s %-8s %-8s\n","z","dl","diff");
  for (int i = 0 ; i < nsn; ++i) {
    printf("%8.5f %8.5f %8.5f\n",zhelio[i],dl[i],dl[i]-realval4[i]);
  }
  
  //Test for single distance function
  om = 0.2; ol = 0.8;
  cout << "Single Lumdist for Omega_m = " << om << " Omega_l = " << ol <<
    endl;
  printf("%-8s %-8s %-8s\n","z","dl","diff");
  for (int i = 0; i < nsn; ++i) {
    lm.getLumDist(sne,dl,om,ol,-1,0.0,lumdist::SINGLE);
    printf("%8.5f %8.5f %8.5f\n",zhelio[i],dl[i],dl[i]-realval1[i]);
  }

  om = 0.2; ol = 0.8;
  lm.getLumDist(sne,dl,om,ol,-1,0.0,lumdist::RK4);
  //lm.lumDistSum_rk4(sne,dl,om,ol);
  cout << "RK4 Results for Omega_m=" << om << " Omega_l=" << ol << endl;
  printf("%-8s %-8s %-8s\n","z","dl","diff");
  for (int i = 0 ; i < nsn; ++i) {
    printf("%8.5f %8.5f %8.5f\n",zhelio[i],dl[i],dl[i]-realval1[i]);
  }

  om = 0.5; ol = 0.5; w = -0.6;
  lm.getLumDist(sne,dl,om,ol,w,0.0,lumdist::RK4);
  //lm.lumDistSumW_rk4(sne,dl,w,om,ol);
  cout << "RK4 Results for Omega_m=" << om << " Omega_l=" << ol << 
    " w=" << w << endl;
  printf("%-8s %-8s %-8s\n","z","dl","diff");
  for (int i = 0 ; i < nsn; ++i) {
    printf("%8.5f %8.5f %8.5f\n",zhelio[i],dl[i],dl[i]-realval4[i]);
  }

  om = 0.5; ol = 0.5; w = -0.6; wa = 0.01;
  lm.getLumDist(sne,dl,om,ol,w,wa,lumdist::RK4);
  //lm.lumDistSumW_rk4(sne,dl,w,om,ol);
  cout << "RK4 Results for Omega_m=" << om << " Omega_DE=" << ol << 
    " w=" << w << " wa=" << wa << endl;
  printf("%-8s %-8s\n","z","dl");
  for (int i = 0 ; i < nsn; ++i) {
    printf("%8.5f %8.5f\n",zhelio[i],dl[i]);
  }

  //Failure tests
  w = -1; om = 0.1; ol = 2.0;
  int res;
  res = lm.getLumDist(sne,dl,om,ol,w,0.0,lumdist::TRAPEZOIDAL);
  //res = lm.lumDistSum(sne,dl,om,ol,w);
  if (res != 0) {
    cout << "Trap rule fails correctly for bouncy universe\n";
  } else {
    cout << "Trap rule does not fail correctly for bouncy universe\n";
  }
  res = lm.getLumDist(sne,dl,om,ol,w,0.0,lumdist::RK4);
  if (res != 0) {
    cout << "RK4 fails correctly for bouncy universe\n";
  } else {
    cout << "RK4 does not fail correctly for bouncy universe\n";
  }

  //Timing tests
  const int nrep = 3000;
  const int nsn_time = 250;
  double z_time[nsn_time];
  vector<double> dl_time(nsn_time);
  const double maxz = 1.1;
  const double dz_time = maxz / (nsn_time - 1);
  for (int i = 0; i < nsn_time; ++i) z_time[i] = i*dz_time;

  sne.resize( nsn_time );
  for (int i = 0; i < nsn_time; ++i) {
    sndata.zhel = z_time[i];
    sndata.zcmb = z_time[i];
    sne[i] = sndata;
  }

  clock_t starttime, endtime;
  
  om = 0.4; ol = 0.9;

  //Single distance
  starttime = clock();
  for (int i = 0; i < nrep; ++i)
    lm.getLumDist (sne, dl_time, om, ol, -1.0, 0.0, lumdist::SINGLE );
  endtime = clock();
  cout << "Single distances, " << nsn_time << " sn " << nrep << " times: " <<
    (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) << endl;

  //Trapezoidal rule
  starttime = clock();
  for (int i = 0; i < nrep; ++i)
    lm.getLumDist (sne, dl_time, om, ol, -1.0, 0.0, lumdist::TRAPEZOIDAL );
  endtime = clock();
  cout << "Trapezoidal Rule, " << nsn_time << " sn " << nrep << " times: " <<
    (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) << endl;

  //RK45
  starttime = clock();
  for (int i = 0; i < nrep; ++i)
    lm.getLumDist (sne, dl_time, om, ol, -1.0, 0.0, lumdist::RK4 );
  endtime = clock();
  cout << "Runge-Kutta, " << nsn_time << " sn " << nrep << " times: " <<
    (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) << endl;

}
