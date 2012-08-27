//Program to test some of the capabilities of snedata

#include<iostream>
#include<cstdio>
#include<cmath>
#include<ctime>

#include <cosfitterexcept.h>
#include <snedata.h>

using namespace std;

int main() {

  SNeData sne("Test");

  clock_t starttime, endtime;

  try {
    starttime = clock();
    sne.readData("sample/knop03_extended_lowe.dat");
    endtime = clock();
    cout << "Read with no covariance: " <<
      (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) << " sec"
	 << endl;
    cout << "Number of SNe: " << sne.size() << endl;

    //Sort by zcmb
    starttime = clock();
    sne.zcmbsort();
    endtime = clock();
    cout << "Sort with no covariance: " <<
      (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) << " sec" 
	 << endl;

    starttime = clock();
    sne.readData("sample/knop03_extended_lowe.dat");
    std::vector<std::string> fnames(1);
    fnames[0] = "sample/sample_covfile.txt";
    sne.readCovData(fnames);
    endtime = clock();
    cout << "Read with covariance: " <<
      (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) << " sec" 
	 << endl;
    cout << "For i=4 (" << sne[4].name << ") mag cov [i,i] term: " <<
      sne.getCovMatrixSNeRef().getMagCovElement(4,4) << 
      " (should be 5.0)" << endl;
    cout << "Width cov element [4,4]: " << 
      sne.getCovMatrixSNeRef().getWidthCovElement(4,4)
	 << " (should be 0.0 since not present)" << std::endl;

    starttime = clock();
    sne.zcmbsort();
    endtime = clock(); 
    cout << "Sort with covariance: " <<
      (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) << " sec" 
	 << endl;
    cout << "After sorting, i=0 (" << sne[0].name << ") variance term: "
	 << sne.getCovMatrixSNeRef().getMagCovElement(0,0) 
	 << " (should be 5.0)" << endl;

    //Forward iteration
    cout << "Forward iterator test:" << endl << "=========================="
	 << endl;;
    for (SNeData::iterator i = sne.begin(); i != sne.end(); ++i)
      printf("SN: %-6s at zcmb: %6.4f\n",i->name.c_str(),i->zcmb);
    cout << endl;
    
    //Reverse iterate and correct mag
    cout << "Reverse iterator test: " << endl << "============================="
	 << endl;
    for (SNeData::reverse_iterator i = sne.rbegin(); i != sne.rend(); ++i)
      printf("SN: %-6s at zcmb: %6.4f raw mag: % 6.3f corr mag: %6.3f\n",
	     i->name.c_str(),i->zcmb,i->mag,i->getCorrectedMag(1.5,4.1));
    cout << endl;
    
    //Now remove some SNe
    vector<int> indxarr(4); 
    indxarr[0] = 1; indxarr[1] = 3; indxarr[2] = 5; indxarr[3] = 7;
    starttime = clock();
    SNeData tmp1 = sne.copy_remove(indxarr);
    endtime = clock();
    cout << "Index removal test time: " <<
      (endtime-starttime)/static_cast<double>(CLOCKS_PER_SEC) << " sec" 
	 << endl;
    
    cout << "Index Removal test: 94S, 92bc, 92P and 90O should be removed" << 
      endl;
    cout << "Printing first 15 entries" << endl;
    cout << "=======================================================" <<
      "============" << endl;
    for (int i = 0; i < 15; ++i) 
      printf("SN: %-6s at zcmb: %6.4f stretch: %6.4f dstretch: %6.4f\n",
	     tmp1[i].name.c_str(),tmp1[i].zcmb,tmp1[i].widthpar,
	     sqrt(tmp1[i].var_widthpar) );
    cout << endl;
    
    //Try it with strings 
    vector<string> namearr(3);
    namearr[0] = "1992al"; namearr[1] = "1992bo";
    namearr[2] = "1994M";
    sne.remove(namearr,true);
    cout << "Index Removal test" << endl;
    cout << "Removing: " << endl;
    for (vector<string>::iterator i = namearr.begin(); i != namearr.end();
	 ++i) cout << *i << endl;
    cout << "Printing first 15 entries" << endl;
    cout << "=======================================================" <<
      "============" << endl;
    for (int i = 0; i < 15; ++i) 
      printf("SN: %-6s at zcmb: %6.4f E(B-V): %7.4f dE(B-V): %6.4f\n",
	     sne[i].name.c_str(),sne[i].zcmb,sne[i].colourpar,
	     sqrt(sne[i].var_colourpar) );
    cout << endl;
  } catch (const CosFitterExcept& ex ) {
    cerr << "Error in snedatatest\n" << endl;
    cerr << ex << endl;
    return 1;
  }
  return 0;
}
