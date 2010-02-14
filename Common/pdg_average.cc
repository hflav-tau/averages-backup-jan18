#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include <math.h>
#include <stddef.h>
int pdg_average(){
  double x[99],ex1[99],ex2[99],ex[99],wt[99]; TString experiment[99]; int year[99];
  const char* name="pdg_average.input";  ifstream ifs(name) ; if (!ifs.good()) {cout << "Cannot open input file '" << name << "'" << endl ; assert(0) ;}
  char buf[1024] ;  Int_t nbuf = 1024 ;  int n=0 ;
  while(ifs.good()) {
    char firstch ; ifs.get(firstch) ;
    if (ifs.eof()) break;
    if (firstch=='#' || firstch=='*') { // Skip this line
      ifs.getline(buf,nbuf) ;
    } else {
      // Put back first char
      ifs.putback(firstch) ;
      // Parse content
      ifs >> x[n] >> ex1[n] >> ex2[n] >> experiment[n] >> year[n];
      ifs.ignore(100,'\n') ;
      n++ ;
    }
  }
  cout << "Read " << n << " lines from file " << name << endl ;
  double num=0,den=0,aver=0,err=0;
  for (int i=0;i<n;i++) {
    ex[i]=sqrt(ex1[i]*ex1[i] + ex2[i]*ex2[i]);
    wt[i]=1./(ex[i]*ex[i]);
    cout << i << " " << x[i] << " +- " << ex1[i] << "(stat) +- " << ex2[i] << "(syst) +- " << ex[i] << " (tot) from " << experiment[i] << " " << year[i] 
	 << endl;
    num+=x[i]*wt[i];
    den+=wt[i];
  }
  aver=num/den;
  err =sqrt(1/den);
  double delta =  3 * sqrt(n) * err;
  cout << "Weighted average of n = "<< n << " measurements = " << aver << " +- " << err << " delta (=3sqrt(n)err) = " << delta << endl;
  //
  double chi2=0; int nchi2=0;
  for (int i=0;i<n;i++){
    double chi2tmp=0;
    if (ex[i] < delta) {
      chi2tmp=(x[i]-aver)/ex[i];
      chi2tmp*=chi2tmp;
      chi2+=chi2tmp;
      nchi2++;
    }  
    cout << i << " " << x[i] << " +- " << ex1[i] << "(stat) +- " << ex2[i] << "(syst) +- " << ex[i] << " (tot) from " << experiment[i] << " " << year[i]
         << " chi2tmp = " << chi2tmp << " chi2 = " << chi2 << " nchi2 = " << nchi2 << endl;
  }
  double scale = (nchi2>1) ? sqrt(chi2/(nchi2-1)) : 1;
  cout << "Total Chi2 for " << nchi2 << " measurements (out of " << n << ") = " << chi2 << " Scale factor = " << scale << endl;
  if (scale > 1 ) {
    cout << "Weighted average (with scale factor) = " << aver << " +- " << err*scale << endl;
  } else {
    cout << "Weighted average (with scale set to 1) = " << aver << " +- " << err << endl;
  }
  return 0;
}
