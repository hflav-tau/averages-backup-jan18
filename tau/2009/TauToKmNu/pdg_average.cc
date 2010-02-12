#include <iostream>
#include <math.h>
#include <stddef.h>
int tau_K_1(){
  double x[99],ex1[99],ex2[99],ex[99],wt[99]; TString experiment[99];
  // tau- -> K- nu
  int n=0;
  //  x[n]=0.692;   ex1[n]=0.006;       ex2[n]=0.010;       experiment[n]="BaBar '08"; n++; 
  x[n]=0.658;   ex1[n]=0.027;       ex2[n]=0.029 ;      experiment[n]="OPAL  '01"; n++; 
  x[n]=0.696;   ex1[n]=0.025;       ex2[n]=0.014 ;      experiment[n]="ALEPH '99"; n++; 
  x[n]=0.85;    ex1[n]=0.0;         ex2[n]=0.18  ;      experiment[n]="DELPHI'94"; n++; 
  x[n]=0.66;    ex1[n]=0.07;        ex2[n]=0.09  ;      experiment[n]="CLEO  '94"; n++;
  //
  double num=0,den=0,aver=0,err=0;
  for (int i=0;i<n;i++) {
    ex[i]=sqrt(ex1[i]*ex1[i] + ex2[i]*ex2[i]);
    wt[i]=1./(ex[i]*ex[i]);
    cout << i << " " << x[i] << " +- " << ex1[i] << "(stat) +- " << ex2[i] << "(syst) +- " << ex[i] << " (tot) from " << experiment[i] << endl;
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
    cout << i << " " << x[i] << " +- " << ex1[i] << "(stat) +- " << ex2[i] << "(syst) +- " << ex[i] << " (tot) from " << experiment[i] 
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
