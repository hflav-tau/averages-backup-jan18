#include <iostream>
#include <math.h>
#include <stddef.h>
int PDG_average(){
  double x[99],ex1[99],ex2[99],ex[99],wt[99]; TString experiment[99];
  // tau- -> K0 pi- nu
  int n=0;
  x[n]=0.840; ex1[n]=0.004; ex2[n]=0.023; experiment[n]="BaBar '08";n++;//Babar '08
  x[n]=0.808; ex1[n]=0.004; ex2[n]=0.026; experiment[n]="Belle '07";n++;//Belle '07
  x[n]=0.933; ex1[n]=0.068; ex2[n]=0.049; experiment[n]="OPAL  '00";n++;//OPAL '00
  x[n]=0.928; ex1[n]=0.045; ex2[n]=0.034; experiment[n]="ALEPH '99";n++;//ALEPH '99
  x[n]=0.855; ex1[n]=0.117; ex2[n]=0.066; experiment[n]="ALEPH '98";n++;//ALEPH '98
  x[n]=0.704; ex1[n]=0.041; ex2[n]=0.072; experiment[n]="CLEO  '96";n++;//CLEO '96
  x[n]=0.950; ex1[n]=0.150; ex2[n]=0.060; experiment[n]="L3    '95";n++;//L3 '95
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
