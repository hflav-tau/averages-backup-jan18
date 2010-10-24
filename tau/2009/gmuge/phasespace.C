#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
using namespace std;
//______________________________________________________________________________
double MyGradientPar(TF1* func, int ipar, const int npar, const double *x, const double *par, const double *epar, double eps){
   if(eps< 1e-10 || eps > 1) {
      Warning("Derivative","parameter esp=%g out of allowed range[1e-10,1], reset to 0.01",eps);
      eps = 0.01;
   }
   //save original parameters
   double par0 = par[ipar];
   // new array
   double parv[npar];
   for (int i=0;i<npar;++i) parv[i] = par[i];
   //
   func->InitArgs(x, parv);
   double f=func->EvalPar(x,parv) ; 
   //
   double f1, f2, g1, g2, h2, d0, d2;
   double h = eps*epar[ipar];
   parv[ipar] = par0 + h;     f1 = func->EvalPar(x,parv);
   parv[ipar] = par0 - h;     f2 = func->EvalPar(x,parv);
   parv[ipar] = par0 + h/2;   g1 = func->EvalPar(x,parv);
   parv[ipar] = par0 - h/2;   g2 = func->EvalPar(x,parv);
   //compute the central differences
   h2    = 1/(2.*h);
   d0    = f1 - f2;
   d2    = 2*(g1 - g2);
   double  grad = h2*(4*d2 - d0)/3.;
//   //
//   cout << "eps = " << eps  << " epar = " << epar[ipar] << " h = " << h << endl;   
//   cout << "par0     = " << par0   << " f  = " << f  << endl;
//   cout << "par0+h   = " << par0+h << " f1 = " << f1 << endl;
//   cout << "par0-h   = " << par0-h << " f2 = " << f2 << endl;
//   cout << "par0+h/2 = " << par0+h << " g1 = " << g1 << endl;
//   cout << "par0-h/2 = " << par0-h << " g2 = " << g2 << endl;
//   //
   return grad;
}
//______________________________________________________________________________
double func(double *x, double * par){
  double value = x[0];
  double ratio1 = (par[0]*par[0]) / (par[2]*par[2]);
  double value1 = 1 - (8 * ratio1) + (8 * ratio1 * ratio1 * ratio1) - (ratio1 * ratio1 * ratio1 * ratio1) - (12 * ratio1 * ratio1 * log(ratio1));
  double ratio2 = (par[1]*par[1]) / (par[2]*par[2]);
  double value2 = 1 - (8 * ratio2) + (8 * ratio2 * ratio2 * ratio2) - (ratio2 * ratio2 * ratio2 * ratio2) - (12 * ratio2 * ratio2 * log(ratio2));
  value = value2/value1;
  return value;
}
//______________________________________________________________________________
int phasespace(){
  int ipar,jpar;
  //--- from PDG 2010                      
  double m_e =  0.510998910; double e_m_e = 0.000000013 ; // MeV 
  double m_mu = 105.6583668 ; double e_m_mu = 0.0000038 ; // MeV 
  //--- from HFAG 2009
  double m_tau =  1776.7673082; double e_m_tau = 0.1507259; // MeV 
  //
  const int npar=3;
  const double par[npar] = {m_e, m_mu, m_tau};
  const double epar[npar] = {e_m_e, e_m_mu, e_m_tau};
  double CovarianceMatrix[npar][npar] = {{epar[0]*epar[0],0,0},
					 {0,epar[1]*epar[1],0},
					 {0,0,epar[2]*epar[2]}};
  //
  const double xdummy[1]={0};
  TF1* func1 = new TF1 ("func1", func, xdummy[0], xdummy[0]+1, npar);
  func1->InitArgs(xdummy,par);
  double func1_val = func1->EvalPar(xdummy,par);
  //
  double deriv[npar]={0,0,0};
  const double eps=1e-2;
  for (ipar=0;ipar<npar;++ipar) {
    deriv[ipar] = MyGradientPar(func1,ipar,npar,xdummy,par,epar,eps);
    cout << "ipar = " << ipar << " par = " << par[ipar] << " epar = " << epar[ipar] << " deriv = " << deriv[ipar] << endl;
  }
  //
  double func1_err2=0;
  for (ipar=0;ipar<npar;++ipar) {
    for (jpar=0;jpar<npar;++jpar) {
      func1_err2+= deriv[ipar] * CovarianceMatrix[ipar][jpar] * deriv[jpar]; 
    }
  }
  double func1_err = TMath::Sqrt(func1_err2);
  //
  cout << "func1_val = " << func1_val << " func1_err = " << func1_err << endl;
  //
  // OLD CALCULATION
  //
  cout << "old calculations without error propagation ... " << endl;
  const double elecmass=0.000510999;
  const double muonmass=0.1056584;
  const double taumass=1.77677;
  // 
  double x1 = (elecmass * elecmass) / (taumass * taumass);
  double fx1 = 
    1 
    - 8 * x1 
    + 8 * x1 * x1 * x1 
    - x1 * x1 * x1 * x1
    - 12 * x1 * x1 * log(x1);
  //
  double x2 = (muonmass * muonmass) / (taumass * taumass);
  double fx2 = 
    1 
    - 8 * x2 
    + 8 * x2 * x2 * x2 
    - x2 * x2 * x2 * x2
    - 12 * x2 * x2 * log(x2);
  //
  cout << x1 << " " << fx1 << " " << x2 << " " << fx2  << " " << fx1/fx2 << " " << fx2/fx1 << endl;
  //
  return 0;
}
