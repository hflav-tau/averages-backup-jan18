#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "TMath.h"
#include "TF1.h"
#include "TMatrixD.h"
//
using namespace std;
//
//--- from PDG
//            MOHR     08 RVUE                2006 CODATA value
//  MOHR 2008         Reviews of Modern Physics 80 (2008) 633 
// Published in Rev.Mod.Phys.80:633-730,2008.  e-Print: arXiv:0801.0028 [physics.atom-ph]
// CODATA Recomended Values of the Fundamental Physical Constants: 2006
// http://www.slac.stanford.edu/spires/find/hep/www?eprint=arXiv:0801.0028                     
const double m_e =  0.510998910; const double e_m_e = 0.000000013 ; // MeV 
const double m_mu = 105.6583668 ; const double e_m_mu = 0.0000038 ; // MeV 
//--- from PDG 
//const double m_tau =  1776.82; const double e_m_tau = 0.16; // MeV 
// -- from HFAG
const double m_tau =  1776.7673082; const double e_m_tau = 0.1507259; // MeV 
//--- from PDG
const double tau_mu  = 2.197034e-6;const double e_tau_mu  = 0.000021e-6; // s
const double tau_tau = 290.6e-15 ; const double e_tau_tau = 1.0e-15; // s
//--- from PDG
const double m_pim = 139.57018;   const double e_m_pim   = 0.00035; // MeV
const double tau_pim = 2.6033e-8; const double e_tau_pim = 0.0005e-8; // s
const double BR_PimToMumNu = 99.98770e-2 ; const double e_BR_PimToMumNu = 0.00004e-2;
//--- from PDG [Average]
const double m_km   = 493.677; const double e_m_km = 0.013; // MeV
const double tau_km = 1.2379e-8; const double e_tau_km = 0.0021e-8; // s
const double BR_KmToMumNu = 63.60e-2; const double e_BR_KmToMumNu = 0.16e-2;
//--- from Marciano:1993sh,Decker:1994ea,Decker:1994dd
const double Delta_TauToPim_over_PimToMu = 1.0016 ; const double e_Delta_TauToPim_over_PimToMu = 0.0014; 
const double Delta_TauToKm_over_KmToMu = 1.0090 ;   const double e_Delta_TauToKm_over_KmToMu = 0.0022; 
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
void calc_average(const int n, const double* value, double** covmat, double& average, double& error, double* wt){
  //
  int i,j;
  //
  TMatrixD covmatrix(n,n);
  for (i=0;i<n;++i) {
    for (j=0;j<n;++j) {
      covmatrix[i][j] = covmat[i][j];
    }
  }
  //
  TMatrixD inverse = covmatrix;
  double determinant;
  inverse.Invert(&determinant);
  //
  for (i=0;i<n;++i) {
    wt[i] = 0;
    for (j=0;j<n;++j) {
      wt[i] += inverse[i][j];
    }
  }
  //
  double sumwt=0;
  for (i=0;i<n;++i) {
    sumwt += wt[i];
  }
  //
  for (i=0;i<n;++i) {
    wt[i] /= sumwt;
  }
  //
  average = 0;
  for (i=0;i<n;++i) {
    average += wt[i] * value[i];
  }
  error = sqrt(1/sumwt);
  //
  double pull[n];
  for (i=0;i<n;++i) {
    pull[i] = value[i] - average;
  }
  //
  double chi2 = 0;
  for (i=0;i<n;++i) {
    for (j=0;j<n;++j) {
      chi2 += pull[i] * inverse[i][j] * pull[j];
    }
  }
  //  cout << "chi2 = " << chi2 << " n = " << n << " Prob = " << TMath::Prob(chi2,n-1) << endl;
}
// ----------------------------------------------------------------------
double func_gtaugmu_h(double *x, double * par){
  double dummy = x[0];
  // E. Gamiz, et. al., http://arxiv.org/pdf/0709.0282v1 Eqn 4.1
  // BR_TauToPiNu = BR_PimToMumNu * pow(M_Tau,3) * Tau_Tau * 1./(2 * pow(M_Mu,2) * M_Pi * Tau_Pi) * pow((1-pow(M_Pi/M_Tau,2))/(1-pow(M_Mu/M_Pi,2)),2) * Delta_Pi ;
  return 
    par[0] * 
    TMath::Power(par[1],3) * 
    par[2] * 
    1/(2 * TMath::Power(par[3],2) * par[4] * par[5]) * 
    TMath::Power((1-TMath::Power(par[4]/par[1],2))/(1-TMath::Power(par[3]/par[4],2)),2) * 
    par[6];
}
// ----------------------------------------------------------------------
void calc_gtaugmu_h(const double BR_TauToHmNu, const double e_BR_TauToHmNu,
		    const double BR_HmToMumNu, const double e_BR_HmToMumNu,
		    const double m_h, const double e_m_h,
		    const double tau_h, const double e_tau_h,
		    const double Delta_TauToHm_over_HmToMu, const double e_Delta_TauToHm_over_HmToMu,
		    double& func1_val, double& func1_err, 
		    double& func1_err_0, double& func1_err_1, double& func1_err_2, double& func1_err_3, 
		    double& func1_err_4, double& func1_err_5, double& func1_err_6, 
		    double& func2_val, double& func2_err, double& func2_err_exp, 
		    double& func2_err_0, double& func2_err_1, double& func2_err_2, double& func2_err_3,
		    double& func2_err_4, double& func2_err_5, double& func2_err_6){
  //
  int ipar,jpar;
  //
  const double eps=1e-2;
  //
  const int npar=7;
  const double par[npar]  = {  BR_HmToMumNu,   m_tau,   tau_tau,   m_mu,   m_h,   tau_h,   Delta_TauToHm_over_HmToMu};
  const double epar[npar] = {e_BR_HmToMumNu, e_m_tau, e_tau_tau, e_m_mu, e_m_h, e_tau_h, e_Delta_TauToHm_over_HmToMu};
  //
  double** CovarianceMatrix = new double*[npar]; for (ipar=0;ipar<npar;++ipar) CovarianceMatrix[ipar] = new double[npar];
  //
  for (ipar=0;ipar<npar;++ipar) {
    for (jpar=0;jpar<npar;++jpar) {
      if (ipar==jpar) {
	CovarianceMatrix[ipar][ipar] = epar[ipar]*epar[ipar];
      } else {
	CovarianceMatrix[ipar][jpar] = 0;
      }
    }
  }
  //
  const double xdummy[1]={0};
  TF1* func1 = new TF1 ("func1", func_gtaugmu_h, xdummy[0], xdummy[0]+1, npar);
  func1->InitArgs(xdummy,par);
  func1_val = func1->EvalPar(xdummy,par);
  //
  double* deriv = new double[npar];
  for (ipar=0;ipar<npar;++ipar) {
    deriv[ipar] = (epar[ipar]) ? MyGradientPar(func1,ipar,npar,xdummy,par,epar,eps) : 0;
  }
  //
 double func1_err2=0;
  for (ipar=0;ipar<npar;++ipar) {
    for (jpar=0;jpar<npar;++jpar) {
      func1_err2 += deriv[ipar] * CovarianceMatrix[ipar][jpar] * deriv[jpar]; 
    }
  }
  func1_err = TMath::Sqrt(func1_err2);
  //
  const int npar2=2;
  const double par2[npar2]  = {  BR_TauToHmNu,func1_val};
  const double epar2[npar2] = {e_BR_TauToHmNu,func1_err};
  double CovarianceMatrix2[npar2][npar2] = {{epar2[0]*epar2[0],0},{0,epar2[1]*epar2[1]}};
  //
  func2_val = TMath::Sqrt(par2[0]/par2[1]);
  //
  double deriv2[2];
  deriv2[0] = (1/(2*TMath::Sqrt(par2[0]/par2[1]))) * (1/par2[1]);
  deriv2[1] = (1/(2*TMath::Sqrt(par2[0]/par2[1]))) * (par2[0]/(par2[1]*par2[1])) * -1;
  //
  double func2_err2=0;
  for (ipar=0;ipar<npar2;++ipar) {
    for (jpar=0;jpar<npar2;++jpar) {
      func2_err2 += deriv2[ipar] * CovarianceMatrix2[ipar][jpar] * deriv2[jpar]; 
    }
  }
  func2_err     = TMath::Sqrt(func2_err2);
  //
  func2_err_exp = deriv2[0] * epar2[0];
  //
  func1_err_0 = deriv[1] * epar[1]; // m(tau)
  func1_err_1 = deriv[3] * epar[3]; // m(mu) 
  func1_err_2 = deriv[4] * epar[4]; // m(h)
  func1_err_3 = deriv[2] * epar[2]; // tau(tau)
  func1_err_4 = deriv[5] * epar[5]; // tau(h)
  func1_err_5 = deriv[0] * epar[0]; // BR(h)
  func1_err_6 = deriv[6] * epar[6]; // Delta
  //
  func2_err_0 = deriv2[1] * func1_err_0; // m(tau)
  func2_err_1 = deriv2[1] * func1_err_1; // m(mu) 
  func2_err_2 = deriv2[1] * func1_err_2; // m(h)
  func2_err_3 = deriv2[1] * func1_err_3; // tau(tau)
  func2_err_4 = deriv2[1] * func1_err_4; // tau(h)
  func2_err_5 = deriv2[1] * func1_err_5; // BR(h)
  func2_err_6 = deriv2[1] * func1_err_6; // Delta
  //
  delete [] CovarianceMatrix;
  delete [] deriv;
  delete func1;
}
// ----------------------------------------------------------------------
void calc_results(const double Bpi, const double e_Bpi, const double BK, const double e_BK, const double corr_Bpi_BK){
  //
  // G(tau- --> pi- nu(tau))
  //
  double Bpi_exp = 0;
  double e_Bpi_exp = 0;
  double e0_Bpi_exp = 0;
  double e1_Bpi_exp = 0;
  double e2_Bpi_exp = 0;
  double e3_Bpi_exp = 0;
  double e4_Bpi_exp = 0;
  double e5_Bpi_exp = 0;
  double e6_Bpi_exp = 0;
  double gtaugmu_pi = 0;
  double e_gtaugmu_pi = 0;
  double e0_gtaugmu_pi = 0;
  double e1_gtaugmu_pi = 0;
  double e2_gtaugmu_pi = 0;
  double e3_gtaugmu_pi = 0;
  double e4_gtaugmu_pi = 0;
  double e5_gtaugmu_pi = 0;
  double e6_gtaugmu_pi = 0;
  double experr_gtaugmu_pi = 0;
  //
  calc_gtaugmu_h(Bpi, e_Bpi, 
		 BR_PimToMumNu, e_BR_PimToMumNu,
		 m_pim, e_m_pim,
		 tau_pim, e_tau_pim,
		 Delta_TauToPim_over_PimToMu, e_Delta_TauToPim_over_PimToMu,
		 Bpi_exp, e_Bpi_exp, e0_Bpi_exp, e1_Bpi_exp, e2_Bpi_exp, e3_Bpi_exp, e4_Bpi_exp, e5_Bpi_exp, e6_Bpi_exp,
		 gtaugmu_pi, e_gtaugmu_pi, experr_gtaugmu_pi,
		 e0_gtaugmu_pi, e1_gtaugmu_pi, e2_gtaugmu_pi, e3_gtaugmu_pi, e4_gtaugmu_pi, e5_gtaugmu_pi, e6_gtaugmu_pi);
  cout << Form("%s = (%6.3f +- %6.3f) %%\n","B(tau- -> pi- nu)", 100.*Bpi, 100.*e_Bpi);
  cout << Form("%s = (%6.3f +- %6.3f [Total] +- %6.3f [mtau] +- %6.3f [mmu] +- %6.3f [mh] +- %6.3f [tautau] +- %6.3f [tauh] +- %6.3f [Bh] +- %6.3f [Delta]) %%\n",
	  "B(tau- -> pi- nu)_univ", 
	  100.*Bpi_exp, 100.*e_Bpi_exp,
	  100.*e0_Bpi_exp, 100.*e1_Bpi_exp, 100.*e2_Bpi_exp, 100.*e3_Bpi_exp, 100.*e4_Bpi_exp, 100.*e5_Bpi_exp, 100.*e6_Bpi_exp);
  cout << Form("%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [mmu] +- %6.4f [mh] +- %6.4f [tautau] +- %6.4f [tauh] +- %6.4f [Bh] +- %6.4f [Delta]\n",
	  "(gtau/gmu)_pi", gtaugmu_pi, e_gtaugmu_pi, experr_gtaugmu_pi,
	  e0_gtaugmu_pi, e1_gtaugmu_pi, e2_gtaugmu_pi, e3_gtaugmu_pi, e4_gtaugmu_pi, e5_gtaugmu_pi, e6_gtaugmu_pi);
  //
  // G(tau- --> K- nu(tau))
  //
  cout << endl;
  double BK_exp = 0;
  double e_BK_exp = 0;
  double e0_BK_exp = 0;
  double e1_BK_exp = 0;
  double e2_BK_exp = 0;
  double e3_BK_exp = 0;
  double e4_BK_exp = 0;
  double e5_BK_exp = 0;
  double e6_BK_exp = 0;
  double gtaugmu_K = 0;
  double e_gtaugmu_K = 0;
  double e0_gtaugmu_K = 0;
  double e1_gtaugmu_K = 0;
  double e2_gtaugmu_K = 0;
  double e3_gtaugmu_K = 0;
  double e4_gtaugmu_K = 0;
  double e5_gtaugmu_K = 0;
  double e6_gtaugmu_K = 0;
  double experr_gtaugmu_K = 0;
  //
  calc_gtaugmu_h(BK, e_BK,
		 BR_KmToMumNu, e_BR_KmToMumNu,
		 m_km, e_m_km,
		 tau_km, e_tau_km,
		 Delta_TauToKm_over_KmToMu, e_Delta_TauToKm_over_KmToMu,
		 BK_exp, e_BK_exp, e0_BK_exp, e1_BK_exp, e2_BK_exp, e3_BK_exp, e4_BK_exp, e5_BK_exp, e6_BK_exp,
		 gtaugmu_K, e_gtaugmu_K, experr_gtaugmu_K,
		 e0_gtaugmu_K, e1_gtaugmu_K, e2_gtaugmu_K, e3_gtaugmu_K, e4_gtaugmu_K, e5_gtaugmu_K, e6_gtaugmu_K);
  cout << Form("%s = (%6.3f +- %6.3f) %%\n","B(tau- -> K- nu)", 100.*BK, 100.*e_BK);
  cout << Form("%s = (%6.3f +- %6.3f [Total] +- %6.3f [mtau] +- %6.3f [mmu] +- %6.3f [mh] +- %6.3f [tautau] +- %6.3f [tauh] +- %6.3f [Bh] +- %6.3f [Delta]) %%\n",
	  "B(tau- -> K- nu)_univ", 
	  100.*BK_exp, 100.*e_BK_exp,
	  100.*e0_BK_exp, 100.*e1_BK_exp, 100.*e2_BK_exp, 100.*e3_BK_exp, 100.*e4_BK_exp, 100.*e5_BK_exp, 100.*e6_BK_exp);
  cout << Form("%s = %6.4f +- %6.4f [Total] +- %6.4f [Btau] +- %6.4f [mtau] +- %6.4f [mmu] +- %6.4f [mh] +- %6.4f [tautau] +- %6.4f [tauh] +- %6.4f [Bh] +- %6.4f [Delta]\n",
	  "(gtau/gmu)_K", gtaugmu_K, e_gtaugmu_K, experr_gtaugmu_K,
	  e0_gtaugmu_K, e1_gtaugmu_K, e2_gtaugmu_K, e3_gtaugmu_K, e4_gtaugmu_K, e5_gtaugmu_K, e6_gtaugmu_K);
  //
  // Average of (gtau/gmu)_pi and (gtau/gmu)_K
  //
  cout << endl;
  cout << Form("Corr between %s and %s = %6.2f %%\n","B_pi","B_K ",100.*corr_Bpi_BK);
  double gtaugmu_pik[2] = {gtaugmu_pi,gtaugmu_K};
  double** cov_gtaumu_pik = new double*[2]; for (int i=0;i<2;++i) cov_gtaumu_pik[i] = new double[2];
  cov_gtaumu_pik[0][0] = e_gtaugmu_pi * e_gtaugmu_pi;
  cov_gtaumu_pik[1][1] = e_gtaugmu_K  * e_gtaugmu_K;
  cov_gtaumu_pik[0][1] = cov_gtaumu_pik[1][0]
    = corr_Bpi_BK * experr_gtaugmu_pi * experr_gtaugmu_K 
    + e0_gtaugmu_pi * e0_gtaugmu_K + e1_gtaugmu_pi * e1_gtaugmu_K + e3_gtaugmu_pi * e3_gtaugmu_K;
  cout << Form("Corr between %s and %s = %6.2f %%\n","(gtau/gmu)_pi","(gtau/gmu)_K ",100.*cov_gtaumu_pik[0][1]/TMath::Sqrt(cov_gtaumu_pik[0][0]*cov_gtaumu_pik[1][1]));
  double gtaugmu_pik_ave = 0;
  double gtaugmu_pik_err = 0;
  double gtaugmu_pik_wt[2] = {0,0};
  calc_average(2,gtaugmu_pik,cov_gtaumu_pik,gtaugmu_pik_ave,gtaugmu_pik_err,gtaugmu_pik_wt);
  cout << Form("%s = %6.4f +- %6.4f has weights = %6.4f [pi], %6.4f [K]\n",
	       "<(gtau/gmu)_pik>",gtaugmu_pik_ave,gtaugmu_pik_err,gtaugmu_pik_wt[0],gtaugmu_pik_wt[1]);
  //
  // Clean Up
  //
  delete [] cov_gtaumu_pik;
}
int main(int argc, char* argv[]){
  const double Bpi = 10.831e-2;
  const double e_Bpi = 0.051e-2;
  const double BK = 0.697e-2;
  const double e_BK = 0.010e-2;
  const double corr_Bpi_BK = -0.49e-2;
  calc_results(Bpi, e_Bpi, BK, e_BK, corr_Bpi_BK);
  return 0;
}
