#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TKey.h"
#include "TLine.h"
#include "TBox.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TSpline.h"
#include "TF1.h"
#include "TH2D.h"
#include "TString.h"
#include "TCut.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPaveText.h"

double calc_sfactor(const int &n, const double* meas, const double *merr, const double &ave, const double& aerr){
  int nchi2tot=0;
  double chi2tot=0, chi2temp=0;
  double delta = 3 * sqrt(n) * aerr;
  for (int i = 0; i< n; ++i) {
    //    if (merr[i] < delta) {
    chi2temp=(meas[i]-ave)/merr[i];
    chi2temp*=chi2temp;
    chi2tot+=chi2temp;
    nchi2tot++;
    //    }
    cout <<  i << " " << ave << " " << meas[i] << " " <<  merr[i] << " " << chi2temp << " " << chi2tot << " " << nchi2tot << endl;
  }
  double sfactor = (nchi2tot>1) ? sqrt(chi2tot/(nchi2tot-1)) : 1;
  cout << sfactor << endl;
  return sfactor;
}
//int main(){
int sfactor(){

  // Inputs are from printouts after " Corrected measurements:" in TauToHmHmHpNu/average.log

  const int n_PIMPIMPIPNU = 3;
  TString expname_PIMPIMPIPNU[n_PIMPIMPIPNU] = {TString("BELLE") ,
						TString("BABAR") ,
						TString("CLEO3") };
  const double  meas_PIMPIMPIPNU[n_PIMPIMPIPNU] = { 0.8420000E-01 , 
					            0.8833700E-01 ,  
					            0.9130000E-01 } ; 
  const double  merr_PIMPIMPIPNU[n_PIMPIMPIPNU] = { 0.2551987E-02 , 
					            0.1269399E-02 ,  
					            0.4627094E-02 } ; 
  
  const int n_PIMKMPIPNU = 6;
  TString expname_PIMKMPIPNU[n_PIMKMPIPNU] = {TString("BELLE") ,
					      TString("BABAR") ,
					      TString("OPAL") ,
					      TString("CLEO3") ,
					      TString("CLEO") ,
					      TString("ALEPH") } ;  
  const double meas_PIMKMPIPNU[n_PIMKMPIPNU] = { 0.3300000E-02 ,
					         0.2725700E-02 , 
					         0.4150000E-02 , 
					         0.3840000E-02 , 
					         0.3460000E-02 , 
					         0.2140000E-02 } ; 
  const double merr_PIMKMPIPNU[n_PIMKMPIPNU] = { 0.1653014E-03 ,
					         0.9417657E-04 , 
					         0.6640030E-03 , 
					         0.4049693E-03 , 
					         0.6053925E-03 , 
					         0.4701064E-03 } ; 

  const int n_PIMKMKPNU = 6;
  TString expname_PIMKMKPNU[n_PIMKMKPNU] = {TString("BELLE") ,
					    TString("BABAR") ,
					    TString("CLEO3") ,
					    TString("OPAL") ,
					    TString("CLEO") ,
					    TString("ALEPH") } ;
  const double meas_PIMKMKPNU[n_PIMKMKPNU] = { 0.1550000E-02 ,
					       0.1346100E-02 ,
					       0.1550000E-02 ,
					       0.8700000E-03 ,
					       0.1450000E-02 ,
					       0.1630000E-02 } ;
  const double merr_PIMKMKPNU[n_PIMKMKPNU] = { 0.5590215E-04 ,
					       0.3776149E-04 ,
					       0.1081665E-03 ,
					       0.6881860E-03 ,
					       0.3087070E-03 ,
					       0.2701851E-03 } ;
  
  const int n_KMKMKPNU = 2;
  TString expname_KMKMKPNU[n_KMKMKPNU] = {TString("BELLE") ,
					  TString("BABAR") } ;
  const double meas_KMKMKPNU[n_KMKMKPNU] = { 0.3290000E-04 ,
					     0.1577700E-04 } ;
  const double merr_KMKMKPNU[n_KMKMKPNU] = { 0.2587015E-05 ,
					     0.1790238E-05 } ;

  const int n_HMHMHPNU = 3;
  TString expname_HMHMHPNU[n_HMHMHPNU] = {TString("DELPHI") ,
					  TString("OPAL") ,
					  TString("CLEO") } ;
  const double meas_HMHMHPNU[n_HMHMHPNU] = { 0.9317000E-01 ,
					     0.9870000E-01 ,
					     0.9510000E-01 } ;
  const double merr_HMHMHPNU[n_HMHMHPNU] = { 0.1217538E-02 ,
					     0.2600000E-02 ,
					     0.2118962E-02 } ;

  // TauToHmHmHpNu
  const double ave_PIMPIMPIPNU = 0.0895163 ;
  const double ave_PIMKMPIPNU  = 0.0029253 ;
  const double ave_PIMKMKPNU   = 0.0014318 ;
  const double ave_KMKMKPNU    = 0.0000216 ;

  const double aerr_PIMPIMPIPNU = 0.0007026 ;
  const double aerr_PIMKMPIPNU  = 0.0000709 ;
  const double aerr_PIMKMKPNU   = 0.0000279 ;
  const double aerr_KMKMKPNU    = 0.0000015 ;
  
  const double Correlation_12 =      0.2864430 ;
  const double Correlation_13 =      0.2247855 ;
  const double Correlation_14 =     -0.0059314 ;
  const double Correlation_23 =      0.0479168 ;
  const double Correlation_24 =      0.0730910 ;
  const double Correlation_34 =      0.0483768 ;

  const double ave_HMHMHPNU = ave_PIMPIMPIPNU + ave_PIMKMPIPNU + ave_PIMKMKPNU + ave_KMKMKPNU;
  const double aerr_HMHMHPNU = sqrt ( aerr_PIMPIMPIPNU * aerr_PIMPIMPIPNU +
				      aerr_PIMKMPIPNU * aerr_PIMKMPIPNU +
				      aerr_PIMKMKPNU * aerr_PIMKMKPNU +
				      aerr_KMKMKPNU * aerr_KMKMKPNU +
				      2 * aerr_PIMPIMPIPNU * aerr_PIMKMPIPNU * Correlation_12 + 
				      2 * aerr_PIMPIMPIPNU * aerr_PIMKMKPNU * Correlation_13 + 
				      2 * aerr_PIMPIMPIPNU * aerr_KMKMKPNU * Correlation_14 + 
				      2 * aerr_PIMKMPIPNU * aerr_PIMKMKPNU * Correlation_23 + 
				      2 * aerr_PIMKMPIPNU * aerr_KMKMKPNU * Correlation_24 + 
				      2 * aerr_PIMKMKPNU * aerr_KMKMKPNU * Correlation_34);
  
  cout << "PIMPIMPIPNU : " << endl ;

  double sfactor_PIMPIMPIPNU = calc_sfactor(n_PIMPIMPIPNU, meas_PIMPIMPIPNU, merr_PIMPIMPIPNU, ave_PIMPIMPIPNU, aerr_PIMPIMPIPNU);

  cout << "PIMKMPIPNU : " << endl ;

  double sfactor_PIMKMPIPNU = calc_sfactor(n_PIMKMPIPNU, meas_PIMKMPIPNU, merr_PIMKMPIPNU, ave_PIMKMPIPNU, aerr_PIMKMPIPNU);
 
  cout << "PIMKMKPNU " << endl ;

  double sfactor_PIMKMKPNU = calc_sfactor(n_PIMKMKPNU, meas_PIMKMKPNU, merr_PIMKMKPNU, ave_PIMKMKPNU, aerr_PIMKMKPNU);

  cout << "KMKMKPNU " << endl ;

  double sfactor_KMKMKPNU = calc_sfactor(n_KMKMKPNU, meas_KMKMKPNU, merr_KMKMKPNU, ave_KMKMKPNU, aerr_KMKMKPNU);

  cout << "HMHMHPNU " << endl ;

  double sfactor_HMHMHPNU = calc_sfactor(n_HMHMHPNU, meas_HMHMHPNU, merr_HMHMHPNU, ave_HMHMHPNU, aerr_HMHMHPNU);

  cout << "Scale Factors :" << endl;
  cout << "PIMPIMPIPNU = " << sfactor_PIMPIMPIPNU << " "
       << "PIMKMPIPNU = " << sfactor_PIMKMPIPNU << " "
       << "PIMKMKPNU = " << sfactor_PIMKMKPNU << " "
       << "KMKMKPNU = " << sfactor_KMKMKPNU << " "
       << "HMHMHPNU = " << sfactor_HMHMHPNU << " "
       << endl << endl;

  return 0;
}

