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
#include <cstdlib>
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
#include "Riostream.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
//
using namespace std;
// ----------------------------------------------------------------------
int main(int argc, char* argv[]){
  int ifile;
  char firstch;
  string dummy;
  //
  // INPUTS:
  //
  string fname_common = "../Common/Parameters.input" ;
  //
  vector<string> fnames;
  fnames.push_back("../TauFit/babar/gamma3by5.input");
  fnames.push_back("../TauFit/babar/gamma9by5.input");
  fnames.push_back("../TauFit/babar/gamma10by5.input");
  fnames.push_back("../TauFit/babar/gamma16.input");
  fnames.push_back("../TauFit/babar/gamma35.input");
  fnames.push_back("../TauFit/babar/gamma40.input");
  fnames.push_back("../TauFit/babar/gamma60.input");
  fnames.push_back("../TauFit/babar/gamma85.input");
  fnames.push_back("../TauFit/babar/gamma93.input");
  fnames.push_back("../TauFit/babar/gamma103.input");
  fnames.push_back("../TauFit/babar/gamma96.input");
  fnames.push_back("../TauFit/babar/gamma136.input");
  fnames.push_back("../TauFit/babar/gamma128.input");
  fnames.push_back("../TauFit/belle/gamma13.input");
  fnames.push_back("../TauFit/belle/gamma35.input");
  fnames.push_back("../TauFit/belle/gamma60.input");
  fnames.push_back("../TauFit/belle/gamma85.input");
  fnames.push_back("../TauFit/belle/gamma93.input");
  fnames.push_back("../TauFit/belle/gamma126.input");
  fnames.push_back("../TauFit/belle/gamma128.input");
  fnames.push_back("../TauFit/belle/gamma130.input");
  fnames.push_back("../TauFit/belle/gamma132.input");
  fnames.push_back("../TauFit/belle/gamma96.input");
  //  fnames.push_back("../TauFit/belle/gamma115.input");
  //
  int ipar;
  string par_name[100];
  double par_val[100], par_poserr[100], par_negerr[100];
  //
  ifstream ifs(fname_common.data());
  if (!ifs.good()) {
    cout << Form("Cannot open input file : %s\n", fname_common.data()); 
    exit(1);
  }
  //
  ipar=-1;
  int nlines=-1;
  while(ifs.good()) {
    ++nlines;
    firstch=' '; ifs.get(firstch) ;  ifs.putback(firstch) ;
    if (firstch=='*'||firstch=='\n') { // Skip such lines
      ifs.ignore(256,'\n') ;
    } else {
      ifs >> dummy;
      if (strcmp(dummy.data(),"PARAMETERS")==0) {
	ifs.ignore(256,'\n') ; 
      } else {
	++ipar;
	par_name[ipar].assign(dummy);
	ifs >> par_val[ipar] >> par_poserr[ipar] >> par_negerr[ipar];
	ifs.ignore(256,'\n') ;
      }
    }
  }
  ifs.close();
  int npar=ipar;
  //
  cout << Form("Read %d lines, %d parameters from : %s\n", nlines,npar,fname_common.data());
  for (ipar=0;ipar<npar;++ipar){
    cout << "ipar = " << ipar << " " << par_name[ipar] << " " << par_val[ipar] << " " << par_poserr[ipar] << " " << par_negerr[ipar] << endl;
  }
  //
  //
  int imeas=-1;
  int jmeas, kpar, kpar_sv, esys, esys_sv;
  vector<string> nsys_name;
  string meas_name[100], meas_exp[100], meas_gamma[100], meas_pub[100];
  string meas_corr_exp[100][100], meas_corr_gamma[100][100], meas_corr_pub[100][100], meas_par_name[100][100], meas_esys_name[100][100];
  double meas_val[100], meas_stat[100], meas_syst[100], meas_corr_val[100][100], meas_val_adj[100], meas_err_adj[100];
  double meas_par_value[100][100], meas_par_poserr[100][100], meas_par_negerr[100][100], meas_esys_val[100][100], meas_esys_val_adj[100][100];
  double meas_cov[100][100], meas_corr[100][100];
  int meas_corr_num[100], meas_corrij[100][100], meas_par_num[100], meas_sys_num[100];
  for (ifile=0;ifile<fnames.size(); ++ifile) {
    vector<string> vdummy;
    ifs.open(fnames[ifile].data(),ifstream::in);
    if (!ifs.good()) {
      cout << Form("Cannot open input file : %s\n", fnames[ifile].data()); 
      exit(1);
    } else {
      while (ifs.good()){
	firstch=' '; ifs.get(firstch); ifs.putback(firstch) ;
	if (firstch=='*' || firstch=='\n') {// Skip such lines
	  ifs.ignore(256,'\n'); 
	  continue;
	}
	ifs >> dummy;
	if (strcmp(dummy.data(),"*")==0||strcmp(dummy.data(),"!")==0) { // Skip rest of the line
	  ifs.ignore(256,'\n'); 
	  continue;
	}
	if (
	    (strcmp(dummy.data(),"statistical")==0) ||
	    (strcmp(dummy.data(),"systematic")==0) ||
	    (strcmp(dummy.data(),"stat")==0) ||
	    (strcmp(dummy.data(),"syst")==0)
	    ) {/*skip*/} else { // Keep the rest
	  vdummy.push_back(dummy);
	}
      }
    }
    ifs.close();
    //    cout << Form("Read from file: %s\n",fnames[ifile].data());
    //
    bool kpar_start=false;
    bool esys_start=false;
    for (vector<string>::iterator vsi=vdummy.begin();vsi!=vdummy.end();++vsi) {
      dummy=*vsi;
      int keyword = 0;
      if (
	  (strcmp(dummy.data(),"BEGIN")==0) ||
	  (strcmp(dummy.data(),"PARAMETERS")==0) ||
	  (strcmp(dummy.data(),"MEASUREMENT")==0) ||
	  (strcmp(dummy.data(),"STAT_CORR_WITH")==0) ||
	  (strcmp(dummy.data(),"DATA")==0) ||
	  (strcmp(dummy.data(),"END")==0)
	  ) {
	keyword = 1;
      }
      if (strcmp(dummy.data(),"BEGIN")==0) {
	++imeas;
	jmeas=-1;
	kpar=-2;
	esys=-2;
	meas_exp[imeas]=*(vsi+1);
	meas_gamma[imeas]=*(vsi+2); meas_gamma[imeas].erase(0,5);
	meas_pub[imeas]=*(vsi+3);
      }
      if (strcmp(dummy.data(),"MEASUREMENT")==0){
	meas_name[imeas]=*(vsi+1);
	meas_val[imeas]=strtod((*(vsi+4)).data(),0);
	meas_stat[imeas]=strtod((*(vsi+5)).data(),0);
	meas_syst[imeas]=strtod((*(vsi+6)).data(),0);
      }
      if (strcmp(dummy.data(),"STAT_CORR_WITH")==0){
	++jmeas;
	meas_corr_exp[imeas][jmeas]=*(vsi+1);
	meas_corr_gamma[imeas][jmeas]=*(vsi+2); meas_corr_gamma[imeas][jmeas].erase(0,5);
	meas_corr_pub[imeas][jmeas]=*(vsi+3);
	meas_corr_val[imeas][jmeas]=strtod((*(vsi+4)).data(),0);
      }
      if (strcmp(dummy.data(),"PARAMETERS")==0){
	kpar_start=true;
      } else {
	if (kpar_start && keyword) kpar_start=false; 
      }
      if (kpar_start) ++kpar;
      if (kpar>-1 && kpar!=kpar_sv) {
	if (kpar%4==0) meas_par_name[imeas][kpar/4] = *vsi;
	if (kpar%4==1) meas_par_value[imeas][kpar/4] =  strtod((*vsi).data(),0);
	if (kpar%4==2) meas_par_poserr[imeas][kpar/4] =  strtod((*vsi).data(),0);
	if (kpar%4==3) meas_par_negerr[imeas][kpar/4] =  strtod((*vsi).data(),0);
      }
      if (strcmp(dummy.data(),"DATA")==0 && strcmp((*(vsi+1)).substr(0,2).data(),"m_")!=0) {
	//	cout << *(vsi+1) << " " << strcmp((*(vsi+1)).substr(0,2).data(),"m_") << endl;
	esys_start=true;
      } else {
	if (esys_start && keyword) esys_start=false;
      }
      if (esys_start) ++esys;
      //
      if (esys>-1 && esys!=esys_sv) {
	if (esys%2==0) {
	  meas_esys_name[imeas][esys/2] = *vsi; 
	  string stemp=*vsi;
	  size_t found=stemp.find("%");
	  if (found!=string::npos) stemp.erase(found);
	  nsys_name.push_back(stemp);
	}
	if (esys%2==1) meas_esys_val_adj[imeas][esys/2] = meas_esys_val[imeas][esys/2] = strtod((*vsi).data(),0); 
      }
      //
      kpar_sv = kpar;
      esys_sv = esys;
      /*
      cout  << " keyword = " << keyword << " dummy = " << dummy // << " substr  = " << dummy.substr(0,2) 
	    << " kpar = " << kpar << " kpar_sv = " << kpar_sv 
	    << " esys = " << esys << " esys_sv = " << esys_sv 
	    << endl;
      */
    }
    //
    meas_corr_num[imeas] = jmeas+1;
    meas_par_num[imeas] = (kpar+1)/4;
    meas_sys_num[imeas] = (esys+1)/2;
    //
  }
  //
  sort(nsys_name.begin(),nsys_name.end());
  vector<string>::iterator new_end=unique(nsys_name.begin(),nsys_name.end());
  int nsys=nsys_name.size();
  cout << "nsys = " << nsys << " :: " ; for (esys=0;esys<nsys;++esys) {cout << esys << " " << nsys_name[esys] << " , "; } cout << endl << endl;
  nsys=new_end-nsys_name.begin();
  cout << "nsys = " << nsys << " :: " ; for (esys=0;esys<nsys;++esys) {cout << esys << " " << nsys_name[esys] << " , "; } cout << endl << endl;
  //
  int nmeas=imeas+1;
  for (imeas=0;imeas<nmeas;++imeas){
    cout << endl;
    cout << Form("Read from file: %s\n",fnames[imeas].data());
    cout << "imeas = " << imeas 
	 << " name = " << meas_name[imeas] 
	 << " exp = " <<  meas_exp[imeas] 
	 << " gamma = " << meas_gamma[imeas] 
	 << " pub = " << meas_pub[imeas]  
	 << " val = " << meas_val[imeas] 
	 << " stat = " << meas_stat[imeas] 
	 << " syst = " << meas_syst[imeas] 
	 << " corr_n = " << meas_corr_num[imeas]
	 << " par_n = " << meas_par_num[imeas] 
	 << " sys_n = " << meas_sys_num[imeas]
	 << endl;
    //
    for (esys=0;esys<meas_sys_num[imeas];++esys) {
      size_t found=meas_esys_name[imeas][esys].find("%");
      if (found!=string::npos) {
	meas_esys_name[imeas][esys].erase(found);
	meas_esys_val[imeas][esys]*=meas_val[imeas]/100.;	
	meas_esys_val_adj[imeas][esys]*=meas_val[imeas]/100.;	
      }
      cout << "esys = " << esys << " name = " << meas_esys_name[imeas][esys] << " val = " << meas_esys_val[imeas][esys] << " adj = " << meas_esys_val_adj[imeas][esys] << endl;
    }
    //
    meas_val_adj[imeas] = meas_val[imeas];
    for (kpar=0;kpar<meas_par_num[imeas];++kpar) {
      int this_esys=-1;
      for (esys=0;esys<meas_sys_num[imeas];++esys) {
	if (strcmp(meas_par_name[imeas][kpar].data(),meas_esys_name[imeas][esys].data())==0) {
	  this_esys=esys;
	  break;
	}
      }
      if (this_esys==-1) continue;
      for (ipar=0;ipar<npar;++ipar){
	if (strcmp(meas_par_name[imeas][kpar].data(),par_name[ipar].data())==0) {
	  meas_val_adj[imeas]+= ((meas_par_value[imeas][kpar] - par_val[ipar])/(0.5*(meas_par_poserr[imeas][kpar] - meas_par_negerr[imeas][kpar]))) * fabs(meas_esys_val[imeas][this_esys]);
	  meas_esys_val_adj[imeas][this_esys] = meas_esys_val[imeas][this_esys] * (0.5*(par_poserr[ipar] - par_negerr[ipar])) / (0.5*(meas_par_poserr[imeas][kpar] - meas_par_negerr[imeas][kpar])) ;
	}
      }
      cout << "kpar = " << kpar << " name = " << meas_par_name[imeas][kpar] << " val = " << meas_par_value[imeas][kpar] << " poserr = " << meas_par_poserr[imeas][kpar] << " negerr = " << meas_par_negerr[imeas][kpar] << " meas_val_adj = " << meas_val_adj[imeas] << " esys = " << this_esys << " esys_val_adj = " << meas_esys_val_adj[imeas][this_esys] << endl;
    }
    meas_err_adj[imeas] = pow(meas_stat[imeas],2);
    for (esys=0;esys<meas_sys_num[imeas];++esys) {
      meas_err_adj[imeas]+= pow(meas_esys_val_adj[imeas][esys],2);
    }
    meas_err_adj[imeas] = sqrt(max(0.,meas_err_adj[imeas]));
    cout << Form("adjusted value = %20.10fg +- %20.10f\n",meas_val_adj[imeas],meas_err_adj[imeas]) ;
    //
    for (jmeas=0;jmeas<meas_corr_num[imeas];++jmeas) {
      for (int ijmeas=0;ijmeas<nmeas;++ijmeas) {
	if ((strcmp(meas_exp[ijmeas].data(),meas_corr_exp[imeas][jmeas].data())==0) &&
	    (strcmp(meas_gamma[ijmeas].data(),meas_corr_gamma[imeas][jmeas].data())==0) &&
	    (strcmp(meas_pub[ijmeas].data(),meas_corr_pub[imeas][jmeas].data())==0)) {
	  meas_corrij[imeas][jmeas]=ijmeas;
	  break;
	}
      }
      cout << "jmeas = " << jmeas << " exp = " << meas_corr_exp[imeas][jmeas] << " gamma = " << meas_corr_gamma[imeas][jmeas] << " pub = " << meas_corr_pub[imeas][jmeas] << " corrij = " << meas_corrij[imeas][jmeas] << " val = " << meas_corr_val[imeas][jmeas] << endl;
    }
    //
    for (jmeas=0;jmeas<imeas;++jmeas) {
      meas_cov[imeas][jmeas]=0;
      for (int ijmeas=0;ijmeas<meas_corr_num[imeas];++ijmeas) {
	if (jmeas==meas_corrij[imeas][ijmeas]) {
	  meas_cov[imeas][jmeas]+=meas_corr_val[imeas][ijmeas]*meas_stat[imeas]*meas_stat[jmeas];
	  cout << "imeas = " << imeas << " stat = " << meas_stat[imeas] << " jmeas = " << jmeas << " stat = " << meas_stat[jmeas] << " ijmeas = " << ijmeas << " corr = " << meas_corr_val[imeas][ijmeas] << " cov = " << meas_cov[imeas][jmeas] << endl;
	}
      }
      for (esys=0;esys<nsys;++esys){	  
	int i_esys=-1;
	for (i_esys=0;i_esys<meas_sys_num[imeas];++i_esys) {
	  if (strcmp(meas_esys_name[imeas][i_esys].data(),nsys_name[esys].data())==0) {
	    //	    cout << meas_esys_name[imeas][i_esys] << endl;
	    break;
	  }
	}
	int j_esys=-1;
	for (j_esys=0;j_esys<meas_sys_num[jmeas];++j_esys) {
	  if (strcmp(meas_esys_name[jmeas][j_esys].data(),nsys_name[esys].data())==0) {
	    //	    cout << meas_esys_name[jmeas][j_esys] << endl;
	    break;
	  }
	}
	if (i_esys<meas_sys_num[imeas] && j_esys<meas_sys_num[jmeas]){
	  cout << " imeas = " << imeas << " name = " << meas_name[imeas] << " i_esys = " << i_esys << " name = " << meas_esys_name[imeas][i_esys]
	       << " jmeas = " << jmeas << " name = " << meas_name[jmeas] << " j_esys = " << j_esys << " name = " << meas_esys_name[jmeas][j_esys] << endl;
	  meas_cov[imeas][jmeas]+=meas_esys_val_adj[imeas][i_esys] * meas_esys_val_adj[jmeas][j_esys];
	}
      }
      meas_corr[imeas][jmeas] = meas_cov[imeas][jmeas] / (meas_err_adj[imeas]*meas_err_adj[jmeas]) ;
      if ( meas_cov[imeas][jmeas] != 0) {
	cout << "meas_cov["<<imeas<<"]["<<jmeas<<"] = " << Form("%20.10f",meas_cov[imeas][jmeas]) << " " 
	     << "meas_cor["<<imeas<<"]["<<jmeas<<"] = " << Form("%20.10f",meas_corr[imeas][jmeas]) << " "
	     << endl;
      }
    }
    //
  }
  for (imeas=0;imeas<nmeas;++imeas){ 
    for (jmeas=0;jmeas<nmeas;++jmeas) {
      if (imeas==jmeas){
	meas_cov[imeas][jmeas] = pow(meas_err_adj[imeas],2);
	meas_corr[imeas][jmeas]=1;
      }
      if (jmeas>imeas) {
	meas_cov[imeas][jmeas] = meas_cov[jmeas][imeas];
	meas_corr[imeas][jmeas]= meas_corr[jmeas][imeas];
      }
    }
  }
  //
  FILE *thisfile = fopen("readBB.out","w");
  fprintf (thisfile,"\nbabar belle measurements \n\n");
  for (imeas=0;imeas<nmeas;++imeas){ 
    string istr= meas_exp[imeas]+".Gamma"+meas_gamma[imeas]+"."+meas_pub[imeas];
    fprintf (thisfile,"    %-40s\n", istr.data());
    fprintf (thisfile,"val       %20.10f\n", meas_val_adj[imeas]);
    fprintf (thisfile,"err       %20.10f\n", meas_err_adj[imeas]);
  }
  fprintf (thisfile,"\nbabar belle correlation\n\n");
  for (imeas=0;imeas<nmeas;++imeas){ 
    string istr= meas_exp[imeas]+".Gamma"+meas_gamma[imeas]+"."+meas_pub[imeas];
    fprintf (thisfile,"                           %-40s\n", istr.data());
    for (jmeas=0;jmeas<nmeas;++jmeas){ 
      string jstr= meas_exp[jmeas]+".Gamma"+meas_gamma[jmeas]+"."+meas_pub[jmeas];
      fprintf (thisfile,"%-30s  %20.10f\n", jstr.data(), meas_corr[imeas][jmeas]);
    }
  }
  fclose(thisfile);
  //
  return 0;
  //
}
