
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <Rtypes.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>
#include <TCut.h>
#include <TString.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include "math.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooExponential.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include <vector>

#include <unistd.h>

void script()
{
  gROOT->SetBatch(true);

  TCanvas *canvas = new TCanvas("aaa","aaa", 800,600);
  //canvas->Divide(4,4);
  //canvas->cd(1);
  gSystem->Load("libRooFit") ;
  //  gSystem->Load("/mnt/belle_data1/HOME/mchrzaszcz/cmtuser/Bender_v22r7/Phys/Tau23Mu/params/libLimitROOT.so");

  gStyle->SetOptFit(111111);
  gSystem->Umask(0002);
  // TFile *_file0 = TFile::Open("stp17reblind_mud_BDTv3e2_1017pb_flatV5.root");
  //TTree *a = (TTree*)(gDirectory->Get("Tau23Mu"));
  

  TChain* a= new TChain("Tau23Mu");
  //  a->Add("/xdat1/lhcb/mchrzaszcz/tau23mu_FINAL/2011DOWN.root/Tau23Mu");
  // a->Add("/xdat1/lhcb/mchrzaszcz/tau23mu_FINAL/2011UP.root/Tau23Mu");
  a->Add("/mnt/data2/mchrzasz/Tau23mu_2012/2011DOWN.root/Tau23Mu");
  a->Add("/mnt/data2/mchrzasz/Tau23mu_2012/2011UP.root/Tau23Mu");
  

  


  UInt_t nev=TChain::kBigNumber;


  //test2
  //   ETA CUT            ///////////////////////////////////////////////
  // double dc_glmc[5] ={ 0.1, 0.4, 0.65, 0.7, 1};  
  //  double dc_pidmc[5] = {-0.2,-0.095,-0.04,0.05,1.}; 
  /////////////////////////////////////////////////////////////////////////
  /*
    double dc_glmc[8] ={0, 0.28,0.32, 0.46, 0.53, 0.65, 0.8, 1.01};
    double dc_pidmc[7] = {0.,0.4, 0.45, 0.54, 0.63, 0.75,1.01};
  */


  vector < RooPlot*> xframes;
 
  f = fopen("results.txt","w");

   limit_interface* li = new limit_interface();
   
   
   //
   // limit_interface* li2 = new limit_interface();

   
   for (int i = 0 ; i < li->the_pidbins() ; ++i) {
     for (int j = 0 ; j < li->the_geobins() ; ++j) {
       float pidmin=li->GetPIDlow(i);
       float pidmax=li->GetPIDhigh(i);
       float MNmax = li->GetGeohigh(j);                                                                  
       float MNmin = li->GetGeolow(j);


       

	 string cut=li->CutString2011();
	 cut+="&&mass_p0p2>450&&mass_p1p2>450";
	 
	 stringstream ss;
	 ss<<"&&min(p0_ANNmuon,min(p1_ANNmuon,p2_ANNmuon)) >"<< pidmin<<"&& min(p0_ANNmuon,min(p1_ANNmuon,p2_ANNmuon))<= "<<pidmax<<"&& MN_BLEND_FLAT >"<<MNmin<<"&& MN_BLEND_FLAT <="<<MNmax;
	 cut+=ss.str();
	 cout<<cut<<endl;
	 
	 
	 sleep(10);
	 a->SetBranchStatus("*", 0);
	 a->SetBranchStatus("cleaningcut", 1);
	 a->SetBranchStatus("TriggerTOSHlt1TrackMuonDecision", 1);
	 a->SetBranchStatus("mas*", 1);
	 a->SetBranchStatus("TriggerTOSL0Muon", 1);
	 a->SetBranchStatus("TriggerTOSHlt2TriMuonTauDecision", 1);
	 a->SetBranchStatus("TriggerTOSHlt2CharmSemilepD2HMuMuDecision", 1);
	 a->SetBranchStatus("*ANNmuon*", 1);
	 a->SetBranchStatus("MN_BLEND_FLAT", 1);
	 a->SetBranchStatus("GeoMVA*", 1);

	 //Form("iDet == %d", idet);
	 // TCut cuts = c1 && c2 && c3 && c4 && c5;

	 TTree *tree = a->CopyTree(cut.c_str());
	 cout<<"NUMBER OF EVENTS "<<tree->GetEntries()<<endl;
	 // TTree *tree = a->CopyTree(c1&&c3&&c4&&c5);
	 
	 


 for (int m = 0 ; m < li->GetNMassbins() ; ++m) {
   
   float lowm = li->GetMasslow(m);
   float highm = li->GetMasshigh(m);

   int obs0=tree->Draw("mass", "mass>1776.82-20. && mass<1776.82-15.");
   int obs1=tree->Draw("mass", "mass>1776.82-15. && mass<1776.82-10.");
   int obs2=tree->Draw("mass", "mass>1776.82-10. && mass<1776.82-5.");
   int obs3=tree->Draw("mass", "mass>1776.82-5. && mass<1776.82");
   int obs4=tree->Draw("mass", "mass>1776.82 && mass<1776.82+5.");
   int obs5=tree->Draw("mass", "mass>1776.82+5. && mass<1776.82+10.");
   int obs6=tree->Draw("mass", "mass>1776.82+10. && mass<1776.82+15.");
   int obs7=tree->Draw("mass", "mass>1776.82+15. && mass<1776.82+20.");
   
   cout<<obs0<<"  "<<obs1<<"  "<<obs2<<endl;
   cout<<"TOTAL: "<<tree->GetEntries()<<endl;

   if(m==0) li->SetPrivate(j,i,(lowm+highm)*.5,obs0,obs0,obs0);
   if(m==1) li->SetPrivate(j,i,(lowm+highm)*.5,obs1,obs1,obs1);
   if(m==2) li->SetPrivate(j,i,(lowm+highm)*.5,obs2,obs2,obs2);
   if(m==3) li->SetPrivate(j,i,(lowm+highm)*.5,obs3,obs3,obs3);
   if(m==4) li->SetPrivate(j,i,(lowm+highm)*.5,obs4,obs4,obs4);
   if(m==5) li->SetPrivate(j,i,(lowm+highm)*.5,obs5,obs5,obs5);
   if(m==6) li->SetPrivate(j,i,(lowm+highm)*.5,obs6,obs6,obs6);
   if(m==7) li->SetPrivate(j,i,(lowm+highm)*.5,obs7,obs7,obs7);

	 



    } // end mass loop

 delete tree;

       }
   }


 

  TFile* fj = new TFile("bck_eta_cut450_OBSERVED.root","recreate");
 fj->WriteTObject(li->GetFilledHistogram("observed"));
//fj->WriteTObject(li2->GetFilledHistogram("expected_lin"));
 fj->WriteTObject(li->GetFilledHistogram_PLUS("also_observed"));
 fj->WriteTObject(li->GetFilledHistogram_MINUS("also_observed2"));
 fj->Close();


 






 fclose(f);

} 

