#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream, std::stringbuf

void create_sys()
{
  TFile *f= new TFile("bck_eta_cut450_notrash_UNBLINDED.root");  
  TH1D *hist= dynamic_cast<TH1D*>(gROOT->FindObject("expected_exp"));
  TH1D *hist_up= dynamic_cast<TH1D*>(gROOT->FindObject("expected_plus"));
  TH1D *hist_down= dynamic_cast<TH1D*>(gROOT->FindObject("expected_minus"));
  int xbins=expected_exp->GetNbinsX();
  int ybins=expected_exp->GetNbinsY();

  int counter=0;

  for(int x=1;x<=xbins;++x)
    {

      stringstream ss;
      ss<<"Sys_plus_"<<counter;
      stringstream ss1;
      ss1<<"Sys_minus_"<<counter;
      
      TH1F *hist_sys_plus=hist->Clone(ss.str().c_str());
      TH1F *hist_sys_minus=hist->Clone(ss1.str().c_str());
      
      

      for(int y=1;y<=ybins;++y)
	{
	  double mean=hist->GetBinContent(x,y);
	  double sigma=hist_up->GetBinContent(x,y)-mean;
	  hist_sys_plus->SetBinContent(x,y,mean+3*sigma);
	  if(mean-3*sigma>0) hist_sys_minus->SetBinContent(x,y,mean-3*sigma);
	  else hist_sys_minus->SetBinContent(x,y,0.);
	  
	}
      stringstream ss2;
      ss2<<"Sys_ETAVETO/Sys_"<<counter<<"_UNBLIND.root";
      TFile *ff= new TFile(ss2.str().c_str(), "RECREATE");
      hist_sys_plus->Write();
      hist_sys_minus->Write();
      ff->Close();


      counter++;



 }







}
