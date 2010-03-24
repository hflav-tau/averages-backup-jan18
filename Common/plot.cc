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
int SetUp(){
 // gROOT->SetBatch(kTRUE);
  //gROOT->SetStyle("Plain");
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPadBorderMode(0); 
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(0.35);
  gStyle->SetErrorX(0.001);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTopMargin(0.15);
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadBottomMargin(0.175);
  gStyle->SetLabelFont(132,"X");
  gStyle->SetLabelFont(132,"Y");
  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");
  gStyle->SetLabelOffset(99.,"Y");
  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(100.,"Y");
  gStyle->SetTitleSize(0.055,"X");
  gStyle->SetTitleSize(0.055,"Y");
  gStyle->SetTitleFont(132,"X");
  gStyle->SetTitleFont(132,"Y");
  gStyle->SetTickLength(1.e-5,"Y");
  gStyle->SetNdivisions(504, "X");
  gStyle->SetNdivisions(505, "Y");
  gStyle->SetCanvasBorderSize(0);
}
//int main(){
int plot(){
  SetUp();
  //
  int nPoints=0;
  float xmin,xmax,units,meas[99],stat[99],syst[99],tot[99],wt[99]; 
  TString title,expname[99];
  const char* filename="plot.input";  
  ifstream ifs(filename) ; if (!ifs.good()) {cout << "Cannot open input file '" << filename << "'" << endl ; assert(0) ;}
  char buffer[200] ; 
  while(ifs.good()) {
    if (ifs.eof()) break;
    char firstch ; ifs.get(firstch) ;
    if (firstch=='#'||firstch=='\n') { // Skip this line
      ifs.ignore(200,'\n') ;
    } else if (firstch=='*') { 
      ifs >> xmin >> xmax >> units;
      ifs.getline(buffer,200,'\n');
      title=TString(buffer).Strip((TString::EStripType)1,' ');//remove kLeading whitespace
    } else {
      // Put back first char
      ifs.putback(firstch) ;
      // Parse content
      ifs >> meas[nPoints] >> stat[nPoints] >> syst[nPoints];
      tot[nPoints]=sqrt(stat[nPoints]*stat[nPoints] + syst[nPoints]*syst[nPoints]);
      ifs.getline(buffer,200,'\n') ;
      expname[nPoints]=TString(buffer).Strip((TString::EStripType)1,' ');//remove kLeading whitespace
      if (expname[nPoints].Length()) nPoints++;
    }
  }
  cout << "Read " << nPoints << " lines from filename " << filename << endl ;
  for (int i=0;i<nPoints;i++) {
    cout << i << " " << meas[i] << " +- " << stat[i] << "(stat) +- " << syst[i] << "(syst) +- " << tot[i] << " (tot) from " << expname[i] << endl;
  }
  //
  TCanvas *canvas=new TCanvas("canvas","canvas",356,0,656,700);
  canvas->Clear();
  //  canvas->SetFillColor(10);
  canvas->SetLeftMargin(0.04);
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.04);
  canvas->SetBottomMargin(0.12);  

  // -- set up frame for plot
  float fXmin = xmin; 
  float fXmax = xmax; 
  float fYmin = 0.0;
  float fYmax = nPoints*1.0 + 1.5;

  cout << "Drawing frame " << fXmin << "  " << fYmin << "  " << fXmax << "  " << fYmax << endl;

  TH1F *h = canvas->DrawFrame(fXmin, fYmin, fXmax, fYmax, "");  h->SetXTitle(title);  h->Draw();
  cout << "h drawn with x-axis title = " << title << endl;
  
  int ci = 999; // new color index
  TColor *color = new TColor(ci, 0.50, 1.00, 0.50); // lightgreen for average

  TBox b1;
  b1.SetFillStyle(1000);
  b1.SetFillColor(ci);
  b1.DrawBox(meas[0]-tot[0], fYmin,meas[0]+tot[0], fYmax);

  Float_t y[99], ey[99], eyl[99], eyh[99];
  int fColor[99], fSymbol[99];
  double fMarkerSize = 1.2;
  int fMarkerStyle = 24;
  int fLineW=1;
  for (int i=0;i<nPoints;++i) {
    y[i]=i+2.0; ey[i]=0.0; eyl[i]=0.0; eyh[i]=0.0;
    fColor[i]  = 1;
    fSymbol[i] = 20;
  }

  float fTxtX= 0.05;
  float xtext=fXmin+fTxtX*(fXmax-fXmin);
  float fTxtY= 0.4;
  float ytext;
  TLatex  tl;
  tl.SetTextSize(0.028);
  for (int i=0;i<nPoints;++i) {
    TGraphErrors *GraphA = new TGraphErrors(1,&meas[i],&y[i],&tot[i],&ey[i]);
    GraphA->SetMarkerColor(fColor[i]);
    GraphA->SetMarkerStyle(fSymbol[i]);
    GraphA->SetMarkerSize(fMarkerSize);
    GraphA->SetLineColor(fColor[i]);
    GraphA->SetLineWidth(fLineW);
    // -- Change if it's the average
    if (expname[i].Contains("Average")) {
      GraphA->SetMarkerColor(2);
      GraphA->SetLineColor(2);
    }
    GraphA->Draw("P");
    //
    TGraph* GraphB = new TGraphAsymmErrors(1,&meas[i],&y[i],&tot[i],&tot[i],&eyl[i],&eyh[i]);
    GraphB->SetMarkerColor(fColor[i]);
    GraphB->SetMarkerStyle(fSymbol[i]);
    GraphB->SetMarkerSize(fMarkerSize);
    GraphB->SetLineColor(fColor[i]);
    GraphB->SetLineWidth(fLineW);
    // -- Change if it's the average
    if (expname[i].Contains("Average")) {
      GraphB->SetMarkerColor(2);
      GraphB->SetLineColor(2);
    }
    GraphB->Draw("P");
    //
    tl.SetTextColor(fColor[i]);
    // -- Change if it's the average
    if (expname[i].Contains("Average")) {
      tl.SetTextColor(2);
    }
    ytext=y[i]*1.0;
    tl.DrawLatex(xtext,ytext,expname[i]);
    if (syst[i]<1e-8) {
      tl.DrawLatex(xtext,ytext-fTxtY,Form("(%6.3f #pm %6.3f)",meas[i],tot[i]));
    } else {
      tl.DrawLatex(xtext,ytext-fTxtY,Form("(%6.3f #pm %6.3f #pm %6.3f)",meas[i],stat[i],syst[i]));
    }
  }

  HFAGTauLabel(2010001,.33,.225,.25,.05,.8);

  canvas->SaveAs("plot.eps");
  //  canvas->SaveAs("plot.gif"); // ROOT's native gif saving looks garbled
  gSystem->Exec("rm -f plot.gif; convert plot.eps plot.gif");
  //canvas->Clear();

  //
  return 0;
}
// ----------------------------------------------------------------------
void HFAGTauLabel(Int_t yearversion=2010001, Double_t xpos= .99, Double_t ypos= 0.825, Double_t xsiz=0.1, Double_t ysiz=0.12, Double_t scale= 1.0) {
  TPaveText *tbox1 = new TPaveText(xpos-xsiz,ypos-ysiz,xpos,ypos,"BRNDC");
  tbox1->SetLineColor(1);
  tbox1->SetLineStyle(1);
  tbox1->SetLineWidth(1);
  tbox1->SetFillColor(1);
  tbox1->SetFillStyle(1001);
  tbox1->SetBorderSize(1);
  tbox1->SetShadowColor(0);
  tbox1->SetTextColor(0);
  tbox1->SetTextFont(76);
  tbox1->SetTextSize(24*scale);
  tbox1->SetTextAlign(22); //center-adjusted and vertically centered
  TString ts1("HFAG-Tau");
  tbox1->AddText(ts1);
  tbox1->Draw();

  TPaveText *tbox2 = new TPaveText(xpos-xsiz,ypos-1.8*ysiz,xpos,ypos-ysiz,"BRNDC");
  tbox2->SetLineColor(1);
  tbox2->SetLineStyle(1);
  tbox2->SetLineWidth(1);
  tbox2->SetFillColor(0);
  tbox2->SetFillStyle(1001);
  tbox2->SetBorderSize(1);
  tbox2->SetShadowColor(0);
  tbox2->SetTextColor(1);
  tbox2->SetTextFont(76);
  tbox2->SetTextSize(20*scale);
  tbox2->SetTextAlign(22); //center-adjusted and vertically centered
  TString ts2;
  if (yearversion==2009001) ts2="Spring 2009";
  if (yearversion==2009002) ts2="Summer 2009";
  if (yearversion==2010001) ts2="Winter 2010";
  tbox2->AddText(ts2);
  tbox2->Draw();

}
