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
int SetUp(){
 // gROOT->SetBatch(kTRUE);
  //gROOT->SetStyle("Plain");
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetPadBorderMode(0); 
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
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
  gStyle->SetEndErrorSize(3);
  gStyle->SetNdivisions(504, "X");
  gStyle->SetNdivisions(505, "Y");
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
  //  tbox1->SetShadowColor(0);
  tbox1->SetTextColor(0);
  tbox1->SetTextFont(76);
  tbox1->SetTextSize(24*scale);
  tbox1->SetTextAlign(22); //center-adjusted and vertically centered
  TString ts1("HFAG-Tau");
  tbox1->AddText(ts1);
  tbox1->Draw();
  //
  TPaveText *tbox2 = new TPaveText(xpos-xsiz,ypos-1.8*ysiz,xpos,ypos-ysiz,"BRNDC");
  tbox2->SetLineColor(1);
  tbox2->SetLineStyle(1);
  tbox2->SetLineWidth(1);
  tbox2->SetFillColor(0);
  tbox2->SetFillStyle(1001);
  tbox2->SetBorderSize(1);
  //  tbox2->SetShadowColor(0);
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
// ----------------------------------------------------------------------
//int main(std::string filename_string = "plot.input"){
int plot(std::string filename_string = "plot.input"){
  SetUp();
  //
  int nPoints=0;
  TString title,expname[99],tempstring;
  Double_t xmin,xmax,precision,meas[99],stath[99],statl[99],stat[99],systh[99],systl[99],syst[99];
  int statasy[99],systasy[99];
  const char* filename = filename_string.c_str();
  ifstream ifs(filename) ; if (!ifs.good()) {cout << "Cannot open input file '" << filename << "'" << endl ; exit(1) ;}
  char buffer[200] ; 
  while(ifs.good()) {
    if (ifs.eof()) break;
    char firstch ; ifs.get(firstch) ;
    if (firstch=='#'||firstch=='\n') { // Skip this line
      ifs.ignore(200,'\n') ;
    } else if (firstch=='*') {  // First line with xmin, xmax, precision and title should have '*' as the first character
      // Parse content : xmin xmax precision
      ifs >> xmin >> xmax >> precision;
      ifs.getline(buffer,200,'\n');
      title=TString(buffer).Strip((TString::EStripType)1,' ');//remove kLeading whitespace
    } else if (firstch=='&') {  // Next lines with meas, error, CL (or -ScaleFactor) HFAG Average Description should have '&' as the first character
      // Parse content : meas, error, CL (or -ScaleFactor)
      ifs >>  meas[nPoints] >> stat[nPoints] >> syst[nPoints];
      statasy[nPoints] = 0; // always false
      systasy[nPoints] = 0; // always false
      ifs.getline(buffer,200,'\n');
      expname[nPoints++]=TString(buffer).Strip((TString::EStripType)1,' ');//remove kLeading whitespace
    } else if (firstch=='%') {  // Next lines with meas, error, ScaleFactor (or 0) PDG Average Description should have '%' as the first character
      // Parse content : meas, error, ScaleFactor (or 0)
      ifs >>  meas[nPoints] >> stat[nPoints] >> syst[nPoints];
      statasy[nPoints] = 0; // always false
      systasy[nPoints] = 0; // always false
      ifs.getline(buffer,200,'\n');
      expname[nPoints++]=TString(buffer).Strip((TString::EStripType)1,' ');//remove kLeading whitespace
    } else { // Actual Inputs of the measurments 
      // Put back first char
      ifs.putback(firstch) ;
      // Parse content : meas stath statl[with negative sign] systh systl[with negative sign]
      ifs >> meas[nPoints] >> stath[nPoints] >> statl[nPoints] >> systh[nPoints] >> systl[nPoints];
      statl[nPoints]*=-1;
      systl[nPoints]*=-1;
      stat[nPoints] = 0.5 * (stath[nPoints] + statl[nPoints]);
      syst[nPoints] = 0.5 * (systh[nPoints] + systl[nPoints]);
      statasy[nPoints] = (stath[nPoints]==statl[nPoints]) ? 0: 1 ; // true if asymmetric
      systasy[nPoints] = (systh[nPoints]==systl[nPoints]) ? 0: 1 ; // true if assymetric
      ifs.getline(buffer,200,'\n') ;
      expname[nPoints]=TString(buffer).Strip((TString::EStripType)1,' ');//remove kLeading whitespace
      if (expname[nPoints].Length()) nPoints++;
    }
  }
  cout << "Read " << nPoints << " lines from filename " << filename << endl ;
  for (int i=0;i<nPoints;i++) {
    if (expname[i].Contains("Average") && expname[i].Contains("HFAG")) {
      cout << i << " " << meas[i] << " +- " << stat[i] << " (error) with CL (or -Scale Factor) = " << syst[i] << " from " << expname[i] << endl;
    } else if (expname[i].Contains("Average") && expname[i].Contains("PDG")) {
      cout << i << " " << meas[i] << " +- " << stat[i] << " (error) with Scale Factor = " << syst[i] << " from " << expname[i] << endl;
    } else {
      cout << i << " " << meas[i] << " +- " << stat[i] << " (stat) +- " << syst[i] << " (syst) from " << expname[i] << endl;
    }
  }
  //
  TString sprecision = Form("%s%3.1f%s","%",precision,"f");
  cout << "sprecision = " << sprecision << endl;
  //
  TCanvas *canvas=new TCanvas("canvas","canvas",200,0,600,600);
  canvas->Clear();
  //  canvas->SetFillColor(10);
  canvas->SetLeftMargin(0.04);
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.04);
  canvas->SetBottomMargin(0.12);  

  // -- set up frame for plot
  Double_t fXmin = xmin; 
  Double_t fXmax = xmax; 
  Double_t fYmin = 0.0;
  Double_t fYmax = nPoints*1.0 + 1.5;

  cout << "Drawing frame " << fXmin << "  " << fYmin << "  " << fXmax << "  " << fYmax << endl;

  TH1F *h = canvas->DrawFrame(fXmin, fYmin, fXmax, fYmax, "");  h->SetXTitle(title);  h->Draw();
  cout << "h drawn with x-axis title = " << title << endl;
  
  int ci = 1756; // new color index
  TColor *color = new TColor(ci, 0.50, 1.00, 0.50); // lightgreen for average

  // draw boxes using the default final average value listed at position 0
  TBox b1;
  b1.SetFillStyle(1000);
  b1.SetFillColor(ci);
  tempstring=Form(sprecision.Data(),meas[0]-stat[0]); Double_t boxl=tempstring.Atof();
  tempstring=Form(sprecision.Data(),meas[0]+stat[0]); Double_t boxh=tempstring.Atof();
  b1.DrawBox(boxl, fYmin, boxh, fYmax);

  Double_t y[99], ey[99], eyl[99], eyh[99];
  int fColor[99], fSymbol[99];
  double fMarkerSize = 1.2;
  int fMarkerStyle = 24;
  int fLineW=1;
  for (int i=0;i<nPoints;++i) {
    y[i]=i+2.0; ey[i]=0.0; eyl[i]=0.0; eyh[i]=0.0;
    fColor[i]  = 1;
    fSymbol[i] = 20;
  }

  Double_t fTxtX= 0.05;
  Double_t xtext=fXmin+fTxtX*(fXmax-fXmin);
  Double_t fTxtY= 0.4;
  Double_t ytext;
  TLatex  tl;
  tl.SetTextSize(0.028);
  for (int i=0;i<nPoints;++i) {
    Double_t totl[1],toth[1];
    if (expname[i].Contains("Average")) {
      totl[0] = stat[i];
      toth[0] = stat[i];
    }else{
      totl[0] = sqrt(statl[i]*statl[i] + systl[i]*systl[i]);
      toth[0] = sqrt(stath[i]*stath[i] + systh[i]*systh[i]);
    }
    //
    tempstring=Form(sprecision.Data(),meas[i]); meas[i]=tempstring.Atof();
    tempstring=Form(sprecision.Data(),totl[0]); totl[0]=tempstring.Atof();
    tempstring=Form(sprecision.Data(),toth[0]); toth[0]=tempstring.Atof();
    TGraph* GraphA = new TGraphAsymmErrors(1,&meas[i],&y[i],&totl[0],&toth[0],&eyl[i],&eyh[i]);
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
    tempstring=Form(sprecision.Data(),statl[0]); statl[0]=tempstring.Atof();
    tempstring=Form(sprecision.Data(),stath[0]); stath[0]=tempstring.Atof();
    TGraph* GraphB = new TGraphAsymmErrors(1,&meas[i],&y[i],&statl[i],&stath[i],&eyl[i],&eyh[i]);
    GraphB->SetMarkerColor(fColor[i]);
    GraphB->SetMarkerStyle(fSymbol[i]);
    GraphB->SetMarkerSize(fMarkerSize);
    GraphB->SetLineColor(fColor[i]);
    GraphB->SetLineWidth(fLineW);
    if (!expname[i].Contains("Average")&&fabs(stat[i])>1e-9) GraphB->Draw("P");
    //
    tl.SetTextColor(fColor[i]);
    // -- Change if it's the average
    if (expname[i].Contains("Average")) {
      tl.SetTextColor(2);
    }
    ytext=y[i]*1.0;
    tl.DrawLatex(xtext,ytext,expname[i]);
    if (expname[i].Contains("Average") && expname[i].Contains("PDG")) {
      if (syst[i]<1e-9) { // zero means no Scale Factor quoted
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #pm %s",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stat[i])));
      } else {
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #pm %s (Scale Factor = %s)",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stat[i]),
					    Form("%4.2f",syst[i])));
      }
    } else if (expname[i].Contains("Average") && expname[i].Contains("HFAG")) {
      if (fabs(syst[i])<1e-9) { // zero means no CL quoted
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #pm %s",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stat[i])));
      } else if (syst[i]<0) { // negative means Scale Factor is quoted
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #pm %s (Scale Factor = %s)",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stat[i]),
					    Form("%4.2f",fabs(syst[i]))));
      } else { // positive means CL is quoted
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #pm %s (CL = %s)",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stat[i]),
					    Form("%.2g",syst[i])));
      }
    } else {
      if (!statasy[i]&&syst[i]<1e-9) { // zero means no syst-error quoted
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #pm %s",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stat[i])));
      } else if (statasy[i]&&syst[i]<1e-9) { // zero means no syst-error quoted
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #splitline{+%s}{ -%s}",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stath[i]),
					    Form(sprecision.Data(),statl[i])));
      } else if ( statasy[i]&&!systasy[i]) {
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #splitline{+%s}{ -%s} #pm %s",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stath[i]),
					    Form(sprecision.Data(),statl[i]),
					    Form(sprecision.Data(),syst[i])));
      } else if (!statasy[i]&& systasy[i]) {
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #pm %s #splitline{+%s}{ -%s}",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stat[i]),
					    Form(sprecision.Data(),systh[i]),
					    Form(sprecision.Data(),systl[i])));
      } else if ( systasy[i]&& systasy[i]) {
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #splitline{+%s}{ -%s} #splitline{+%s}{ -%s}",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stath[i]),
					    Form(sprecision.Data(),statl[i]),
					    Form(sprecision.Data(),systh[i]),
					    Form(sprecision.Data(),systl[i])));
      } else {
	tl.DrawLatex(xtext,ytext-fTxtY,Form("%s #pm %s #pm %s",
					    Form(sprecision.Data(),meas[i]),
					    Form(sprecision.Data(),stat[i]),
					    Form(sprecision.Data(),syst[i])));
      }
    }
  }

  HFAGTauLabel(2010001,.3,.225,.2,.05,.8);

  size_t extPos = filename_string.rfind('.');
  std::string outFilename;
  if (extPos != std::string::npos) {
    // Erase the current extension.
    outFilename.assign(filename_string, 0, extPos);
    // Add the new extension.
    outFilename.append(".eps");
  }

  canvas->SaveAs(outFilename.c_str());
  //  canvas->SaveAs("plot.gif"); // ROOT's native gif saving looks garbled
  //  gSystem->Exec("rm -f plot.gif; convert plot.eps plot.gif");
  //canvas->Clear();

  //
  return 0;
}
