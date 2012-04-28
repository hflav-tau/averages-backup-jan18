//
// aluAverPlot.cc
// root -b -q "aluAverPlot.cc+(\"plot.input\")"
//
// plot measurements and averages listed in ASCII .input file
// inspired by S.Banerjee code for HFAG-tau plots
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TLine.h"
#include "TBox.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TColor.h"
#include "TPaveText.h"

#include "HFAGTauLabel.cc"

void SetUp(){
 // gROOT->SetBatch(kTRUE);

  gROOT->SetStyle("Pub");

  //--- canvas size, normal plots
  // gStyle->SetCanvasDefW(700);
  // gStyle->SetCanvasDefH(680);

  //--- canvas size, average plots
  gStyle->SetCanvasDefW(500);
  gStyle->SetCanvasDefH(Int_t(0.8*500));

  //--- normal plots
  // gStyle->SetPadLeftMargin(0.20);
  // gStyle->SetPadRightMargin(0.06);
  // gStyle->SetPadTopMargin(0.08);
  // gStyle->SetPadBottomMargin(0.15);

  //--- average plots
  gStyle->SetPadLeftMargin(0.05);
  gStyle->SetPadRightMargin(0.60);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetLineWidth(3);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetEndErrorSize(8);

  gStyle->SetFuncColor(kBlue);
  gStyle->SetLineColor(kBlack);
  gStyle->SetStatColor(kWhite);

  gStyle->SetTitleOffset(1.1, "x");
  gStyle->SetTitleOffset(1.6, "y");

  gStyle->SetLabelOffset(0.01, "xyz");
  gStyle->SetTickLength(0.06, "x");
  gStyle->SetNdivisions(505, "x");

  //--- do not print y axis
  gStyle->SetAxisColor(kWhite, "y");
  gStyle->SetLabelSize(0, "y");
  gStyle->SetLabelColor(kWhite, "y");
  gStyle->SetTickLength(0, "y");

  const int plotFont(42);
  gStyle->SetLabelFont(plotFont, "xyz");
  gStyle->SetTitleFont(plotFont, "xyz");
  gStyle->SetTextFont(plotFont);
  gStyle->SetStatFont(plotFont);  

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.5);
  gStyle->SetMarkerColor(kBlack);

  gStyle->SetStatH(0.10);
  gStyle->SetStatW(0.16);
  gStyle->SetFitFormat(".4g");
  gStyle->SetStatFormat(".4g");

  // gStyle->SetOptTitle(0);
  gStyle->SetOptStat(000000);
  gStyle->SetOptFit(000000);

  //--- average specific
  gStyle->SetTitleTextColor(kBlack);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetCanvasBorderSize(0);

  gROOT->ForceStyle();
}

// ////////////////////////////////////////////////////////////////////////////

void aluAverPlot(std::string filename = "plot.input");

void aluAverPlot(std::string filename){
  SetUp();
  //
  Int_t nPoints(0);
  TString title,expname[99],tempstring;
  Double_t precision;

  // -- set up frame for plot
  Double_t fXmin;
  Double_t fXmax;
  Double_t fYmin(0);
  Double_t fYmax(6);

  Double_t meas[99];
  Double_t stath[99], statl[99], stat[99];
  Double_t systh[99], systl[99], syst[99];
  Double_t errorh[99], errorl[99], error[99];
  Double_t meascl[99];
  Double_t sfact[99];
  TString measinfo[99];
  ifstream ifs(filename.c_str()) ; if (!ifs.good()) {cout << "Cannot open input file '" << filename << "'" << endl ; exit(1) ;}
  char buffer[200] ; 
  int ref_meas(0);

  while(ifs.good()) {
    if (ifs.eof()) break;
    char first_ch; ifs.get(first_ch) ;

    if (first_ch=='#'||first_ch=='\n') { // Skip this line
      ifs.ignore(200,'\n') ;

    } else if (first_ch=='i') {
      //--- plot limits and title: i <xmin> <xmax> <format> <title>
      ifs >> fXmin >> fXmax >> precision;
      ifs.getline(buffer, 200, '\n');
      title = TString(buffer).Strip((TString::EStripType)1, ' '); // remove kLeading whitespace
      
    } else if (first_ch=='h' || first_ch=='p') {
      //--- other average: h <mean> <sigma> <s-factor> <CL> <label>
      //--- pdg average: p <mean> <sigma> <s-factor> <CL> <label> -- other average
      if (first_ch=='p') ref_meas = nPoints;
      //--- get averages: mean, sigma, sfact, CL
      measinfo[nPoints] = first_ch;
      ifs >>  meas[nPoints] >> stath[nPoints] >> sfact[nPoints] >> meascl[nPoints];
      statl[nPoints] = stath[nPoints];
      systh[nPoints] = systl[nPoints] = 0;
      ifs.getline(buffer,200,'\n');
      expname[nPoints++]=TString(buffer).Strip((TString::EStripType)1, ' '); // remove kLeading whitespace

    } else if (first_ch=='m') {
      //
      // measurements
      // - m   symm errors
      // - ma  asymm errors
      // - mt  "this work" symm errors
      // - mat "this work" asymm errors
      //
      measinfo[nPoints] = first_ch;
      bool asymm_input = false;
      char ch2 ; ifs.get(ch2);

      if (ch2 == 'a') {
	asymm_input = true;
	char ch3 ; ifs.get(ch3);
	if (ch3=='t') {
	  measinfo[nPoints] = 't';
	} else {
	  ifs.putback(ch3);
	}
      } else if (ch2 == 't') {
	measinfo[nPoints] = 't';
      } else {
	ifs.putback(ch2);
      }

      //--- parse input
      if (asymm_input) {
	ifs >> meas[nPoints] >> stath[nPoints] >> statl[nPoints] >> systh[nPoints] >> systl[nPoints];
	stath[nPoints] = fabs(stath[nPoints]);
	statl[nPoints] = fabs(statl[nPoints]);
	systh[nPoints] = fabs(systh[nPoints]);
	systl[nPoints] = fabs(systl[nPoints]);
      } else {
	ifs >> meas[nPoints] >> stath[nPoints] >> systh[nPoints];
	statl[nPoints] = stath[nPoints];
	systl[nPoints] = systh[nPoints];
      }

      ifs.getline(buffer,200,'\n') ;
      expname[nPoints]=TString(buffer).Strip((TString::EStripType)1, ' '); // remove kLeading whitespace
      if (expname[nPoints].Length()) nPoints++;

    } else { // discard
      ifs.ignore(200,'\n') ;
    }
  }
  cout << "Read " << nPoints << " lines from filename " << filename.c_str() << endl ;

  int nexp=0;
  double num=0,den=0,aver=0,err=0,wt;
  Double_t meas_xmin(meas[0]);
  Double_t meas_xmax(meas[0]);
  Double_t meas_width_min(stath[0]+statl[0]+systh[0]+systl[0]);
  Double_t meas_width_xmin(meas[0]);
  Double_t meas_width_xmax(meas[0]);
  for (int i=0;i<nPoints;++i) {

    stat[i] = 0.5*(stath[i]+statl[i]);
    syst[i] = 0.5*(systh[i]+systl[i]);
    
    errorh[i] = sqrt(pow(stath[i], 2) + pow(systh[i], 2));
    errorl[i] = sqrt(pow(statl[i], 2) + pow(systl[i], 2));
    error[i] = 0.5*(errorh[i]+errorl[i]);

    if (measinfo[i] == "h" || measinfo[i] == "p") {
      if (measinfo[i] == "h") {
	cout << i << " hfag ave  ";
      } else {
	cout << i << " pdg  ave  ";
      }
      cout << meas[i] << " +- " << stath[i] << " sfact " << sfact[i] << " CL " << meascl[i]
	   << " " << expname[i] << endl;
    } else if (measinfo[i] == "m" || measinfo[i] == "t") {
      if (measinfo[i] == "m") {
	cout << i << " meas      " << meas[i];
      } else {
	cout << i << " this work " << meas[i];
      }
      if (stath[i] == statl[i]) {
	cout << " +- " << stath[i];
      } else {
	cout << " +" << stath[i] << "-" << statl[i];
      }
      if (systh[i] == systl[i]) {
	cout << " +- " << systh[i];
      } else {
	cout << " +" << systh[i] << "-" << systl[i];
      }
      cout << " " << expname[i] << endl;

      wt=1./(stat[i]*stat[i]+syst[i]*syst[i]);
      num+=meas[i]*wt;
      den+=wt;
      nexp++;
    } else {
      cout << "+++ error, expname = " << expname[i] << endl;
    }

    //--- compute min/max x position of drawn bars
    Double_t xmintmp = meas[i] - errorl[i];
    Double_t xmaxtmp = meas[i] + errorh[i];
    if (xmintmp < meas_xmin) meas_xmin = xmintmp;
    if (xmaxtmp > meas_xmax) meas_xmax = xmaxtmp;
    
    Double_t meas_width_min_tmp(errorh[i]+errorl[i]);
    // cout <<  meas_width_min_tmp << " " << meas_width_min << endl;
    if (meas_width_min_tmp < meas_width_min) {
      meas_width_min = meas_width_min_tmp;
      meas_width_xmin = meas[i] - errorl[i];
      meas_width_xmax = meas[i] + errorh[i];
    }
  }

  if (den>0) {
    aver=num/den;
    err =sqrt(1/den);
  } else {
    aver = 0;
    err = 1;
  }

  //--- if xmin >= xmax then use auto-computed x axis limits
  if (fXmin >= fXmax) {
    //--- set x axis limits where drawing arrives plus 10% of the total width
    fXmin = meas_xmin - 0.10 * (meas_xmax - meas_xmin);
    fXmax = meas_xmax + 0.10 * (meas_xmax - meas_xmin);

    //--- avoid excessive x axis limits when errors are large
    const Double_t max_width_mult(8);
    if (fXmin < meas_width_xmin - max_width_mult*meas_width_min) {
      fXmin = meas_width_xmin - max_width_mult*meas_width_min;
    }
    if (fXmax > meas_width_xmax + max_width_mult*meas_width_min) {
      fXmax = meas_width_xmax + max_width_mult*meas_width_min;
    }
    // cout << meas_width_min << " " << meas_width_xmax << " " << meas_width_xmax + max_width_mult*meas_width_min << endl;
  }
  
  //--- if more than 6 measurements, increase y axis limit
  fYmax = (nPoints < 6 ? 6 : nPoints) + 0.6;
  
  cout << "x min/max: " << fXmin << " " << fXmax
       << " meas min/max: " << meas_xmax << " " << meas_xmin << endl;
  cout << "meas num: " << nexp << ", average: " << aver << " +- " << err << endl;

  TString sprecision = Form("%s%3.1f%s", "%", precision, "f");
  cout << "sprecision = " << sprecision << endl;

  //
  // Canvas
  //
  TCanvas* canvas;
  {
    Int_t canvas_height;

    canvas_height = Int_t(0.8*gStyle->GetCanvasDefW() + (nPoints > 8 ? 75*(nPoints-8) : 0));
    canvas_height = Int_t(float(canvas_height) * float(1 - 0.02 - 0.06)/float(1 -0.02 - 0.14));
    gStyle->SetCanvasDefH(canvas_height);
    
    if (title != "") {
      canvas_height = Int_t(float(canvas_height) * float(1 - 0.02 - 0.06)/float(1 -0.02 - 0.14));
    }
    canvas = new TCanvas("canvas", "canvas");
  }
  
  //--- adjust pad margin if no title
  if (title == "") {
    canvas->SetBottomMargin(0.06);
  } else {
    canvas->SetBottomMargin(0.15);
  }
  // canvas->SetGrid(1, 1);

  cout << "Drawing frame: " << fXmin << ", " << fYmin << ", " << fXmax << ", " << fYmax
       << ", x axis title=\"" << title << "\"" << endl;

  TH1F *hpx = canvas->DrawFrame(fXmin, fYmin, fXmax, fYmax, "");
  hpx->GetXaxis()->CenterTitle(kTRUE);
  hpx->SetXTitle(title);

  //--- draw band for average measurement
  TBox average_band;
  average_band.SetFillStyle(1000);
  // average_band.SetFillColor(kGreen-9);  
  // average_band.SetFillColor(kOrange-9);  
  average_band.SetFillColor(kYellow-9);  
  average_band.DrawBox(meas[ref_meas]-statl[ref_meas], fYmin,
		       meas[ref_meas]+stath[ref_meas], fYmax);

  //--- redraw axis above the average band
  canvas->RedrawAxis();

  Double_t y[nPoints], ey[nPoints], eyl[nPoints], eyh[nPoints];
  int fColor[nPoints], fSymbol[nPoints];

  for (int i=0; i<nPoints; ++i) {
    //--- vertical position
    y[i] = double(nPoints-i) -1 + 0.8;
    ey[i] = 0.0;
    eyl[i] = 0.0;
    eyh[i] = 0.0;
    fColor[i]  = kBlack;
    fSymbol[i] = 20;
  }

  //--- x label position at 107% of plot width
  const Double_t fTxtX(0.07);
  Double_t xtext(fXmax + fTxtX*(fXmax-fXmin));

  //--- y label base position at position dependent on plot height
  Double_t fTxtY(0.40 + 0.20*(fYmax-6.6)/6);
  Double_t ytext;

  //--- result label
  TLatex  tlexp;
  tlexp.SetTextSize(gStyle->GetLabelSize());

  //--- result value
  TLatex  tlval;
  tlval.SetTextSize(gStyle->GetLabelSize()*0.9);

  Int_t color;
  for (int i=0; i<nPoints; ++i) {
    color = kBlack;
    if (measinfo[i] == "p") {
      color = kBlue;
    } else if (measinfo[i] == "t") {
      color = kRed;
    }

    //--- total errors
    TGraph* GraphA = new TGraphAsymmErrors(1, &meas[i], &y[i], &errorl[i], &errorh[i], &eyl[i], &eyh[i]);
    GraphA->SetLineColor(color);
    GraphA->SetMarkerColor(color);
    GraphA->SetMarkerStyle(fSymbol[i]);
    GraphA->Draw("p");
    
    //--- stat errors
    if (systl[i] != 0 || systh[i] != 0) {
      TGraph* GraphB = new TGraphAsymmErrors(1, &meas[i], &y[i], &statl[i], &stath[i], &eyl[i], &eyh[i]);
      GraphB->SetLineColor(color);
      GraphB->Draw("]");
    }

    tlexp.SetTextColor(color);
    tlval.SetTextColor(color);

    ytext = y[i];
    tlexp.DrawLatex(xtext, ytext, expname[i]);
    // cout << "ytext " << ytext << " " << i << endl;
    
    TString valtxt(Form(sprecision.Data(),meas[i]));
    
    if (stath[i] == statl[i]) {
      valtxt += " #pm ";
      valtxt += Form(sprecision.Data(),stath[i]);
    } else {
      valtxt += " #splitline{+";
      valtxt += Form(sprecision.Data(),stath[i]);
      valtxt += "}{- ";
      valtxt += Form(sprecision.Data(),statl[i]);
      valtxt += "}";
    }
    
    if (systh[i] != 0 || systl[i] != 0) {
      if (systh[i] == systl[i]) {
	valtxt += " #pm ";
	valtxt += Form(sprecision.Data(),systh[i]);
      } else {
	valtxt += " #splitline{+";
	valtxt += Form(sprecision.Data(),systh[i]);
	valtxt += "}{- ";
	valtxt += Form(sprecision.Data(),systl[i]);
	valtxt += "}";
      }
    }
    
    if (measinfo[i] == "h" || measinfo[i] == "p") {
      if (sfact[i] >0) { // S-factor
	valtxt += Form(" (S = %3.1f)", sfact[i]);
      } else if (meascl[i] >0) {
	valtxt += Form(" (CL = %4.1f%%)", 100*meascl[i]);
      }
    }

    tlval.DrawLatex(xtext, ytext-fTxtY, valtxt);
    // cout << "value " << valtxt << endl;
  }
  
  HFAGTauLabel("Winter 2012", 0.5, 0.05);

  std::string basefname(filename);
  size_t extPos = basefname.rfind('.');
  if (extPos != std::string::npos) {
    //--- Erase the current extension.
    basefname.assign(filename, 0, extPos);
  }
  
  canvas->SaveAs(std::string(basefname + ".eps").c_str());
  // canvas->SaveAs(std::string(basefname + ".gif").c_str());
  canvas->SaveAs(std::string(basefname + ".png").c_str());
  // canvas->Clear();

  return;
}
