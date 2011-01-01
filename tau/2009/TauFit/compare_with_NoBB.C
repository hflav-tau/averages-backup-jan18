#include <fstream>
compare_with_NoBB(){
  cout << "SetUp Plain " << endl;
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
 // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //
  gStyle->SetNdivisions(7,"X");
  gStyle->SetNdivisions(6,"Y");
  gStyle->SetHatchesSpacing(1);
  //
  gStyle->SetTitleFont(12,"Y"); gStyle->SetTitleSize(.05,"Y"); gStyle->SetTitleOffset(1,"Y");
  gStyle->SetTitleFont(12,"X"); gStyle->SetTitleSize(.05,"X"); gStyle->SetTitleOffset(1,"X");
  //
  ifstream inFile;
  inFile.open("compare_with_NoBB.txt");
  TH1D *h0 = new TH1D("h0","BaBar;Sigma Difference;Number per 0.5",12,-3,3);
  TH1D *h1 = new TH1D("h1","Belle;Sigma Difference;Number per 0.5",12,-3,3);
  TH1D *h01= new TH1D("h01","Belle;Sigma Difference;Number per 0.5",12,-3,3);
  h0->GetXaxis()->CenterTitle(1);  h0->GetYaxis()->CenterTitle(1);  h0->SetFillColor(19);  h0->SetFillStyle(4000);
  h1->GetXaxis()->CenterTitle(1);  h1->GetYaxis()->CenterTitle(1);  h1->SetFillColor(10);  h1->SetFillStyle(4000);
  h01->GetXaxis()->CenterTitle(1); h01->GetYaxis()->CenterTitle(1); h01->SetFillColor(19); h01->SetFillStyle(4000);
  int ibb;
  double x;
  char buffer[256]; 
  while (inFile.good()) { 
    inFile >> ibb >> x;
    inFile.getline(buffer,256,'\n'); // ignore the first line
    if (ibb==0) {
      if (x<0){
	cout << "BaBar : " << x << " " << h0->FindBin(x) << " " << h0->GetBinLowEdge(h0->FindBin(x)) << " : " <<  h0->GetBinLowEdge(h0->FindBin(x))+0.5 << endl;
      } else {
	cout << "BaBar :  " << x << " " << h0->FindBin(x) << " " << h0->GetBinLowEdge(h0->FindBin(x)) << " : " <<  h0->GetBinLowEdge(h0->FindBin(x))+0.5 << endl;
      }
      h0->Fill(x);
      h01->Fill(x);
    } else {
      if (x<0) {
	cout << "Belle : " << x << " " << h0->FindBin(x) << " " << h0->GetBinLowEdge(h0->FindBin(x)) << " : " <<  h0->GetBinLowEdge(h0->FindBin(x))+0.5 << endl;
      } else {
	cout << "Belle :  " << x << " " << h0->FindBin(x) << " " << h0->GetBinLowEdge(h0->FindBin(x)) << " : " <<  h0->GetBinLowEdge(h0->FindBin(x))+0.5 << endl;
      }
      h1->Fill(x);
      h01->Fill(x);
    }
  }
  inFile.close();
  //
  TCanvas *c1 = new TCanvas("c1","c1",400,400);
  c1->cd();
  gPad->SetTopMargin(0.025);
  gPad->SetBottomMargin(0.1);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.025);
  h01->SetMaximum(6);
  h01->Draw();
  h1->Draw("same");
  gPad->RedrawAxis();
  c1->Update();
  //
  TLegend legend(0.75, 0.75, 0.9, 0.9); // x1, y1, x2, y2
  legend.SetEntrySeparation(.05);
  legend.SetMargin(.4);
  legend.SetTextFont(12);
  legend.SetTextSize(0.04);
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  legend.SetFillStyle(4000); // transparent
  legend.AddEntry(h01, "BaBar","f");
  legend.AddEntry(h1, "Belle","f");
  legend.Draw();
  c1->Update();
  c1->SaveAs("compare_with_NoBB.eps");
  //
  cout << "h0 (BaBar): Mean, RMS, Int = " << h0->GetMean() << " " << h0->GetRMS() << " " << h0->Integral() << endl;
  cout << "h1 (Belle): Mean, RMS, Int = " << h1->GetMean() << " " << h1->GetRMS() << " " << h1->Integral() << endl;
  cout << "h01 (Both): Mean, RMS, Int = " << h01->GetMean()<< " " << h01->GetRMS()<< " " << h01->Integral()<< endl;
}
