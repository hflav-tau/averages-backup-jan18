//
// lepton universality plot (weak charged current couplings)
//
// root -b -q tauLifeSMpred.cc++
//

void tauLifeSMpred()
{
  // TLine(0,0,0,0).Draw();
  // gPad->Update();
  // gPad->Clear();

  // new TCanvas("TauLifeCmp", "", 0, 0, 600, 400);
  new TCanvas("TauLifeCmp", "");
  // gPad->Update();

  // gPad->SetRightMargin(gStyle->GetPadRightMargin()*700/560);
  // gPad->SetLeftMargin(gStyle->GetPadLeftMargin()*1.1*700/560);

  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetRightMargin(gStyle->GetPadRightMargin()*0.2);
  gPad->SetLeftMargin(gStyle->GetPadLeftMargin()*1.2);

  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetTitleOffset(1.5, "y");

  Stat_t xmin(289);
  Stat_t xmax(292.0);
  Stat_t ymin(0.177);
  Stat_t ymax(0.179);

  TH2F* frame = new TH2F("frame", "",
			 100, xmin, xmax,
			 100, ymin, ymax);

  (frame->GetXaxis())->SetTitle("#tau_{#tau} (fs)");
  (frame->GetYaxis())->SetTitle("B_{#tau l}");

  // (frame->GetXaxis())->SetNdivisions(510);
  // (frame->GetYaxis())->SetNdivisions(510);
  frame->SetStats(0);
  frame->Draw();
  
  //--- PDG 2015
  Stat_t massTau(1776.86);
  Stat_t massTauPlus(0.12);
  Stat_t massTauMinus(0.12);

  //--- CODATA 2010 and PDG 2014 with 2015 update
  Stat_t massMu(105.658369);
  Stat_t massMuErr(0.0000035);

  //--- HFAG 2014
  Stat_t phspf_mebymtau(0.99999933833);
  Stat_t phspf_mmubymtau(0.9725588);
  Stat_t phspf_mebymmu(0.999812949174);

  Stat_t delta_mu_gamma(1 - 42.4e-4);
  Stat_t delta_tau_gamma(1 - 43.2e-4);

  Stat_t delta_mu_W(1.00000103622);
  Stat_t delta_tau_W(1.00029304);

  Stat_t tauEWcorrMu(phspf_mebymtau * delta_tau_W * delta_tau_gamma);
  Stat_t tauEWcorrTau(phspf_mebymmu * delta_mu_W * delta_mu_gamma);

  //
  // HFAG 2014
  // Be_lept                                     0.178493 0.000320
  // average of Be and Be_from_Bmu
  //
  Stat_t brTauMuEl(0.178493);
  Stat_t brTauMuElErr(0.000320);

  //--- PDG 2014 with 2015 update
  Stat_t tauMu(2.1969811e-6);
  Stat_t tauMuErr(0.0000022e-6);

  //--- PDG 2014 with 2015 update
  Stat_t tauTau(290.3e-15);
  Stat_t tauTauErr(0.5e-15);  

  //                   5
  //   B_tau_e     mTau  r_tau tauTau
  // ---------- = --------------------
  //   B_mu_e         5
  //               mMu   r_mu  tauMu
  
  cout
    << "BR tau -> e univ. impr. " << brTauMuEl << " +- " << brTauMuElErr
    << endl;

  Stat_t slope;
  Stat_t slopeMin;
  Stat_t slopeMax;

  slope = (pow(massTau,5) * tauEWcorrTau) / (pow(massMu,5) * tauEWcorrMu * tauMu);
  slopeMin = (pow(massTau-massTauMinus,5) * tauEWcorrTau) / (pow(massMu+tauMuErr,5) * tauEWcorrMu * (tauMu+tauMuErr));
  slopeMax = (pow(massTau+massTauMinus,5) * tauEWcorrTau) / (pow(massMu-tauMuErr,5) * tauEWcorrMu * (tauMu-tauMuErr));
  
  slope /= 1e15;
  slopeMin /= 1e15;
  slopeMax /= 1e15;

  cout
    << "slope " << slope << " " << slopeMin << " " << slopeMax
    << endl;

  Double_t tauMassBandX[4] = {ymin/slopeMax, ymin/slopeMin, xmax, xmax};
  Double_t tauMassBandY[4] = {ymin, ymin, xmax*slopeMin, xmax*slopeMax};
  TPolyLine *pline = new TPolyLine(4, tauMassBandX, tauMassBandY);
  pline->SetFillColor(5);
  pline->SetLineWidth(0);
  pline->Draw("f");

  frame->Draw("axissame");

  Float_t x[1]  = {tauTau*1e15};
  Float_t y[1]  = {brTauMuEl};
  Float_t ex[1] = {tauTauErr*1e15};
  Float_t ey[1] = {brTauMuElErr};
  gr = new TGraphErrors(1, x, y, ex, ey);
  gr->SetLineWidth(3);
  gr->SetLineColor(2);
  gr->Draw("+");

  // gStyle->SetTextSize(0.7*gStyle->GetTextSize());

  TLatex text;
  Float_t xPos, yPos;

  Font_t oldFont = text.GetTextFont();
  Float_t textSize = text.GetTextSize();

  text.SetTextSize(0.9*textSize);

  xPos = 290+0.05;
  yPos = 290*slope;
  text.SetText(xPos,  yPos, "#leftarrow SM lepton universality, 1#sigma");
  text.DrawClone();

  gPad->Update();
  (gPad->GetCanvas())->SaveAs("tauLifeSMpred.pdf");
}
