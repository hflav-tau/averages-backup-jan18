///////////////////////////////////////////////////////////////////////////////
//
// TauLFV_UL_plot.cc
//
// produce tau LFV upper limits plot
// execute as follows:
//
// > root -l -b -q TauLFV_UL_plot.cc
//
// This will create a .eps file with the plot.
//
///////////////////////////////////////////////////////////////////////////////

enum DECAY {
  EGAMMA = 1,
  MUGAMMA = 2,
  EPI0 = 3,
  MUPI0 = 4,
  EETA = 5,
  MUETA = 6,
  EETAP = 7,
  MUETAP = 8,
  EKS0 = 9,
  MUKS0 = 10,
  EF0 = 11, 
  MUF0 = 12,
  ERHO = 13,  
  MURHO = 14,
  EKSTAR = 15,
  MUKSTAR = 16,
  EAKSTAR = 17,
  MUAKSTAR = 18,
  EPHI = 19,
  MUPHI = 20,
  EOMEGA = 21,
  MUOMEGA = 22,
  EEE = 23, //order of charge -+-     
  MEE = 24, 
  EMM = 25,
  MMM = 26, 
  EME = 27,   
  MEM = 28, 
  EPIPI = 29,
  MUPIPI = 30,
  EPIK = 31,
  MUPIK = 32,
  EKPI = 33,
  MUKPI = 34,
  EKK = 35,
  MUKK = 36,
  EKS0KS0 = 37,
  MUKS0KS0 = 38,
  PIEPI = 39,
  PIMUPI = 40,
  PIEK = 41,
  PIMUK = 42,
  KEK = 43,
  KMUK = 44,
  PILAM = 45,
  PILAMBAR = 46,
  KLAM = 47,
  KLAMBAR = 48,
  NDECAY = 49
};

void SetUp();

void HFAGTauLabel(Int_t yearversion=2012001, Double_t xpos= .99, Double_t ypos= 0.825, Double_t scale= 1.0);

void TauLFV_UL_plot(Int_t when=2012001)
{
  SetUp();
  unsigned int ibin;
  
  double Belle[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) Belle[ibin] = 0;
  
  double BaBar[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) BaBar[ibin] = 0;
  
  double CLEO[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) CLEO[ibin] = 0;
  
  double histbins[NDECAY];
  for (ibin = 0; ibin < NDECAY; ++ibin) histbins[ibin] = ibin+1;
  
  if (when==2009001){
    fillBelle_2009001(Belle);
    fillBaBar_2009001(BaBar);
  }else if (when==2009002) {
    fillBelle_2009002(Belle);
    fillBaBar_2009002(BaBar);
  }else if (when==2010001) {
    fillBelle_2010001(Belle);
    fillBaBar_2010001(BaBar);
  }else if (when==2012001) {
    fillBelle_2012001(Belle);
    fillBaBar_2010001(BaBar);
  }
  fillCLEO(CLEO);
  
  TH1F * hBelle = new TH1F("hBelle","",NDECAY+1, -0.5, float(NDECAY)+0.5);
  TH1F * hBaBar = new TH1F("hBaBar","",NDECAY+1, -0.5, float(NDECAY)+0.5);
  TH1F * hCLEO  = new TH1F("hCLEO","", NDECAY+1, -0.5, float(NDECAY)+0.5);
  
  hBelle->FillN(NDECAY, histbins, Belle);
  hBaBar->FillN(NDECAY, histbins, BaBar);
  hCLEO-> FillN(NDECAY, histbins, CLEO);
  
  hBelle->Scale(1.e-8);
  hBaBar->Scale(1.e-8);
  hCLEO->Scale(1.e-6);
  
  setLabels(hBelle);
  
  hBelle->GetYaxis()->SetTitle("90% C.L. upper limits for LFV #tau decays");
  hBelle->LabelsOption("v");
  hBelle->GetYaxis()->SetLabelSize(36);

  hBelle->SetMarkerStyle(kOpenCircle); hBelle->SetMarkerColor(kRed); hBelle->SetMarkerSize(1.0);
  hBaBar->SetMarkerStyle(kOpenDiamond); hBaBar->SetMarkerColor(kBlue); hBaBar->SetMarkerSize(1.3);
  hCLEO->SetMarkerStyle(kOpenSquare);  hCLEO->SetMarkerColor(kGreen);  hCLEO->SetMarkerSize(0.8);
  
  TCanvas *c1 = new TCanvas("c1","",1500,800); 
  c1->SetBottomMargin(.2);
  c1->SetTopMargin(.02);
  c1->SetLeftMargin(.10);
  c1->SetRightMargin(.12);
  c1->SetLogy(1) ;
  
  double ul_min=4.e-5;
  hBelle->SetMaximum(ul_min);
  double ul_max=4.e-9;
  hBelle->SetMinimum(ul_max);
  hBelle->Draw("p");
  c1->SetGridy();
  hBaBar->Draw("p,same");
  hCLEO->Draw("p,same");
  c1->Update();
  
  double y_latex=2.e-5;
  
  TLine l0(0+.5,ul_min,0+.5,ul_max); l0.SetLineStyle(2); l0.Draw();

  TLine l1(MUGAMMA+.5,ul_min,MUGAMMA+.5,ul_max); l1.SetLineStyle(2); l1.Draw();
  TLatex t1(((MUGAMMA-EGAMMA)*1.0/2.0)+EGAMMA*1.0,y_latex,"l#gamma"); t1.SetTextAlign(21); t1.SetTextFont(12); t1.Draw();
  
  TLine l2(MUKS0+.5,ul_min,MUKS0+.5,ul_max); l2.SetLineStyle(2); l2.Draw();
  TLatex t2(((MUKS0-EPI0)*1.0/2.0)+EPI0*1.0,y_latex,"lP^{0}"); t2.SetTextAlign(21); t2.SetTextFont(12); t2.Draw();
  
  TLine l3(MUF0+.5,ul_min,MUF0+.5,ul_max); l3.SetLineStyle(2); l3.Draw();
  TLatex t3(((MUF0-EF0)*1.0/2.0)+EF0*1.0,y_latex,"lS^{0}"); t3.SetTextAlign(21); t3.SetTextFont(12); t3.Draw();
  
  TLine l4(MUOMEGA+.5,ul_min,MUOMEGA+.5,ul_max); l4.SetLineStyle(2); l4.Draw();
  TLatex t4(((MUOMEGA-ERHO)*1.0/2.0)+ERHO*1.0,y_latex,"lV^{0}"); t4.SetTextAlign(21); t4.SetTextFont(12); t4.Draw();
  
  TLine l5(MEM+.5,ul_min,MEM+.5,ul_max); l5.SetLineStyle(2); l5.Draw();
  TLatex t5(((MEM-EEE)*1.0/2.0)+EEE*1.0,y_latex,"lll"); t5.SetTextAlign(21); t5.SetTextFont(12); t5.Draw();
  
  TLine l6(KMUK+.5,ul_min,KMUK+.5,ul_max); l6.SetLineStyle(2); l6.Draw();
  TLatex t6(((KMUK-EPIPI)*1.0/2.0)+EPIPI*1.0,y_latex,"lhh"); t6.SetTextAlign(21); t6.SetTextFont(12); t6.Draw();
  
  TLine l7(KLAMBAR+.5,ul_min,KLAMBAR+.5,ul_max); l7.SetLineStyle(2); l7.Draw();
  TLatex t7(((KLAMBAR-PILAM)*1.0/2.0)+PILAM*1.0,y_latex,"#Lambdah"); t7.SetTextAlign(21); t7.SetTextFont(12); t7.Draw();
  
  TLegend *leg = new TLegend(0.88,0.3,0.95,.5);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4000);
  leg->SetFillColor(0);
  leg->SetTextSize(.05);
  leg->SetMargin(.5);
  leg->AddEntry(hCLEO,"CLEO","p");
  leg->AddEntry(hBaBar,"BaBar","p");
  leg->AddEntry(hBelle,"Belle","p");
  leg->Draw();
  c1->Update();
  //
  HFAGTauLabel(when,.99,.825,1.0);
  c1->Update();  
  
  c1->SaveAs(Form("TauLFV_UL_%d.eps",when));
}

void fillBelle_2009001(double * array)
{
  // at 10^-8
  array[-1+EGAMMA]  = 12;
  array[-1+MUGAMMA] = 4.5 + .2; // added .2 for clarity in display
  
  array[-1+EPI0]    = 8.0; 
  array[-1+MUPI0]   = 12;
  
  array[-1+EETA]    = 9.2;
  array[-1+MUETA]   = 6.5;
  
  array[-1+EETAP]   = 16;
  array[-1+MUETAP]  = 13;
  
  array[-1+EKS0]    = 5.6;//2.6;
  array[-1+MUKS0]   = 4.8;//2.3;
  
  array[-1+EF0]     = 3.2;
  array[-1+MUF0]    = 3.4;
  
  array[-1+ERHO]    = 6.3;
  array[-1+MURHO]   = 6.8;
  array[-1+EKSTAR]  = 7.8;
  array[-1+MUKSTAR] = 5.9;
  array[-1+EAKSTAR] = 7.7;
  array[-1+MUAKSTAR]= 10;
  array[-1+EPHI]    = 7.3;
  array[-1+MUPHI]   = 13;
  array[-1+EOMEGA]  = 18;
  array[-1+MUOMEGA] = 8.9;
  
  array[-1+EEE]     = 3.6;//2.7;
  array[-1+MEE]     = 2.7;//1.8;
  array[-1+EMM]     = 4.1;//2.7; 
  array[-1+MMM]     = 3.2;//2.1; 
  array[-1+EME]     = 2.0;//1.5; 
  array[-1+MEM]     = 2.3;//1.7; 
  
  array[-1+EPIPI]   = 4.4;
  array[-1+MUPIPI]  = 3.3;
  array[-1+EPIK]    = 5.8;
  array[-1+MUPIK]   = 16;
  array[-1+EKPI]    = 5.2;
  array[-1+MUKPI]   = 10;
  array[-1+EKK]     = 5.4;
  array[-1+MUKK]    = 6.8;
  //  array[-1+EKS0KS0] = 7.1;
  //  array[-1+MUKS0KS0]= 8.0;
  array[-1+PIEPI]   = 8.8;
  array[-1+PIMUPI]  = 3.7;
  array[-1+PIEK]    = 6.7;
  array[-1+PIMUK]   = 9.4;
  array[-1+KEK]     = 6.0;
  array[-1+KMUK]    = 9.6;
  
  array[-1+PILAM]   = 7.2;
  array[-1+PILAMBAR]= 14;
}

void fillBelle_2009002(double * array)
{
  // at 10^-8
  array[-1+EGAMMA]  = 12;
  array[-1+MUGAMMA] = 4.5 + .3; // added .3 for clarity in display
  
  array[-1+EPI0]    = 8.0; 
  array[-1+MUPI0]   = 12;
  
  array[-1+EETA]    = 9.2;
  array[-1+MUETA]   = 6.5;
  
  array[-1+EETAP]   = 16;
  array[-1+MUETAP]  = 13;
  
  array[-1+EKS0]    = 2.6;
  array[-1+MUKS0]   = 2.3;
  
  array[-1+EF0]     = 3.2;
  array[-1+MUF0]    = 3.4;
  
  array[-1+ERHO]    = 6.3;
  array[-1+MURHO]   = 6.8;
  array[-1+EKSTAR]  = 7.8;
  array[-1+MUKSTAR] = 5.9;
  array[-1+EAKSTAR] = 7.7;
  array[-1+MUAKSTAR]= 10;
  array[-1+EPHI]    = 7.3;
  array[-1+MUPHI]   = 13;
  array[-1+EOMEGA]  = 18;
  array[-1+MUOMEGA] = 8.9;
  
  array[-1+EEE]     = 2.7;
  array[-1+MEE]     = 1.8;
  array[-1+EMM]     = 2.7; 
  array[-1+MMM]     = 2.1; 
  array[-1+EME]     = 1.5; 
  array[-1+MEM]     = 1.7; 
  
  array[-1+EPIPI]   = 4.4;
  array[-1+MUPIPI]  = 3.3;
  array[-1+EPIK]    = 5.8;
  array[-1+MUPIK]   = 16;
  array[-1+EKPI]    = 5.2;
  array[-1+MUKPI]   = 10;
  array[-1+EKK]     = 5.4;
  array[-1+MUKK]    = 6.8;
  array[-1+EKS0KS0] = 7.1;
  array[-1+MUKS0KS0]= 8.0;
  array[-1+PIEPI]   = 8.8;
  array[-1+PIMUPI]  = 3.7;
  array[-1+PIEK]    = 6.7;
  array[-1+PIMUK]   = 9.4;
  array[-1+KEK]     = 6.0;
  array[-1+KMUK]    = 9.6;
  
  array[-1+PILAM]   = 7.2;
  array[-1+PILAMBAR]= 14; 
}

void fillBelle_2010001(double * array)
{
  // at 10^-8
  array[-1+EGAMMA]  = 12;
  array[-1+MUGAMMA] = 4.5 + .3; // added .3 for clarity in display
  
  array[-1+EPI0]    = 2.2; // 8.0; 
  array[-1+MUPI0]   = 2.7; // 12;
  
  array[-1+EETA]    = 4.4; // 9.2;
  array[-1+MUETA]   = 2.3; // 6.5;
  
  array[-1+EETAP]   = 3.6; // 16;
  array[-1+MUETAP]  = 3.8; // 13;
  
  array[-1+EKS0]    = 2.6;
  array[-1+MUKS0]   = 2.3;
  
  array[-1+EF0]     = 3.2;
  array[-1+MUF0]    = 3.4;
  
  array[-1+ERHO]    = 1.8; // 6.3;
  array[-1+MURHO]   = 1.2; // 6.8;
  array[-1+EKSTAR]  = 3.2; // 7.8;
  array[-1+MUKSTAR] = 7.2; // 5.9;
  array[-1+EAKSTAR] = 3.4; // 7.7;
  array[-1+MUAKSTAR]= 7.0 - 0.2; // 10; // substracted .2 for clarity in display
  array[-1+EPHI]    = 3.1 - 0.2; // 7.3; // substracted .2 for clarity in display
  array[-1+MUPHI]   = 8.4; // 13;
  array[-1+EOMEGA]  = 4.8; // 18;
  array[-1+MUOMEGA] = 4.7; // 8.9;
  
  array[-1+EEE]     = 2.7;
  array[-1+MEE]     = 1.8;
  array[-1+EMM]     = 2.7; 
  array[-1+MMM]     = 2.1; 
  array[-1+EME]     = 1.5; 
  array[-1+MEM]     = 1.7; 
  
  array[-1+EPIPI]   = 4.4;
  array[-1+MUPIPI]  = 3.3;
  array[-1+EPIK]    = 5.8;
  array[-1+MUPIK]   = 16;
  array[-1+EKPI]    = 5.2;
  array[-1+MUKPI]   = 10;
  array[-1+EKK]     = 5.4;
  array[-1+MUKK]    = 6.8;
  array[-1+EKS0KS0] = 7.1;
  array[-1+MUKS0KS0]= 8.0;
  array[-1+PIEPI]   = 8.8;
  array[-1+PIMUPI]  = 3.7;
  array[-1+PIEK]    = 6.7;
  array[-1+PIMUK]   = 9.4;
  array[-1+KEK]     = 6.0;
  array[-1+KMUK]    = 9.6;
  
  array[-1+PILAM]   = 7.2;
  array[-1+PILAMBAR]= 14;
  
}

void fillBelle_2012001(double * array)
{
  // at 10^-8
  array[-1+EGAMMA]  = 12;
  array[-1+MUGAMMA] = 4.5 + .3; // added .3 for clarity in display
  //
  // Belle 901/fb data // previous 401/fb
  array[-1+EPI0]    = 2.2; // 8.0; 
  array[-1+MUPI0]   = 2.7; // 12;
  
  array[-1+EETA]    = 4.4; // 9.2;
  array[-1+MUETA]   = 2.3; // 6.5;
  
  array[-1+EETAP]   = 3.6; // 16;
  array[-1+MUETAP]  = 3.8; // 13;
  
  array[-1+EKS0]    = 2.6;
  array[-1+MUKS0]   = 2.3;
  
  array[-1+EF0]     = 3.2;
  array[-1+MUF0]    = 3.4;
  //
  // Belle 854 fb-1  //previous 543 fb-1
  array[-1+ERHO]    = 1.8; // 6.3;
  array[-1+MURHO]   = 1.2; // 6.8;
  array[-1+EKSTAR]  = 3.2; // 7.8;
  array[-1+MUKSTAR] = 7.2; // 5.9;
  array[-1+EAKSTAR] = 3.4; // 7.7;
  array[-1+MUAKSTAR]= 7.0 - 0.2; // 10; // substracted .2 for clarity in display
  array[-1+EPHI]    = 3.1 - 0.2; // 7.3; // substracted .2 for clarity in display
  array[-1+MUPHI]   = 8.4; // 13;
  array[-1+EOMEGA]  = 4.8; // 18;
  array[-1+MUOMEGA] = 4.7; // 8.9;
  
  array[-1+EEE]     = 2.7;
  array[-1+MEE]     = 1.8;
  array[-1+EMM]     = 2.7; 
  array[-1+MMM]     = 2.1; 
  array[-1+EME]     = 1.5; 
  array[-1+MEM]     = 1.7; 
  
  // 671 fb-1    //Previous   671 fb-1.
  array[-1+EPIPI]   = 2.3;  //4.4;
  array[-1+MUPIPI]  = 2.1;  //3.3;
  array[-1+EPIK]    = 3.7;  //5.8;
  array[-1+MUPIK]   = 8.6;   //16;
  array[-1+EKPI]    = 3.1;   //5.2;
  array[-1+MUKPI]   = 4.5;  //10;
  array[-1+EKK]     = 3.4;  //5.4;
  array[-1+MUKK]    = 4.4;  // 6.8;
  array[-1+EKS0KS0] = 7.1;
  array[-1+MUKS0KS0]= 8.0;
  array[-1+PIEPI]   = 2.0; //8.8;
  array[-1+PIMUPI]  = 3.9; //3.7;
  array[-1+PIEK]    = 3.2; //6.7;
  array[-1+PIMUK]   = 4.8;  //9.4;
  array[-1+KEK]     = 3.3;  //6.0;
  array[-1+KMUK]    = 4.7;  //9.6;
  
  //  906 fb-1 //Previous 154 fb-1
  array[-1+PILAM]   = 3.0; //7.2;
  array[-1+PILAMBAR]= 2.8; //14;
  array[-1+KLAM]    = 4.2;   
  array[-1+KLAMBAR] = 3.1; 
  
}

void fillBaBar_2009001(double * array)
{
  // at 10^-8
  array[-1+EGAMMA]  = 11;//3.3;
  array[-1+MUGAMMA] = 6.8;//4.4 - 0.2;
  
  array[-1+EPI0]    = 13; 
  array[-1+MUPI0]   = 11;
  
  array[-1+EETA]    = 16;
  array[-1+MUETA]   = 15;
  
  array[-1+EETAP]   = 24;
  array[-1+MUETAP]  = 14;
  
  array[-1+EKS0]    = 3.3;
  array[-1+MUKS0]   = 4.0;
  
  
  
  
  array[-1+ERHO]    = 4.6;
  array[-1+MURHO]   = 2.6;
  array[-1+EKSTAR]  = 5.9;
  array[-1+MUKSTAR] = 17;
  array[-1+EAKSTAR] = 4.6;
  array[-1+MUAKSTAR]= 7.3;
  array[-1+EPHI]    = 3.1;
  array[-1+MUPHI]   = 19;
  array[-1+EOMEGA]  = 11;
  array[-1+MUOMEGA] = 10;
  
  array[-1+EEE]     = 4.3;//2.9;
  array[-1+MEE]     = 8.0;//2.2;
  array[-1+EMM]     = 3.7;//3.2; 
  array[-1+MMM]     = 5.3;//3.3; 
  array[-1+EME]     = 5.8;//1.8; 
  array[-1+MEM]     = 5.6;//2.6; 
  
  array[-1+EPIPI]   = 12;
  array[-1+MUPIPI]  = 29;
  array[-1+EPIK]    = 32;
  array[-1+MUPIK]   = 26;
  array[-1+EKPI]    = 17;
  array[-1+MUKPI]   = 32;
  array[-1+EKK]     = 14;
  array[-1+MUKK]    = 25;
  
  
  array[-1+PIEPI]   = 27;
  array[-1+PIMUPI]  = 7.0;
  array[-1+PIEK]    = 18;
  array[-1+PIMUK]   = 22;
  array[-1+KEK]     = 15;
  array[-1+KMUK]    = 48;
  
  array[-1+PILAM]   = 5.8;
  array[-1+PILAMBAR]= 5.9;
  array[-1+KLAM]    = 15;
  array[-1+KLAMBAR] = 7.2;
}

void fillBaBar_2009002(double * array)
{
  // at 10^-8
  array[-1+EGAMMA]  = 3.3;
  array[-1+MUGAMMA] = 4.4 - 0.3; // substracted .3 for clarity in display
  
  array[-1+EPI0]    = 13; 
  array[-1+MUPI0]   = 11;
  
  array[-1+EETA]    = 16;
  array[-1+MUETA]   = 15;
  
  array[-1+EETAP]   = 24;
  array[-1+MUETAP]  = 14;
  
  array[-1+EKS0]    = 3.3;
  array[-1+MUKS0]   = 4.0;
  
  
  
  
  array[-1+ERHO]    = 4.6;
  array[-1+MURHO]   = 2.6;
  array[-1+EKSTAR]  = 5.9;
  array[-1+MUKSTAR] = 17;
  array[-1+EAKSTAR] = 4.6;
  array[-1+MUAKSTAR]= 7.3;
  array[-1+EPHI]    = 3.1;
  array[-1+MUPHI]   = 19;
  array[-1+EOMEGA]  = 11;
  array[-1+MUOMEGA] = 10;
  
  array[-1+EEE]     = 2.9;
  array[-1+MEE]     = 2.2;
  array[-1+EMM]     = 3.2; 
  array[-1+MMM]     = 3.3; 
  array[-1+EME]     = 1.8; 
  array[-1+MEM]     = 2.6; 
  
  array[-1+EPIPI]   = 12;
  array[-1+MUPIPI]  = 29;
  array[-1+EPIK]    = 32;
  array[-1+MUPIK]   = 26;
  array[-1+EKPI]    = 17;
  array[-1+MUKPI]   = 32;
  array[-1+EKK]     = 14;
  array[-1+MUKK]    = 25;
  
  
  array[-1+PIEPI]   = 27;
  array[-1+PIMUPI]  = 7.0;
  array[-1+PIEK]    = 18;
  array[-1+PIMUK]   = 22;
  array[-1+KEK]     = 15;
  array[-1+KMUK]    = 48;
  
  array[-1+PILAM]   = 5.8;
  array[-1+PILAMBAR]= 5.9;
  array[-1+KLAM]    = 15;
  array[-1+KLAMBAR] = 7.2;
}

void fillBaBar_2010001(double * array)
{
  // at 10^-8
  array[-1+EGAMMA]  = 3.3;
  array[-1+MUGAMMA] = 4.4 - 0.3; // substracted .3 for clarity in display
  
  array[-1+EPI0]    = 13; 
  array[-1+MUPI0]   = 11;
  
  array[-1+EETA]    = 16;
  array[-1+MUETA]   = 15;
  
  array[-1+EETAP]   = 24;
  array[-1+MUETAP]  = 14;
  
  array[-1+EKS0]    = 3.3;
  array[-1+MUKS0]   = 4.0;
  
  
  
  
  array[-1+ERHO]    = 4.6;
  array[-1+MURHO]   = 2.6;
  array[-1+EKSTAR]  = 5.9;
  array[-1+MUKSTAR] = 17;
  array[-1+EAKSTAR] = 4.6;
  array[-1+MUAKSTAR]= 7.3 + 0.2; //added .2 for clarity in display
  array[-1+EPHI]    = 3.1 + 0.2; //added .2 for clarity in display
  array[-1+MUPHI]   = 19;
  array[-1+EOMEGA]  = 11;
  array[-1+MUOMEGA] = 10;
  
  array[-1+EEE]     = 2.9;
  array[-1+MEE]     = 2.2;
  array[-1+EMM]     = 3.2; 
  array[-1+MMM]     = 3.3; 
  array[-1+EME]     = 1.8; 
  array[-1+MEM]     = 2.6; 
  
  array[-1+EPIPI]   = 12;
  array[-1+MUPIPI]  = 29;
  array[-1+EPIK]    = 32;
  array[-1+MUPIK]   = 26;
  array[-1+EKPI]    = 17;
  array[-1+MUKPI]   = 32;
  array[-1+EKK]     = 14;
  array[-1+MUKK]    = 25;
  
  
  array[-1+PIEPI]   = 27;
  array[-1+PIMUPI]  = 7.0;
  array[-1+PIEK]    = 18;
  array[-1+PIMUK]   = 22;
  array[-1+KEK]     = 15;
  array[-1+KMUK]    = 48;
  
  array[-1+PILAM]   = 5.8;
  array[-1+PILAMBAR]= 5.9;
  array[-1+KLAM]    = 15;
  array[-1+KLAMBAR] = 7.2;
}

void fillCLEO(double * array)
{
  // at 10^-6
  array[-1+EGAMMA]  = 2.7;
  array[-1+MUGAMMA] = 1.1;
  
  array[-1+EPI0]    = 3.7; 
  array[-1+MUPI0]   = 4.0;
  
  array[-1+EETA]    = 8.2;
  array[-1+MUETA]   = 9.6;
  
  array[-1+EETAP]   = 0;
  array[-1+MUETAP]  = 0;
  
  array[-1+EKS0]    = 0.91;
  array[-1+MUKS0]   = 0.95;
  
  array[-1+EF0]     = 0;
  array[-1+MUF0]    = 0;
  
  array[-1+ERHO]    = 2.0;
  array[-1+MURHO]   = 6.3;
  array[-1+EKSTAR]  = 5.1;
  array[-1+MUKSTAR] = 7.5;
  array[-1+EAKSTAR] = 7.4;
  array[-1+MUAKSTAR]= 7.5;
  array[-1+EPHI]    = 6.9;
  array[-1+MUPHI]   = 7.0;
  array[-1+EOMEGA]  = 0;
  array[-1+MUOMEGA] = 0;
  
  array[-1+EEE]     = 2.9;
  array[-1+MEE]     = 1.7;
  array[-1+EMM]     = 1.8; 
  array[-1+MMM]     = 1.9; 
  array[-1+EME]     = 1.5; 
  array[-1+MEM]     = 1.5; 
  
  array[-1+EPIPI]   = 2.2;
  array[-1+MUPIPI]  = 8.2;
  array[-1+EPIK]    = 6.4;
  array[-1+MUPIK]   = 7.5;
  array[-1+EKPI]    = 3.8;
  array[-1+MUKPI]   = 7.4;
  array[-1+EKK]     = 6.0;
  array[-1+MUKK]    = 15;
  array[-1+EKS0KS0] = 2.2;
  array[-1+MUKS0KS0]= 3.4;
  array[-1+PIEPI]   = 1.9;
  array[-1+PIMUPI]  = 3.4;
  array[-1+PIEK]    = 2.1;
  array[-1+PIMUK]   = 7.0;
  array[-1+KEK]     = 3.8;
  array[-1+KMUK]    = 6.0;
}

void setLabels(TH1* hist)
{
  TAxis * axis = hist->GetXaxis();
  axis->SetBinLabel(EGAMMA  ,"e^{-} #gamma");
  axis->SetBinLabel(MUGAMMA ,"#mu^{-} #gamma");
  axis->SetBinLabel(EPI0    ,"e^{-} #pi^{0}");
  axis->SetBinLabel(MUPI0   ,"#mu^{-} #pi^{0}");
  axis->SetBinLabel(EETA    ,"e^{-} #eta");
  axis->SetBinLabel(MUETA   ,"#mu^{-} #eta");
  axis->SetBinLabel(EETAP   ,"e^{-} #eta'");
  axis->SetBinLabel(MUETAP  ,"#mu^{-} #eta'");
  axis->SetBinLabel(EKS0    ,"e^{-} K_{S}^{0}");
  axis->SetBinLabel(MUKS0   ,"#mu^{-} K_{S}^{0}");
  axis->SetBinLabel(EF0     ,"e^{-} f_{0}");
  axis->SetBinLabel(MUF0    ,"#mu^{-} f_{0}");
  axis->SetBinLabel(ERHO    ,"e^{-} #rho_{0}");
  axis->SetBinLabel(MURHO   ,"#mu^{-} #rho_{0}");
  axis->SetBinLabel(EKSTAR  ,"e^{-} K*");
  axis->SetBinLabel(MUKSTAR ,"#mu^{-} K*");
  axis->SetBinLabel(EAKSTAR ,"e^{-} #bar{K*}");
  axis->SetBinLabel(MUAKSTAR,"#mu^{-} #bar{K*}");
  axis->SetBinLabel(EPHI    ,"e^{-} #phi");
  axis->SetBinLabel(MUPHI   ,"#mu^{-} #phi");
  axis->SetBinLabel(EOMEGA  ,"e^{-} #omega");
  axis->SetBinLabel(MUOMEGA ,"#mu^{-} #omega");
  axis->SetBinLabel(EEE,     "e^{-} e^{+} e^{-}");
  axis->SetBinLabel(MEE,     "#mu^{-} e^{+} e^{-}");
  axis->SetBinLabel(EMM,     "e^{-} #mu^{+} #mu^{-}");
  axis->SetBinLabel(MMM,     "#mu^{-} #mu^{+} #mu^{-}");
  axis->SetBinLabel(EME,     "e^{-} #mu^{+} e^{-}");
  axis->SetBinLabel(MEM,     "#mu^{-} e^{+} #mu^{-}");
  axis->SetBinLabel(EPIPI   ,"e^{-} #pi^{+} #pi^{-}");
  axis->SetBinLabel(MUPIPI  ,"#mu^{-} #pi^{+} #pi^{-}");
  axis->SetBinLabel(EPIK    ,"e^{-} #pi^{+} K^{-}");
  axis->SetBinLabel(MUPIK   ,"#mu^{-} #pi^{+} K^{-}");
  axis->SetBinLabel(EKPI    ,"e^{-} K^{+} #pi^{-}");
  axis->SetBinLabel(MUKPI   ,"#mu^{-} K^{+} #pi^{-}");
  axis->SetBinLabel(EKK     ,"e^{-} K^{+} K^{-}");
  axis->SetBinLabel(MUKK    ,"#mu^{-} K^{+} K^{-}");
  axis->SetBinLabel(EKS0KS0 ,"e^{-} K_{S}^{0} K_{S}^{0}");
  axis->SetBinLabel(MUKS0KS0,"#mu^{-} K_{S}^{0} K_{S}^{0}");
  axis->SetBinLabel(PIEPI   ,"#pi^{-} e^{+} #pi^{-}");
  axis->SetBinLabel(PIMUPI  ,"#pi^{-} #mu^{+} #pi^{-}");
  axis->SetBinLabel(PIEK    ,"#pi^{-} e^{+} K^{-}");
  axis->SetBinLabel(PIMUK   ,"#pi^{-} #mu^{+} K^{-}");
  axis->SetBinLabel(KEK     ,"K^{-} e^{+} K^{-}");
  axis->SetBinLabel(KMUK    ,"K^{-} #mu^{+} K^{-}");
  axis->SetBinLabel(PILAM   ,"#pi^{-} #Lambda");
  axis->SetBinLabel(PILAMBAR,"#pi^{-} #bar{#Lambda}");
  axis->SetBinLabel(KLAM   , "K^{-} #Lambda");
  axis->SetBinLabel(KLAMBAR, "K^{-} #bar{#Lambda}");
}

void SetUp()
{
  //  cout << "SetUp Plain " << endl;
  gROOT->SetStyle("Plain");
  //  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle();
  gStyle->SetStatStyle(0);
  gStyle->SetOptStat(00000000);
  
  gStyle->SetLabelFont(64,"X");
  gStyle->SetLabelSize(42,"X");
  gStyle->SetLabelOffset(.006,"X");
  
  gStyle->SetTitleFont(64,"X");
  gStyle->SetTitleSize(42,"X");
  gStyle->SetTitleOffset(.5,"X");
  
  gStyle->SetLabelFont(64,"Y");
  gStyle->SetLabelSize(42,"Y");
  gStyle->SetLabelOffset(.001,"Y");
  
  gStyle->SetTitleFont(64,"Y");
  gStyle->SetTitleSize(42,"Y");
  gStyle->SetTitleOffset(.8,"Y");
  
  gStyle->SetTickLength(.03,"X");
  gStyle->SetTickLength(.02,"Y");
  gStyle->SetNdivisions(5,"Y");
  
  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(0);//not for x-axis
  //gStyle->SetTickLength(0,"x");
  //..Get rid of X error bars
  gStyle->SetErrorX(0.001);
}

void HFAGTauLabel(Int_t yearversion, Double_t xpos, Double_t ypos, Double_t scale)
{
  TPaveText *tbox1 = new TPaveText(xpos-0.1,ypos-.12,xpos,ypos,"BRNDC");
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
  
  TPaveText *tbox2 = new TPaveText(xpos-0.1,ypos-0.18,xpos,ypos-0.12,"BRNDC");
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
  if (yearversion==2010001) ts2="Summer 2010";
  if (yearversion==2012001) ts2="Winter 2012";
  tbox2->AddText(ts2);
  tbox2->Draw();
}
