{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
  gStyle->SetStatStyle(0);
  gStyle->SetHistMinimumZero();

  int i;
  const Int_t nx = 34;
  char *os_X[nx]   = {"Gamma3","Gamma5","Gamma9","Gamma10","Gamma14",
		      "Gamma16","Gamma20","Gamma23","Gamma27","Gamma28",
		      "Gamma30","Gamma35","Gamma37","Gamma40","Gamma42",
		      "Gamma47","Gamma48","Gamma62","Gamma70","Gamma77",
		      "Gamma78","Gamma85","Gamma89","Gamma93","Gamma94",
		      "Gamma103","Gamma104","Gamma126","Gamma128","Gamma130",
		      "Gamma132","Gamma150","Gamma152", "Gamma999"};
  
  // COMBOS
  // Unitarity Error = 1.0e-5
  float val_0[nx] = {0.1738393, 0.1781079, 0.1080913, 0.0069759, 0.2548368, 0.0043066, 0.0924671, 0.0006041, 0.0104080, 0.0004384, 0.0011184, 0.0082626, 0.00
15714, 0.0034940, 0.0015844, 0.0002420, 0.0011339, 0.0897836, 0.0271264, 0.0010112, 0.0002225, 0.0029298, 0.0007532, 0.0014347, 0.0000602, 0.0008259, 0.00017
94, 0.0013891, 0.0001613, 0.0000483, 0.0000937, 0.0198781, 0.0040655, 0.0025550};
  float err_0[nx] = {0.0004576, 0.0004887, 0.0005309, 0.0000987, 0.0009022, 0.0001491, 0.0009885, 0.0002227, 0.0007088, 0.0002162, 0.0003906, 0.0001832, 0.00
01595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005185, 0.0007056, 0.0003522, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.00002
57, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010241};

  TH1F *h1b = new TH1F("h1b","FITTED VALUES: COMBOS (BLUE) vs ALUCOMB (RED)",nx,0,nx);
  h1b->SetFillColor(4);
  h1b->SetBarWidth(0.45);
  h1b->SetBarOffset(0.0);
  h1b->SetStats(0);
  h1b->SetMinimum(3e-5);
  h1b->SetMaximum(3e-1);
  for (i=1; i<=nx; i++) {
    h1b->Fill(os_X[i-1], val_0[i-1]);
    h1b->GetXaxis()->SetBinLabel(i,os_X[i-1]);
  }
  TH1F *h1e = new TH1F("h1e","FITTED ERRORS: COMBOS (BLUE) vs ALUCOMB (RED)",nx,0,nx);
  h1e->SetFillColor(4);
  h1e->SetBarWidth(0.45);
  h1e->SetBarOffset(0.0);
  h1e->SetStats(0);
  h1e->SetMinimum(0.9e-5);
  h1e->SetMaximum(1.1e-3);
  for (i=1; i<=nx; i++) {
    h1e->Fill(os_X[i-1], err_0[i-1]);
    h1e->GetXaxis()->SetBinLabel(i,os_X[i-1]);
  }

  // ALUCOMB
  //                  0             1             2             3             4             5             6             7             8             9             
  float val_1[nx] = { 0.1738590845, 0.1782541809, 0.1080968657, 0.0069767951, 0.2548308682, 0.0043065871, 0.0924654495, 0.0006043352, 0.0104063050, 0.0004385660, 
		      //10          11            12            13            14            15            16            17            18            19          
                      0.0011171835, 0.0082625667, 0.0015714111, 0.0034940424, 0.0015844691, 2.420462e-04, 0.0011340548, 8.978043e-02, 0.0271247054, 0.0010109089, 
		      //20          21            22            23            24            25            26            27            28            29
		      2.225196e-04, 2.930080e-03, 0.0007533600, 1.434775e-03, 6.022149e-05, 8.260009e-04, 1.793745e-04, 0.0013889683, 1.613101e-04, 4.832282e-05, 
		      //30          31            32            33           
		      9.373595e-05, 0.0198787276, 0.0040656142, 0.002543951};//,  0.2591374553, 0.0900845710}

  float err_1[nx] = { 0.0004072637, 0.0004808258, 0.0005282905, 0.0000982477, 0.0009005460, 0.0001490496, 0.0009883572, 0.0002227264, 0.0007084371, 0.0002161436, 
		      0.0003903924, 0.0001831599, 0.0001595133, 0.0001479371, 0.0001956770, 5.198715e-05, 0.0002462339, 5.183849e-04, 0.0007055538, 0.0003521454, 
		      4.999808e-05, 6.913499e-05, 0.0001160973, 2.749040e-05, 1.799956e-05, 3.069295e-05, 2.568275e-05, 0.0000725263, 1.015869e-05, 1.160769e-05,
		      1.491316e-05, 0.0006385443, 0.0004178982, 0.001017934 };//0.0009014570 0.0005176572


  
  TH1F *h2b = new TH1F("h2b","h2b",nx,0,nx);
  h2b->SetFillColor(2);
  h2b->SetBarWidth(0.45);
  h2b->SetBarOffset(0.5);
  h2b->SetStats(0);
  for (i=1;i<=nx;i++) h2b->Fill(os_X[i-1], val_1[i-1]);
  
  TH1F *h2e = new TH1F("h2e","h2e",nx,0,nx);
  h2e->SetFillColor(2);
  h2e->SetBarWidth(0.45);
  h2e->SetBarOffset(0.5);
  h2e->SetStats(0);
  for (i=1;i<=nx;i++) h2e->Fill(os_X[i-1], err_1[i-1]);

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetFillColor(30);

  c1->SetLogx();
  c1->SetGrid();
  
  h1b->Draw("hbar0");
  h2b->Draw("hbar0,same");

  c1->Update();
  c1->SaveAs("compare-val_combos_alucomb.eps");
  
  TCanvas *c1 = new TCanvas("c2","c2",600,600);
  c2->SetFillColor(30);
  c2->SetLogx();
  c1->SetGrid();
  
  h1e->Draw("hbar0");
  h2e->Draw("hbar0,same");

  c2->Update();
  c2->SaveAs("compare-err_combos_alucomb.eps");
  
}
