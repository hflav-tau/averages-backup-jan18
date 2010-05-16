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
		      "Gamma132","Gamma150","Gamma152","Gamma999"};
  
  // Unitarity Error = 1.0e-0
  float val_0[nx] = {0.1738591, 0.1780827, 0.1080966, 0.0069768, 0.2548309, 0.0043065, 0.0924652, 0.0006044, 0.0104062, 0.0004386, 0.0011171, 0.0082625, 0.0015714, 0.0034940, 0.0015845, 0.0002420, 0.0011341, 0.0897825, 0.0271242, 0.0010109, 0.0002225, 0.0029298, 0.0007533, 0.0014347, 0.0000602, 0.0008260, 0.0001794, 0.0013890, 0.0001613, 0.0000483, 0.0000937, 0.0198786, 0.0040656, 0.0025439};
  float err_0[nx] = {0.0004073, 0.0004040, 0.0005282, 0.0000982, 0.0009005, 0.0001491, 0.0009884, 0.0002227, 0.0007084, 0.0002161, 0.0003904, 0.0001832, 0.0001595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005183, 0.0007056, 0.0003521, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.0000257, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010179};
  // Unitarity Error = 1.0e-1
  float val_1[nx] = {0.1738591, 0.1780827, 0.1080966, 0.0069768, 0.2548309, 0.0043065, 0.0924652, 0.0006044, 0.0104062, 0.0004386, 0.0011171, 0.0082625, 0.0015714, 0.0034940, 0.0015845, 0.0002420, 0.0011341, 0.0897825, 0.0271242, 0.0010109, 0.0002225, 0.0029298, 0.0007533, 0.0014347, 0.0000602, 0.0008260, 0.0001794, 0.0013890, 0.0001613, 0.0000483, 0.0000937, 0.0198786, 0.0040656, 0.0025439};
  float err_1[nx] = {0.0004073, 0.0004040, 0.0005282, 0.0000982, 0.0009005, 0.0001491, 0.0009884, 0.0002227, 0.0007084, 0.0002161, 0.0003904, 0.0001832, 0.0001595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005183, 0.0007056, 0.0003521, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.0000257, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010179};
  // Unitarity Error = 1.0e-2
  float val_2[nx] = {0.1738591, 0.1780827, 0.1080966, 0.0069768, 0.2548309, 0.0043065, 0.0924652, 0.0006044, 0.0104062, 0.0004386, 0.0011171, 0.0082625, 0.0015714, 0.0034940, 0.0015845, 0.0002420, 0.0011341, 0.0897825, 0.0271242, 0.0010109, 0.0002225, 0.0029298, 0.0007533, 0.0014347, 0.0000602, 0.0008260, 0.0001794, 0.0013890, 0.0001613, 0.0000483, 0.0000937, 0.0198786, 0.0040656, 0.0025439};
  float err_2[nx] = {0.0004072, 0.0004040, 0.0005282, 0.0000982, 0.0009005, 0.0001491, 0.0009884, 0.0002227, 0.0007084, 0.0002161, 0.0003904, 0.0001832, 0.0001595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005183, 0.0007056, 0.0003521, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.0000257, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010179};
  // Unitarity Error = 1.0e-3
  float val_3[nx] = {0.1738604, 0.1780810, 0.1080969, 0.0069768, 0.2548305, 0.0043065, 0.0924651, 0.0006044, 0.0104060, 0.0004386, 0.0011171, 0.0082625, 0.0015714, 0.0034940, 0.0015845, 0.0002420, 0.0011341, 0.0897824, 0.0271242, 0.0010109, 0.0002225, 0.0029298, 0.0007534, 0.0014347, 0.0000602, 0.0008260, 0.0001794, 0.0013890, 0.0001613, 0.0000483, 0.0000937, 0.0198786, 0.0040656, 0.0025432};
  float err_3[nx] = {0.0004037, 0.0003978, 0.0005281, 0.0000982, 0.0009004, 0.0001491, 0.0009883, 0.0002227, 0.0007084, 0.0002161, 0.0003904, 0.0001832, 0.0001595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005183, 0.0007056, 0.0003521, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.0000257, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010175};
  // Unitarity Error = 1.0e-4
  float val_4[nx] = {0.1738355, 0.1781138, 0.1080906, 0.0069757, 0.2548370, 0.0043067, 0.0924673, 0.0006041, 0.0104087, 0.0004384, 0.0011185, 0.0082626, 0.0015714, 0.0034941, 0.0015844, 0.0002420, 0.0011339, 0.0897842, 0.0271255, 0.0010113, 0.0002225, 0.0029298, 0.0007532, 0.0014347, 0.0000602, 0.0008260, 0.0001794, 0.0013889, 0.0001613, 0.0000483, 0.0000936, 0.0198787, 0.0040657, 0.0025565};
  float err_4[nx] = {0.0004666, 0.0005033, 0.0005314, 0.0000988, 0.0009025, 0.0001491, 0.0009886, 0.0002227, 0.0007089, 0.0002162, 0.0003906, 0.0001832, 0.0001595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005186, 0.0007057, 0.0003522, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.0000257, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010253};
  // Unitarity Error = 1.0e-5
  float val_5[nx] = {0.1738393, 0.1781079, 0.1080913, 0.0069759, 0.2548368, 0.0043066, 0.0924671, 0.0006041, 0.0104080, 0.0004384, 0.0011184, 0.0082626, 0.0015714, 0.0034940, 0.0015844, 0.0002420, 0.0011339, 0.0897836, 0.0271264, 0.0010112, 0.0002225, 0.0029298, 0.0007532, 0.0014347, 0.0000602, 0.0008259, 0.0001794, 0.0013891, 0.0001613, 0.0000483, 0.0000937, 0.0198781, 0.0040655, 0.0025550};
  float err_5[nx] = {0.0004576, 0.0004887, 0.0005309, 0.0000987, 0.0009022, 0.0001491, 0.0009885, 0.0002227, 0.0007088, 0.0002162, 0.0003906, 0.0001832, 0.0001595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005185, 0.0007056, 0.0003522, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.0000257, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010241};
  // Unitarity Error = 1.0e-6
  float val_6[nx] = {0.1738449, 0.1781108, 0.1081373, 0.0069764, 0.2548755, 0.0043061, 0.0924471, 0.0006053, 0.0104107, 0.0004379, 0.0011133, 0.0082670, 0.0015693, 0.0034947, 0.0015860, 0.0002424, 0.0011337, 0.0897695, 0.0270360, 0.0010181, 0.0002226, 0.0029303, 0.0007558, 0.0014348, 0.0000603, 0.0008141, 0.0001805, 0.0013942, 0.0001617, 0.0000481, 0.0000937, 0.0199321, 0.0040672, 0.0025229};
  float err_6[nx] = {0.0004575, 0.0004886, 0.0005309, 0.0000987, 0.0009022, 0.0001491, 0.0009885, 0.0002227, 0.0007088, 0.0002162, 0.0003906, 0.0001832, 0.0001595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005185, 0.0007056, 0.0003522, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.0000257, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010241};
  // Unitarity Error = 1.0e-7
  float val_7[nx] = {0.1730281, 0.1809366, 0.1056545, 0.0070420, 0.2602372, 0.0041224, 0.0850751, 0.0007058, 0.0106586, 0.0004895, 0.0016718, 0.0080965, 0.0013976, 0.0034159, 0.0018471, 0.0003028, 0.0009824, 0.0874583, 0.0366541, 0.0015489, 0.0001829, 0.0028193, 0.0003967, 0.0013673, 0.0000614, -0.0004407, 0.0002965, 0.0004799, 0.0001339, 0.0000787, 0.0000853, 0.0140710, 0.0028901, 0.0062525};
  float err_7[nx] = {0.0004575, 0.0004886, 0.0005309, 0.0000987, 0.0009022, 0.0001491, 0.0009885, 0.0002227, 0.0007088, 0.0002162, 0.0003906, 0.0001832, 0.0001595, 0.0001479, 0.0001957, 0.0000520, 0.0002462, 0.0005185, 0.0007056, 0.0003522, 0.0000500, 0.0000691, 0.0001161, 0.0000275, 0.0000180, 0.0000307, 0.0000257, 0.0000725, 0.0000102, 0.0000116, 0.0000149, 0.0006385, 0.0004179, 0.0010241};
  
  float val_1_re[nx];
  float val_2_re[nx];
  float val_3_re[nx];
  float val_4_re[nx];
  float val_5_re[nx];
  float val_6_re[nx];
  float val_7_re[nx];

  float err_1_re[nx];
  float err_2_re[nx];
  float err_3_re[nx];
  float err_4_re[nx];
  float err_5_re[nx];
  float err_6_re[nx];
  float err_7_re[nx];

  for (i=1; i <=nx; i++){
    val_1_re[i-1] = 100.*(val_1[i-1] - val_4[i-1] ) / val_4[i-1];
    err_1_re[i-1] = 100.*(err_1[i-1] - err_4[i-1] ) / err_4[i-1];
    val_2_re[i-1] = 100.*(val_2[i-1] - val_4[i-1] ) / val_4[i-1];
    err_2_re[i-1] = 100.*(err_2[i-1] - err_4[i-1] ) / err_4[i-1];
    val_3_re[i-1] = 100.*(val_3[i-1] - val_4[i-1] ) / val_4[i-1];
    err_3_re[i-1] = 100.*(err_3[i-1] - err_4[i-1] ) / err_4[i-1];
    val_4_re[i-1] = 100.*(val_4[i-1] - val_4[i-1] ) / val_4[i-1];
    err_4_re[i-1] = 100.*(err_4[i-1] - err_4[i-1] ) / err_4[i-1];
    val_5_re[i-1] = 100.*(val_5[i-1] - val_4[i-1] ) / val_4[i-1];
    err_5_re[i-1] = 100.*(err_5[i-1] - err_4[i-1] ) / err_4[i-1];
    val_6_re[i-1] = 100.*(val_6[i-1] - val_4[i-1] ) / val_4[i-1];
    err_6_re[i-1] = 100.*(err_6[i-1] - err_4[i-1] ) / err_4[i-1];
    val_7_re[i-1] = 100.*(val_7[i-1] - val_4[i-1] ) / val_4[i-1];
    err_7_re[i-1] = 100.*(err_7[i-1] - err_4[i-1] ) / err_4[i-1];

  }
  
  int todraw=6;
  int istart=1;
  double BARWIDTH=.9/(todraw*1.);
  double OFFSET[7]={0,0,0,0,0,0,0};
  double thisoffset=0;
  for (i=0; i<7; i++){
    if (i > istart-1) {
      thisoffset+=BARWIDTH;
      OFFSET[i]=thisoffset;
    }
    cout << "i = " << i << " OFFSET = " << OFFSET[i] << endl;
  }

  //                 1         2      3        4     5       6      7
  int USECOLOR[7] = {kMagenta, kCyan, kOrange, kBlack, kRed, kTeal, kYellow};

  TH1F *h1b = new TH1F("h1b",Form("Difference of fitted Values by varying sum error from 1e-%d to 1e-%d Normalized to 1e-4 [expressed in %]",istart,todraw),nx,0,nx);
  h1b->SetFillColor(USECOLOR[0]);
  h1b->SetBarWidth(BARWIDTH);
  h1b->SetBarOffset(OFFSET[0]);
  h1b->SetStats(0);
  h1b->SetMinimum(-10);
  h1b->SetMaximum( 10);
  for (i=1; i<=nx; i++) {
    h1b->Fill(os_X[i-1], val_1_re[i-1]);
    h1b->GetXaxis()->SetBinLabel(i,os_X[i-1]);
  }
  
  TH1F *h2b = new TH1F("h2b","h2b",nx,0,nx);
  h2b->SetFillColor(USECOLOR[1]);
  h2b->SetBarWidth(BARWIDTH);
  h2b->SetBarOffset(OFFSET[1]);
  h2b->SetStats(0);
  for (i=1;i<=nx;i++) h2b->Fill(os_X[i-1], val_2_re[i-1]);
  
  TH1F *h3b = new TH1F("h3b","h3b",nx,0,nx);
  h3b->SetFillColor(USECOLOR[2]);
  h3b->SetBarWidth(BARWIDTH);
  h3b->SetBarOffset(OFFSET[2]);
  h3b->SetStats(0);
  for (i=1;i<=nx;i++) h3b->Fill(os_X[i-1], val_3_re[i-1]);
  
  TH1F *h4b = new TH1F("h4b","h4b",nx,0,nx);
  h4b->SetFillColor(USECOLOR[3]);
  h4b->SetBarWidth(BARWIDTH);
  h4b->SetBarOffset(OFFSET[3]);
  h4b->SetStats(0);
  for (i=1;i<=nx;i++) h4b->Fill(os_X[i-1], val_4_re[i-1]);
  
  TH1F *h5b = new TH1F("h5b","h5b",nx,0,nx);
  h5b->SetFillColor(USECOLOR[4]);
  h5b->SetBarWidth(BARWIDTH);
  h5b->SetBarOffset(OFFSET[4]);
  h5b->SetStats(0);
  for (i=1;i<=nx;i++) h5b->Fill(os_X[i-1], val_5_re[i-1]);
  
  TH1F *h6b = new TH1F("h6b","h6b",nx,0,nx);
  h6b->SetFillColor(USECOLOR[5]);
  h6b->SetBarWidth(BARWIDTH);
  h6b->SetBarOffset(OFFSET[5]);
  h6b->SetStats(0);
  for (i=1;i<=nx;i++) h6b->Fill(os_X[i-1], val_6_re[i-1]);
  
  TH1F *h7b = new TH1F("h7b","h7b",nx,0,nx);
  h7b->SetFillColor(USECOLOR[6]);
  h7b->SetBarWidth(BARWIDTH);
  h7b->SetBarOffset(OFFSET[6]);
  h7b->SetStats(0);
  for (i=1;i<=nx;i++) h7b->Fill(os_X[i-1], val_7_re[i-1]);
  
  TH1F *h1e = new TH1F("h1e",Form("Difference of fitted Errors by varying sum error from 1e-%d to 1e-%d Normalized to 1e-4 [expressed in %]",istart,todraw),nx,0,nx);
  h1e->SetFillColor(USECOLOR[0]);
  h1e->SetBarWidth(BARWIDTH);
  h1e->SetBarOffset(OFFSET[0]);
  h1e->SetStats(0);
  h1e->SetMinimum(-25);
  h1e->SetMaximum( 25);
  for (i=1; i<=nx; i++) {
    h1e->Fill(os_X[i-1], err_1_re[i-1]);
    h1e->GetXaxis()->SetBinLabel(i,os_X[i-1]);
  }
  
  TH1F *h2e = new TH1F("h2e","h2e",nx,0,nx);
  h2e->SetFillColor(USECOLOR[1]);
  h2e->SetBarWidth(BARWIDTH);
  h2e->SetBarOffset(OFFSET[1]);
  h2e->SetStats(0);
  for (i=1;i<=nx;i++) h2e->Fill(os_X[i-1], err_2_re[i-1]);
  
  TH1F *h3e = new TH1F("h3e","h3e",nx,0,nx);
  h3e->SetFillColor(USECOLOR[2]);
  h3e->SetBarWidth(BARWIDTH);
  h3e->SetBarOffset(OFFSET[2]);
  h3e->SetStats(0);
  for (i=1;i<=nx;i++) h3e->Fill(os_X[i-1], err_3_re[i-1]);
  
  TH1F *h4e = new TH1F("h4e","h4e",nx,0,nx);
  h4e->SetFillColor(USECOLOR[3]);
  h4e->SetBarWidth(BARWIDTH);
  h4e->SetBarOffset(OFFSET[3]);
  h4e->SetStats(0);
  for (i=1;i<=nx;i++) h4e->Fill(os_X[i-1], err_4_re[i-1]);
  
  TH1F *h5e = new TH1F("h5e","h5e",nx,0,nx);
  h5e->SetFillColor(USECOLOR[4]);
  h5e->SetBarWidth(BARWIDTH);
  h5e->SetBarOffset(OFFSET[4]);
  h5e->SetStats(0);
  for (i=1;i<=nx;i++) h5e->Fill(os_X[i-1], err_5_re[i-1]);
  
  TH1F *h6e = new TH1F("h6e","h6e",nx,0,nx);
  h6e->SetFillColor(USECOLOR[5]);
  h6e->SetBarWidth(BARWIDTH);
  h6e->SetBarOffset(OFFSET[5]);
  h6e->SetStats(0);
  for (i=1;i<=nx;i++) h6e->Fill(os_X[i-1], err_6_re[i-1]);
  
  TH1F *h7e = new TH1F("h7e","h7e",nx,0,nx);
  h7e->SetFillColor(USECOLOR[6]);
  h7e->SetBarWidth(BARWIDTH);
  h7e->SetBarOffset(OFFSET[6]);
  h7e->SetStats(0);
  for (i=1;i<=nx;i++) h7e->Fill(os_X[i-1], err_7_re[i-1]);
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->SetGrid();
  //  c1->SetLogx(1);
  h1b->Draw("hbar");
  h2b->Draw("hbar,same");
  h3b->Draw("hbar,same");
  h4b->Draw("hbar,same");
  h5b->Draw("hbar,same");
  h6b->Draw("hbar,same");
  //  h7b->Draw("hbar,same");
  c1->SaveAs("compare-val_unitarity_error.eps");

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  c2->SetGrid();
  //  c2->SetLogx(1);
  h1e->Draw("hbar");
  h2e->Draw("hbar,same");
  h3e->Draw("hbar,same");
  h4e->Draw("hbar,same");
  h5e->Draw("hbar,same");
  h6e->Draw("hbar,same");
  //  h7e->Draw("hbar,same");
  c2->SaveAs("compare-err_unitarity_error.eps");

  return 0;
}
