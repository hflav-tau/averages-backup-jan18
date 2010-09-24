#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "TROOT.h"
#include "TString.h"
#include "TMath.h"

int aleph_hcorr(){
  int i;
  //                               0       1       2        3        4        5        6        7        8         9       10      11       12
  //                               e      mu       h     hpi0    h2pi0    h3pi0    h4pi0       3h    3hpi0    3h2pi0    3h3pi0     5h    5hpi0
  TString channels[13] =  {       "e",   "mu",    "h",  "hpi0", "h2pi0", "h3pi0", "h4pi0",    "3h", "3hpi0", "3h2pi0", "3h3pi0",  "5h", "5hpi0"};
  int    measnumber[13] = {       72,     66,     115,   120,       29,      30,       5,      51,      55,       33,       0,     94,     100};

  //  Table 14
  //                               e      mu      pi    pipi0   pi2pi0   pi3pi0   pi4pi0      3pi   3pipi0   3pi2pi0   3pi3pi0     5pi  5pipi0
  double measval_pi[13] =   {  17.837, 17.319, 10.828,  25.471,  9.239,   0.977,   0.112,   9.041,   4.590,    0.392,    0.013,  0.072,  0.014};
  //                               e      mu      pi    pipi0   pi2pi0   pi3pi0   pi4pi0      3pi   3pipi0   3pi2pi0   3pi3pi0     5pi  5pipi0
  double staterr_pi[13] =   {   0.072,  0.070,  0.070,  0.097,   0.086,   0.069,   0.037,   0.060,   0.057,    0.030,     0.00,  0.009,  0.007};
  //                               e      mu      pi    pipi0   pi2pi0   pi3pi0   pi4pi0      3pi   3pipi0   3pi2pi0   3pi3pi0     5pi  5pipi0
  double systerr_pi[13] =   {   0.036,  0.032,  0.078,  0.085,   0.090,   0.058,   0.035,   0.076,   0.064,    0.035,    0.010,  0.012,  0.006};

  for (i=0;i<13;++i) {
    measval_pi[i]*=1e-2;
    staterr_pi[i]*=1e-2;
    systerr_pi[i]*=1e-2;
  }

  double measval_h[13];
  double staterr_h[13];
  double systerr_h[13];
  
  // e [Fig 40 : same as Table 14]
  measval_h[0] = measval_pi[0];
  staterr_h[0] = staterr_pi[0];
  systerr_h[0] = systerr_pi[0];

  // mu [Fig 41 : same as Table 14]
  measval_h[1] = measval_pi[1];
  staterr_h[1] = staterr_pi[1];
  systerr_h[1] = systerr_pi[1];

  // h [Fig 42]
  measval_h[2] = 11.524e-2;
  staterr_h[2] = 0.070e-2;
  systerr_h[2] = 0.078e-2;

  // hpi0 [Fig 43]
  measval_h[3] = 25.924e-2;
  staterr_h[3] = 0.097e-2;
  systerr_h[3] = 0.085e-2;

  // h2pi0 [Fig 44]
  measval_h[4] = 9.295e-2;
  staterr_h[4] = 0.084e-2;
  systerr_h[4] = 0.088e-2;

  // h4pi0 [same as Table 14]
  measval_h[6] = measval_pi[6];
  staterr_h[6] = staterr_pi[6];
  systerr_h[6] = systerr_pi[6];

  // h3,4pi0 [Fig 45] = (1.194 +- .080 +- 0.069) % => sqrt ( 0.08**2 - 0.037**2 ) = 0.0709295, sqrt ( 0.069**2 - 0.035**2 ) = 0.0594643
  double h34pi0_val = 1.194e-2;
  double h34pi0_ste = 0.080e-2;
  double h34pi0_sye = 0.069e-2;

  // h3pi0 [from above]
  measval_h[5] = h34pi0_val - measval_h[6];
  staterr_h[5] = TMath::Sqrt(TMath::Power(h34pi0_ste,2) - TMath::Power(staterr_h[6],2));
  systerr_h[5] = TMath::Sqrt(TMath::Power(h34pi0_sye,2) - TMath::Power(systerr_h[6],2));

  // 3h [Fig 46]
  measval_h[7] = 9.469e-2;
  staterr_h[7] = 0.062e-2;
  systerr_h[7] = 0.073e-2;

  // 3hpi0 [Fig 47]
  measval_h[8] = 4.726e-2;
  staterr_h[8] = 0.059e-2;
  systerr_h[8] = 0.049e-2;

  // 3h2pi0 [same as Table 14]
  measval_h[9] = measval_pi[9];
  staterr_h[9] = staterr_pi[9];
  systerr_h[9] = systerr_pi[9];

  // 3h3pi0 [same as Table 14]
  measval_h[10] = measval_pi[10];
  staterr_h[10] = staterr_pi[10];
  systerr_h[10] = systerr_pi[10];

  // 5h [Fig 49 : same as Table 14]
  measval_h[11] = measval_pi[11];
  staterr_h[11] = staterr_pi[11];
  systerr_h[11] = systerr_pi[11];

  // 5hpi0 [Fig 50 : same as Table 14]
  measval_h[12] = measval_pi[12];
  staterr_h[12] = staterr_pi[12];
  systerr_h[12] = systerr_pi[12];// NOTE: PDG thinks this number is .009e-3, instead of .006e-3

  // PDG Notes:
  measval_h[8]+= 0.008e-2; // 3pipi0 
  measval_h[9]+= 0.043e-2; // 3pi2pi0
  measval_h[12]+= 0.7e-4;  // 5pipi0

  //Table 15
  //Correlation matrix of the statistical errors on the branchingfractions
  //                               e    mu      h     hpi0   h2pi0   h3pi0   h4pi0      3h   3hpi0  3h2pi0  3h3pi0      5h   5hpi0
  double statcorr[13][13] = { { 1.00, -0.21, -0.15,  -0.25,  -0.09,  -0.01,   0.00,  -0.15,  -0.10,   0.03,  -0.06,   0.00,   0.01},    // e     
			      { 0.00,  1.00, -0.13,  -0.21,  -0.07,  -0.06,   0.00,  -0.09,  -0.07,   0.00,  -0.02,   0.00,  -0.04},   	// mu	
			      { 0.00,  0.00,  1.00,  -0.31,  -0.02,   0.01,  -0.06,  -0.12,  -0.06,  -0.02,   0.01,  -0.01,   0.02},   	// h	
			      { 0.00,  0.00,  0.00,   1.00,  -0.40,   0.05,   0.00,  -0.11,  -0.06,  -0.02,   0.00,  -0.04,  -0.04},   	// hpi0	
			      { 0.00,  0.00,  0.00,   0.00,   1.00,  -0.51,   0.26,  -0.09,   0.01,  -0.07,   0.06,  -0.01,   0.03},   	// h2pi0	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   1.00,  -0.75,   0.01,  -0.03,   0.05,  -0.02,  -0.01,   0.01},   	// h3pi0	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   1.00,  -0.02,  -0.02,  -0.03,   0.01,   0.02,  -0.03},   	// h4pi0	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.33,   0.08,  -0.05,  -0.04,   0.00},   	// 3h	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.45,   0.19,  -0.02,  -0.02},   	// 3h1pi0
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.65,   0.03,   0.02},   	// 3h2pi0
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.01,  -0.04},   	// 3h3pi0
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.24},   	// 5h	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00} }; 	// 5hpi0 
  
  //Table 16
  //Correlation matrix of the systematic errors on the branchingfractions
  //                               e    mu      h     hpi0   h2pi0   h3pi0   h4pi0      3h   3hpi0  3h2pi0  3h3pi0      5h    5hpi0
  double systcorr[13][13] = { { 1.00, -0.17, -0.01,   0.02,   0.01,   0.03,  -0.08,  -0.17,  -0.22,  -0.05,   0.02,   0.00,   0.00},     // e     
			      { 0.00,  1.00,  0.05,   0.09,  -0.03,   0.02,  -0.13,  -0.11,  -0.24,  -0.06,   0.01,   0.03,  -0.04},	 // mu	
			      { 0.00,  0.00,  1.00,   0.36,  -0.29,  -0.32,  -0.42,   0.34,  -0.40,  -0.40,  -0.07,   0.16,  -0.09},	 // h	
			      { 0.00,  0.00,  0.00,   1.00,  -0.35,  -0.02,  -0.33,   0.01,  -0.54,  -0.26,   0.02,   0.11,  -0.06},	 // hpi0	
			      { 0.00,  0.00,  0.00,   0.00,   1.00,  -0.01,   0.13,  -0.24,   0.07,   0.13,   0.06,  -0.13,   0.03},	 // h2pi0	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   1.00,  -0.13,  -0.29,  -0.02,   0.15,   0.09,  -0.06,   0.04},	 // h3pi0	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   1.00,  -0.14,   0.34,   0.27,   0.00,  -0.12,  -0.05},	 // h4pi0	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.03,  -0.16,  -0.11,   0.17,  -0.06},	 // 3h	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00,   0.04,  -0.03,  -0.07,  -0.01},	 // 3h1pi0
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.14,  -0.09,   0.07},	 // 3h2pi0
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.02,   0.02},	 // 3h3pi0
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00,  -0.26},	 // 5h	
			      { 0.00,  0.00,  0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   1.00} };	 // 5hpi0 

  double toterr_h[13];
  for (i=0;i<13;++i) toterr_h[i] = sqrt(staterr_h[i]*staterr_h[i] + systerr_h[i]*systerr_h[i]);

  for (i=1;i<150;++i) {
    for (int ii=0;ii<13;++ii) {
      if (measnumber[ii]==i) {
	cout << "*" << "  ";
	cout << Form("%4d",measnumber[ii]) << "  " ;
	cout << Form("%10s",channels[ii].Data()) << "  " ;
	cout << Form("%10.5e",measval_h[ii]) << "  " ;
	//	cout << Form("%10.5e",staterr_h[ii]) << "  " ;
	//	cout << Form("%10.5e",systerr_h[ii]) << "  " ;
	cout << Form("%10.5e",toterr_h[ii])  << endl;
      }
    }
  }

  int j;
  
  for (i=0;i<13;++i) {
    for (j=0;j<13;++j) {
      if (i>j) {
	statcorr[i][j] = statcorr[j][i];
	systcorr[i][j] = systcorr[j][i];
      }
    }
  }
  for (i=1;i<150;++i) {
    for (j=1;j<150;++j) {
      for (int ii=0;ii<13;++ii) {
	if (measnumber[ii]==i) {
	  for (int jj=0;jj<13;++jj) {
	    if (measnumber[jj]==j) {
	      if (measnumber[ii]<measnumber[jj]) {
		cout << "%" << "  ";
		cout << Form("%4d",measnumber[ii]) << " " ;//<< channels[i].Data() << " " ;
		cout << Form("%4d",measnumber[jj]) << " " ;//<< channels[j].Data() << " " ;
		cout << Form("%7.3f",100.*(staterr_h[ii]*statcorr[ii][jj]*staterr_h[jj] + systerr_h[ii]*systcorr[ii][jj]*systerr_h[jj])/(toterr_h[ii]*toterr_h[jj]));
		cout << endl;
	      }
	    }
	  }
	}
      }
    }
  }
  //
  return 0;
}
