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
#include "Riostream.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"

using namespace std;
// ----------------------------------------------------------------------
int main() {
  gSystem->Load("libMatrix");
  //
  int ipar,iimeas,jjmeas;
  // from s035-fit-no-babar-belle.fit
  const int nbase=31;//33;
  vector<int> basegamma;
  vector<int> baseparm;
  double    baseseed[nbase],basefitvalue[nbase],basefiterror[nbase],baserescalederror[nbase],basescalefactor[nbase];
  string    basetitle[nbase];
  //
  int ibase = 0;
  // INPUT PARAMETERS
  // GAMMA PARMETER SEED TITLE //BASE = NQUAN (in COMBOS)
  // RESULTS FOR PARAMETERS
  // ----- -------- ---- ----- 
  basegamma.push_back(  3); baseparm.push_back(  1) ; baseseed[ibase] = 1.740000E+01 ; basetitle[ibase] = "tau- --> mu- nubar(mu) nu(tau)"; basefitvalue[ibase]=0.173556 ; basefiterror[ibase]=0.000463 ; baserescalederror[ibase]=0.000465 ; basescalefactor[ibase]=1.01 ;   ++ibase; // 1
  basegamma.push_back(  5); baseparm.push_back(  2) ; baseseed[ibase] = 1.780000E+01 ; basetitle[ibase] = "tau- --> e- nubar(e) nu(tau)";   basefitvalue[ibase]=0.178413 ; basefiterror[ibase]=0.000479 ; baserescalederror[ibase]=0.000483 ; basescalefactor[ibase]=1.01 ;   ++ibase; // 2
  basegamma.push_back(  9); baseparm.push_back( 12) ; baseseed[ibase] = 1.140000E+01 ; basetitle[ibase] = "tau- --> pi- nu(tau)";           basefitvalue[ibase]=0.109002 ; basefiterror[ibase]=0.000617 ; baserescalederror[ibase]=0.000660 ; basescalefactor[ibase]=1.07 ;   ++ibase; // 3
  basegamma.push_back( 10); baseparm.push_back(  7) ; baseseed[ibase] = 7.000000E-01 ; basetitle[ibase] = "tau- --> K- nu(tau)" ;           basefitvalue[ibase]=0.006910 ; basefiterror[ibase]=0.000219 ; baserescalederror[ibase]=0.000228 ; basescalefactor[ibase]=1.04 ;   ++ibase; // 4
  basegamma.push_back( 14); baseparm.push_back( 16) ; baseseed[ibase] = 2.500000E+01 ; basetitle[ibase] = "tau- --> pi- pi0 nu(tau)";       basefitvalue[ibase]=0.254980 ; basefiterror[ibase]=0.000924 ; baserescalederror[ibase]=0.001023 ; basescalefactor[ibase]=1.11 ;   ++ibase; // 5
  basegamma.push_back( 16); baseparm.push_back(182) ; baseseed[ibase] = 4.500000E-01 ; basetitle[ibase] = "tau- --> K- pi0 nu(tau)";        basefitvalue[ibase]=0.004525 ; basefiterror[ibase]=0.000265 ; baserescalederror[ibase]=0.000267 ; basescalefactor[ibase]=1.01 ;   ++ibase; // 6
  basegamma.push_back( 20); baseparm.push_back(201) ; baseseed[ibase] = 9.000000E+00 ; basetitle[ibase] = "tau- --> pi- 2pi0 nu(tau) (ex.K0)";basefitvalue[ibase]=0.092494 ; basefiterror[ibase]=0.000975 ; baserescalederror[ibase]=0.001240 ; basescalefactor[ibase]=1.27 ; ++ibase; // 7
  basegamma.push_back( 23); baseparm.push_back(115) ; baseseed[ibase] = 6.000000E-02 ; basetitle[ibase] = "tau- --> K- 2pi0 nu(tau) (ex.K0)"; basefitvalue[ibase]=0.000581 ; basefiterror[ibase]=0.000225 ; baserescalederror[ibase]=0.000232 ; basescalefactor[ibase]=1.03 ; ++ibase; // 8
  basegamma.push_back( 27); baseparm.push_back(203) ; baseseed[ibase] = 1.000000E+00 ; basetitle[ibase] = "tau- --> pi- 3pi0 nu(tau) (ex.K0)";basefitvalue[ibase]=0.010403 ; basefiterror[ibase]=0.000709 ; baserescalederror[ibase]=0.000760 ; basescalefactor[ibase]=1.07 ; ++ibase; // 9
  basegamma.push_back( 28); baseparm.push_back(116) ; baseseed[ibase] = 4.000000E-01 ; basetitle[ibase] = "tau- --> K- 3pi0 nu(tau) (ex.K0, eta)";basefitvalue[ibase]=0.000417 ; basefiterror[ibase]=0.000219 ; baserescalederror[ibase]=0.000220 ; basescalefactor[ibase]=1.01 ; ++ibase; // 10
  basegamma.push_back( 30); baseparm.push_back(110) ; baseseed[ibase] = 1.000000E-01 ; basetitle[ibase] = "tau- --> h- 4pi0 nu(tau) (ex.K0,eta)"; basefitvalue[ibase]=0.001024 ; basefiterror[ibase]=0.000392 ; baserescalederror[ibase]=0.000398 ; basescalefactor[ibase]=1.01 ; ++ibase; // 11
  basegamma.push_back( 35); baseparm.push_back(117) ; baseseed[ibase] = 9.600000E-01 ; basetitle[ibase] = "tau- --> pi- Kbar0 nu(tau)";     basefitvalue[ibase]=0.008963 ; basefiterror[ibase]=0.000367 ; baserescalederror[ibase]=0.000409 ; basescalefactor[ibase]=1.11 ;   ++ibase; // 12
  basegamma.push_back( 37); baseparm.push_back( 62) ; baseseed[ibase] = 1.600000E-01 ; basetitle[ibase] = "tau- --> K- K0 nu(tau)";         basefitvalue[ibase]=0.001530 ; basefiterror[ibase]=0.000161 ; baserescalederror[ibase]=0.000163 ; basescalefactor[ibase]=1.01 ;   ++ibase; // 13
  basegamma.push_back( 40); baseparm.push_back(118) ; baseseed[ibase] = 4.000000E-01 ; basetitle[ibase] = "tau- --> pi- Kbar0 pi0 nu(tau)"; basefitvalue[ibase]=0.003778 ; basefiterror[ibase]=0.000368 ; baserescalederror[ibase]=0.000374 ; basescalefactor[ibase]=1.02 ;   ++ibase; // 14
  basegamma.push_back( 42); baseparm.push_back(119) ; baseseed[ibase] = 2.000000E-01 ; basetitle[ibase] = "tau- --> K- K0 pi0 nu(tau)";     basefitvalue[ibase]=0.001541 ; basefiterror[ibase]=0.000200 ; baserescalederror[ibase]=0.000200 ; basescalefactor[ibase]=1.00 ;   ++ibase; // 15
  basegamma.push_back( 47); baseparm.push_back(214) ; baseseed[ibase] = 1.000000E-01 ; basetitle[ibase] = "tau- --> pi- K(S)0 K(S)0 nu(tau)";basefitvalue[ibase]=0.000241 ; basefiterror[ibase]=0.000052 ; baserescalederror[ibase]=0.000052 ; basescalefactor[ibase]=1.00 ;  ++ibase; // 16
  basegamma.push_back( 48); baseparm.push_back(240) ; baseseed[ibase] = 1.000000E-01 ; basetitle[ibase] = "tau- --> pi- K(S)0 K(L)0 nu(tau)";basefitvalue[ibase]=0.001121 ; basefiterror[ibase]=0.000245 ; baserescalederror[ibase]=0.000301 ; basescalefactor[ibase]=1.23 ;  ++ibase; // 17
  basegamma.push_back( 62); baseparm.push_back(259) ; baseseed[ibase] = 9.100000E+00 ; basetitle[ibase] = "tau- --> pi- pi+ pi- nu(tau) (ex.K0,omega)";basefitvalue[ibase]=0.089866 ; basefiterror[ibase]=0.000600 ; baserescalederror[ibase]=0.000757 ; basescalefactor[ibase]=1.26 ; ++ibase; // 18
  basegamma.push_back( 70); baseparm.push_back(263) ; baseseed[ibase] = 2.500000E+00 ; basetitle[ibase] = "tau- --> pi- pi+ pi- pi0 nu(tau) (ex.K0,omega)"; basefitvalue[ibase]=0.026914 ; basefiterror[ibase]=0.000707 ; baserescalederror[ibase]=0.000826 ; basescalefactor[ibase]=1.17 ; ++ibase; // 19
  basegamma.push_back( 77); baseparm.push_back(216) ; baseseed[ibase] = 1.000000E-01 ; basetitle[ibase] = "tau- --> h- h- h+ 2pi0 nu(tau) (ex.K0,omega,eta)"; basefitvalue[ibase]=0.000906 ; basefiterror[ibase]=0.000357 ; baserescalederror[ibase]=0.000368 ; basescalefactor[ibase]=1.03 ;++ibase; // 20
  basegamma.push_back( 78); baseparm.push_back(204) ; baseseed[ibase] = 1.300000E-01 ; basetitle[ibase] = "tau- --> h- h- h+ 3pi0 nu(tau)"; basefitvalue[ibase]=0.000223 ; basefiterror[ibase]=0.000050 ; baserescalederror[ibase]=0.000050 ; basescalefactor[ibase]=1.00 ;   ++ibase; // 21
  basegamma.push_back( 85); baseparm.push_back(260) ; baseseed[ibase] = 3.000000E-01 ; basetitle[ibase] = "tau- --> K- pi+ pi- nu(tau) (ex.K0)";basefitvalue[ibase]=0.003335 ; basefiterror[ibase]=0.000223 ; baserescalederror[ibase]=0.000352 ; basescalefactor[ibase]=1.58 ;  ++ibase; // 22
  basegamma.push_back( 89); baseparm.push_back(285) ; baseseed[ibase] = 6.000000E-02 ; basetitle[ibase] = "tau- --> K- pi+ pi- pi0 nu(tau) (ex.K0,eta)"; basefitvalue[ibase]=0.000730 ; basefiterror[ibase]=0.000117 ; baserescalederror[ibase]=0.000122 ; basescalefactor[ibase]=1.04 ; ++ibase; // 23
  basegamma.push_back( 93); baseparm.push_back(  5) ; baseseed[ibase] = 1.600000E-01 ; basetitle[ibase] = "tau- --> K- K+ pi- nu(tau)";     basefitvalue[ibase]=0.001531 ; basefiterror[ibase]=0.000070 ; baserescalederror[ibase]=0.000097 ; basescalefactor[ibase]=1.38 ;  ++ibase; // 24
  basegamma.push_back( 94); baseparm.push_back(247) ; baseseed[ibase] = 4.000000E-02 ; basetitle[ibase] = "tau- --> K- K+ pi- pi0 nu(tau)"; basefitvalue[ibase]=0.000061 ; basefiterror[ibase]=0.000018 ; baserescalederror[ibase]=0.000020 ; basescalefactor[ibase]=1.10 ;  ++ibase; // 25
  basegamma.push_back(104); baseparm.push_back(  4) ; baseseed[ibase] = 2.000000E-02 ; basetitle[ibase] = "tau- --> 3h- 2h+ pi0 nu(tau) (ex.K0)"; basefitvalue[ibase]=0.000181 ; basefiterror[ibase]=0.000026 ; baserescalederror[ibase]=0.000026 ; basescalefactor[ibase]=1.01 ;  ++ibase; // 26
  basegamma.push_back(126); baseparm.push_back( 58) ; baseseed[ibase] = 1.700000E-01 ; basetitle[ibase] = "tau- --> eta pi- pi0 nu(tau)"; basefitvalue[ibase]=0.001774 ; basefiterror[ibase]=0.000235 ; baserescalederror[ibase]=0.000236 ; basescalefactor[ibase]=1.00 ; ++ibase; // 27
  basegamma.push_back(128); baseparm.push_back(109) ; baseseed[ibase] = 3.000000E-02 ; basetitle[ibase] = "tau- --> eta K- nu(tau)";      basefitvalue[ibase]=0.000268 ; basefiterror[ibase]=0.000063 ; baserescalederror[ibase]=0.000063 ; basescalefactor[ibase]=1.00 ; ++ibase; // 28
  basegamma.push_back(150); baseparm.push_back(  8) ; baseseed[ibase] = 2.000000E+00 ; basetitle[ibase] = "tau- --> h- omega nu(tau)";    basefitvalue[ibase]=0.019860 ; basefiterror[ibase]=0.000638 ; baserescalederror[ibase]=0.000788 ; basescalefactor[ibase]=1.24 ; ++ibase; // 29
  basegamma.push_back(152); baseparm.push_back(113) ; baseseed[ibase] = 4.000000E-01 ; basetitle[ibase] = "tau- --> h- omega pi0 nu(tau)";basefitvalue[ibase]=0.004063 ; basefiterror[ibase]=0.000418 ; baserescalederror[ibase]=0.000429 ; basescalefactor[ibase]=1.03 ; ++ibase; // 30
  basegamma.push_back(103); baseparm.push_back(  3) ; baseseed[ibase] = 8.000000E-02 ; basetitle[ibase] = "tau- --> 3h- 2h+ nu(tau) (ex.K0)"; basefitvalue[ibase]=0.000810 ; basefiterror[ibase]=0.000053 ; baserescalederror[ibase]=0.000055 ; basescalefactor[ibase]=1.05 ;  ++ibase; // 31
  // FOR CROSS-CHECKING PDG FIT THE FOLLOWING SHOULD NOT BE USED
  //  basegamma.push_back(130); baseparm.push_back(266) ; baseseed[ibase] = 4.600000E-03 ; basetitle[ibase] = "tau- --> eta K- pi0 nu(tau)"; ++ibase;//32
  //  basegamma.push_back(132); baseparm.push_back(267) ; baseseed[ibase] = 8.800000E-03 ; basetitle[ibase] = "tau- --> eta Kbar0 pi- nu(tau)"; ++ibase;//33
  //
  double totalseed = 0;
  for (ibase=0;ibase<nbase;++ibase) {
    baseseed[ibase]*=1.e-2; // percent to fraction
    totalseed+=baseseed[ibase];
  }
  //  cout << "Totalseed = " << totalseed << " is re-normalized to 1 by adjusting each baseseed value ... " << endl << endl;
  //  for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]/=totalseed;
  //
  // From previous iteration :
  ibase=-1;
  baseseed[++ibase] = /*           M_GAMMA3    */  0.1735575 ; //+/-      0.0004626 +/-      0.0000000 Tot Err:     0.0004626 Check Sys:     0.0000000
  baseseed[++ibase] = /*           M_GAMMA5    */  0.1784136 ; //+/-      0.0004788 +/-      0.0000000 Tot Err:     0.0004788 Check Sys:     0.0000000
  baseseed[++ibase] = /*           M_GAMMA9    */  0.1090047 ; //+/-      0.0006175 +/-      0.0000000 Tot Err:     0.0006175 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA10    */  0.0069093 ; //+/-      0.0002188 +/-      0.0000000 Tot Err:     0.0002188 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA14    */  0.2549936 ; //+/-      0.0009233 +/-      0.0000000 Tot Err:     0.0009233 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA16    */  0.0045244 ; //+/-      0.0002645 +/-      0.0000000 Tot Err:     0.0002645 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA20    */  0.0924794 ; //+/-      0.0009745 +/-      0.0000000 Tot Err:     0.0009745 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA23    */  0.0005801 ; //+/-      0.0002251 +/-      0.0000000 Tot Err:     0.0002251 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA27    */  0.0104186 ; //+/-      0.0007060 +/-      0.0000000 Tot Err:     0.0007060 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA28    */  0.0004181 ; //+/-      0.0002186 +/-      0.0000000 Tot Err:     0.0002186 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA30    */  0.0010186 ; //+/-      0.0003920 +/-      0.0000000 Tot Err:     0.0003920 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA35    */  0.0089618 ; //+/-      0.0003671 +/-      0.0000000 Tot Err:     0.0003671 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA37    */  0.0015298 ; //+/-      0.0001606 +/-      0.0000000 Tot Err:     0.0001606 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA40    */  0.0037774 ; //+/-      0.0003676 +/-      0.0000000 Tot Err:     0.0003676 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA42    */  0.0015407 ; //+/-      0.0001995 +/-      0.0000000 Tot Err:     0.0001995 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA47    */  0.0002411 ; //+/-      0.0000516 +/-      0.0000000 Tot Err:     0.0000516 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA48    */  0.0011208 ; //+/-      0.0002454 +/-      0.0000000 Tot Err:     0.0002454 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA62    */  0.0898792 ; //+/-      0.0006003 +/-      0.0000000 Tot Err:     0.0006003 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA70    */  0.0268568 ; //+/-      0.0006932 +/-      0.0000000 Tot Err:     0.0006932 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA77    */  0.0009086 ; //+/-      0.0003595 +/-      0.0000000 Tot Err:     0.0003595 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA78    */  0.0002233 ; //+/-      0.0000499 +/-      0.0000000 Tot Err:     0.0000499 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA85    */  0.0033354 ; //+/-      0.0002232 +/-      0.0000000 Tot Err:     0.0002232 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA89    */  0.0007347 ; //+/-      0.0001157 +/-      0.0000000 Tot Err:     0.0001157 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA93    */  0.0015313 ; //+/-      0.0000699 +/-      0.0000000 Tot Err:     0.0000699 Check Sys:     0.0000000
  baseseed[++ibase] = /*          M_GAMMA94    */  0.0000609 ; //+/-      0.0000183 +/-      0.0000000 Tot Err:     0.0000183 Check Sys:     0.0000000
  baseseed[++ibase] = /*         M_GAMMA104    */  0.0001807 ; //+/-      0.0000260 +/-      0.0000000 Tot Err:     0.0000260 Check Sys:     0.0000000
  baseseed[++ibase] = /*         M_GAMMA126    */  0.0017735 ; //+/-      0.0002354 +/-      0.0000000 Tot Err:     0.0002354 Check Sys:     0.0000000
  baseseed[++ibase] = /*         M_GAMMA128    */  0.0002679 ; //+/-      0.0000632 +/-      0.0000000 Tot Err:     0.0000632 Check Sys:     0.0000000
  baseseed[++ibase] = /*         M_GAMMA150    */  0.0198809 ; //+/-      0.0006285 +/-      0.0000000 Tot Err:     0.0006285 Check Sys:     0.0000000
  baseseed[++ibase] = /*         M_GAMMA152    */  0.0040670 ; //+/-      0.0004175 +/-      0.0000000 Tot Err:     0.0004175 Check Sys:     0.0000000
  baseseed[++ibase] = 1- 0.9991898;
  // INPUT NODES
  // NODE# NODE    EQ ADJ. UN. COUNT PAR_CODE PARAM.  SUM COEFFICIENT COEFF-EXTRA  
  // ----- ------- -- -------- ----- -------- ------- --- ----------- ------------ 
  const int      nnode=95;
  vector<string> nodename;
  double         nodeunit[nnode];
  int            node_num_npar[nnode];
  vector<int>    node_num_parm[nnode];
  vector<double> node_num_coef[nnode];
  int            node_den_npar[nnode];
  vector<int>    node_den_parm[nnode];
  vector<double> node_den_coef[nnode];
  string         nodetitle[nnode];
  int inode=0;
  //
  //1
  nodename.push_back("S035B20"); nodeunit[inode]=1.0E-04; nodetitle[inode]="G(eta K- nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //2
  nodename.push_back("S035B21"); nodeunit[inode]=1.;      nodetitle[inode]="G(h- 2pi0 nu(tau) (ex.K0))/G(h- pi0 nu(tau))";
  node_num_npar[inode]=2; node_den_npar[inode]=2;
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(16 );       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(182);       node_den_coef[inode].push_back(1.        );               
  ++inode;
  //3  
  nodename.push_back("S035B22"); nodeunit[inode]=1.;      nodetitle[inode]="G(h- 3pi0 nu(tau))/G(h- pi0 nu(tau))";
  node_num_npar[inode]=5; node_den_npar[inode]=2;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2200E-01);              
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(16 );       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(182);       node_den_coef[inode].push_back(1.        );               
  ++inode;
  //4  
  nodename.push_back("S035B23"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(h- 4pi0 nu(tau) (ex.K0,eta))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //5
  nodename.push_back("S035B25"); nodeunit[inode]=1.;         nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex.K0))/G(h- h- h+ >=0 neutrals >=0K(L)0 nu(tau) )";
  node_num_npar[inode]=3; node_den_npar[inode]=18;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);              
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(2.8500E-01);               
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(9.1010E-01);               
  node_den_parm[inode].push_back(117);       node_den_coef[inode].push_back(3.4310E-01);               
  node_den_parm[inode].push_back(118);       node_den_coef[inode].push_back(3.4310E-01);               
  node_den_parm[inode].push_back(119);       node_den_coef[inode].push_back(3.4310E-01);               
  node_den_parm[inode].push_back(204);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(214);       node_den_coef[inode].push_back(4.3070E-01);               
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(240);       node_den_coef[inode].push_back(6.8610E-01);               
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(260);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(5  );       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(2.8500E-01);               
  node_den_parm[inode].push_back(62 );       node_den_coef[inode].push_back(3.4310E-01);               
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(9.1010E-01);               
  ++inode;
  //6
  nodename.push_back("S035B26"); nodeunit[inode]=1.;         nodetitle[inode]="G(h- omega pi0 nu(tau))/G(h- h- h+ >=0 neutrals >=0K(L)0 nu(tau) )";
  node_num_npar[inode]=1; node_den_npar[inode]=18;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(2.8500E-01);               
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(9.1010E-01);               
  node_den_parm[inode].push_back(117);       node_den_coef[inode].push_back(3.4310E-01);               
  node_den_parm[inode].push_back(118);       node_den_coef[inode].push_back(3.4310E-01);               
  node_den_parm[inode].push_back(119);       node_den_coef[inode].push_back(3.4310E-01);               
  node_den_parm[inode].push_back(204);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(214);       node_den_coef[inode].push_back(4.3070E-01);               
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(240);       node_den_coef[inode].push_back(6.8610E-01);               
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(260);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(5  );       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(2.8500E-01);               
  node_den_parm[inode].push_back(62 );       node_den_coef[inode].push_back(3.4310E-01);               
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(9.1010E-01);               
  ++inode;
  //7
  nodename.push_back("S035B27"); nodeunit[inode]=1.;         nodetitle[inode]="G(h- omega pi0 nu(tau))/G(h- h- h+ 2pi0 nu(tau) (ex.K0))";
  node_num_npar[inode]=1; node_den_npar[inode]=3;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(8.8800E-01);               
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(2.2600E-01);               
  ++inode;
  //8
  nodename.push_back("S035B29"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(K- pi0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //9
  nodename.push_back("S035B30"); nodeunit[inode]=1.0E-04;    nodetitle[inode]="G(K- 2pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //10
  nodename.push_back("S035B31"); nodeunit[inode]=1.0E-04;    nodetitle[inode]="G(K- 3pi0 nu(tau) (ex.K0, eta))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //11
  nodename.push_back("S035B32"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(pi- Kbar0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //12
  nodename.push_back("S035B33"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(pi- Kbar0 pi0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //13
  nodename.push_back("S035B34"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(K- K0 pi0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //14
  nodename.push_back("S035B37"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(K- K+ pi- >=0 neut.  nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //15  
  nodename.push_back("S035B43"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(K(S)0 (particles)- nu(tau))/G(total)";
  node_num_npar[inode]=6; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(5.0000E-01);              
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(5.0000E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(5.0000E-01);               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(5.0000E-01);               
  ++inode;
  //16
  nodename.push_back("S035B45"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(( 5pi )- nu(tau))/G(total)";
  node_num_npar[inode]=6; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(5.5300E-01);               
  ++inode;
  //17  
  nodename.push_back("S035B51"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(pi- K0 Kbar0 nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(2.0000E+00);              
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //18
  nodename.push_back("S035B53"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(h- h- h+ pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=6; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.7000E-02);               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //19
  nodename.push_back("S035B54"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(h- h- h+ pi0 nu(tau) (ex. K0, omega))/G(total)";
  node_num_npar[inode]=4; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //20
  nodename.push_back("S035B55"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(pi- 2pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //21
  nodename.push_back("S035B56"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(pi- 3pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //22
  nodename.push_back("S035B57"); nodeunit[inode]=1.0E-04;    nodetitle[inode]="G(h- h- h+ 3pi0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //23
  nodename.push_back("S035B58"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(h- omega pi0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //24
  nodename.push_back("S035B59"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=3; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);              
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  ++inode;
  //25
  nodename.push_back("S035B62"); nodeunit[inode]=1.0E-02;    nodetitle[inode]="G(h- h- h+ nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=4; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);               
  ++inode;
  //26
  nodename.push_back("S035B63"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G( h- h- h+ >=0 neutrals  nu(tau)  (ex. K(S)0 --> pi+ pi-)(``3-prong''))/G(total)";
  node_num_npar[inode]=12; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.8500E-01);              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(9.1010E-01);               
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.8500E-01);               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.1010E-01);               
  ++inode;
  //27 
  nodename.push_back("S035B64"); nodeunit[inode]=1.;        nodetitle[inode]="G(h- h- h+ nu(tau) (ex.K0))/G( h- h- h+ >=0 neutrals  nu(tau)  (ex. K(S)0 --> pi+ pi-)(``3-prong''))";
  node_num_npar[inode]=4; node_den_npar[inode]=12;
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );            
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);               
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(2.8500E-01);               
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(9.1010E-01);               
  node_den_parm[inode].push_back(204);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(260);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(5  );       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(2.8500E-01);               
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(9.1010E-01);               
  ++inode;
  //28 
  nodename.push_back("S035B67"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- Kbar0 nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //29 
  nodename.push_back("S035B68"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- Kbar0 pi0 nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //30 
  nodename.push_back("S035B69"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(pi- K(S)0 K(S)0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //31 
  nodename.push_back("S035B71"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- h- h+ nu(tau) (ex.K0,omega))/G(total)";
  node_num_npar[inode]=3; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //32 
  nodename.push_back("S035B72"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex.K0,omega,eta))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //33 
  nodename.push_back("S035B73"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //34 
  nodename.push_back("S035B74"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- 2pi0 nu(tau))/G(total)";
  node_num_npar[inode]=4; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.5700E-01);               
  ++inode;
  //35 
  nodename.push_back("S035B75"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(particle->=0 neutrals >=0K0 nu(tau)  (``1-prong''))/G(total)";
  node_num_npar[inode]=21; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(7.1500E-01);               
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(9.0000E-02);               
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(2  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(2.0000E+00);               
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(7.0800E-01);               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.0000E-02);               
  ++inode;
  //36 
  nodename.push_back("S035B76"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- h- h+ pi0 nu(tau))/G(total)";
  node_num_npar[inode]=8; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.7000E-02);               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //37 
  nodename.push_back("S035B77"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau))/G(total)";
  node_num_npar[inode]=4; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);              
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(4.3070E-01);               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  ++inode;
  //38 
  nodename.push_back("S035B78"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- h- h+ >=1 pi0 nu(tau) (ex. K0))/G(total)";
  node_num_npar[inode]=9; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(9.1010E-01);               
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //39 
  nodename.push_back("S035B79"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- 4pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.1900E-01);               
  ++inode;
  //40 
  nodename.push_back("S035B97"); nodeunit[inode]=1.;        nodetitle[inode]="G(h- nu(tau))/G(e- nubar(e) nu(tau))";
  node_num_npar[inode]=2; node_den_npar[inode]=1;
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(1.        );               
  ++inode;
  //41 
  nodename.push_back("S035C01"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- >= 1pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=9; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2500E-01);              
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.2500E-01);               
  ++inode;
  //42 
  nodename.push_back("S035C02"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- >= 3pi0 nu(tau) (ex. K0))/G(total)";
  node_num_npar[inode]=5; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2500E-01);              
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.2500E-01);               
  ++inode;
  //43 
  nodename.push_back("S035C03"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- h- h+ >= 2pi0 nu(tau) (ex. K0))/G(total)";
  node_num_npar[inode]=4; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);              
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  ++inode;
  //44 
  nodename.push_back("S035C1"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(pi- K(S)0 K(L)0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //45 
  nodename.push_back("S035C18"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(pi- pi+ pi- nu(tau))/G(total)";
  node_num_npar[inode]=3; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(3.4310E-01);              
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);               
  ++inode;
  //46 
  nodename.push_back("S035C19"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(pi- pi+ pi- nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);               
  ++inode;
  //47 
  nodename.push_back("S035C20"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(pi- pi+ pi- nu(tau) (ex.K0,omega))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //48 
  nodename.push_back("S035C21"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- pi+ pi- nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //49 
  nodename.push_back("S035C22"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(pi- pi+ pi- pi0 nu(tau))/G(total)";
  node_num_npar[inode]=4; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.7000E-02);              
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //50 
  nodename.push_back("S035C23"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(pi- pi+ pi- pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=3; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.7000E-02);              
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //51 
  nodename.push_back("S035C24"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(pi- pi+ pi- pi0 nu(tau) (ex.K0,omega))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //52 
  nodename.push_back("S035C25"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(K- pi+ pi- pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //53 
  nodename.push_back("S035C31"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- h+ pi- nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //54 
  nodename.push_back("S035C32"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- h+ pi- nu(tau) (ex.K0))/G(pi- pi+ pi- nu(tau) (ex.K0))";
  node_num_npar[inode]=2; node_den_npar[inode]=2;
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(1.7000E-02);               
  ++inode;
  //55 
  nodename.push_back("S035C33"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(K- h+ pi- pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=3; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //56 
  nodename.push_back("S035C34"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- h+ pi- pi0 nu(tau) (ex.K0))/G(pi- pi+ pi- pi0 nu(tau) (ex.K0))";
  node_num_npar[inode]=3; node_den_npar[inode]=3;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(1.7000E-02);               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //57 
  nodename.push_back("S035C35"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- K+ pi- nu(tau))/G(pi- pi+ pi- nu(tau) (ex.K0))";
  node_num_npar[inode]=1; node_den_npar[inode]=2;
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(1.7000E-02);               
  ++inode;
  //58 
  nodename.push_back("S035C36"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- K+ pi- pi0 nu(tau))/G(pi- pi+ pi- pi0 nu(tau) (ex.K0))";
  node_num_npar[inode]=1; node_den_npar[inode]=3;
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(1.7000E-02);               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //59 
  nodename.push_back("S035C38"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- K0 >=0pi0 nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //60 
  nodename.push_back("S035C40"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- pi+ pi- >=0pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=3; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //61 
  nodename.push_back("S035C47"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(X- (S=-1) nu(tau))/G(total)";
  node_num_npar[inode]=9; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //62 
  nodename.push_back("S035C54"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(K- pi+ pi- pi0 nu(tau) (ex.K0,eta))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //63 
  nodename.push_back("S035C6"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- pi+ pi- nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  ++inode;
  //64 
  nodename.push_back("S035C7"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(K- pi+ pi- pi0 nu(tau))/G(total)";
  node_num_npar[inode]=3; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //65 
  nodename.push_back("S035C8"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(K- K+ pi- pi0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );             
  ++inode;
  //66 
  nodename.push_back("S035R1"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(mu- nubar(mu) nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //67 
  nodename.push_back("S035R14"); nodeunit[inode]=1.;        nodetitle[inode]="G(h- omega nu(tau))/G(h- h- h+ pi0 nu(tau) (ex.K0))";
  node_num_npar[inode]=1; node_den_npar[inode]=6;
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(2.2600E-01);               
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(1.7000E-02);               
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //68 
  nodename.push_back("S035R15"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- omega  >=0 neutrals  nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //69 
  nodename.push_back("S035R2"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(e- nubar(e) nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(2  );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //70 
  nodename.push_back("S035R20"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- 2pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //71 
  nodename.push_back("S035R21"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- 3pi0 nu(tau))/G(total)";
  node_num_npar[inode]=5; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2200E-01);              
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //72 
  nodename.push_back("S035R23"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- omega nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //73 
  nodename.push_back("S035R24"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(particle->=0 neutrals >=0K(L)0 nu(tau) )/G(total)";
  node_num_npar[inode]=21; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(7.1500E-01);               
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(9.0000E-02);               
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(6.5690E-01);               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(6.5690E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(6.5690E-01);               
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(2  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.0985E+00);               
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(3.1390E-01);               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(7.0800E-01);               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(6.5690E-01);               
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.0000E-02);               
  ++inode;
  //74 
  nodename.push_back("S035R26"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- >=0pi0 >=0K0  >=0gamma  nu(tau))/G(total)";
  node_num_npar[inode]=7; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(7.1500E-01);              
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //75 
  nodename.push_back("S035R27"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- >=1 (pi0 or K0 or gamma)  nu(tau))/G(total)";
  node_num_npar[inode]=6; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(7.1500E-01);              
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //76 
  nodename.push_back("S035R28"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- h- h+ nu(tau))/G(total)";
  node_num_npar[inode]=6; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(3.4310E-01);              
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);               
  ++inode;
  //77 
  nodename.push_back("S035R30"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- h- h+ >=1 neutrals  nu(tau))/G(total)";
  node_num_npar[inode]=13; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.8500E-01);              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(9.1010E-01);               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(4.3070E-01);               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(6.8610E-01);               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.8500E-01);               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);               
  ++inode;
  //78 
  nodename.push_back("S035R31"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- h- h+ >=0 neutrals >=0K(L)0 nu(tau) )/G(total)";
  node_num_npar[inode]=18; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.8500E-01);              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(9.1010E-01);               
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(4.3070E-01);               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(6.8610E-01);               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.8500E-01);               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.1010E-01);               
  ++inode;
  //79 
  nodename.push_back("S035R32"); nodeunit[inode]=1.0E-03;   nodetitle[inode]="G(eta pi- pi0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //80 
  nodename.push_back("S035R33"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G( 3h- 2h+ >=0 neutrals  nu(tau)  (ex. K(S)0 --> pi- pi+)(``5-prong''))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(4  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //81 
  nodename.push_back("S035R34"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- h+ h- >=0 neutrals  nu(tau))/G(total)";
  node_num_npar[inode]=7; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.8500E-01);              
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  ++inode;
  //82 
  nodename.push_back("S035R38"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(3h- 2h+ nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //83 
  nodename.push_back("S035R39"); nodeunit[inode]=1.0E-04;   nodetitle[inode]="G(3h- 2h+ pi0 nu(tau) (ex.K0))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(4  );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //84 
  nodename.push_back("S035R40"); nodeunit[inode]=1.0E-03;   nodetitle[inode]="G(K- K+ pi- nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //85 
  nodename.push_back("S035R41"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- pi+ pi- >=0 neutrals  nu(tau))/G(total)";
  node_num_npar[inode]=5; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.8500E-01);              
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  ++inode;
  //86 
  nodename.push_back("S035R42"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- >=1 neutrals nu(tau))/G(total)";
  node_num_npar[inode]=16; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(7.1500E-01);              
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(9.0000E-02);               
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(9.8500E-02);               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(7.0800E-01);               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.0000E-02);               
  ++inode;
  //87 
  nodename.push_back("S035R43"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- >=0K(L)0  nu(tau))/G(total)";
  node_num_npar[inode]=5; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(5.0000E-01);              
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(5.0000E-01);               
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //88 
  nodename.push_back("S035R44"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- >= 2pi0 nu(tau))/G(total)";
  node_num_npar[inode]=12; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2200E-01);              
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(9.8500E-02);               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.1900E-01);               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.5700E-01);               
  ++inode;
  //89 
  nodename.push_back("S035R46"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- K0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //90 
  nodename.push_back("S035R5"); nodeunit[inode]=1.;        nodetitle[inode]="G(mu- nubar(mu) nu(tau))/G(e- nubar(e) nu(tau))";
  node_num_npar[inode]=1; node_den_npar[inode]=1;
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(1.        );               
  ++inode;
  //91 
  nodename.push_back("S035R6"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(pi- nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //92 
  nodename.push_back("S035R7"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(K- nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //93 
  nodename.push_back("S035R8"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(pi- pi0 nu(tau))/G(total)";
  node_num_npar[inode]=1; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );              
  ++inode;
  //94 
  nodename.push_back("S035R84"); nodeunit[inode]=1.0E-02;   nodetitle[inode]="G(h- pi0 nu(tau))/G(total)";
  node_num_npar[inode]=2; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  ++inode;
  //95 
  nodename.push_back("S035R97"); nodeunit[inode]=1.0E-02; nodetitle[inode]="G(h- >= 3pi0 nu(tau))/G(total)";
  node_num_npar[inode]=8; node_den_npar[inode]=0;
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2200E-01);              
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(9.8500E-02);               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.1900E-01);               
  ++inode;
  //
  // Input R [node] as a linearized sum of variables P_i [base quantities]
  //
  // General Case:
  // R(P_i) = (sum_j=1,n alpha_j P_j) / (sum_k=1,m beta_k P_k)
  // dR/dP_i= (1/sum_k=1,m beta_i P_i) * alpha_j * delta_ij - (R/sum_k=1,m beta_k P_k) * beta_k * delta_ik
  //
  // Special Case 1:
  // R(P_i) = (sum_j=1,n alpha_j P_j) 
  // dR/dP_i= alpha_i 
  //
  // Special Case 2: 
  // R(P_i) = (sum_j=1,n alpha_j P_j) / (sum_k=1,m beta_k P_k)
  // Case 2a: if P_i only in numerator: dR/dP_i = (1/sum_k=1,m beta_k P_k) * alpha_i
  // Case 2b: if P_i only in denominator: dR/dP_i = - (R/sum_k=1,m beta_k P_k) * beta_i
  // Case 2c: if P_i both in numerator and denominator: dR/dP_i = (1/sum_k=1,m beta_i P_i) * alpha_i - (R/sum_k=1,m beta_k P_k) * beta_i
  //
  vector<int> node_parm[nnode]; // vector of parameters for each node
  vector<int> node_quan[nnode]; // vector of quantities for each node
  vector<int> node_gamma[nnode]; // vector of gamma for each node
  double node_num[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double node_den[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  vector<double> node_part[nnode]; // vector of partial derivatives w.r.t quantities for each node
  int first_quan[nnode];  // first quantity measured for each node
  for (inode=0;inode<nnode;++inode){
    node_parm[inode].insert(node_parm[inode].end(),node_num_parm[inode].begin(),node_num_parm[inode].end());
    node_parm[inode].insert(node_parm[inode].end(),node_den_parm[inode].begin(),node_den_parm[inode].end());
    sort(node_parm[inode].begin(), node_parm[inode].end());//sort must be called before unique
    vector<int>::iterator new_end=unique(node_parm[inode].begin(), node_parm[inode].end());
    node_parm[inode].erase(new_end, node_parm[inode].end()); // <--
    //
    double numerator, denominator, partial;
    //
    numerator=0; // sum_j=1,n alpha_j P_j
    if (node_num_npar[inode]>0) {
      for (ipar=0;ipar<node_num_npar[inode];++ipar){
	int parm=node_num_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	numerator+=(node_num_coef[inode].at(ipar))*(baseseed[quan-1]); 
	//	numerator+=(node_num_coef[inode].at(ipar))*(basefitvalue[quan-1]);
      }
    }
    node_num[inode]=numerator; // <--
    //
    denominator=0; // sum_k=1,m beta_k P_k
    if (node_den_npar[inode]>0) {
      for (ipar=0;ipar<node_den_npar[inode];++ipar){
	int parm=node_den_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	denominator+=(node_den_coef[inode].at(ipar))*(baseseed[quan-1]); 
	//	denominator+=(node_den_coef[inode].at(ipar))*(basefitvalue[quan-1]);
      }
    }
    node_den[inode]=denominator; // <--
    //
    first_quan[inode]=99; // <--
    for (ipar=0;ipar<node_parm[inode].size();++ipar){
      int parm=node_parm[inode].at(ipar);
      vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
      int quan=ibase-baseparm.begin()+1;
      node_quan[inode].push_back(quan); // <--
      int gamma = basegamma.at(quan-1); 
      node_gamma[inode].push_back(gamma); // <--
      if (quan<first_quan[inode]) first_quan[inode]=quan; // <--
      //
      vector<int>::iterator it_num=find(node_num_parm[inode].begin(),node_num_parm[inode].end(),parm); bool is_in_num = it_num != node_num_parm[inode].end();
      vector<int>::iterator it_den=find(node_den_parm[inode].begin(),node_den_parm[inode].end(),parm); bool is_in_den = it_den != node_den_parm[inode].end();
      partial=0;
      if (node_den_npar[inode]==0 && is_in_num) {//case 1
	partial = node_num_coef[inode].at(it_num - node_num_parm[inode].begin());
      } else if ( is_in_num && !is_in_den) { // case 2a
	partial = (1./denominator) * (node_num_coef[inode].at(it_num - node_num_parm[inode].begin()));
      } else if (!is_in_num &&  is_in_den) { // case 2b
	partial = -1. * (numerator/(denominator*denominator)) * (node_den_coef[inode].at(it_den - node_den_parm[inode].begin()));
      } else if ( is_in_num &&  is_in_den) { // case 2c
	partial = (1./denominator) * (node_num_coef[inode].at(it_num - node_num_parm[inode].begin())) -1. * (numerator/(denominator*denominator)) * (node_den_coef[inode].at(it_den - node_den_parm[inode].begin()));
      }
      node_part[inode].push_back(partial); // <--
      //      if (inode==14) cout << " inode = " << inode << " ipar = " << ipar << " parm = " << parm << " quan = " << quan << " gamma = " << gamma << " partial = " << partial 
      //			  << " it_num - node_num_parm[inode].begin() = " << it_num - node_num_parm[inode].begin() 
      //			  << endl;
    }
  }
  //
  int nmeas=0;
  int imeas[200], measweak[200], imeas1, imeas2;
  string gammaname[200], measnodename[200], author[200], year[200], meastitle[200];
  double fitvalue[200], fiterror[200], scalederror[200], measvalue[200], measerror[200], chisquare[200], pullav[200], scale[200], corrtemp, corrmat[200][200];
  for (imeas1=0;imeas1<200;imeas1++) for (imeas2=0;imeas2<200;imeas2++) corrmat[imeas1][imeas2] = 0;
  //
  ifstream ifs("s035-fit-no-babar-belle.data") ;
  if (!ifs.good()) {cout << "Cannot open input file " << endl ; exit(1) ;}
  char buffer[256]; 
  while(ifs.good()) {
    if (ifs.eof()) break;
    char firstch(' ') ; ifs.get(firstch) ;
    if (firstch=='#'||firstch=='\n') { // Skip this linex
      ifs.ignore(256,'\n') ;
    } else if (firstch=='*') {  // measurement line
      ifs >> imeas[nmeas] >> gammaname[nmeas] >> measnodename[nmeas]
	  >> fitvalue[nmeas] >> fiterror[nmeas] >> scalederror[nmeas]
	  >> measvalue[nmeas] >> measerror[nmeas] >> measweak[nmeas] >> chisquare[nmeas] >> pullav[nmeas] >> scale[nmeas]
	  >> author[nmeas] >> year[nmeas];
      ifs.getline(buffer,256,'\n');
      string stringbuffer=string(buffer);
      meastitle[nmeas]="";
      bool first=false;
      for (string::iterator si=stringbuffer.begin();si!=stringbuffer.end();++si){
	if (13 == (int)*si) break; // ^M is 13
	if (' ' != *si) first=true;
	if (first) meastitle[nmeas]+=*si;
      }
//cout << imeas[nmeas] << " " << gammaname[nmeas] << " " << measnodename[nmeas] << " " 
//	   << fitvalue[nmeas] << " " << fiterror[nmeas] << " " << scalederror[nmeas] << " " 
//	   << measvalue[nmeas] << " " << measerror[nmeas] << " " << measweak[nmeas] << " " << chisquare[nmeas] << " " << pullav[nmeas] << " " << scale[nmeas] <<  " " 
//	   << author[nmeas] << " " << year[nmeas] << " " << meastitle[nmeas] << endl;
      ++nmeas;
    } else if (firstch=='%') {  // correlation line
      ifs >> imeas1 >> imeas2 >> corrtemp;       ifs.ignore(256,'\n') ;
      if (corrmat[imeas1-1][imeas2-1] != 0) {cout << "WATCH OUT" << endl; exit(1);}
      corrmat[imeas1-1][imeas2-1] = corrmat[imeas2-1][imeas1-1] = corrtemp*1e-2;
      //      cout << "Correlation between imeas1 = " << imeas[imeas1-1] << " and imeas2 = " << imeas[imeas2-1] << " is = " << corrmat[imeas1-1][imeas2-1] << endl;
    }
  }
  //
  //  cout << "nmeas = " << nmeas << endl;
  FILE *measfile[2];
  for (int p=0;p<2;++p){
    if (p==0) measfile[0]=fopen("combos_pdginput_measurements.input","w");
    if (p==1) measfile[1]=fopen("alucomb_pdginput_measurements.input","w");
    for (int i=0;i<nmeas;++i){
      //      cout << "i = " << i << endl;
      fprintf (measfile[p], "* IMEAS = %d \n",     imeas[i]);
      fprintf (measfile[p], "* GAMMANAME = %s \n", gammaname[i].data());
      fprintf (measfile[p], "* DECAYNAME = %s \n", meastitle[i].data());
      fprintf (measfile[p], "* IWEAK = %d \n",   measweak[i]);
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      fprintf (measfile[p], "* NODENAME = %s found at inode+1 = %d with NODENAME[inode] = %s has %d, %d, %d base quantities in numerator, denominator and both [excluding overlap] \n", measnodename[i].data(), inode+1, nodename[inode].data(),node_num_npar[inode], node_den_npar[inode],node_parm[inode].size());
      for (ipar=0;ipar<node_num_parm[inode].size();++ipar){
	int parm=node_num_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	fprintf (measfile[p],"*                numerator of inode+1 = %d has gamma = %d parm = %d quan = %d title = %s seed = %f coef = %f\n",inode+1, basegamma[quan-1], parm, quan, basetitle[quan-1].data(), baseseed[quan-1], node_num_coef[inode].at(ipar));
      }
      for (ipar=0;ipar<node_den_parm[inode].size();++ipar){
	int parm=node_den_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	fprintf (measfile[p],"*              denominator of inode+1 = %d has gamma = %d parm = %d quan = %d title = %s seed = %f coef = %f\n",inode+1, basegamma[quan-1], parm, quan, basetitle[quan-1].data(), baseseed[quan-1], node_den_coef[inode].at(ipar));
      }
      fprintf (measfile[p], "*  first quantity measured by inode+1 = %d has gamma = %d parm = %d quan = %d title = %s\n",inode+1,basegamma[first_quan[inode]-1],baseparm[first_quan[inode]-1],first_quan[inode],basetitle[first_quan[inode]-1].data());
      //
      fprintf (measfile[p], "\nBEGIN %s Gamma%s published.%s \n\n", author[i].data(), gammaname[i].data(), year[i].data());
      if (p==0) {//COMBOS
	if (// SPECIAL CASE [because these nodes contain Gamma103 [used to express unitarity constraint]
	    (inode+1)==80 || // NODE = 79 NAME = S035R33 GAMMA = 102
	    (inode+1)==82) { // NODE = 81 NAME = S035R38 GAMMA = 103
	  fprintf (measfile[p], "MEASUREMENT  m_Gamma%d statistical systematic \n",3);
	  fprintf (measfile[p], "DATA         m_Gamma%d statistical systematic \n",3);
	}else{
	  fprintf (measfile[p], "MEASUREMENT  m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
	  fprintf (measfile[p], "DATA         m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
	}
      }else if (p==1) {//ALUCOMB
	fprintf (measfile[p], "MEASUREMENT  m_Gamma%s statistical systematic \n",gammaname[i].data());
	fprintf (measfile[p], "DATA         m_Gamma%s statistical systematic \n",gammaname[i].data());
      }
      fprintf (measfile[p], "             %10.5g %10.5g  0 \n",measvalue[i],measerror[i]);
      bool firstcorr=true;
      for (int j=0;j<nmeas;++j) {
	if (corrmat[i][j]!=0) {
	  if (firstcorr) {fprintf (measfile[p], " \n"); firstcorr=false;}
	  fprintf (measfile[p], "STAT_CORR_WITH %s Gamma%s published.%s %f ! IMEAS = %d \n",author[j].data(), gammaname[j].data(), year[j].data(), corrmat[i][j], imeas[j]);
	}
      }
      fprintf (measfile[p], " \nEND \n\n");
    }
    fclose(measfile[p]);
  }
  //
  FILE *avefile[2];
  for (int p=0;p<2;++p){
    if (p==0) avefile[0]=fopen("combos_average_pdginput_no_babar_belle.input","w");
    if (p==1) avefile[1]=fopen("alucomb_average_pdginput_no_babar_belle.input","w");
    if (p==0) fprintf (avefile[p], "INCLUDE combos_pdginput_measurements.input \n\n"); 
    if (p==1) fprintf (avefile[p], "INCLUDE alucomb_pdginput_measurements.input \n\n"); 
    fprintf (avefile[p], "BEGIN   PDG-BABAR-BELLE all_methods \n\n");
    fprintf (avefile[p], "COMBINE * * * \n\n");
    for (ibase=0;ibase<nbase;++ibase){
      if (p==0&&ibase==(nbase-1)){/* skip */}else{
	fprintf (avefile[p], "MEASUREMENT m_Gamma%d statistical systematic   ! NQUAN = %d \n",basegamma.at(ibase),ibase+1);  
      }
    }
    if (p==0){
      fprintf (avefile[p], "\nCALL DUMP_MASTER_INC \n\n");
      fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_PRT -1.0 0. \n\n");
      fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_INV 0 0    \n\n");
    }
    //
    int lastnode=-1;
    int usum=0;//number of [unique] measurements to be expressed as linearized sum of base quantities
    for (int i=0;i<nmeas;++i){
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      if ((node_num_npar[inode]+node_den_npar[inode])>1) {
	if (p==1&&inode!=lastnode){//new node
	  ++usum;
	  lastnode=inode;
	  fprintf (avefile[p], "MEASUREMENT m_Gamma%s statistical systematic   ! NQUAN = %d \n",gammaname[i].data(),ibase+usum);  
	}
      }
    }
    //
    int isum=0;//number of          measurements to be expressed as linearized sum of base quantities
    lastnode=-1;
    for (int i=0;i<nmeas;++i){
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      if ((node_num_npar[inode]+node_den_npar[inode])>1) {
	//
	if (p==0&&(((inode+1)==80)||((inode+1)==82))) continue; // SPECIAL CASE [because these are derived nodes containing Gamma103 ]
	++isum; // translate C index to Fortran index
	//
	if (inode!=lastnode){//new node
	  lastnode=inode;
	}else{
	  if (p==1) continue;
	}
	//
	// PRINT NODE DEFINITION
	//
	fprintf (avefile[p], "\n* Gamma%s = ",gammaname[i].data());
	for (ipar=0; ipar < node_num_parm[inode].size(); ++ipar) {
	  if (ipar==0) { fprintf (avefile[p], "(") ; } else {fprintf (avefile[p], " + ") ;}
	  int parm=node_num_parm[inode].at(ipar);
	  vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	  int quan=ibase-baseparm.begin()+1;
	  fprintf (avefile[p], "%f*Gamma%d",node_num_coef[inode].at(ipar), basegamma[quan-1]);
	  if (ipar==node_num_parm[inode].size()-1) fprintf (avefile[p], ")");
	}
	if (node_den_parm[inode].size()==0) fprintf (avefile[p], "\n") ; 
	for (ipar=0; ipar < node_den_parm[inode].size(); ++ipar) {
	  if (ipar==0) { fprintf (avefile[p], " / (") ; } else {fprintf (avefile[p], " + ") ;}
	  int parm=node_den_parm[inode].at(ipar);
	  vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	  int quan=ibase-baseparm.begin()+1;
	  fprintf (avefile[p], "%f*Gamma%d",node_den_coef[inode].at(ipar), basegamma[quan-1]);
	  if (ipar==node_den_parm[inode].size()-1) fprintf (avefile[p], ")\n");
	}
	//
	double offset = -node_num[inode]; // - [ f(x0,y0) - df/dx(x=x0) x0 - df/dy(y=y0) y0 - ...]
	if (node_den_npar[inode]>0) offset /= node_den[inode];
	//	  cout << inode << " " << node_num[inode] << " " << node_den[inode] << " " << offset << " " << measvalue[i] << endl;
	for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	  int parm=node_parm[inode].at(ipar);
	  int quan=node_quan[inode].at(ipar);
	  double partial=node_part[inode].at(ipar);
	  offset += partial*baseseed[quan-1]; 
	  //	  offset += partial*(basefitvalue[quan-1]);
	}
	if (p==0) { // COMBOS
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d    %2d %d\n",isum,i+1,node_parm[inode].size()); 
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_AD %f 1.0\n",isum,offset); 
	} else if (p==1) { // ALUCOMB
	  fprintf (avefile[p], "CONSTRAINT Gamma%s.c %f Gamma%s -1", gammaname[i].data(), offset, gammaname[i].data());
	}
	for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	  int parm=node_parm[inode].at(ipar);
	  int quan=node_quan[inode].at(ipar);
	  double partial=node_part[inode].at(ipar);
	  if (p==0) { // COMBOS
	    fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_%2.2d %2d %f ! Gamma%d \n",isum,ipar+1,quan,partial,basegamma.at(quan-1));
	  }else if (p==1) { // ALUCOMB
	    fprintf(avefile[p], " Gamma%d %f", basegamma[quan-1], partial);
	  }
	}
	if (p==1) fprintf(avefile[p], "\n");
      }
    }
    if (p==0) {
      for (int i=0;i<nmeas;++i){
	vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
	inode=it-nodename.begin();
	//
	if ((inode+1)==80) { // SPECIAL CASE : NODE = 79 NAME = S035R33 GAMMA = 102 :: Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104)
	  fprintf(avefile[p], "\n*Gamma102 = 1 - Gamma3   - Gamma5   - Gamma9   - Gamma10  - Gamma14  - Gamma16\n");
	  fprintf(avefile[p], "*             - Gamma20  - Gamma23  - Gamma27  - Gamma28  - Gamma30 - Gamma35\n");
	  fprintf(avefile[p], "*             - Gamma37 - Gamma40  - Gamma42  - Gamma47  - Gamma48  - Gamma62\n");
	  fprintf(avefile[p], "*             - Gamma70 - Gamma77  - Gamma78  - Gamma85  - Gamma89  - Gamma93\n");
	  fprintf(avefile[p], "*             - Gamma94 - Gamma126 - Gamma128 - Gamma150 - Gamma152\n");
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d    %2d 29 \n",++isum,i+1); 
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_AD -1 +1 \n",isum); // becomes a measurement of -1+Gamma102; thats why the coefficients below have - sign 
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_01  1 -1 ! Gamma3  \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_02  2 -1 ! Gamma5  \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_03  3 -1 ! Gamma9  \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_04  4 -1 ! Gamma10 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_05  5 -1 ! Gamma14 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_06  6 -1 ! Gamma16 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_07  7 -1 ! Gamma20 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_08  8 -1 ! Gamma23 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_09  9 -1 ! Gamma27 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_10 10 -1 ! Gamma28 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_11 11 -1 ! Gamma30 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_12 12 -1 ! Gamma35 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_13 13 -1 ! Gamma37 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_14 14 -1 ! Gamma40 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_15 15 -1 ! Gamma42 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_16 16 -1 ! Gamma47 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_17 17 -1 ! Gamma48 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_18 18 -1 ! Gamma62 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_19 19 -1 ! Gamma70 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_20 20 -1 ! Gamma77 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_21 21 -1 ! Gamma78 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma85 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma89 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_24 24 -1 ! Gamma93 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_25 25 -1 ! Gamma94 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_26 27 -1 ! Gamma126\n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_27 28 -1 ! Gamma128\n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_28 29 -1 ! Gamma150\n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_29 30 -1 ! Gamma152\n",isum);
	}
	//
	if ((inode+1)==82) { // SPECIAL CASE : NODE = 81 NAME = S035R38 GAMMA = 103 
	  fprintf(avefile[p], "\n*Gamma103 = 1 - Gamma3   - Gamma5   - Gamma9   - Gamma10  - Gamma14  - Gamma16\n");
	  fprintf(avefile[p], "*             - Gamma20  - Gamma23  - Gamma27  - Gamma28  - Gamma30 - Gamma35\n");
	  fprintf(avefile[p], "*             - Gamma37 - Gamma40  - Gamma42  - Gamma47  - Gamma48  - Gamma62\n");
	  fprintf(avefile[p], "*             - Gamma70 - Gamma77  - Gamma78  - Gamma85  - Gamma89  - Gamma93\n");
	  fprintf(avefile[p], "*             - Gamma94 - Gamma104 - Gamma126 - Gamma128 - Gamma150 - Gamma152\n");
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d    %2d 30 \n",++isum,i+1); 
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_AD -1 +1 \n",isum); // becomes a measurement of -1+Gamma103; thats why the coefficients below have - sign 
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_01  1 -1 ! Gamma3  \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_02  2 -1 ! Gamma5  \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_03  3 -1 ! Gamma9  \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_04  4 -1 ! Gamma10 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_05  5 -1 ! Gamma14 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_06  6 -1 ! Gamma16 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_07  7 -1 ! Gamma20 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_08  8 -1 ! Gamma23 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_09  9 -1 ! Gamma27 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_10 10 -1 ! Gamma28 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_11 11 -1 ! Gamma30 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_12 12 -1 ! Gamma35 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_13 13 -1 ! Gamma37 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_14 14 -1 ! Gamma40 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_15 15 -1 ! Gamma42 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_16 16 -1 ! Gamma47 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_17 17 -1 ! Gamma48 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_18 18 -1 ! Gamma62 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_19 19 -1 ! Gamma70 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_20 20 -1 ! Gamma77 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_21 21 -1 ! Gamma78 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma85 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma89 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_24 24 -1 ! Gamma93 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_25 25 -1 ! Gamma94 \n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_26 26 -1 ! Gamma104\n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_27 27 -1 ! Gamma126\n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_28 28 -1 ! Gamma128\n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_29 29 -1 ! Gamma150\n",isum);
	  fprintf (avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_30 30 -1 ! Gamma152\n",isum);
	}
	//
      }
      if (p==0) fprintf (avefile[p], "\nSPARAMETER CHI2_N_SYM_NSUM  %d 0 \n",isum); 
    }
    if (p==1) {
      fprintf(avefile[p], "\n* unitarity constraint (sum of basic modes, possibly adding also dummy)\n");
      fprintf(avefile[p], "CONSTRAINT GammaAll 1\n");
      fprintf(avefile[p], "  Gamma3   1 Gamma5   1 Gamma9   1 Gamma10  1 Gamma14  1 Gamma16  1\n");
      fprintf(avefile[p], "  Gamma20  1 Gamma23  1 Gamma27  1 Gamma28  1 Gamma30  1 Gamma35  1\n");
      fprintf(avefile[p], "  Gamma37  1 Gamma40  1 Gamma42  1 Gamma47  1 Gamma48  1 Gamma62  1\n");
      fprintf(avefile[p], "  Gamma70  1 Gamma77  1 Gamma78  1 Gamma85  1 Gamma89  1 Gamma93  1\n");
      fprintf(avefile[p], "  Gamma94  1 Gamma103 1 Gamma104 1 Gamma126 1 Gamma128 1 Gamma150 1 Gamma152 1\n");
      fprintf(avefile[p], "* Gamma998 1\n");
    }
    fprintf(avefile[p], "\nCALL CHI2_N_SYM\n");
    fprintf(avefile[p], "\nEND\n");
    fclose(avefile[p]);
  }
  //
  vector<int> basemodes_used_in_measured_derivednodes;
  for (int i=0;i<nmeas;++i){
    vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
    inode=it-nodename.begin();
    if ((node_num_npar[inode]+node_den_npar[inode])>1) {
      basemodes_used_in_measured_derivednodes.insert(basemodes_used_in_measured_derivednodes.end(),node_quan[inode].begin(),node_quan[inode].end());
      sort(basemodes_used_in_measured_derivednodes.begin(),basemodes_used_in_measured_derivednodes.end());
      vector<int>::iterator new_end=unique(basemodes_used_in_measured_derivednodes.begin(),basemodes_used_in_measured_derivednodes.end());
      basemodes_used_in_measured_derivednodes.erase(new_end,basemodes_used_in_measured_derivednodes.end());
    }
  }
  cout << "basemodes_used_in_measured_derivednodes.size() = " << basemodes_used_in_measured_derivednodes.size() << endl;
  for (ibase=0;ibase<basemodes_used_in_measured_derivednodes.size();++ibase) {
    int quan = basemodes_used_in_measured_derivednodes.at(ibase) ;
    cout << "basetitle = " << basetitle[quan-1] << " with quan = " << quan << " gamma = " << basegamma[quan-1]  << " appears in ";
    int iquan_occurance=0;
    int inode_occurance=0;
    int lastnode=-1;
    for (int i=0;i<nmeas;++i){
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      for (int ii=0;ii<node_quan[inode].size();++ii){
	if (quan == node_quan[inode].at(ii)) { 
	  ++iquan_occurance;
	  if (inode!=lastnode){
	    lastnode=inode;
	    ++inode_occurance;
	  }
	}
      }
    }
    cout << iquan_occurance << " measurements of " << inode_occurance << " nodes. \nThese " << iquan_occurance << " measurements are :" ;
    iquan_occurance = 0;
    vector<int> vector_meas_with_correlation;
    for (int i=0;i<nmeas;++i){
      bool meas_with_correlation = false;
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      for (int ii=0;ii<node_quan[inode].size();++ii){
	if (quan == node_quan[inode].at(ii)) { 
	  cout << " (" << ++iquan_occurance << ") IMEAS = " << i+1 << " NODE = " << inode << " NAME = " << nodename[inode].data() << " GAMMA = " << gammaname[i].data() ;
	  for (int j=0;j<nmeas;++j){
	    if (corrmat[i][j] != 0) {
	      meas_with_correlation = true;
	      break;
	    }
	  }
	}
      }
      if (meas_with_correlation) vector_meas_with_correlation.push_back(i+1);
    }
    cout << endl;
    if (vector_meas_with_correlation.size()!=0) {
      cout << "vector_meas_with_correlation has size = "<< vector_meas_with_correlation.size() << " : ";
      for (int i=0; i < vector_meas_with_correlation.size(); ++i) cout << vector_meas_with_correlation[i] << " " ; cout << endl;
    }else{
      cout << "vector_meas_with_correlation has size = 0"<< endl;
    }
    cout << endl;
  }
  //
  // Cross-check understanding of chisquare and pullaverage in the input file
  //
  int i,j;
  TMatrixD ErrorMatrix(nmeas,nmeas);
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) {
	ErrorMatrix(i,j)=measerror[i]*measerror[i];
      }else{
	ErrorMatrix(i,j)=corrmat[i][j]*measerror[i]*measerror[j];
      }
    }
  }
  Double_t determinant;
  TMatrixD InvErrorMatrix = ErrorMatrix;
  InvErrorMatrix.Invert(&determinant);  

  TMatrixD Unitary(InvErrorMatrix,TMatrixD::kMult,ErrorMatrix);
  TMatrixDDiag diagonal(Unitary); diagonal = 0.0;
  const Double_t Unitary_max_offdiag = (Unitary.Abs()).Max();
  //  cout << "  Maximum off-diagonal = " << Unitary_max_offdiag << endl;
  //  cout << "  Determinant          = " << determinant <<endl;
  //
  double chisquare_re[200];
  double chisquare_tot=0;
  for (i=0;i<nmeas;++i){
    chisquare_re[i]=0.0;
    for (j=0;j<nmeas;++j){
      chisquare_re[i]+=(measvalue[i]-fitvalue[i])*InvErrorMatrix[i][j]*(measvalue[j]-fitvalue[j]);
    }
    chisquare_tot+=chisquare_re[i];
  }
  cout << "chisquare_tot = " << chisquare_tot << endl;
  //
  vector<int> icorrj[200];
  for (i=0;i<nmeas;++i){
    for (j=0;j<nmeas;++j){
      if (InvErrorMatrix[i][j]!=0.0) {
	icorrj[i].push_back(j);
      }
    }
  }
  double chisquare_tmp[200];
  for (i=0; i<nmeas; ++i){
    if ( icorrj[i].size() > 1 ) {
      chisquare_tmp[i]=0;
      //      cout << "meas = i+1 = " << i+1 << " icorrj[i].size = " << icorrj[i].size() << " ; ";
      //      cout << "icorrj[i][j]+1 = " ;
      for (j=0; j < icorrj[i].size(); ++j) {
	//	cout << icorrj[i][j]+1 << " "; 
	chisquare_tmp[i]+=chisquare_re[icorrj[i][j]];
      }
      //      cout << endl;
      chisquare_tmp[i]/=icorrj[i].size();
    }
  }
  //
  for (i=0; i<nmeas; ++i){
    if ( icorrj[i].size() > 1 ) {
      chisquare_re[i]=chisquare_tmp[i];
    }
  }
  //
  int ncorrij=0;
  int ifirstj[200];
  for (i=0;i<nmeas;++i) {
    if (icorrj[i].size() > 1) {
      bool isnew=true;
      for (j=0;j<i;++j){
	if (icorrj[j].size() > 1) {
	  if (find(icorrj[j].begin(),icorrj[j].end(),i)!=icorrj[j].end()) {
	    isnew=false;
	    break;
	  }
	}
      }
      if (isnew) {
	ifirstj[ncorrij] = i;
	++ncorrij;
      }
      //      cout << "i,isnew,ncorrij = " << i << " " << isnew << " " << ncorrij << endl;
    }
  }
  cout << "ncorrij = " << ncorrij << endl;
  vector<int> veccorrij[200];
  for (i=0;i<ncorrij;++i){
    cout << "i = " << i << " ifirstj[i]+1 = " << ifirstj[i]+1 << " : ";
    veccorrij[i].insert(veccorrij[i].end(),icorrj[ifirstj[i]].begin(),icorrj[ifirstj[i]].end());
    cout << "veccorrij[i].size() : " << veccorrij[i].size() << " : veccorrij[i][j]+1 = " ;
    for (j=0;j<veccorrij[i].size();++j) cout << veccorrij[i][j]+1 << " ";
    cout << endl;
  }
  //
  int lastnode_all=-1;
  int newnode_all=-1;
  int measnode_all[200];
  double pullsq_all[200];
  for (i=0;i<nmeas;++i) {
    if (icorrj[i].size()>1) continue;
    vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
    inode=it-nodename.begin();
    if (inode!=lastnode_all) {
      lastnode_all=inode;
      ++newnode_all;
    }
    measnode_all[i] = newnode_all;
    pullsq_all[i] = ((measvalue[i]-fitvalue[i])*(measvalue[i]-fitvalue[i]))/(measerror[i]*measerror[i] - fiterror[i]*fiterror[i]);
  }
  //
  for (i=0;i<ncorrij;++i) {
    ++newnode_all;
    for (j=0;j<veccorrij[i].size();++j) {
      int meas=veccorrij[i][j];
      measnode_all[meas] = newnode_all;
      pullsq_all[meas] = ((measvalue[meas]-fitvalue[meas])*(measvalue[meas]-fitvalue[meas])) / (measerror[meas]*measerror[meas] - fiterror[meas]*fiterror[meas]);
    }
  }
  int n_pernode[200];
  for (inode=0;inode<=newnode_all;++inode){
    n_pernode[inode]=0;
    for (i=0;i<nmeas;++i){
      if (measnode_all[i]==inode){
	if (icorrj[i].size()==1) {
	  n_pernode[inode]+=1;
	}else{
	  n_pernode[inode]=1;
	}
      }
    }
  }
  //
  double nsig_fit[200];
  bool weak_re[200];
  for (i=0;i<nmeas;++i){
    if (n_pernode[measnode_all[i]] > 0 ) {
      nsig_fit[i] = measerror[i]/((sqrt(n_pernode[measnode_all[i]]))*fiterror[i]);
    }else{
      nsig_fit[i]=0;
    }
    weak_re[i] = (icorrj[i].size()==1) && (nsig_fit[i] > 3.);
  }
  //
  double pullsq_pernode[200];
  double pullav_pernode[200];
  for (inode=0;inode<=newnode_all;++inode){
    pullsq_pernode[inode]=0;
    for (i=0;i<nmeas;++i){
      if (measnode_all[i]==inode){
	if (weak_re[i]==0) {
	  pullsq_pernode[inode]+=pullsq_all[i];
	}
	//	cout << "i+1 = " << i+1 << " pullav = " << pullav[i] << " inode = " << inode << " n_pernode[inode] = " << n_pernode[inode] << " pullsq_pernode[inode] = " << pullsq_pernode[inode] << " weak_re = " << weak_re[i] << endl;
      }
    }
    pullav_pernode[inode]=(n_pernode[inode]>0) ? sqrt(pullsq_pernode[inode]/n_pernode[inode]) : 0;
    //    cout << "inode = " << inode << " n_pernode = " << n_pernode[inode] << " pullsq_pernode = " << pullsq_pernode[inode] << " pullav_pernode = " << pullav_pernode[inode] << endl;
  }
  double pullav_re[200];
  for (i=0;i<nmeas;++i){
    if (weak_re[i]==0) {
      pullav_re[i]=pullav_pernode[measnode_all[i]];
    }else{
      pullav_re[i]=0;
    }
    if (icorrj[i].size()==1 && fitvalue[i] < 3. * fiterror[i] ) pullav_re[i]*=sqrt(3./(fitvalue[i]/fiterror[i]));
    cout << "i+1 = " << i+1 << " node = " << measnode_all[i] << " corr = " << icorrj[i].size() << " n = " << n_pernode[measnode_all[i]] << " meas  = " << measvalue[i] << " +- " << measerror[i] << " fit = " << fitvalue[i] << " +- " << fiterror[i] << " pullsq = " << pullsq_all[i] << " pullav = " << pullav_re[i] << " nsig_fit = " << nsig_fit[i] << " weak_re = " << weak_re[i]  << " weak = " << measweak[i]    << endl;
  }
  //
  chisquare_tot=0;
  for (i=0; i<nmeas; ++i){
    chisquare_tot+=chisquare_re[i];
    TString spull=Form("%4.2f",pullav_re[i]);
    double dpull=atof(spull.Data());
    bool same=(dpull==pullav[i]);
    //
    //    if (!same) {
    cout << "Summary: i+1 = " << i+1 << " icorrj[i] = " << icorrj[i].size() << " chi2_re = " << chisquare_re[i] << " chi2 = " << chisquare[i] << " chi2sum = " << chisquare_tot << " pull_re = " << pullav_re[i] << " pull = " << pullav[i] << " same = " << same << " weak = " << measweak[i] <<  " weak_re = " << weak_re[i] << " fit/err = " << fitvalue[i]/fiterror[i] << " meas/err = " << measvalue[i]/measerror[i] << endl;
      //    }
  } 
  //
  FILE *scaled_measfile[2];
  for (int p=0;p<1;++p){
    if (p==0) scaled_measfile[0]=fopen("combos_pdginput_measurements_scaled.input","w");
    if (p==1) scaled_measfile[1]=fopen("alucomb_pdginput_measurements_scaled.input","w");
    iimeas=0;
    for (int i=0;i<nmeas;++i){
      if (weak_re[i]==1) continue;
      //      cout << "i = " << i << endl;
      fprintf (scaled_measfile[p], "* IMEAS = %d \n",    ++iimeas);
      fprintf (scaled_measfile[p], "* GAMMANAME = %s \n", gammaname[i].data());
      fprintf (scaled_measfile[p], "* DECAYNAME = %s \n", meastitle[i].data());
      fprintf (scaled_measfile[p], "* IWEAK = %d \n",   measweak[i]);
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      fprintf (scaled_measfile[p], "* NODENAME = %s found at inode+1 = %d with NODENAME[inode] = %s has %d, %d, %d base quantities in numerator, denominator and both [excluding overlap] \n", measnodename[i].data(), inode+1, nodename[inode].data(),node_num_npar[inode], node_den_npar[inode],node_parm[inode].size());
      for (ipar=0;ipar<node_num_parm[inode].size();++ipar){
	int parm=node_num_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	fprintf (scaled_measfile[p],"*                numerator of inode+1 = %d has gamma = %d parm = %d quan = %d title = %s seed = %f coef = %f\n",inode+1, basegamma[quan-1], parm, quan, basetitle[quan-1].data(), baseseed[quan-1], node_num_coef[inode].at(ipar));
      }
      for (ipar=0;ipar<node_den_parm[inode].size();++ipar){
	int parm=node_den_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	fprintf (scaled_measfile[p],"*              denominator of inode+1 = %d has gamma = %d parm = %d quan = %d title = %s seed = %f coef = %f\n",inode+1, basegamma[quan-1], parm, quan, basetitle[quan-1].data(), baseseed[quan-1], node_den_coef[inode].at(ipar));
      }
      fprintf (scaled_measfile[p], "*  first quantity measured by inode+1 = %d has gamma = %d parm = %d quan = %d title = %s\n",inode+1,basegamma[first_quan[inode]-1],baseparm[first_quan[inode]-1],first_quan[inode],basetitle[first_quan[inode]-1].data());
      //
      fprintf (scaled_measfile[p], "\nBEGIN %s Gamma%s published.%s \n\n", author[i].data(), gammaname[i].data(), year[i].data());
      if (p==0) {//COMBOS
	if (// SPECIAL CASE [because these nodes contain Gamma103 [used to express unitarity constraint]
	    (inode+1)==80 || // NODE = 79 NAME = S035R33 GAMMA = 102
	    (inode+1)==82) { // NODE = 81 NAME = S035R38 GAMMA = 103
	  fprintf (scaled_measfile[p], "MEASUREMENT  m_Gamma%d statistical systematic \n",3);
	  fprintf (scaled_measfile[p], "DATA         m_Gamma%d statistical systematic \n",3);
	}else{
	  fprintf (scaled_measfile[p], "MEASUREMENT  m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
	  fprintf (scaled_measfile[p], "DATA         m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
	}
      }else if (p==1) {//ALUCOMB
	fprintf (scaled_measfile[p], "MEASUREMENT  m_Gamma%s statistical systematic \n",gammaname[i].data());
	fprintf (scaled_measfile[p], "DATA         m_Gamma%s statistical systematic \n",gammaname[i].data());
      }
      if (pullav_re[i]>1&&icorrj[i].size()>1) {
	fprintf (scaled_measfile[p], "             %10.5g %10.5g  0 \n",measvalue[i],measerror[i]*(pullav_re[i]/sqrt(icorrj[i].size() )));
      }else if (pullav_re[i]>1&&icorrj[i].size()==1) {
	fprintf (scaled_measfile[p], "             %10.5g %10.5g  0 \n",measvalue[i],measerror[i]*(pullav_re[i]));
      }else{
	fprintf (scaled_measfile[p], "             %10.5g %10.5g  0 \n",measvalue[i],measerror[i]);
      }
      bool firstcorr=true;
      jjmeas=0;
      for (int j=0;j<nmeas;++j) {
	if (weak_re[j]==1) continue;
	++jjmeas;
	if (corrmat[i][j]!=0) {
	  if (firstcorr) {fprintf (scaled_measfile[p], " \n"); firstcorr=false;}
	  fprintf (scaled_measfile[p], "STAT_CORR_WITH %s Gamma%s published.%s %f ! IMEAS = %d \n",author[j].data(), gammaname[j].data(), year[j].data(), corrmat[i][j], jjmeas);
	}
      }
      fprintf (scaled_measfile[p], " \nEND \n\n");
    }
    fclose(scaled_measfile[p]);
  }
  //
  //
  FILE *scaled_avefile[2];
  for (int p=0;p<1;++p){
    if (p==0) scaled_avefile[0]=fopen("combos_average_pdginput_no_babar_belle_scaled.input","w");
    if (p==1) scaled_avefile[1]=fopen("alucomb_average_pdginput_no_babar_belle_scaled.input","w");
    if (p==0) fprintf (scaled_avefile[p], "INCLUDE combos_pdginput_measurements_scaled.input \n\n"); 
    if (p==1) fprintf (scaled_avefile[p], "INCLUDE alucomb_pdginput_measurements_scaled.input \n\n"); 
    fprintf (scaled_avefile[p], "BEGIN   PDG-BABAR-BELLE all_methods \n\n");
    fprintf (scaled_avefile[p], "COMBINE * * * \n\n");
    for (ibase=0;ibase<nbase;++ibase){
      if (p==0&&ibase==(nbase-1)){/* skip */}else{
	fprintf (scaled_avefile[p], "MEASUREMENT m_Gamma%d statistical systematic   ! NQUAN = %d \n",basegamma.at(ibase),ibase+1);  
      }
    }
    if (p==0){
      fprintf (scaled_avefile[p], "\nCALL DUMP_MASTER_INC \n\n");
      fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_PRT -1.0 0. \n\n");
      fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_INV 0 0    \n\n");
    }
    //
    int lastnode=-1;
    int usum=0;//number of [unique] measurements to be expressed as linearized sum of base quantities
    for (int i=0;i<nmeas;++i){
      if (weak_re[i]==1) continue;
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      if ((node_num_npar[inode]+node_den_npar[inode])>1) {
	if (p==1&&inode!=lastnode){//new node
	  ++usum;
	  lastnode=inode;
	  fprintf (scaled_avefile[p], "MEASUREMENT m_Gamma%s statistical systematic   ! NQUAN = %d \n",gammaname[i].data(),ibase+usum);  
	}
      }
    }
    //
    int isum=0;//number of          measurements to be expressed as linearized sum of base quantities
    lastnode=-1;
    iimeas=0;
    for (int i=0;i<nmeas;++i){
      if (weak_re[i]==1) continue;
      ++iimeas;
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      if ((node_num_npar[inode]+node_den_npar[inode])>1) {
	//
	if (p==0&&(((inode+1)==80)||((inode+1)==82))) continue; // SPECIAL CASE [because these are derived nodes containing Gamma103 ]
	++isum; // translate C index to Fortran index
	//
	if (inode!=lastnode){//new node
	  lastnode=inode;
	}else{
	  if (p==1) continue;
	}
	//
	// PRINT NODE DEFINITION
	//
	fprintf (scaled_avefile[p], "\n* Gamma%s = ",gammaname[i].data());
	for (ipar=0; ipar < node_num_parm[inode].size(); ++ipar) {
	  if (ipar==0) { fprintf (scaled_avefile[p], "(") ; } else {fprintf (scaled_avefile[p], " + ") ;}
	  int parm=node_num_parm[inode].at(ipar);
	  vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	  int quan=ibase-baseparm.begin()+1;
	  fprintf (scaled_avefile[p], "%f*Gamma%d",node_num_coef[inode].at(ipar), basegamma[quan-1]);
	  if (ipar==node_num_parm[inode].size()-1) fprintf (scaled_avefile[p], ")");
	}
	if (node_den_parm[inode].size()==0) fprintf (scaled_avefile[p], "\n") ; 
	for (ipar=0; ipar < node_den_parm[inode].size(); ++ipar) {
	  if (ipar==0) { fprintf (scaled_avefile[p], " / (") ; } else {fprintf (scaled_avefile[p], " + ") ;}
	  int parm=node_den_parm[inode].at(ipar);
	  vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	  int quan=ibase-baseparm.begin()+1;
	  fprintf (scaled_avefile[p], "%f*Gamma%d",node_den_coef[inode].at(ipar), basegamma[quan-1]);
	  if (ipar==node_den_parm[inode].size()-1) fprintf (scaled_avefile[p], ")\n");
	}
	//
	double offset = -node_num[inode]; // - [ f(x0,y0) - df/dx(x=x0) x0 - df/dy(y=y0) y0 - ...]
	if (node_den_npar[inode]>0) offset /= node_den[inode];
	//	  cout << inode << " " << node_num[inode] << " " << node_den[inode] << " " << offset << " " << measvalue[i] << endl;
	for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	  int parm=node_parm[inode].at(ipar);
	  int quan=node_quan[inode].at(ipar);
	  double partial=node_part[inode].at(ipar);
	  offset += partial*baseseed[quan-1]; 
	  //	  offset += partial*(basefitvalue[quan-1]);
	}
	if (p==0) { // COMBOS
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d    %2d %d\n",isum,iimeas,node_parm[inode].size()); 
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_AD %f 1.0\n",isum,offset); 
	} else if (p==1) { // ALUCOMB
	  fprintf (scaled_avefile[p], "CONSTRAINT Gamma%s.c %f Gamma%s -1", gammaname[i].data(), offset, gammaname[i].data());
	}
	for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	  int parm=node_parm[inode].at(ipar);
	  int quan=node_quan[inode].at(ipar);
	  double partial=node_part[inode].at(ipar);
	  if (p==0) { // COMBOS
	    fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_%2.2d %2d %f ! Gamma%d \n",isum,ipar+1,quan,partial,basegamma.at(quan-1));
	  }else if (p==1) { // ALUCOMB
	    fprintf(scaled_avefile[p], " Gamma%d %f", basegamma[quan-1], partial);
	  }
	}
	if (p==1) fprintf(scaled_avefile[p], "\n");
      }
    }
    if (p==0) {
      iimeas=0;
      for (int i=0;i<nmeas;++i){
	if (weak_re[i]==1) continue;
	++iimeas;
	vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
	inode=it-nodename.begin();
	//
	if ((inode+1)==80) { // SPECIAL CASE : NODE = 79 NAME = S035R33 GAMMA = 102 :: Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104)
	  fprintf(scaled_avefile[p], "\n*Gamma102 = 1 - Gamma3   - Gamma5   - Gamma9   - Gamma10  - Gamma14  - Gamma16\n");
	  fprintf(scaled_avefile[p], "*             - Gamma20  - Gamma23  - Gamma27  - Gamma28  - Gamma30 - Gamma35\n");
	  fprintf(scaled_avefile[p], "*             - Gamma37 - Gamma40  - Gamma42  - Gamma47  - Gamma48  - Gamma62\n");
	  fprintf(scaled_avefile[p], "*             - Gamma70 - Gamma77  - Gamma78  - Gamma85  - Gamma89  - Gamma93\n");
	  fprintf(scaled_avefile[p], "*             - Gamma94 - Gamma126 - Gamma128 - Gamma150 - Gamma152\n");
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d    %2d 29 \n",++isum,iimeas); 
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_AD -1 +1 \n",isum); // becomes a measurement of -1+Gamma102; thats why the coefficients below have - sign 
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_01  1 -1 ! Gamma3  \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_02  2 -1 ! Gamma5  \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_03  3 -1 ! Gamma9  \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_04  4 -1 ! Gamma10 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_05  5 -1 ! Gamma14 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_06  6 -1 ! Gamma16 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_07  7 -1 ! Gamma20 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_08  8 -1 ! Gamma23 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_09  9 -1 ! Gamma27 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_10 10 -1 ! Gamma28 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_11 11 -1 ! Gamma30 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_12 12 -1 ! Gamma35 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_13 13 -1 ! Gamma37 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_14 14 -1 ! Gamma40 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_15 15 -1 ! Gamma42 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_16 16 -1 ! Gamma47 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_17 17 -1 ! Gamma48 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_18 18 -1 ! Gamma62 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_19 19 -1 ! Gamma70 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_20 20 -1 ! Gamma77 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_21 21 -1 ! Gamma78 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma85 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma89 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_24 24 -1 ! Gamma93 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_25 25 -1 ! Gamma94 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_26 27 -1 ! Gamma126\n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_27 28 -1 ! Gamma128\n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_28 29 -1 ! Gamma150\n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_29 30 -1 ! Gamma152\n",isum);
	}
	//
	if ((inode+1)==82) { // SPECIAL CASE : NODE = 81 NAME = S035R38 GAMMA = 103 
	  fprintf(scaled_avefile[p], "\n*Gamma103 = 1 - Gamma3   - Gamma5   - Gamma9   - Gamma10  - Gamma14  - Gamma16\n");
	  fprintf(scaled_avefile[p], "*             - Gamma20  - Gamma23  - Gamma27  - Gamma28  - Gamma30 - Gamma35\n");
	  fprintf(scaled_avefile[p], "*             - Gamma37 - Gamma40  - Gamma42  - Gamma47  - Gamma48  - Gamma62\n");
	  fprintf(scaled_avefile[p], "*             - Gamma70 - Gamma77  - Gamma78  - Gamma85  - Gamma89  - Gamma93\n");
	  fprintf(scaled_avefile[p], "*             - Gamma94 - Gamma104 - Gamma126 - Gamma128 - Gamma150 - Gamma152\n");
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d    %2d 30 \n",++isum,iimeas); 
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_AD -1 +1 \n",isum); // becomes a measurement of -1+Gamma103; thats why the coefficients below have - sign 
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_01  1 -1 ! Gamma3  \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_02  2 -1 ! Gamma5  \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_03  3 -1 ! Gamma9  \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_04  4 -1 ! Gamma10 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_05  5 -1 ! Gamma14 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_06  6 -1 ! Gamma16 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_07  7 -1 ! Gamma20 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_08  8 -1 ! Gamma23 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_09  9 -1 ! Gamma27 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_10 10 -1 ! Gamma28 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_11 11 -1 ! Gamma30 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_12 12 -1 ! Gamma35 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_13 13 -1 ! Gamma37 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_14 14 -1 ! Gamma40 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_15 15 -1 ! Gamma42 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_16 16 -1 ! Gamma47 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_17 17 -1 ! Gamma48 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_18 18 -1 ! Gamma62 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_19 19 -1 ! Gamma70 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_20 20 -1 ! Gamma77 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_21 21 -1 ! Gamma78 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma85 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma89 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_24 24 -1 ! Gamma93 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_25 25 -1 ! Gamma94 \n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_26 26 -1 ! Gamma104\n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_27 27 -1 ! Gamma126\n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_28 28 -1 ! Gamma128\n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_29 29 -1 ! Gamma150\n",isum);
	  fprintf (scaled_avefile[p], "SPARAMETER CHI2_N_SYM_%2.2d_30 30 -1 ! Gamma152\n",isum);
	}
	//
      }
      if (p==0) fprintf (scaled_avefile[p], "\nSPARAMETER CHI2_N_SYM_NSUM  %d 0 \n",isum); 
    }
    if (p==1) {
      fprintf(scaled_avefile[p], "\n* unitarity constraint (sum of basic modes, possibly adding also dummy)\n");
      fprintf(scaled_avefile[p], "CONSTRAINT GammaAll 1\n");
      fprintf(scaled_avefile[p], "  Gamma3   1 Gamma5   1 Gamma9   1 Gamma10  1 Gamma14  1 Gamma16  1\n");
      fprintf(scaled_avefile[p], "  Gamma20  1 Gamma23  1 Gamma27  1 Gamma28  1 Gamma30  1 Gamma35  1\n");
      fprintf(scaled_avefile[p], "  Gamma37  1 Gamma40  1 Gamma42  1 Gamma47  1 Gamma48  1 Gamma62  1\n");
      fprintf(scaled_avefile[p], "  Gamma70  1 Gamma77  1 Gamma78  1 Gamma85  1 Gamma89  1 Gamma93  1\n");
      fprintf(scaled_avefile[p], "  Gamma94  1 Gamma103 1 Gamma104 1 Gamma126 1 Gamma128 1 Gamma150 1 Gamma152 1\n");
      fprintf(scaled_avefile[p], "* Gamma998 1\n");
    }
    fprintf(scaled_avefile[p], "\nCALL CHI2_N_SYM\n");
    fprintf(scaled_avefile[p], "\nEND\n");
    fclose(scaled_avefile[p]);
  }
  //
}
