#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
// ----------------------------------------------------------------------
int main() {
  int ipar;
  // from s035-fit-no-babar-belle.fit
  const int nbase=31;//33;
  vector<int> basegamma;
  vector<int> baseparm;
  double    baseseed[nbase];
  string    basetitle[nbase];
  int ibase = 0;
  // INPUT PARAMETERS
  // GAMMA PARMETER SEED TITLE //BASE = NQUAN (in COMBOS)
  // ----- -------- ---- ----- 
  basegamma.push_back(  3); baseparm.push_back(  1) ; baseseed[ibase] = 1.740000E+01 ; basetitle[ibase] = "tau- --> mu- nubar(mu) nu(tau)"; ++ibase;//1
  basegamma.push_back(  5); baseparm.push_back(  2) ; baseseed[ibase] = 1.780000E+01 ; basetitle[ibase] = "tau- --> e- nubar(e) nu(tau)"; ++ibase;//2
  basegamma.push_back(  9); baseparm.push_back( 12) ; baseseed[ibase] = 1.140000E+01 ; basetitle[ibase] = "tau- --> pi- nu(tau)"; ++ibase;//3
  basegamma.push_back( 10); baseparm.push_back(  7) ; baseseed[ibase] = 7.000000E-01 ; basetitle[ibase] = "tau- --> K- nu(tau)"; ++ibase;//4
  basegamma.push_back( 14); baseparm.push_back( 16) ; baseseed[ibase] = 2.500000E+01 ; basetitle[ibase] = "tau- --> pi- pi0 nu(tau)"; ++ibase;//5
  basegamma.push_back( 16); baseparm.push_back(182) ; baseseed[ibase] = 4.500000E-01 ; basetitle[ibase] = "tau- --> K- pi0 nu(tau)"; ++ibase;//6
  basegamma.push_back( 20); baseparm.push_back(201) ; baseseed[ibase] = 9.000000E+00 ; basetitle[ibase] = "tau- --> pi- 2pi0 nu(tau) (ex.K0)"; ++ibase;//7
  basegamma.push_back( 23); baseparm.push_back(115) ; baseseed[ibase] = 6.000000E-02 ; basetitle[ibase] = "tau- --> K- 2pi0 nu(tau) (ex.K0)"; ++ibase;//8
  basegamma.push_back( 27); baseparm.push_back(203) ; baseseed[ibase] = 1.000000E+00 ; basetitle[ibase] = "tau- --> pi- 3pi0 nu(tau) (ex.K0)"; ++ibase;//9
  basegamma.push_back( 28); baseparm.push_back(116) ; baseseed[ibase] = 4.000000E-01 ; basetitle[ibase] = "tau- --> K- 3pi0 nu(tau) (ex.K0, eta)"; ++ibase;//10
  basegamma.push_back( 30); baseparm.push_back(110) ; baseseed[ibase] = 1.000000E-01 ; basetitle[ibase] = "tau- --> h- 4pi0 nu(tau) (ex.K0,eta)"; ++ibase;//11
  basegamma.push_back( 35); baseparm.push_back(117) ; baseseed[ibase] = 9.600000E-01 ; basetitle[ibase] = "tau- --> pi- Kbar0 nu(tau)"; ++ibase;//12
  basegamma.push_back( 37); baseparm.push_back( 62) ; baseseed[ibase] = 1.600000E-01 ; basetitle[ibase] = "tau- --> K- K0 nu(tau)"; ++ibase;//13
  basegamma.push_back( 40); baseparm.push_back(118) ; baseseed[ibase] = 4.000000E-01 ; basetitle[ibase] = "tau- --> pi- Kbar0 pi0 nu(tau)"; ++ibase;//14
  basegamma.push_back( 42); baseparm.push_back(119) ; baseseed[ibase] = 2.000000E-01 ; basetitle[ibase] = "tau- --> K- K0 pi0 nu(tau)"; ++ibase;//15
  basegamma.push_back( 47); baseparm.push_back(214) ; baseseed[ibase] = 1.000000E-01 ; basetitle[ibase] = "tau- --> pi- K(S)0 K(S)0 nu(tau)"; ++ibase;//16
  basegamma.push_back( 48); baseparm.push_back(240) ; baseseed[ibase] = 1.000000E-01 ; basetitle[ibase] = "tau- --> pi- K(S)0 K(L)0 nu(tau)"; ++ibase;//17
  basegamma.push_back( 62); baseparm.push_back(259) ; baseseed[ibase] = 9.100000E+00 ; basetitle[ibase] = "tau- --> pi- pi+ pi- nu(tau) (ex.K0,omega)"; ++ibase;//18
  basegamma.push_back( 70); baseparm.push_back(263) ; baseseed[ibase] = 2.500000E+00 ; basetitle[ibase] = "tau- --> pi- pi+ pi- pi0 nu(tau) (ex.K0,omega)"; ++ibase;//19
  basegamma.push_back( 77); baseparm.push_back(216) ; baseseed[ibase] = 1.000000E-01 ; basetitle[ibase] = "tau- --> h- h- h+ 2pi0 nu(tau) (ex.K0,omega,eta)"; ++ibase;//20
  basegamma.push_back( 78); baseparm.push_back(204) ; baseseed[ibase] = 1.300000E-01 ; basetitle[ibase] = "tau- --> h- h- h+ 3pi0 nu(tau)"; ++ibase;//21
  basegamma.push_back( 85); baseparm.push_back(260) ; baseseed[ibase] = 3.000000E-01 ; basetitle[ibase] = "tau- --> K- pi+ pi- nu(tau) (ex.K0)"; ++ibase;//22
  basegamma.push_back( 89); baseparm.push_back(285) ; baseseed[ibase] = 6.000000E-02 ; basetitle[ibase] = "tau- --> K- pi+ pi- pi0 nu(tau) (ex.K0,eta)"; ++ibase;//23
  basegamma.push_back( 93); baseparm.push_back(  5) ; baseseed[ibase] = 1.600000E-01 ; basetitle[ibase] = "tau- --> K- K+ pi- nu(tau)"; ++ibase;//24
  basegamma.push_back( 94); baseparm.push_back(247) ; baseseed[ibase] = 4.000000E-02 ; basetitle[ibase] = "tau- --> K- K+ pi- pi0 nu(tau)"; ++ibase;//25
  basegamma.push_back(103); baseparm.push_back(  3) ; baseseed[ibase] = 8.000000E-02 ; basetitle[ibase] = "tau- --> 3h- 2h+ nu(tau) (ex.K0)"; ++ibase;//26
  basegamma.push_back(104); baseparm.push_back(  4) ; baseseed[ibase] = 2.000000E-02 ; basetitle[ibase] = "tau- --> 3h- 2h+ pi0 nu(tau) (ex.K0)"; ++ibase;//27
  basegamma.push_back(126); baseparm.push_back( 58) ; baseseed[ibase] = 1.700000E-01 ; basetitle[ibase] = "tau- --> eta pi- pi0 nu(tau)"; ++ibase;//28
  basegamma.push_back(128); baseparm.push_back(109) ; baseseed[ibase] = 3.000000E-02 ; basetitle[ibase] = "tau- --> eta K- nu(tau)"; ++ibase;//29
  basegamma.push_back(150); baseparm.push_back(  8) ; baseseed[ibase] = 2.000000E+00 ; basetitle[ibase] = "tau- --> h- omega nu(tau)"; ++ibase;//30
  basegamma.push_back(152); baseparm.push_back(113) ; baseseed[ibase] = 4.000000E-01 ; basetitle[ibase] = "tau- --> h- omega pi0 nu(tau)"; ++ibase;//31
  // FOR CROSS-CHECKING PDG FIT THE FOLLOWING SHOULD NOT BE USED
  //  basegamma.push_back(130); baseparm.push_back(266) ; baseseed[ibase] = 4.600000E-03 ; basetitle[ibase] = "tau- --> eta K- pi0 nu(tau)"; ++ibase;//32
  //  basegamma.push_back(132); baseparm.push_back(267) ; baseseed[ibase] = 8.800000E-03 ; basetitle[ibase] = "tau- --> eta Kbar0 pi- nu(tau)"; ++ibase;//33

  ibase=0;
  double    basefitvalue[nbase],basefiterror[nbase],baserescalederror[nbase],basescalefactor[nbase];
  // RESULTS FOR PARAMETERS
  basefitvalue[ibase]=0.173556 ; basefiterror[ibase]=0.000463 ; baserescalederror[ibase]=0.000465 ; basescalefactor[ibase]=1.01 ; ++ibase; 
  basefitvalue[ibase]=0.178413 ; basefiterror[ibase]=0.000479 ; baserescalederror[ibase]=0.000483 ; basescalefactor[ibase]=1.01 ; ++ibase; 
  basefitvalue[ibase]=0.109002 ; basefiterror[ibase]=0.000617 ; baserescalederror[ibase]=0.000660 ; basescalefactor[ibase]=1.07 ; ++ibase; 
  basefitvalue[ibase]=0.006910 ; basefiterror[ibase]=0.000219 ; baserescalederror[ibase]=0.000228 ; basescalefactor[ibase]=1.04 ; ++ibase; 
  basefitvalue[ibase]=0.254980 ; basefiterror[ibase]=0.000924 ; baserescalederror[ibase]=0.001023 ; basescalefactor[ibase]=1.11 ; ++ibase; 
  basefitvalue[ibase]=0.004525 ; basefiterror[ibase]=0.000265 ; baserescalederror[ibase]=0.000267 ; basescalefactor[ibase]=1.01 ; ++ibase; 
  basefitvalue[ibase]=0.092494 ; basefiterror[ibase]=0.000975 ; baserescalederror[ibase]=0.001240 ; basescalefactor[ibase]=1.27 ; ++ibase; 
  basefitvalue[ibase]=0.000581 ; basefiterror[ibase]=0.000225 ; baserescalederror[ibase]=0.000232 ; basescalefactor[ibase]=1.03 ; ++ibase; 
  basefitvalue[ibase]=0.010403 ; basefiterror[ibase]=0.000709 ; baserescalederror[ibase]=0.000760 ; basescalefactor[ibase]=1.07 ; ++ibase; 
  basefitvalue[ibase]=0.000417 ; basefiterror[ibase]=0.000219 ; baserescalederror[ibase]=0.000220 ; basescalefactor[ibase]=1.01 ; ++ibase; 
  basefitvalue[ibase]=0.001024 ; basefiterror[ibase]=0.000392 ; baserescalederror[ibase]=0.000398 ; basescalefactor[ibase]=1.01 ; ++ibase; 
  basefitvalue[ibase]=0.008963 ; basefiterror[ibase]=0.000367 ; baserescalederror[ibase]=0.000409 ; basescalefactor[ibase]=1.11 ; ++ibase; 
  basefitvalue[ibase]=0.001530 ; basefiterror[ibase]=0.000161 ; baserescalederror[ibase]=0.000163 ; basescalefactor[ibase]=1.01 ; ++ibase; 
  basefitvalue[ibase]=0.003778 ; basefiterror[ibase]=0.000368 ; baserescalederror[ibase]=0.000374 ; basescalefactor[ibase]=1.02 ; ++ibase; 
  basefitvalue[ibase]=0.001541 ; basefiterror[ibase]=0.000200 ; baserescalederror[ibase]=0.000200 ; basescalefactor[ibase]=1.00 ; ++ibase; 
  basefitvalue[ibase]=0.000241 ; basefiterror[ibase]=0.000052 ; baserescalederror[ibase]=0.000052 ; basescalefactor[ibase]=1.00 ; ++ibase; 
  basefitvalue[ibase]=0.001121 ; basefiterror[ibase]=0.000245 ; baserescalederror[ibase]=0.000301 ; basescalefactor[ibase]=1.23 ; ++ibase; 
  basefitvalue[ibase]=0.089866 ; basefiterror[ibase]=0.000600 ; baserescalederror[ibase]=0.000757 ; basescalefactor[ibase]=1.26 ; ++ibase; 
  basefitvalue[ibase]=0.026914 ; basefiterror[ibase]=0.000707 ; baserescalederror[ibase]=0.000826 ; basescalefactor[ibase]=1.17 ; ++ibase; 
  basefitvalue[ibase]=0.000906 ; basefiterror[ibase]=0.000357 ; baserescalederror[ibase]=0.000368 ; basescalefactor[ibase]=1.03 ; ++ibase; 
  basefitvalue[ibase]=0.000223 ; basefiterror[ibase]=0.000050 ; baserescalederror[ibase]=0.000050 ; basescalefactor[ibase]=1.00 ; ++ibase; 
  basefitvalue[ibase]=0.003335 ; basefiterror[ibase]=0.000223 ; baserescalederror[ibase]=0.000352 ; basescalefactor[ibase]=1.58 ; ++ibase; 
  basefitvalue[ibase]=0.000730 ; basefiterror[ibase]=0.000117 ; baserescalederror[ibase]=0.000122 ; basescalefactor[ibase]=1.04 ; ++ibase; 
  basefitvalue[ibase]=0.001531 ; basefiterror[ibase]=0.000070 ; baserescalederror[ibase]=0.000097 ; basescalefactor[ibase]=1.38 ; ++ibase; 
  basefitvalue[ibase]=0.000061 ; basefiterror[ibase]=0.000018 ; baserescalederror[ibase]=0.000020 ; basescalefactor[ibase]=1.10 ; ++ibase; 
  basefitvalue[ibase]=0.000810 ; basefiterror[ibase]=0.000053 ; baserescalederror[ibase]=0.000055 ; basescalefactor[ibase]=1.05 ; ++ibase; 
  basefitvalue[ibase]=0.000181 ; basefiterror[ibase]=0.000026 ; baserescalederror[ibase]=0.000026 ; basescalefactor[ibase]=1.01 ; ++ibase; 
  basefitvalue[ibase]=0.001774 ; basefiterror[ibase]=0.000235 ; baserescalederror[ibase]=0.000236 ; basescalefactor[ibase]=1.00 ; ++ibase; 
  basefitvalue[ibase]=0.000268 ; basefiterror[ibase]=0.000063 ; baserescalederror[ibase]=0.000063 ; basescalefactor[ibase]=1.00 ; ++ibase; 
  basefitvalue[ibase]=0.019860 ; basefiterror[ibase]=0.000638 ; baserescalederror[ibase]=0.000788 ; basescalefactor[ibase]=1.24 ; ++ibase; 
  basefitvalue[ibase]=0.004063 ; basefiterror[ibase]=0.000418 ; baserescalederror[ibase]=0.000429 ; basescalefactor[ibase]=1.03 ; ++ibase; 
  //
  //
  double totalseed = 0;
  for (ibase=0;ibase<nbase;++ibase) {
    baseseed[ibase]*=1.e-2; // percent to fraction
    totalseed+=baseseed[ibase];
  }
  cout << "Totalseed = " << totalseed << " is re-normalized to 1 by adjusting each baseseed value ... " << endl << endl;
  for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]/=totalseed;
  //
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
	numerator+=(node_num_coef[inode].at(ipar))*(baseseed[quan-1]); //*(basefitvalue[quan-1]);
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
	denominator+=(node_den_coef[inode].at(ipar))*(baseseed[quan-1]); //*(basefitvalue[quan-1]);
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
    if (firstch=='#'||firstch=='\n') { // Skip this line
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
      fprintf (measfile[p], "\nBEGIN %s Gamma%s published %s \n\n", author[i].data(), gammaname[i].data(), year[i].data());
      if (p==0) {//COMBOS
	fprintf (measfile[p], "MEASUREMENT  m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
	fprintf (measfile[p], "DATA         m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
      }else if (p==1) {//ALUCOMB
	fprintf (measfile[p], "MEASUREMENT  m_Gamma%s statistical systematic \n",gammaname[i].data());
	fprintf (measfile[p], "DATA         m_Gamma%s statistical systematic \n",gammaname[i].data());
      }
      fprintf (measfile[p], "             %10.5g %10.5g  0 \n",measvalue[i],measerror[i]);
      bool firstcorr=true;
      for (int j=0;j<nmeas;++j) {
	if (corrmat[i][j]!=0) {
	  if (firstcorr) {fprintf (measfile[p], " \n"); firstcorr=false;}
	  fprintf (measfile[p], "STAT_CORR_WITH %s Gamma%s published %f ! IMEAS = %d \n",author[j].data(), gammaname[j].data(),corrmat[i][j],imeas[j]);
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
    fprintf (avefile[p], "INCLUDE pdginput_measurements.input \n\n"); 
    fprintf (avefile[p], "BEGIN   PDG-BABAR-BELLE all_methods \n\n");
    fprintf (avefile[p], "COMBINE * * * \n\n");
    for (ibase=0;ibase<nbase;++ibase){
      fprintf (avefile[p], "MEASUREMENT m_Gamma%d statistical systematic   ! NQUAN = %d \n",basegamma.at(ibase),ibase+1);  
    }
    if (p==0){
      fprintf (avefile[p], "\n*CALL DUMP_MASTER_INC \n\n");
      fprintf (avefile[p], "PARAMETER CHI2_N_SYM_PRT -1.0 0. 0. \n\n");
      fprintf (avefile[p], "PARAMETER CHI2_N_SYM_INV 0 0 0 \n\n");
    }
    //
    int lastnode=-1;
    int isum=0;//number of          measurements to be expressed as linearized sum of base quantities
    int usum=0;//number of [unique] measurements to be expressed as linearized sum of base quantities
    for (int i=0;i<nmeas;++i){
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      if ((node_num_npar[inode]+node_den_npar[inode])>1) {
	++isum;
	if (p==1&&inode!=lastnode){//new node
	  ++usum;
	  lastnode=inode;
	  fprintf (avefile[p], "MEASUREMENT m_Gamma%s statistical systematic   ! NQUAN = %d \n",gammaname[i].data(),ibase+usum);  
	}
      }
    }
    if (p==0) fprintf (avefile[p], "PARAMETER CHI2_N_SYM_NSUM  %d 0 0 \n",isum); 
    //
    isum=0;      
    lastnode=-1;
    for (int i=0;i<nmeas;++i){
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[i]);      
      inode=it-nodename.begin();
      if ((node_num_npar[inode]+node_den_npar[inode])>1) {
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
	if (p==0) { // COMBOS
	  fprintf (avefile[p], "PARAMETER CHI2_N_SYM_%2.2d    %2d %d -1 \n",isum,i+1,node_parm[inode].size()); 
	} else if (p==1) { // ALUCOMB
	  double offset = -node_num[inode];
	  if (node_den_npar[inode]>0) offset /= node_den[inode];
	  //	  cout << inode << " " << node_num[inode] << " " << node_den[inode] << " " << offset << " " << measvalue[i] << endl;
	  for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	    int parm=node_parm[inode].at(ipar);
	    int quan=node_quan[inode].at(ipar);
	    double partial=node_part[inode].at(ipar);
	    offset += partial*baseseed[quan-1];
	  }
	  fprintf (avefile[p], "CONSTRAINT Gamma%s.c %f Gamma%s -1", gammaname[i].data(), offset, gammaname[i].data());
	}
	for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	  int parm=node_parm[inode].at(ipar);
	  int quan=node_quan[inode].at(ipar);
	  double partial=node_part[inode].at(ipar);
	  if (p==0) { // COMBOS
	    if (partial>0) {
	      fprintf (avefile[p], "PARAMETER CHI2_N_SYM_%2.2d_%2.2d %2d %f 0 \n",isum,ipar+1,quan,partial);
	    }else{
	      fprintf (avefile[p], "PARAMETER CHI2_N_SYM_%2.2d_%2.2d %2d 0 %f \n",isum,ipar+1,quan,partial);
	    }
	  }else if (p==1) { // ALUCOMB
	    fprintf(avefile[p], " Gamma%d %f", basegamma[quan-1], partial);
	  }
	}
	if (p==1) fprintf(avefile[p], "\n");
      }
    }
    fclose(avefile[p]);
  }
}
