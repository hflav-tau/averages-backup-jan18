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
//
using namespace std;
// ----------------------------------------------------------------------
void print_node_def(int nnode, char** a_nodename, string* nodetitle,  vector<string> nodegammaname,
		    int* node_num_npar, vector<int> * node_num_parm, vector<double> * node_num_coef,
		    int* node_den_npar, vector<int> * node_den_parm, vector<double> * node_den_coef, 
		    vector<int> baseparm, vector<int> basegamma){
  //
  int p, inode, ipar;
  FILE *nodefile[2];
  nodefile[0]=fopen("all_node_def.txt","w");
  nodefile[1]=fopen("derived_node_def.txt","w");
  for (p=0;p<2;++p) {
    for (inode=0;inode<nnode;++inode) {
      if (p==1&&(node_num_npar[inode]+node_den_npar[inode])==1) continue;
      fprintf (nodefile[p], "\n* %s : %s \n%s = ",a_nodename[inode],nodetitle[inode].data(), nodegammaname[inode].data());
      for (ipar=0; ipar < node_num_parm[inode].size(); ++ipar) {
        if (ipar==0) { fprintf (nodefile[p], "(") ; } else {fprintf (nodefile[p], " + ") ;}
        int parm=node_num_parm[inode].at(ipar);
        vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
        int quan=ibase-baseparm.begin()+1;
        fprintf (nodefile[p], "%20.10f*Gamma%d",node_num_coef[inode].at(ipar), basegamma[quan-1]);
        if (ipar==node_num_parm[inode].size()-1) fprintf (nodefile[p], ")");
      }
      if (node_den_parm[inode].size()==0) fprintf (nodefile[p], "\n") ; 
      for (ipar=0; ipar < node_den_parm[inode].size(); ++ipar) {
        if (ipar==0) { fprintf (nodefile[p], " / (") ; } else {fprintf (nodefile[p], " + ") ;}
        int parm=node_den_parm[inode].at(ipar);
        vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
        int quan=ibase-baseparm.begin()+1;
        fprintf (nodefile[p], "%20.10f*Gamma%d",node_den_coef[inode].at(ipar), basegamma[quan-1]);
        if (ipar==node_den_parm[inode].size()-1) fprintf (nodefile[p], ")\n");
      }
    }
    fclose(nodefile[p]);
  }
}
// ----------------------------------------------------------------------
void define_nodes(int nnode, vector<int> * node_num_parm, vector<int> * node_den_parm, vector<int> baseparm,
		  vector<int> * node_parm, // output
		  vector<int> * node_quan, // output
		  int * first_quan // output
		  ){
  int inode, ipar;
  for (inode=0;inode<nnode;++inode){
    node_parm[inode].insert(node_parm[inode].end(),node_num_parm[inode].begin(),node_num_parm[inode].end());
    node_parm[inode].insert(node_parm[inode].end(),node_den_parm[inode].begin(),node_den_parm[inode].end());
    sort(node_parm[inode].begin(), node_parm[inode].end());//sort must be called before unique
    vector<int>::iterator new_end=unique(node_parm[inode].begin(), node_parm[inode].end());
    node_parm[inode].erase(new_end, node_parm[inode].end()); // <--
    first_quan[inode]=99; // <--
    for (ipar=0;ipar<node_parm[inode].size();++ipar){
      int parm=node_parm[inode].at(ipar);
      vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
      int quan=ibase-baseparm.begin()+1;
      node_quan[inode].push_back(quan); // <--
      if (quan<first_quan[inode]) first_quan[inode]=quan; // <--
    }
  }
}
// ----------------------------------------------------------------------
void get_num_den_part(int nnode, vector<int> * node_parm,
		      int* node_num_npar, vector<int> * node_num_parm, vector<double> * node_num_coef,
		      int* node_den_npar, vector<int> * node_den_parm, vector<double> * node_den_coef,
		      vector<int> baseparm, double* baseseed,
		      double * node_num, // output
		      double * node_den, // output
		      double * node_val, // output
		      vector<double> * node_part// output
		      ){
  int ipar;
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
  for (int inode=0;inode<nnode;++inode){
    //
    double numerator=0; // sum_j=1,n alpha_j P_j
    if (node_num_npar[inode]>0) {
      for (ipar=0;ipar<node_num_npar[inode];++ipar){
	int parm=node_num_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	numerator+=(node_num_coef[inode].at(ipar))*(baseseed[quan-1]); 
      }
    }
    node_num[inode]=numerator; // <--
    //
    double denominator=0; // sum_k=1,m beta_k P_k
    if (node_den_npar[inode]>0) {
      for (ipar=0;ipar<node_den_npar[inode];++ipar){
	int parm=node_den_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	denominator+=(node_den_coef[inode].at(ipar))*(baseseed[quan-1]); 
      }
    }
    node_den[inode]=denominator; // <--
    //
    node_val[inode]=(denominator!=0) ? numerator/denominator : numerator ; // <--
    //
    for (ipar=0;ipar<node_parm[inode].size();++ipar){
      int parm=node_parm[inode].at(ipar);
      //
      vector<int>::iterator it_num=find(node_num_parm[inode].begin(),node_num_parm[inode].end(),parm); bool is_in_num = it_num != node_num_parm[inode].end();
      vector<int>::iterator it_den=find(node_den_parm[inode].begin(),node_den_parm[inode].end(),parm); bool is_in_den = it_den != node_den_parm[inode].end();
      double partial=0;
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
      //      if (inode==14) cout << " inode = " << inode << " ipar = " << ipar << " parm = " << parm 
      //			  << " numerator = " << numerator << " denominator = " << denominator << " partial = " << partial 
      //      			  << " it_num - node_num_parm[inode].begin() = " << it_num - node_num_parm[inode].begin() 
      //      			  << endl;
    }
  }
}
// ----------------------------------------------------------------------
void combine(int uconstrain, 
	     int nmeas, int* weak, int weakcompare, int& nmeas_noweak,
	     int* measnode, double* measvalue, TMatrixD InvMeasErrorMatrix,
	     int* node_num_npar, int* node_den_npar, vector<int> * node_quan, vector<int> * node_parm, 
	     int nbase, vector<int> basegamma,  
	     double* baseseed, double* node_num, double* node_den, vector<double> * node_part,
	     double* basevalue_fit, // output
	     double* baseerror_fit, // output
	     double** basecov_fit,  // output
	     double** basecorr_fit, // output
	     double&_chisquared     // output
	     ){
  //
  // measuredquantity = f(wxyz) = f(0) + f'w (w-w0) + f'x (x-x0) + f'y (y-y0) + f'z (z-z0)
  //                            = [f(0) - f'w w0 - f'x x0 - f'y y0 - f'z z0] + f'w w + f'x x + f'y y + f'z z

  // [measuredvalue - measuredquantity] = [measuredvalue - (f(0) - f'w w0 - f'x x0 - f'y y0 - f'z z0) - (f'w w + f'x x + f'y y + f'z z) ] 
  //                                    = Xvector + Delta^T Vvector
  // Xvector = vector [dimension = nmeas] of adjusted measured value. In this case, adjust measurement = M to M - f(0) + f'w w0 + f'x x0 + f'y y0 + f'z z0
  // Delta^T = matrix of dimension (nmeas,nbase). In this case it has dimension 1x4 = (-f'w -f'x -f'y -f'z)
  // Vvector = vector [dimension = nbase] of base quantity. In this case Vvector^T is 1x4 = (w,x,y,z)
  //
  int i,j,imeas,jmeas,inode,ibase,jbase,ipar;
  //
  nmeas_noweak=0; 
  for (imeas=0;imeas<nmeas;++imeas) {
    if (weak[imeas]==weakcompare) continue; 
    ++nmeas_noweak;
  }
  //
  TMatrixD InvMeasErrorMatrix_noweak(nmeas_noweak,nmeas_noweak);
  i=-1;
  for (imeas=0;imeas<nmeas;++imeas) {
    if (weak[imeas]==weakcompare) continue; 
    ++i;
    j=-1;
    for (jmeas=0;jmeas<nmeas;++jmeas) {
      if (weak[jmeas]==weakcompare) continue; 
      ++j;
      InvMeasErrorMatrix_noweak[i][j]=InvMeasErrorMatrix[imeas][jmeas];
    }
  }
  //
  int nbase_u = (uconstrain) ? nbase-1 : nbase;
  TMatrixD Xvector(nmeas_noweak,1);
  TMatrixD Delta(nbase_u,nmeas_noweak);
  TMatrixD DeltaT(nmeas_noweak,nbase_u);
  TMatrixD SDInvMatrix(nbase_u,nbase_u);
  TMatrixD SDMatrix(nbase_u,nbase_u);
  Double_t sd_det;
  TMatrixD Vvector(nbase_u,1);
  TMatrixD XPlusDTV(nmeas_noweak,1);
  TMatrixD XPlusDTVT(1,nmeas_noweak);
  TMatrixD ChiSquared(1,1);
  //
  i=-1;
  for (imeas=0;imeas<nmeas;++imeas){
    if (weak[imeas]==weakcompare) continue;
    ++i;
    //
    for (ibase=0;ibase<nbase_u;++ibase) Delta[ibase][i] = 0;
    //
    inode=measnode[imeas];
    if (uconstrain &&      // SPECIAL CASE [because these nodes contain Gamma103 [used to express unitarity constraint]
	((inode+1)==80 ||  // NODE = 79 NAME = S035R33 GAMMA = 102, because Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104) 
	 (inode+1)==82)) { // NODE = 81 NAME = S035R38 GAMMA = 103
      Xvector[i][0] = measvalue[imeas] - 1 ;
      for (ibase=0;ibase<nbase_u;++ibase) {
	if ((inode+1)==80 && basegamma[ibase]==104) {}else{Delta[ibase][i] = 1;}
      }
    }else{
      Xvector[i][0] = measvalue[imeas];
      double offset = -node_num[inode]; if (node_den_npar[inode]>0) offset /= node_den[inode];
      if ((node_num_npar[inode]+node_den_npar[inode])>1 || inode==102) { // derived node [inode = 102 is special case : Gamma96 = (1/1.699387) * Gamma801]
	Xvector[i][0]+=offset;
      }
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int quan=node_quan[inode].at(ipar);
	double partial=node_part[inode].at(ipar);
	if ((node_num_npar[inode]+node_den_npar[inode])>1 || inode==102) { // derived node [inode = 102 is special case : Gamma96 = (1/1.699387) * Gamma801]
	  Xvector[i][0] += partial*baseseed[quan-1]; 
	}
	Delta[quan-1][i] = -partial;
      }
    }
  }
  //
  DeltaT.Transpose(Delta);
  SDInvMatrix = Delta * InvMeasErrorMatrix_noweak * DeltaT;
  SDMatrix = SDInvMatrix;
  SDMatrix.Invert(&sd_det);
  Vvector = SDMatrix * Delta * InvMeasErrorMatrix_noweak * Xvector; Vvector*=-1;
  XPlusDTV = Xvector + DeltaT*Vvector;
  XPlusDTVT.Transpose(XPlusDTV);
  ChiSquared = XPlusDTVT * InvMeasErrorMatrix_noweak * XPlusDTV;
  _chisquared = ChiSquared[0][0];
  //
  for (ibase=0;ibase<nbase_u;++ibase) {
    basevalue_fit[ibase] = Vvector[ibase][0] ;
    baseerror_fit[ibase] =  TMath::Sqrt(SDMatrix[ibase][ibase]);
  }
  if (uconstrain) {
    basevalue_fit[nbase-1] = 1;
    double error_squared = 0;
    for (ibase=0;ibase<nbase_u;++ibase) {
      basevalue_fit[nbase-1] -= Vvector[ibase][0] ;
      for (jbase=0;jbase<nbase_u;++jbase) {
	error_squared += SDMatrix[ibase][jbase];
      }
      baseerror_fit[nbase-1] = TMath::Sqrt(error_squared);
    }
  }
  //
  for (ibase=0;ibase<nbase_u;++ibase) {
    for (jbase=0;jbase<nbase_u;++jbase) {
      basecov_fit [ibase][jbase] = SDMatrix[ibase][jbase];
      basecorr_fit[ibase][jbase] = (ibase==jbase) ? 1 : basecov_fit[ibase][jbase]/(baseerror_fit[ibase]*baseerror_fit[jbase]);
    }
  }
  if (uconstrain) {
    for (ibase=0;ibase<nbase_u;++ibase){
      basecov_fit[ibase][nbase-1] = 0;
      basecov_fit[nbase-1][ibase] = 0;
      for (jbase=0;jbase<nbase_u;++jbase) {
	basecov_fit[ibase][nbase-1] -= basecov_fit[ibase][jbase];
	basecov_fit[nbase-1][ibase] -= basecov_fit[jbase][ibase];
      }
    }
    basecov_fit [nbase-1][nbase-1] = baseerror_fit[nbase-1]*baseerror_fit[nbase-1];
    basecorr_fit[nbase-1][nbase-1] = 1;
    for (ibase=0;ibase<nbase_u;++ibase){
      basecorr_fit[ibase][nbase-1] = basecov_fit[ibase][nbase-1] / (baseerror_fit[ibase]*baseerror_fit[nbase-1]);
      basecorr_fit[nbase-1][ibase] = basecov_fit[ibase][nbase-1] / (baseerror_fit[ibase]*baseerror_fit[nbase-1]);
    }
  }
}
// ----------------------------------------------------------------------
void process_measurements(int nmeas, int* measnode, double* measerror, double** corrmat,
			  TMatrixD & MeasErrorMatrix,  // output
			  TMatrixD & InvMeasErrorMatrix, //output
			  vector<int> *icorrj, //output
			  int &ncorrij, //output
			  int *ifirstj, //output
			  vector<int> *veccorrij, //output
			  int & newnode_all, //output
			  int* nodegroup_all, //output
			  int* ncycle, // output
			  int* n_nodegroup, // output
			  int* ngroup //output
			  ){
  //
  int i,j,inode;
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) {
	MeasErrorMatrix(i,j)=measerror[i]*measerror[i];
      }else{
	MeasErrorMatrix(i,j)=corrmat[i][j]*measerror[i]*measerror[j]; // by construction diagonal elements of corrmat has not been filled
      }
    }
  }
  //
  InvMeasErrorMatrix = MeasErrorMatrix;
  Double_t determinant;
  InvMeasErrorMatrix.Invert(&determinant);  
  //  cout << Form("Determinant of MeasErrorMatrix = %10.5g\n", determinant);
  //
  for (i=0;i<nmeas;++i){
    for (j=0;j<nmeas;++j){
      if (InvMeasErrorMatrix[i][j]!=0.0) {
	icorrj[i].push_back(j);
      }
    }
  }
  //
  ncorrij=0;
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
  //
  for (i=0;i<ncorrij;++i){
    veccorrij[i].insert(veccorrij[i].end(),icorrj[ifirstj[i]].begin(),icorrj[ifirstj[i]].end());
  }
  //
  newnode_all=-1;
  vector<int> vector_measnodes_all;
  for (i=0;i<nmeas;++i) {
    if (icorrj[i].size()>1) continue;
    ncycle[i] = -1; // cycle number for uncorrelated measurements
    inode=measnode[i];
    vector<int>::iterator itt=find(vector_measnodes_all.begin(),vector_measnodes_all.end(),inode);
    bool is_newnode_all = itt == vector_measnodes_all.end();
    if (is_newnode_all) {
      ++newnode_all;
    }
    nodegroup_all[i] = newnode_all; // node-group number for uncorrelated measurements
    vector_measnodes_all.push_back(inode);
  }
  for (i=0;i<ncorrij;++i) {
    ++newnode_all;
    for (j=0;j<veccorrij[i].size();++j) {
      int meas=veccorrij[i][j];
      ncycle[meas] = i; // cycle number for each correlated measurement
      nodegroup_all[meas] = newnode_all; // node-group number for correlated measurements
    }
  }
  for (inode=0;inode<=newnode_all;++inode){
    n_nodegroup[inode]=0;
    for (i=0;i<nmeas;++i){
      if (nodegroup_all[i]==inode){
	n_nodegroup[inode]+=1;
      }
    }
  }
  for (inode=0;inode<=newnode_all;++inode){
    for (i=0;i<nmeas;++i){
      if (nodegroup_all[i]==inode) {
	ngroup[i] = n_nodegroup[nodegroup_all[i]];
      }
    }
  }
}
// ----------------------------------------------------------------------
void print_measinfo(FILE* thisfile,
		    int nmeas, int* measnode, double* measvalue, double* measerror, double** corrmat,
		    string* expname, string* meastitle, string* measgammaname,
		    char** a_nodename, int* node_num_npar, int* node_den_npar, vector<int> * node_quan,
		    int nbase, vector<int> baseparm, vector<int> basegamma, string* basenode, string* basetitle,
		    int newnode_all, int* nodegroup_all, int* ncycle, int ncorrij, int* ifirstj, vector<int> *icorrj, vector<int> *veccorrij, int* ngroup){
  //
  int i, inode, ibase, j;
  //
  fprintf (thisfile, "nmeas = %d \n", nmeas);
  //
  vector<int> basequan_used_in_measured_basenodes;
  vector<int> basequan_used_in_measured_derivednodes;
  for (i=0;i<nmeas;++i){
    inode=measnode[i];
    if ((node_num_npar[inode]+node_den_npar[inode])>1 || inode==102) { // derived node [inode = 102 is special case : Gamma96 = (1/1.699387) * Gamma801]
      basequan_used_in_measured_derivednodes.insert(basequan_used_in_measured_derivednodes.end(),node_quan[inode].begin(),node_quan[inode].end());
      sort(basequan_used_in_measured_derivednodes.begin(),basequan_used_in_measured_derivednodes.end());
      vector<int>::iterator new_end=unique(basequan_used_in_measured_derivednodes.begin(),basequan_used_in_measured_derivednodes.end());
      basequan_used_in_measured_derivednodes.erase(new_end,basequan_used_in_measured_derivednodes.end());
    }else{ // base node
      basequan_used_in_measured_basenodes.insert(basequan_used_in_measured_basenodes.end(),node_quan[inode].begin(),node_quan[inode].end());
      sort(basequan_used_in_measured_basenodes.begin(),basequan_used_in_measured_basenodes.end());
      vector<int>::iterator new_end=unique(basequan_used_in_measured_basenodes.begin(),basequan_used_in_measured_basenodes.end());
      basequan_used_in_measured_basenodes.erase(new_end,basequan_used_in_measured_basenodes.end());
    }
  }
  fprintf (thisfile, "basequan_used_in_measured_basenodes.size() = %d ", basequan_used_in_measured_basenodes.size());
  fprintf (thisfile, "[IF THIS IS NOT = %d, IT MEANS SOME BASE QUANTITIES ARE NOT MEASURED DIRECTLY BUT ONLY VIA MEASUREMENTS OF DERIVED NODES]\n", nbase); 
  for (ibase=0;ibase<nbase;++ibase) {
    vector<int>::iterator it=find(basequan_used_in_measured_basenodes.begin(),basequan_used_in_measured_basenodes.end(),ibase+1); 
    bool not_present = it == basequan_used_in_measured_basenodes.end();
    if (not_present) fprintf (thisfile, "BASE = %d GAMMA = %d PARM = %d NODE = %s TITLE = %s not measured directly \n",
			      ibase+1, basegamma[ibase], baseparm[ibase], basenode[ibase].data(), basetitle[ibase].data());
  }
  fprintf (thisfile, "\n");
  fprintf (thisfile, "basequan_used_in_measured_derivednodes.size() = %d \n", basequan_used_in_measured_derivednodes.size());
  for (ibase=0;ibase<basequan_used_in_measured_derivednodes.size();++ibase) {
    int quan = basequan_used_in_measured_derivednodes.at(ibase) ;
    fprintf (thisfile, "basetitle = %s with quan = %d gamma = %d appears in ",basetitle[quan-1].data(), quan, basegamma[quan-1]);
    int iquan_occurance=0;
    int inode_occurance=0;
    vector<int> vector_measnodes;
    for (i=0;i<nmeas;++i){
      inode=measnode[i];
      vector<int>::iterator itt=find(vector_measnodes.begin(),vector_measnodes.end(),inode);
      bool is_newnode = itt == vector_measnodes.end();
      for (int ii=0;ii<node_quan[inode].size();++ii){
	if (quan == node_quan[inode].at(ii)) { 
	  ++iquan_occurance;
	  if (is_newnode) {
	    ++inode_occurance;
	  }
	}
      }
      vector_measnodes.push_back(inode);
    }
    fprintf (thisfile, "%d measurements of %d nodes. \nThese %d measurements are :\n", iquan_occurance, inode_occurance, iquan_occurance);
    iquan_occurance = 0;
    vector<int> vector_of_measurements_correlated_with_other_measurements;
    for (i=0;i<nmeas;++i){
      bool meas_with_correlation = false;
      inode=measnode[i];
      for (int ii=0;ii<node_quan[inode].size();++ii){
	if (quan == node_quan[inode].at(ii)) { 
	  ++iquan_occurance;
	  fprintf (thisfile, " (%d) IMEAS = %d NODE = %d NAME = %s GAMMA = %s TITLE = %s EXP = %s\n",
		   iquan_occurance, i+1, inode, a_nodename[inode], measgammaname[i].data(), meastitle[i].data(), expname[i].data());
	  for (j=0;j<nmeas;++j){
	    if (corrmat[i][j] != 0) {
	      meas_with_correlation = true;
	      break;
	    }
	  }
	}
      }
      if (meas_with_correlation) vector_of_measurements_correlated_with_other_measurements.push_back(i+1); // <- NOTE: +1
    }
    if (vector_of_measurements_correlated_with_other_measurements.size()!=0) {
      fprintf (thisfile, "vector_of_measurements_correlated_with_other_measurements has size = %d : ",vector_of_measurements_correlated_with_other_measurements.size());
      for (int i=0; i < vector_of_measurements_correlated_with_other_measurements.size(); ++i) {
	int imeas = vector_of_measurements_correlated_with_other_measurements[i] -1; // <- NOTE: -1
	fprintf (thisfile, "%d(%s) ",imeas+1, expname[imeas].data()); 
      }
      fprintf (thisfile, "\n");
    }else{
      fprintf (thisfile, "vector_of_measurements_correlated_with_other_measurements has size = 0\n");
    }
    fprintf (thisfile, "\n");
  }
  //
  fprintf (thisfile, "ncorrij = %d \n",ncorrij);
  for (i=0;i<ncorrij;++i){
    fprintf (thisfile, "i = %d ifirstj[i]+1 = %3d ",i,ifirstj[i]+1);
    fprintf (thisfile, "veccorrij[i].size() = %2d veccorrij[i][j]+1 = ",veccorrij[i].size());
    for (j=0;j<veccorrij[i].size();++j) fprintf (thisfile, "%d ",veccorrij[i][j]+1);
    fprintf (thisfile, "GammaName = "); 
    for (j=0;j<veccorrij[i].size();++j) fprintf (thisfile, "%s ",measgammaname[veccorrij[i][j]].data());
    fprintf (thisfile,"\n");
  }
  //
  for (inode=0;inode<=newnode_all;++inode){
    for (i=0;i<nmeas;++i){
      if (nodegroup_all[i]==inode) {
	fprintf (thisfile, "i+1 = %3d nodegroup = %2d ngroup = %2d icorrj = %2d ncycle = %2d meas = %10.5e +- %10.5e %6s %s \n",
		 i+1,nodegroup_all[i],ngroup[i],icorrj[i].size(), ncycle[i], measvalue[i], measerror[i], expname[i].data(), meastitle[i].data());
      }
    }
  }
}
// ----------------------------------------------------------------------
void print_measfile(FILE* thisfile, int p, int uconstrain,
		    int nmeas, int* weak, int weakcompare, string* measgammaname, int* measnode, 
		    string* expname, string* author, string* year, string* meastitle,
		    double* measvalue, double* measerror, double**corrmat,
		    char** a_nodename, vector<int> * node_parm,
		    int* node_num_npar, vector<int> * node_num_parm, vector<double> * node_num_coef,
		    int* node_den_npar, vector<int> * node_den_parm, vector<double> * node_den_coef, 
		    int nbase, vector<int> baseparm, vector<int> basegamma, string* basetitle, int* first_quan, double* baseseed){
  int inode,ipar;
  int iimeas=0;
  for (int i=0;i<nmeas;++i){
    if (weak[i]==weakcompare) continue;
    ++iimeas;
    //
    fprintf (thisfile, "* IMEAS = %d \n",     iimeas);
    fprintf (thisfile, "* GAMMANAME = %s \n", measgammaname[i].data());
    fprintf (thisfile, "* DECAYNAME = %s \n", meastitle[i].data());
    //
    inode=measnode[i];
    //    cout << "i = " << i << " inode = " << inode << " nodename = " << a_nodename[inode]
    //	 << " node_num_npar[inode] = " << node_num_npar[inode] << " node_num_parm[inode].size() = " << node_num_parm[inode].size()  
    //	 << " node_den_npar[inode] = " << node_den_npar[inode] << " node_den_parm[inode].size() = " << node_den_parm[inode].size()  
    //	 << " node_parm[inode].size() = " << node_parm[inode].size() 
    //	 << endl;
    //
    fprintf (thisfile, 
	     "* NODENAME = %s found at inode+1 = %d has %d, %d, %d base quantities in numerator, denominator and both [excluding overlap] \n", 
	     a_nodename[inode], inode+1, node_num_npar[inode], node_den_npar[inode],node_parm[inode].size());
    //
    for (ipar=0;ipar<node_num_parm[inode].size();++ipar){
      int parm=node_num_parm[inode].at(ipar);
      vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
      int quan=ibase-baseparm.begin()+1;
      fprintf (thisfile,"*                numerator of inode+1 = %d has gamma = %d parm = %d quan = %d title = %s seed = %20.10f coef = %20.10f\n",
	       inode+1, basegamma[quan-1], parm, quan, basetitle[quan-1].data(), baseseed[quan-1], node_num_coef[inode].at(ipar));
    }
    for (ipar=0;ipar<node_den_parm[inode].size();++ipar){
      int parm=node_den_parm[inode].at(ipar);
      vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
      int quan=ibase-baseparm.begin()+1;
      fprintf (thisfile,"*              denominator of inode+1 = %d has gamma = %d parm = %d quan = %d title = %s seed = %20.10f coef = %20.10f\n",
	       inode+1, basegamma[quan-1], parm, quan, basetitle[quan-1].data(), baseseed[quan-1], node_den_coef[inode].at(ipar));
    }
    fprintf (thisfile, "*  first quantity measured by inode+1 = %d has gamma = %d parm = %d quan = %d title = %s\n",
	     inode+1,basegamma[first_quan[inode]-1],baseparm[first_quan[inode]-1],first_quan[inode],basetitle[first_quan[inode]-1].data());
    //
    fprintf (thisfile, "\nBEGIN %s Gamma%s pub.%s.%s \n\n", expname[i].data(), measgammaname[i].data(), author[i].data(), year[i].data());
    if (p==0) {//COMBOS
      if (uconstrain &&      // SPECIAL CASE [because these nodes contain Gamma103 [used to express unitarity constraint]
	  ((inode+1)==80 ||  // NODE = 79 NAME = S035R33 GAMMA = 102, because Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104) 
	   (inode+1)==82)) { // NODE = 81 NAME = S035R38 GAMMA = 103
	fprintf (thisfile, "MEASUREMENT  m_Gamma%d statistical systematic \n",3);
	fprintf (thisfile, "DATA         m_Gamma%d statistical systematic \n",3);
      }else{
	fprintf (thisfile, "MEASUREMENT  m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
	fprintf (thisfile, "DATA         m_Gamma%d statistical systematic \n",basegamma[first_quan[inode]-1]);
      }
    }else if (p==1) {//ALUCOMB
      fprintf (thisfile, "MEASUREMENT  m_Gamma%s statistical systematic \n",measgammaname[i].data());
      fprintf (thisfile, "DATA         m_Gamma%s statistical systematic \n",measgammaname[i].data());
    }
    fprintf (thisfile, "             %20.10f %20.10f  0 \n",measvalue[i],measerror[i]);
    bool firstcorr=true;
    int jjmeas=0;
    for (int j=0;j<nmeas;++j) {
      if (weak[j]==weakcompare) continue;
      ++jjmeas;
      if (corrmat[i][j]!=0) {
	if (firstcorr) {fprintf (thisfile, " \n"); firstcorr=false;}
	fprintf (thisfile, "STAT_CORR_WITH %s Gamma%s pub.%s.%s %20.10f ! IMEAS = %d \n", 
		 expname[j].data(), measgammaname[j].data(), author[j].data(), year[j].data(), corrmat[i][j], jjmeas);
      }
    }
    fprintf (thisfile, " \nEND \n\n");
  }
}
// ----------------------------------------------------------------------
void print_avefile(FILE* thisfile, int p, int uconstrain,
		   int nmeas, int* weak, int weakcompare, string* measgammaname, int* measnode,
		   vector<int> * node_parm, vector<int> * node_quan, 
		   int* node_num_npar, vector<int> * node_num_parm, vector<double> * node_num_coef,
		   int* node_den_npar, vector<int> * node_den_parm, vector<double> * node_den_coef, 
		   int nbase, vector<int> baseparm, vector<int> basegamma, string* basetitle, int* first_quan, 
		   double* baseseed, double* node_num, double* node_den,  vector<double> * node_part){
  int i,inode,ipar,iimeas,ibase;
  fprintf (thisfile, "BEGIN   PDG+BABAR+BELLE all_methods \n\n");
  fprintf (thisfile, "COMBINE * * * \n\n");
  for (ibase=0;ibase<nbase;++ibase){
    if (p==0&&uconstrain&&ibase==(nbase-1)){/* skip */}else{
      fprintf (thisfile, "MEASUREMENT m_Gamma%d statistical systematic   ! NQUAN = %d \n",basegamma.at(ibase),ibase+1);  
    }
  }
  if (p==0){
    fprintf (thisfile, "\nCALL DUMP_MASTER_INC \n\n");
    fprintf (thisfile, "SPARAMETER CHI2_N_SYM_PRT -1.0 0. \n\n");
    fprintf (thisfile, "SPARAMETER CHI2_N_SYM_INV 0 0    \n\n");
  }
  //
  vector<int> vector_measnodes;
  int usum=0;//number of [unique] measurements to be expressed as linearized sum of base quantities
  for (i=0;i<nmeas;++i){
    if (weak[i]==weakcompare) continue;
    inode=measnode[i];
    vector<int>::iterator itt=find(vector_measnodes.begin(),vector_measnodes.end(),inode);
    bool is_newnode = itt == vector_measnodes.end();
    if ((node_num_npar[inode]+node_den_npar[inode])>1 || inode==102) { // derived node [inode = 102 is special case : Gamma96 = (1/1.699387) * Gamma801]
      if (p==1&&is_newnode){//new node
	++usum;
	fprintf (thisfile, "MEASUREMENT m_Gamma%s statistical systematic   ! NQUAN = %d \n",measgammaname[i].data(),ibase+usum);  
      }
    }
    vector_measnodes.push_back(inode);      
  }
  //
  if (p==1&&!uconstrain){
    fprintf (thisfile, "*--- dummy node to remove unitarity constraint\nMEASUREMENT m_Gamma998 statistical systematic\n");
  }
  if (p==1) {
    fprintf (thisfile, "*--- Gamma110, only present in constraint\nMEASUREMENT m_Gamma110 statistical systematic\n");
  }
  //
  int isum=0;//number of          measurements to be expressed as linearized sum of base quantities
  vector_measnodes.erase(vector_measnodes.begin(),vector_measnodes.end());
  iimeas=0;
  for (int i=0;i<nmeas;++i){
    if (weak[i]==weakcompare) continue;
    ++iimeas;
    inode=measnode[i];
    vector<int>::iterator itt=find(vector_measnodes.begin(),vector_measnodes.end(),inode);
    bool is_newnode = itt == vector_measnodes.end();
    if (p==1 && !is_newnode) continue; // ALUCOMB needs it only once
    if ((node_num_npar[inode]+node_den_npar[inode])>1 || inode==102) { // derived node [inode = 102 is special case : Gamma96 = (1/1.699387) * Gamma801]
      //
      if (p==0&&(((inode+1)==80)||((inode+1)==82))) continue; // SPECIAL CASE [because these are derived nodes containing Gamma103 ]
      ++isum; // translate C index to Fortran index
      //
      // PRINT NODE DEFINITION
      //
      fprintf (thisfile, "\n* Gamma%s = ",measgammaname[i].data());
      for (ipar=0; ipar < node_num_parm[inode].size(); ++ipar) {
	if (ipar==0) { fprintf (thisfile, "(") ; } else {fprintf (thisfile, " + ") ;}
	int parm=node_num_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	fprintf (thisfile, "%20.10f*Gamma%d",node_num_coef[inode].at(ipar), basegamma[quan-1]);
	if (ipar==node_num_parm[inode].size()-1) fprintf (thisfile, ")");
      }
      if (node_den_parm[inode].size()==0) fprintf (thisfile, "\n") ; 
      for (ipar=0; ipar < node_den_parm[inode].size(); ++ipar) {
	if (ipar==0) { fprintf (thisfile, " / (") ; } else {fprintf (thisfile, " + ") ;}
	int parm=node_den_parm[inode].at(ipar);
	vector<int>::iterator ibase=find(baseparm.begin(),baseparm.end(),parm);
	int quan=ibase-baseparm.begin()+1;
	fprintf (thisfile, "%20.10f*Gamma%d",node_den_coef[inode].at(ipar), basegamma[quan-1]);
	if (ipar==node_den_parm[inode].size()-1) fprintf (thisfile, ")\n");
      }
      //
      double offset = -node_num[inode]; // - [ f(x0,y0) - df/dx(x=x0) x0 - df/dy(y=y0) y0 - ...]
      if (node_den_npar[inode]>0) offset /= node_den[inode];
      //	  cout << inode << " " << node_num[inode] << " " << node_den[inode] << " " << offset << " " << measvalue[i] << endl;
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int quan=node_quan[inode].at(ipar);
	double partial=node_part[inode].at(ipar);
	offset += partial*baseseed[quan-1]; 
      }
      if (p==0) { // COMBOS
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d    %2d %d\n",isum,iimeas,node_parm[inode].size()); 
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_AD %20.10f 1.0\n",isum,offset); 
      } else if (p==1) { // ALUCOMB
	fprintf (thisfile, "CONSTRAINT Gamma%s.c %20.10f Gamma%s -1", measgammaname[i].data(), offset, measgammaname[i].data());
      }
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int quan=node_quan[inode].at(ipar);
	double partial=node_part[inode].at(ipar);
	if (p==0) { // COMBOS
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_%2.2d %2d %20.10f ! Gamma%d \n",isum,ipar+1,quan,partial,basegamma.at(quan-1));
	}else if (p==1) { // ALUCOMB
	  fprintf (thisfile, " Gamma%d %20.10f", basegamma[quan-1], partial);
	}
      }
      if (p==1) fprintf (thisfile, "\n");
    }
    vector_measnodes.push_back(inode);      
  }
  if (p==0) {
    iimeas=0;
    for (int i=0;i<nmeas;++i){
      if (weak[i]==weakcompare) continue;
      ++iimeas;
      inode=measnode[i];
      //
      if ((inode+1)==80) { // SPECIAL CASE : NODE = 79 NAME = S035R33 GAMMA = 102 :: Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104)
	if (!uconstrain) {
	  fprintf (thisfile, "\n*Gamma102 = (1.000000*Gamma103 + 1.000000*Gamma104)\n");
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d    %2d %d\n",++isum,iimeas,node_parm[inode].size()); 
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_AD %20.10f 1.0\n",isum,0); 
	  for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	    int quan=node_quan[inode].at(ipar);
	    double partial=node_part[inode].at(ipar);
	    fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_%2.2d %2d %20.10f ! Gamma%d \n",isum,ipar+1,quan,partial,basegamma.at(quan-1));
	  }
	}else{
	  fprintf (thisfile, "\n*Gamma102 = 1 - Gamma3   - Gamma5   - Gamma9   - Gamma10  - Gamma14  - Gamma16\n");
	  fprintf (thisfile, "*             - Gamma20  - Gamma23  - Gamma27  - Gamma28  - Gamma30 - Gamma35\n");
	  fprintf (thisfile, "*             - Gamma37 - Gamma40  - Gamma42  - Gamma47  - Gamma48  - Gamma62\n");
	  fprintf (thisfile, "*             - Gamma70 - Gamma77  - Gamma78  - Gamma85  - Gamma89  - Gamma93\n");
	  fprintf (thisfile, "*             - Gamma94 - Gamma126 - Gamma128 - Gamma800 - Gamma151 - Gamma152\n");
	  fprintf (thisfile, "*             - Gamma130 - Gamma132 - Gamma44 - Gamma53 - Gamma801\n");
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d    %2d  %2d \n",++isum,iimeas,nbase-2); 
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_AD -1 +1 \n",isum); // becomes a measurement of -1+Gamma102; thats why the coefficients below have - sign 
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_01  1 -1 ! Gamma3  \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_02  2 -1 ! Gamma5  \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_03  3 -1 ! Gamma9  \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_04  4 -1 ! Gamma10 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_05  5 -1 ! Gamma14 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_06  6 -1 ! Gamma16 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_07  7 -1 ! Gamma20 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_08  8 -1 ! Gamma23 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_09  9 -1 ! Gamma27 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_10 10 -1 ! Gamma28 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_11 11 -1 ! Gamma30 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_12 12 -1 ! Gamma35 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_13 13 -1 ! Gamma37 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_14 14 -1 ! Gamma40 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_15 15 -1 ! Gamma42 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_16 16 -1 ! Gamma47 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_17 17 -1 ! Gamma48 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_18 18 -1 ! Gamma62 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_19 19 -1 ! Gamma70 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_20 20 -1 ! Gamma77 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_21 21 -1 ! Gamma78 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma85 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma89 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_24 24 -1 ! Gamma93 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_25 25 -1 ! Gamma94 \n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_26 27 -1 ! Gamma126\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_27 28 -1 ! Gamma128\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_28 29 -1 ! Gamma800\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_29 30 -1 ! Gamma151\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_30 31 -1 ! Gamma152\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_31 32 -1 ! Gamma130\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_32 33 -1 ! Gamma132\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_33 34 -1 ! Gamma44\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_34 35 -1 ! Gamma53\n",isum);
	  fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_35 36 -1 ! Gamma801\n",isum);
	}
      }
      //
      if (uconstrain && (inode+1)==82) { // SPECIAL CASE : NODE = 81 NAME = S035R38 GAMMA = 103 
	fprintf (thisfile, "\n*Gamma103 = 1 - Gamma3   - Gamma5   - Gamma9   - Gamma10  - Gamma14  - Gamma16\n");
	fprintf (thisfile, "*             - Gamma20  - Gamma23  - Gamma27  - Gamma28  - Gamma30 - Gamma35\n");
	fprintf (thisfile, "*             - Gamma37 - Gamma40  - Gamma42  - Gamma47  - Gamma48  - Gamma62\n");
	fprintf (thisfile, "*             - Gamma70 - Gamma77  - Gamma78  - Gamma85  - Gamma89  - Gamma93\n");
	fprintf (thisfile, "*             - Gamma94 - Gamma104 - Gamma126 - Gamma128 - Gamma800 - Gamma151 - Gamma152\n");
	fprintf (thisfile, "*             - Gamma130 - Gamma132 - Gamma44 - Gamma53 - Gamma801\n");
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d    %2d %2d \n",++isum,iimeas,nbase-1); 
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_AD -1 +1 \n",isum); // becomes a measurement of -1+Gamma103; thats why the coefficients below have - sign 
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_01  1 -1 ! Gamma3  \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_02  2 -1 ! Gamma5  \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_03  3 -1 ! Gamma9  \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_04  4 -1 ! Gamma10 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_05  5 -1 ! Gamma14 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_06  6 -1 ! Gamma16 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_07  7 -1 ! Gamma20 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_08  8 -1 ! Gamma23 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_09  9 -1 ! Gamma27 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_10 10 -1 ! Gamma28 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_11 11 -1 ! Gamma30 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_12 12 -1 ! Gamma35 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_13 13 -1 ! Gamma37 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_14 14 -1 ! Gamma40 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_15 15 -1 ! Gamma42 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_16 16 -1 ! Gamma47 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_17 17 -1 ! Gamma48 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_18 18 -1 ! Gamma62 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_19 19 -1 ! Gamma70 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_20 20 -1 ! Gamma77 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_21 21 -1 ! Gamma78 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_22 22 -1 ! Gamma85 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_23 23 -1 ! Gamma89 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_24 24 -1 ! Gamma93 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_25 25 -1 ! Gamma94 \n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_26 26 -1 ! Gamma104\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_27 27 -1 ! Gamma126\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_28 28 -1 ! Gamma128\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_29 29 -1 ! Gamma800\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_30 30 -1 ! Gamma151\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_31 31 -1 ! Gamma152\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_32 32 -1 ! Gamma130\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_33 33 -1 ! Gamma132\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_34 34 -1 ! Gamma44\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_35 35 -1 ! Gamma53\n",isum);
	fprintf (thisfile, "SPARAMETER CHI2_N_SYM_%2.2d_36 36 -1 ! Gamma801\n",isum);
      }
    }
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_NSUM  %d 0 \n",isum); 
  }
  if (p==1) {
    if (uconstrain) {
      fprintf (thisfile, "\n* unitarity constraint applied (sum of base nodes without dummy node)\n");
    }else{
      fprintf (thisfile, "\n* unitarity constraint NOT applied (sum of base nodes with dummy node)\n");
    }
    fprintf (thisfile, "CONSTRAINT GammaAll 1\n");
    fprintf (thisfile, "  Gamma3   1 Gamma5   1 Gamma9   1 Gamma10  1 Gamma14  1 Gamma16  1\n");
    fprintf (thisfile, "  Gamma20  1 Gamma23  1 Gamma27  1 Gamma28  1 Gamma30  1 Gamma35  1\n");
    fprintf (thisfile, "  Gamma37  1 Gamma40  1 Gamma42  1 Gamma47  1 Gamma48  1 Gamma62  1\n");
    fprintf (thisfile, "  Gamma70  1 Gamma77  1 Gamma78  1 Gamma85  1 Gamma89  1 Gamma93  1\n");
    fprintf (thisfile, "  Gamma94  1 Gamma103 1 Gamma104 1 Gamma126 1 Gamma128 1 Gamma800 1 Gamma151 1 Gamma152 1\n");
    fprintf (thisfile, "  Gamma130 1 Gamma132 1 Gamma44  1 Gamma53  1 Gamma801 1\n");
    if (!uconstrain) fprintf (thisfile, "  Gamma998 1\n");
  }
  if (p==0){
    //*    4    10     7   R7     7.000000E-03   0.006910   0.000219   0.000228   1.04 tau- --> K- nu(tau)                                 
    //*    6    16   182   B29    4.500000E-03   0.004525   0.000265   0.000267   1.01 tau- --> K- pi0 nu(tau)                             
    //*    8    23   115   B30    6.000000E-04   0.000581   0.000225   0.000232   1.03 tau- --> K- 2pi0 nu(tau) (ex.K0)                    
    //*   10    28   116   B31    4.000000E-03   0.000417   0.000219   0.000220   1.01 tau- --> K- 3pi0 nu(tau) (ex.K0, eta)               
    //*   12    35   117   B32    9.600000E-03   0.008963   0.000367   0.000409   1.11 tau- --> pi- Kbar0 nu(tau)                          
    //*   14    40   118   B33    4.000000E-03   0.003778   0.000368   0.000374   1.02 tau- --> pi- Kbar0 pi0 nu(tau)                      
    //*   22    85   260   C21    3.000000E-03   0.003335   0.000223   0.000352   1.58 tau- --> K- pi+ pi- nu(tau) (ex.K0)                 
    //*   23    89   285   C54    6.000000E-04   0.000730   0.000117   0.000122   1.04 tau- --> K- pi+ pi- pi0 nu(tau) (ex.K0,eta)         
    //*   28   128   109   B20    3.000000E-04   0.000268   0.000063   0.000063   1.00 tau- --> eta K- nu(tau)                             
    //*   30   151   295   C61    4.100000E-04   0.000410   0.000092   0.000092   1.00 tau- --> K-  omega nu(tau)                                            
    //*   32   130   266   C27    1.770000E-04   0.000177   0.000090   0.000090   1.00 tau- --> eta K- pi0 nu(tau)
    //*   33   132   267   C28    2.200000E-04   0.000220   0.000073   0.000073   1.00 tau- --> eta pi- K0bar nu(tau)
    //*   34    44   238   B98    2.600000E-04   0.000260   0.000240   0.000240   1.00 tau- --> pi- K0bar 2pi0 nu(tau)
    //*   35    53   244   C5     2.300000E-04   0.000230   0.000203   0.000203   1.00 tau- --> K0bar h+ h- h- nu(tau)
    //*   36   801   801   Z02    2.000000E-05   0.000020   0.000010   0.000010   1.00 tau- --> K- phi nu(tau) [phi->KK]                           
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_PSUM   1  0 ! print sum of strange decay nodes\n");
    fprintf (thisfile, "\n* Print Gamma(tau -> X-(S=1) nu)");
    fprintf (thisfile, "\n*Gamma110 = Gamma10  + Gamma16   + Gamma23   + Gamma28  + Gamma35  + Gamma40 + Gamma85 + Gamma89 + Gamma128\n");
    fprintf (thisfile, "*         + Gamma151 + Gamma130  + Gamma132  + Gamma44  + Gamma53  + Gamma801\n");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1     15 0 ");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_01  4   1 ! Gamma10");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_02  6   1 ! Gamma16");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_03  8   1 ! Gamma23");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_04  10  1 ! Gamma28");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_05  12  1 ! Gamma35");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_06  14  1 ! Gamma40");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_07  22  1 ! Gamma85");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_08  23  1 ! Gamma89");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_09  28  1 ! Gamma128");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_10  30  1 ! Gamma151");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_11  32  1 ! Gamma130");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_12  33  1 ! Gamma132");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_13  34  1 ! Gamma44");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_14  35  1 ! Gamma53");
    fprintf (thisfile, "\nSPARAMETER CHI2_N_SYM_P1_15  36  1 ! Gamma801");
    fprintf (thisfile, "\n");
  }
  if (p==1) {
    fprintf (thisfile, "\n* --- compute Gamma(tau -> Xs nu)/G(total)\n");
    fprintf (thisfile, "COMBOFQUANT Gamma110\n");
    fprintf (thisfile, " 1 Gamma10  1 Gamma16  1 Gamma23  1 Gamma28  1 Gamma35  1 Gamma40  1 Gamma85  1 Gamma89  1 Gamma128\n");
    fprintf (thisfile, " 1 Gamma151 1 Gamma130 1 Gamma132 1 Gamma44  1 Gamma53  1 Gamma801\n");
  }
  fprintf (thisfile, "\nCALL CHI2_N_SYM\n");
  fprintf (thisfile, "\nEND\n");
}
// ----------------------------------------------------------------------
int main(int argc, char* argv[]){
  //Argument variables
  Int_t uconstrain   = (argc>1) ? atoi(argv[1]) : 1; // 1: unitarity constrained; 0 : not constrained
  Int_t doalephhcorr = (argc>2) ? atoi(argv[2]) : 1; // 1: do aleph hcorr; 0: dont
  //
  string sconstrain = (uconstrain) ? "constrained" : "unconstrained";
  string salephhcorr = (doalephhcorr) ? "_aleph_hcorr" : "";
  //
  gSystem->Load("libMatrix");
  //
  int ipar,iimeas,jjmeas,jpar;
  int i,j;
  int inode,jnode;
  int ibase,jbase;
  int im, jm;
  //
  const int weakcompare=1;
  //
  // READ BASE PARAMETERS
  //
  const int nbase=37;
  int nbase_u = (uconstrain) ? nbase-1 : nbase;
  vector<int> basequan;
  vector<int> basegamma;
  vector<int> baseparm;
  double      baseseed[nbase];
  double      baseseed_orig[nbase];
  double      basefitval_orig[nbase];
  double      basefiterr_orig[nbase];
  double      basescalerr_orig[nbase];
  double      basescalfac_orig[nbase];
  string      basenode[nbase];
  string      basetitle[nbase];
  //
  ibase=0;
  ifstream ifsbase("base_def.txt") ;
  if (!ifsbase.good()) {cout << "Cannot open input file : base_def.txt" << endl ; exit(1) ;}
  char buffer[256]; 
  while(ifsbase.good()) {
    if (ifsbase.eof()) break;
    char firstch(' ') ; ifsbase.get(firstch) ;
    if (firstch=='#'||firstch=='\n') { // Skip such lines
      ifsbase.ignore(256,'\n') ;
    } else if (firstch=='*') {  // measurement line
      int dummy_int;       
      ifsbase >> dummy_int ; basequan.push_back(dummy_int);
      ifsbase >> dummy_int ; basegamma.push_back(dummy_int);
      ifsbase >> dummy_int ; baseparm.push_back(dummy_int);
      ifsbase >> basenode[ibase];
      ifsbase >> baseseed_orig[ibase];
      ifsbase >> basefitval_orig[ibase];
      ifsbase >> basefiterr_orig[ibase];
      ifsbase >> basescalerr_orig[ibase];
      ifsbase >> basescalfac_orig[ibase];
      ifsbase.getline(buffer,256,'\n');
      string stringbuffer=string(buffer);
      basetitle[ibase]="";
      bool first=false;
      for (string::iterator si=stringbuffer.begin();si!=stringbuffer.end();++si){
	if (13 == (int)*si) break; // ^M is 13
	if (' '!=*si) first=true;
	if (first) {
	  if (' ' == *si && ' ' == *(si+1) ) break;
	  basetitle[ibase]+=*si;
	}
      }
      baseseed[ibase] = baseseed_orig[ibase];
      //      baseseed[ibase] = basefitval_orig[ibase];
      // cout << basequan[ibase] << " " << basegamma[ibase] << " " << baseparm[ibase] << " " << basenode[ibase] << " " << baseseed[ibase] << " " << basetitle[ibase] << endl;
      ++ibase;
    }
  }
  //
  // READ INPUT NODES
  // 
  const int      nnode=108; 
  vector<string> nodegammaname;
  vector<string> nodename;
  vector<int>    node_num_parm[nnode];
  vector<double> node_num_coef[nnode];
  vector<int>    node_den_parm[nnode];
  vector<double> node_den_coef[nnode];
  string         nodetitle[nnode];
  //
  inode=0;
  nodegammaname.push_back("Gamma128");     
  nodename.push_back("S035B20");    
  nodetitle[inode]="G(eta K- nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(1.        );              
  ++inode;  //1
  nodegammaname.push_back("Gamma19by13");  
  nodename.push_back("S035B21"); 
  nodetitle[inode]="G(h- 2pi0 nu(tau) (ex.K0))/G(h- pi0 nu(tau))";
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(16 );       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(182);       node_den_coef[inode].push_back(1.        );               
  ++inode;  //2
  nodegammaname.push_back("Gamma26by13"); 
  nodename.push_back("S035B22"); 
  nodetitle[inode]="G(h- 3pi0 nu(tau))/G(h- pi0 nu(tau))";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2200E-01);              
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(16 );       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(182);       node_den_coef[inode].push_back(1.        );               
  ++inode; //3
  nodegammaname.push_back("Gamma30"); 
  nodename.push_back("S035B23"); 
  nodetitle[inode]="G(h- 4pi0 nu(tau) (ex.K0,eta))/G(total)";
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );              
  ++inode; //4
  nodegammaname.push_back("Gamma76by54"); 
  nodename.push_back("S035B25"); 
  nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex.K0))/G(h- h- h+ >=0 neutrals >=0K(L)0 nu(tau) )";
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
  //node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(9.1010E-01);            // h-  omega
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(9.1010E-01);              // pi- omega
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(9.1010E-01);              // K-  omega
  ++inode;//5
  nodegammaname.push_back("Gamma152by54"); 
  nodename.push_back("S035B26"); 
  nodetitle[inode]="G(h- omega pi0 nu(tau))/G(h- h- h+ >=0 neutrals >=0K(L)0 nu(tau) )";
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
  //node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(9.1010E-01);            // h-  omega
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(9.1010E-01);              // pi- omega
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(9.1010E-01);              // K-  omega
  ++inode;//6
  nodegammaname.push_back("Gamma152by76"); 
  nodename.push_back("S035B27"); 
  nodetitle[inode]="G(h- omega pi0 nu(tau))/G(h- h- h+ 2pi0 nu(tau) (ex.K0))";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(8.8800E-01);               
  node_den_parm[inode].push_back(216);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(58 );       node_den_coef[inode].push_back(2.2600E-01);               
  ++inode;//7
  nodegammaname.push_back("Gamma16"); 
  nodename.push_back("S035B29"); 
  nodetitle[inode]="G(K- pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );              
  ++inode;//8
  nodegammaname.push_back("Gamma23"); 
  nodename.push_back("S035B30"); 
  nodetitle[inode]="G(K- 2pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );              
  ++inode;//9
  nodegammaname.push_back("Gamma28"); 
  nodename.push_back("S035B31"); 
  nodetitle[inode]="G(K- 3pi0 nu(tau) (ex.K0, eta))/G(total)";
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );              
  ++inode;//10
  nodegammaname.push_back("Gamma35"); 
  nodename.push_back("S035B32"); 
  nodetitle[inode]="G(pi- Kbar0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );              
  ++inode;//11
  nodegammaname.push_back("Gamma40"); 
  nodename.push_back("S035B33"); 
  nodetitle[inode]="G(pi- Kbar0 pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );              
  ++inode;//12
  nodegammaname.push_back("Gamma42"); 
  nodename.push_back("S035B34"); 
  nodetitle[inode]="G(K- K0 pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );              
  ++inode;//13
  nodegammaname.push_back("Gamma92"); 
  nodename.push_back("S035B37"); 
  nodetitle[inode]="G(K- K+ pi- >=0 neut.  nu(tau))/G(total)";
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  ++inode;//14
  nodegammaname.push_back("Gamma33"); 
  nodename.push_back("S035B43"); 
  nodetitle[inode]="G(K(S)0 (particles)- nu(tau))/G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(5.0000E-01);              
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(5.0000E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(5.0000E-01);               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(5.0000E-01);               
  ++inode;//15
  nodegammaname.push_back("Gamma106"); 
  nodename.push_back("S035B45"); 
  nodetitle[inode]="G(( 5pi )- nu(tau))/G(total)";
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(5.5300E-01);               
  ++inode;//16
  nodegammaname.push_back("Gamma46"); 
  nodename.push_back("S035B51"); 
  nodetitle[inode]="G(pi- K0 Kbar0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(2.0000E+00);              
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(1.        );               
  ++inode;//17
  nodegammaname.push_back("Gamma66"); 
  nodename.push_back("S035B53"); 
  nodetitle[inode]="G(h- h- h+ pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);      
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.7000E-02);           
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );           
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );           
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );           
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);  // h-  omega             
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(8.8800E-01);      // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(8.8800E-01);      // K-  omega
  ++inode;//18
  nodegammaname.push_back("Gamma67"); 
  nodename.push_back("S035B54"); 
  nodetitle[inode]="G(h- h- h+ pi0 nu(tau) (ex. K0, omega))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;//19
  nodegammaname.push_back("Gamma20"); 
  nodename.push_back("S035B55"); 
  nodetitle[inode]="G(pi- 2pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );              
  ++inode;//20
  nodegammaname.push_back("Gamma27"); 
  nodename.push_back("S035B56"); 
  nodetitle[inode]="G(pi- 3pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );              
  ++inode;//21
  nodegammaname.push_back("Gamma78"); 
  nodename.push_back("S035B57"); 
  nodetitle[inode]="G(h- h- h+ 3pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );              
  ++inode;//22
  nodegammaname.push_back("Gamma152"); 
  nodename.push_back("S035B58"); 
  nodetitle[inode]="G(h- omega pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.        );              
  ++inode;//23
  nodegammaname.push_back("Gamma76"); 
  nodename.push_back("S035B59"); 
  nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);              
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  ++inode;//24
  nodegammaname.push_back("Gamma57"); 
  nodename.push_back("S035B62"); 
  nodetitle[inode]="G(h- h- h+ nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );        
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );            
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);   // h-  omega            
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.7000E-02);       // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.7000E-02);       // K-  omega
  ++inode;//25
  nodegammaname.push_back("Gamma55"); 
  nodename.push_back("S035B63"); 
  nodetitle[inode]="G( h- h- h+ >=0 neutrals  nu(tau)  (ex. K(S)0 --> pi+ pi-)(``3-prong''))/G(total)";
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
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.1010E-01);  // h-  omega             
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(9.1010E-01);      // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(9.1010E-01);      // K-  omega
  ++inode;//26
  nodegammaname.push_back("Gamma57by55"); 
  nodename.push_back("S035B64"); 
  nodetitle[inode]="G(h- h- h+ nu(tau) (ex.K0))/G( h- h- h+ >=0 neutrals  nu(tau)  (ex. K(S)0 --> pi+ pi-)(``3-prong''))";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );            
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);  // h-  omega             
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.7000E-02);      // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.7000E-02);      // K-  omega
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
  //  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(9.1010E-01);  // h-  omega
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(9.1010E-01);      // pi- omega         
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(9.1010E-01);      // K-  omega
  ++inode;//27
  nodegammaname.push_back("Gamma34"); 
  nodename.push_back("S035B67"); 
  nodetitle[inode]="G(h- Kbar0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  ++inode;//28
  nodegammaname.push_back("Gamma39"); 
  nodename.push_back("S035B68"); 
  nodetitle[inode]="G(h- Kbar0 pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );               
  ++inode;//29
  nodegammaname.push_back("Gamma47"); 
  nodename.push_back("S035B69"); 
  nodetitle[inode]="G(pi- K(S)0 K(S)0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );              
  ++inode;//30
  nodegammaname.push_back("Gamma58"); 
  nodename.push_back("S035B71"); 
  nodetitle[inode]="G(h- h- h+ nu(tau) (ex.K0,omega))/G(total)";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  ++inode;//31
  nodegammaname.push_back("Gamma77");
  nodename.push_back("S035B72"); 
  nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau) (ex.K0,omega,eta))/G(total)";
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );              
  ++inode;//32
  nodegammaname.push_back("Gamma8"); 
  nodename.push_back("S035B73"); 
  nodetitle[inode]="G(h- nu(tau))/G(total)";
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  ++inode;//33
  nodegammaname.push_back("Gamma18");
  nodename.push_back("S035B74"); 
  nodetitle[inode]="G(h- 2pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.5700E-01);               
  ++inode;//34
  nodegammaname.push_back("Gamma1");
  nodename.push_back("S035B75"); 
  nodetitle[inode]="G(particle->=0 neutrals >=0K0 nu(tau)  (``1-prong''))/G(total)";
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
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.0000E-02);       // h-  omega        
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(9.0000E-02);           // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(9.0000E-02);           // K-  omega
  ++inode;//35
  nodegammaname.push_back("Gamma65");
  nodename.push_back("S035B76"); 
  nodetitle[inode]="G(h- h- h+ pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.7000E-02);               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);               // h-  omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(8.8800E-01);                   // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(8.8800E-01);                   // K-  omega
  ++inode;//36
  nodegammaname.push_back("Gamma75");
  nodename.push_back("S035B77");
  nodetitle[inode]="G(h- h- h+ 2pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);              
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(4.3070E-01);               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  ++inode;//37
  nodegammaname.push_back("Gamma64");
  nodename.push_back("S035B78"); 
  nodetitle[inode]="G(h- h- h+ >=1 pi0 nu(tau) (ex. K0))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(9.1010E-01);               
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01); // h-  omega              
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(8.8800E-01);     // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(8.8800E-01);     // K-  omega
  ++inode;//38
  nodegammaname.push_back("Gamma29"); 
  nodename.push_back("S035B79"); 
  nodetitle[inode]="G(h- 4pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.1900E-01);               
  ++inode;//39
  nodegammaname.push_back("Gamma8by5");
  nodename.push_back("S035B97"); 
  nodetitle[inode]="G(h- nu(tau))/G(e- nubar(e) nu(tau))";
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(1.        );               
  ++inode;//40
  nodegammaname.push_back("Gamma12");
  nodename.push_back("S035C01"); 
  nodetitle[inode]="G(h- >= 1pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2500E-01);              
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.2500E-01);               
  ++inode;//41
  nodegammaname.push_back("Gamma25"); 
  nodename.push_back("S035C02"); 
  nodetitle[inode]="G(h- >= 3pi0 nu(tau) (ex. K0))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2500E-01);              
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.2500E-01);               
  ++inode;//42
  nodegammaname.push_back("Gamma74");
  nodename.push_back("S035C03");
  nodetitle[inode]="G(h- h- h+ >= 2pi0 nu(tau) (ex. K0))/G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(8.8800E-01);              
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(2.2600E-01);               
  ++inode;//43
  nodegammaname.push_back("Gamma48");
  nodename.push_back("S035C1");
  nodetitle[inode]="G(pi- K(S)0 K(L)0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(1.        );              
  ++inode;//44
  nodegammaname.push_back("Gamma59");
  nodename.push_back("S035C18"); 
  nodetitle[inode]="G(pi- pi+ pi- nu(tau))/G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(3.4310E-01);              
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );               
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);  // h-  omega             
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.7000E-02);      // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.7000E-02);      // K-  omega
  ++inode;//45
  nodegammaname.push_back("Gamma60");
  nodename.push_back("S035C19");
  nodetitle[inode]="G(pi- pi+ pi- nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );              
  // node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);  // h-  omega             
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.7000E-02);     // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.7000E-02);     // K-  omega
  ++inode;//46
  nodegammaname.push_back("Gamma62");
  nodename.push_back("S035C20");
  nodetitle[inode]="G(pi- pi+ pi- nu(tau) (ex.K0,omega))/G(total)";
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );              
  ++inode;//47
  nodegammaname.push_back("Gamma85");
  nodename.push_back("S035C21");
  nodetitle[inode]="G(K- pi+ pi- nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );              
  ++inode;//48
  nodegammaname.push_back("Gamma68");
  nodename.push_back("S035C22");
  nodetitle[inode]="G(pi- pi+ pi- pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.7000E-02);              
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);           // h-  omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(8.8800E-01);               // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(8.8800E-01);               // K-  omega
  ++inode;//49
  nodegammaname.push_back("Gamma69");
  nodename.push_back("S035C23");
  nodetitle[inode]="G(pi- pi+ pi- pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.7000E-02);              
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);           // h-  omega    
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(8.8800E-01);               // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(8.8800E-01);               // K-  omega
  ++inode;//50
  nodegammaname.push_back("Gamma70");
  nodename.push_back("S035C24");
  nodetitle[inode]="G(pi- pi+ pi- pi0 nu(tau) (ex.K0,omega))/G(total)";
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );              
  ++inode;//51
  nodegammaname.push_back("Gamma88"); 
  nodename.push_back("S035C25");
  nodetitle[inode]="G(K- pi+ pi- pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;//52
  nodegammaname.push_back("Gamma80"); 
  nodename.push_back("S035C31"); 
  nodetitle[inode]="G(K- h+ pi- nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  ++inode;//53
  nodegammaname.push_back("Gamma80by60"); 
  nodename.push_back("S035C32"); 
  nodetitle[inode]="G(K- h+ pi- nu(tau) (ex.K0))/G(pi- pi+ pi- nu(tau) (ex.K0))";
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(1.        );               
  //  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(1.7000E-02);               // h-  omega
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(1.7000E-02);                   // pi- omega
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(1.7000E-02);                   // K-  omega
  ++inode;//54
  nodegammaname.push_back("Gamma81"); 
  nodename.push_back("S035C33"); 
  nodetitle[inode]="G(K- h+ pi- pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;//55
  nodegammaname.push_back("Gamma81by69");
  nodename.push_back("S035C34");
  nodetitle[inode]="G(K- h+ pi- pi0 nu(tau) (ex.K0))/G(pi- pi+ pi- pi0 nu(tau) (ex.K0))";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(1.7000E-02);               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  //  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(8.8800E-01);               // h-  omega
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(8.8800E-01);                   // pi- omega
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(8.8800E-01);                   // K-  omega
  ++inode;//56
  nodegammaname.push_back("Gamma93by60");
  nodename.push_back("S035C35");
  nodetitle[inode]="G(K- K+ pi- nu(tau))/G(pi- pi+ pi- nu(tau) (ex.K0))";
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(259);       node_den_coef[inode].push_back(1.        );               
  //  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(1.7000E-02);          // h-  omega     
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(1.7000E-02);              // pi- omega
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(1.7000E-02);              // K-  omega
  ++inode;//57
  nodegammaname.push_back("Gamma94by69");
  nodename.push_back("S035C36");
  nodetitle[inode]="G(K- K+ pi- pi0 nu(tau))/G(pi- pi+ pi- pi0 nu(tau) (ex.K0))";
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(1.7000E-02);               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  //  node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(8.8800E-01);           // h-  omega    
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(8.8800E-01);               // pi- omega
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(8.8800E-01);               // K-  omega
  ++inode;//58
  nodegammaname.push_back("Gamma38");
  nodename.push_back("S035C38");
  nodetitle[inode]="G(K- K0 >=0pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  ++inode;//59
  nodegammaname.push_back("Gamma83");
  nodename.push_back("S035C40"); 
  nodetitle[inode]="G(K- pi+ pi- >=0pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;//60
  nodegammaname.push_back("Gamma110");
  nodename.push_back("S035C47"); 
  nodetitle[inode]="G(X- (S=-1) nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  ++inode;//61
  nodegammaname.push_back("Gamma89");
  nodename.push_back("S035C54"); 
  nodetitle[inode]="G(K- pi+ pi- pi0 nu(tau) (ex.K0,eta))/G(total)";
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );              
  ++inode;//62
  nodegammaname.push_back("Gamma84");
  nodename.push_back("S035C6"); 
  nodetitle[inode]="G(K- pi+ pi- nu(tau))/G(total)";
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  ++inode;//63
  nodegammaname.push_back("Gamma87");
  nodename.push_back("S035C7");
  nodetitle[inode]="G(K- pi+ pi- pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.2600E-01);              
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  ++inode;//64
  nodegammaname.push_back("Gamma94"); 
  nodename.push_back("S035C8");
  nodetitle[inode]="G(K- K+ pi- pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );             
  ++inode;//65
  nodegammaname.push_back("Gamma3"); 
  nodename.push_back("S035R1");
  nodetitle[inode]="G(mu- nubar(mu) nu(tau))/G(total)";
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(1.        );              
  ++inode;//66
  nodegammaname.push_back("Gamma150by66"); 
  nodename.push_back("S035R14");
  nodetitle[inode]="G(h- omega nu(tau))/G(h- h- h+ pi0 nu(tau) (ex.K0))";
  //node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.        );            // h- omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.        );              // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.        );              // K- omega
  node_den_parm[inode].push_back(109);       node_den_coef[inode].push_back(2.2600E-01);               
  node_den_parm[inode].push_back(113);       node_den_coef[inode].push_back(1.7000E-02);               
  node_den_parm[inode].push_back(247);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(263);       node_den_coef[inode].push_back(1.        );               
  node_den_parm[inode].push_back(285);       node_den_coef[inode].push_back(1.        );               
  //node_den_parm[inode].push_back(8  );       node_den_coef[inode].push_back(8.8800E-01);            // h- omega
  node_den_parm[inode].push_back(800);       node_den_coef[inode].push_back(8.8800E-01);              // pi- omega
  node_den_parm[inode].push_back(295);       node_den_coef[inode].push_back(8.8800E-01);              // K- omega
  ++inode;//67
  nodegammaname.push_back("Gamma149"); 
  nodename.push_back("S035R15"); 
  nodetitle[inode]="G(h- omega  >=0 neutrals  nu(tau))/G(total)";
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.        );              
  //node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.        );            // h- omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.        );              // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.        );              // K- omega
  ++inode;//68
  nodegammaname.push_back("Gamma5"); 
  nodename.push_back("S035R2");
  nodetitle[inode]="G(e- nubar(e) nu(tau))/G(total)";
  node_num_parm[inode].push_back(2  );       node_num_coef[inode].push_back(1.        );              
  ++inode;//69
  nodegammaname.push_back("Gamma19"); 
  nodename.push_back("S035R20");
  nodetitle[inode]="G(h- 2pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               
  ++inode;//70
  nodegammaname.push_back("Gamma26");
  nodename.push_back("S035R21");
  nodetitle[inode]="G(h- 3pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2200E-01);              
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  ++inode;//71
  nodegammaname.push_back("Gamma150");
  nodename.push_back("S035R23");
  nodetitle[inode]="G(h- omega nu(tau))/G(total)";
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.        );              // h-  omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.        );                  // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.        );                  // K-  omega
  ++inode;//72
  nodegammaname.push_back("Gamma2");
  nodename.push_back("S035R24");
  nodetitle[inode]="G(particle->=0 neutrals >=0K(L)0 nu(tau) )/G(total)";
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
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.0000E-02);               // h-  omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(9.0000E-02);                   // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(9.0000E-02);                   // K-  omega
  ++inode;//73
  nodegammaname.push_back("Gamma31");
  nodename.push_back("S035R26");
  nodetitle[inode]="G(K- >=0pi0 >=0K0  >=0gamma  nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(7.1500E-01);              
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  ++inode;//74
  nodegammaname.push_back("Gamma32");
  nodename.push_back("S035R27");
  nodetitle[inode]="G(K- >=1 (pi0 or K0 or gamma)  nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(7.1500E-01);              
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );               
  ++inode;//75
  nodegammaname.push_back("Gamma56");
  nodename.push_back("S035R28");
  nodetitle[inode]="G(h- h- h+ nu(tau))/G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(3.4310E-01);              
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(1.7000E-02);               // h-  omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.7000E-02);                   // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.7000E-02);                   // K-  omega
  ++inode;//76
  nodegammaname.push_back("Gamma63");
  nodename.push_back("S035R30");
  nodetitle[inode]="G(h- h- h+ >=1 neutrals  nu(tau))/G(total)";
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
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(8.8800E-01);               // h-  omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(8.8800E-01);                   // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(8.8800E-01);                   // K-  omega
  ++inode;//77
  nodegammaname.push_back("Gamma54");
  nodename.push_back("S035R31");
  nodetitle[inode]="G(h- h- h+ >=0 neutrals >=0K(L)0 nu(tau) )/G(total)";
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
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.1010E-01);      // h-  omega         
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(9.1010E-01);          // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(9.1010E-01);          // K-  omega
  ++inode;//78
  nodegammaname.push_back("Gamma126");
  nodename.push_back("S035R32");
  nodetitle[inode]="G(eta pi- pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(1.        );              
  ++inode;//79
  nodegammaname.push_back("Gamma102");
  nodename.push_back("S035R33");
  nodetitle[inode]="G( 3h- 2h+ >=0 neutrals  nu(tau)  (ex. K(S)0 --> pi- pi+)(``5-prong''))/G(total)";
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(4  );       node_num_coef[inode].push_back(1.        );               
  ++inode;//80
  nodegammaname.push_back("Gamma79");
  nodename.push_back("S035R34");
  nodetitle[inode]="G(K- h+ h- >=0 neutrals  nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.8500E-01);              
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  ++inode;//81
  nodegammaname.push_back("Gamma103");
  nodename.push_back("S035R38");
  nodetitle[inode]="G(3h- 2h+ nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(3  );       node_num_coef[inode].push_back(1.        );              
  ++inode;//82
  nodegammaname.push_back("Gamma104");
  nodename.push_back("S035R39");
  nodetitle[inode]="G(3h- 2h+ pi0 nu(tau) (ex.K0))/G(total)";
  node_num_parm[inode].push_back(4  );       node_num_coef[inode].push_back(1.        );              
  ++inode;//83
  nodegammaname.push_back("Gamma93");
  nodename.push_back("S035R40");
  nodetitle[inode]="G(K- K+ pi- nu(tau))/G(total)";
  node_num_parm[inode].push_back(5  );       node_num_coef[inode].push_back(1.        );              
  ++inode;//84
  nodegammaname.push_back("Gamma82");
  nodename.push_back("S035R41");
  nodetitle[inode]="G(K- pi+ pi- >=0 neutrals  nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(2.8500E-01);              
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(3.4310E-01);               
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(3.4310E-01);               
  ++inode;//85
  nodegammaname.push_back("Gamma11");
  nodename.push_back("S035R42");
  nodetitle[inode]="G(h- >=1 neutrals nu(tau))/G(total)";
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
  //  node_num_parm[inode].push_back(8  );       node_num_coef[inode].push_back(9.0000E-02);               // h-  omega
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(9.0000E-02);                   // pi- omega
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(9.0000E-02);                   // K-  omega
  ++inode;//86
  nodegammaname.push_back("Gamma7");
  nodename.push_back("S035R43");
  nodetitle[inode]="G(h- >=0K(L)0  nu(tau))/G(total)";
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(5.0000E-01);              
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(5.0000E-01);               
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );               
  ++inode;//87
  nodegammaname.push_back("Gamma17");
  nodename.push_back("S035R44");
  nodetitle[inode]="G(h- >= 2pi0 nu(tau))/G(total)";
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
  ++inode;//88
  nodegammaname.push_back("Gamma37");
  nodename.push_back("S035R46");
  nodetitle[inode]="G(K- K0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(62 );       node_num_coef[inode].push_back(1.        );              
  ++inode;//89
  nodegammaname.push_back("Gamma3by5");
  nodename.push_back("S035R5");
  nodetitle[inode]="G(mu- nubar(mu) nu(tau))/G(e- nubar(e) nu(tau))";
  node_num_parm[inode].push_back(1  );       node_num_coef[inode].push_back(1.        );              
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(1.        );               
  ++inode;//90
  nodegammaname.push_back("Gamma9");
  nodename.push_back("S035R6");
  nodetitle[inode]="G(pi- nu(tau))/G(total)";
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );              
  ++inode;//91
  nodegammaname.push_back("Gamma10");
  nodename.push_back("S035R7");
  nodetitle[inode]="G(K- nu(tau))/G(total)";
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );              
  ++inode;//92
  nodegammaname.push_back("Gamma14");
  nodename.push_back("S035R8");
  nodetitle[inode]="G(pi- pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );              
  ++inode;//93
  nodegammaname.push_back("Gamma13"); 
  nodename.push_back("S035R84"); 
  nodetitle[inode]="G(h- pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(16 );       node_num_coef[inode].push_back(1.        );              
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               
  ++inode;//94
  nodegammaname.push_back("Gamma24");
  nodename.push_back("S035R97");
  nodetitle[inode]="G(h- >= 3pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(3.2200E-01);              
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.5700E-01);               
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(9.8500E-02);               
  node_num_parm[inode].push_back(58 );       node_num_coef[inode].push_back(3.1900E-01);               
  ++inode;//95
  nodegammaname.push_back("Gamma130");
  nodename.push_back("S035C27");
  nodetitle[inode]="G(eta K- pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(1.        );
  ++inode;//96
  nodegammaname.push_back("Gamma132");
  nodename.push_back("S035C28");
  nodetitle[inode]="G(eta pi- K0bar nu(tau))/G(total)";
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(1.        );
  ++inode;//97
  nodegammaname.push_back("Gamma43");
  nodename.push_back("S035C37");
  nodetitle[inode]="G(pi- K0bar >= 1pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(1.        );
  ++inode;//98
  nodegammaname.push_back("Gamma44");
  nodename.push_back("S035B98");
  nodetitle[inode]="G(pi- K0bar 2pi0 nu(tau))/G(total)";
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(1.        );
  ++inode;//99
  nodegammaname.push_back("Gamma53");
  nodename.push_back("S035C5");
  nodetitle[inode]="G(K0bar h+ h- h- nu(tau))/G(total)";
  node_num_parm[inode].push_back(244);       node_num_coef[inode].push_back(1.        );
  ++inode;//100
  nodegammaname.push_back("Gamma800");
  nodename.push_back("S035Z01");
  nodetitle[inode]="G(pi- omega nu(tau))/G(total)";
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.        );
  ++inode;//101
  nodegammaname.push_back("Gamma151");
  nodename.push_back("S035C61");
  nodetitle[inode]="G(K- omega nu(tau))/G(total)";
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.        );
  ++inode;//102
  nodegammaname.push_back("Gamma96");
  nodename.push_back("S035C9");
  nodetitle[inode]="G(K- K+ K- nu(tau))/G(total)";
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(0.588447); // 1/1.699387
  ++inode;//103
  nodegammaname.push_back("Gamma801");
  nodename.push_back("S035Z02");
  nodetitle[inode]="G(K- phi nu(tau) [phi->KK])/G(total)";
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(1.        );
  ++inode;//104
  nodegammaname.push_back("Gamma9by5");
  nodename.push_back("S035Z03");
  nodetitle[inode]="G(pi- nu(tau))//G(e- nubar(e) nu(tau))";
  node_num_parm[inode].push_back(12 );       node_num_coef[inode].push_back(1.        );
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(1.        );               
  ++inode;//105
  nodegammaname.push_back("Gamma10by5");
  nodename.push_back("S035Z04");
  nodetitle[inode]="G(K- nu(tau))//G(e- nubar(e) nu(tau))";
  node_num_parm[inode].push_back(7  );       node_num_coef[inode].push_back(1.        );
  node_den_parm[inode].push_back(2  );       node_den_coef[inode].push_back(1.        );               
  ++inode;//106
  nodegammaname.push_back("GammaAll"); // sum of 37 base nodes
  nodename.push_back("S035S01");
  nodetitle[inode]="G(total)";
  node_num_parm[inode].push_back(  1);       node_num_coef[inode].push_back(1.        );               //1
  node_num_parm[inode].push_back(  2);       node_num_coef[inode].push_back(1.        );               //2
  node_num_parm[inode].push_back( 12);       node_num_coef[inode].push_back(1.        );               //3
  node_num_parm[inode].push_back(  7);       node_num_coef[inode].push_back(1.        );               //4
  node_num_parm[inode].push_back( 16);       node_num_coef[inode].push_back(1.        );               //5
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               //6
  node_num_parm[inode].push_back(201);       node_num_coef[inode].push_back(1.        );               //7
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               //8
  node_num_parm[inode].push_back(203);       node_num_coef[inode].push_back(1.        );               //9
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               //10
  node_num_parm[inode].push_back(110);       node_num_coef[inode].push_back(1.        );               //11
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );               //12
  node_num_parm[inode].push_back( 62);       node_num_coef[inode].push_back(1.        );               //13
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );               //14
  node_num_parm[inode].push_back(119);       node_num_coef[inode].push_back(1.        );               //15
  node_num_parm[inode].push_back(214);       node_num_coef[inode].push_back(1.        );               //16
  node_num_parm[inode].push_back(240);       node_num_coef[inode].push_back(1.        );               //17
  node_num_parm[inode].push_back(259);       node_num_coef[inode].push_back(1.        );               //18
  node_num_parm[inode].push_back(263);       node_num_coef[inode].push_back(1.        );               //19
  node_num_parm[inode].push_back(216);       node_num_coef[inode].push_back(1.        );               //20
  node_num_parm[inode].push_back(204);       node_num_coef[inode].push_back(1.        );               //21
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               //22
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               //23
  node_num_parm[inode].push_back(  5);       node_num_coef[inode].push_back(1.        );               //24
  node_num_parm[inode].push_back(247);       node_num_coef[inode].push_back(1.        );               //25
  node_num_parm[inode].push_back(  4);       node_num_coef[inode].push_back(1.        );               //26
  node_num_parm[inode].push_back( 58);       node_num_coef[inode].push_back(1.        );               //27
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(1.        );               //28
  node_num_parm[inode].push_back(800);       node_num_coef[inode].push_back(1.        );               //29
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.        );               //30
  node_num_parm[inode].push_back(113);       node_num_coef[inode].push_back(1.        );               //31
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(1.        );               //32
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(1.        );               //33
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(1.        );               //34
  node_num_parm[inode].push_back(244);       node_num_coef[inode].push_back(1.        );               //35
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(1.        );               //36
  node_num_parm[inode].push_back(  3);       node_num_coef[inode].push_back(1.        );               //37
  ++inode;//107
  nodegammaname.push_back("GammaXs"); // sum of 15 strange base nodes
  nodename.push_back("S035S02");
  nodetitle[inode]="G(strange)/G(total)";
  node_num_parm[inode].push_back(  7);       node_num_coef[inode].push_back(1.        );               //1       
  node_num_parm[inode].push_back(182);       node_num_coef[inode].push_back(1.        );               //2    
  node_num_parm[inode].push_back(115);       node_num_coef[inode].push_back(1.        );               //3 
  node_num_parm[inode].push_back(116);       node_num_coef[inode].push_back(1.        );               //4
  node_num_parm[inode].push_back(117);       node_num_coef[inode].push_back(1.        );               //5
  node_num_parm[inode].push_back(118);       node_num_coef[inode].push_back(1.        );               //6
  node_num_parm[inode].push_back(260);       node_num_coef[inode].push_back(1.        );               //7
  node_num_parm[inode].push_back(285);       node_num_coef[inode].push_back(1.        );               //8
  node_num_parm[inode].push_back(109);       node_num_coef[inode].push_back(1.        );               //9
  node_num_parm[inode].push_back(295);       node_num_coef[inode].push_back(1.        );               //10
  node_num_parm[inode].push_back(266);       node_num_coef[inode].push_back(1.        );               //11
  node_num_parm[inode].push_back(267);       node_num_coef[inode].push_back(1.        );               //12
  node_num_parm[inode].push_back(238);       node_num_coef[inode].push_back(1.        );               //13
  node_num_parm[inode].push_back(244);       node_num_coef[inode].push_back(1.        );               //14
  node_num_parm[inode].push_back(801);       node_num_coef[inode].push_back(1.        );               //15
  ++inode;//108
  //
  // Count number of parameters in numerator and denominator for each node
  //
  int            node_num_npar[nnode];
  int            node_den_npar[nnode];
  for (inode=0;inode<nnode;++inode){
    node_num_npar[inode]=node_num_parm[inode].size();
    node_den_npar[inode]=node_den_parm[inode].size();
  }
  //
  // Define array of a_nodename to simplify passing the strings to a function
  // 
  char** a_nodename = new char*[nnode];  
  for (inode=0;inode<nnode;++inode) {
    a_nodename[inode] = new char[10];
    a_nodename[inode] = (char*)(nodename[inode].data()); 
    //    cout << "inode = " << inode << " nodename = " << nodename[inode] << " a_nodename = " << a_nodename[inode] << endl;
  }
  // PRINT NODE DEFINITION 
  //
  print_node_def(nnode, a_nodename, nodetitle, nodegammaname, 
		 node_num_npar, node_num_parm, node_num_coef, 
		 node_den_npar, node_den_parm, node_den_coef, 
		 baseparm, basegamma);
  //
  // DEFINE NODE DEPENDENCE ON BASE
  //
  vector<int> node_parm[nnode]; // vector of parameters for each node
  vector<int> node_quan[nnode]; // vector of quantities for each node
  int first_quan[nnode];  // first quantity measured for each node
  define_nodes(nnode, node_num_parm, node_den_parm, baseparm,
	       node_parm, // output
	       node_quan, // output
	       first_quan // output
	       );
  //
  // READ MEASUREMENT FILE
  //
  int nmeas=0;
  int nmeas_noweak;
  int measnum[200], measnode[200], imeas1, imeas2;
  string measgammaname[200], measnodename[200], expname[200], author[200], year[200], meastitle[200];
  double measvalue[200], measerror[200], corrtemp;
  double** corrmat = new double*[200]; for (i=0;i<200;++i) corrmat[i] = new double[200];
  for (imeas1=0;imeas1<200;imeas1++) for (imeas2=0;imeas2<200;imeas2++) corrmat[imeas1][imeas2] = 0;
  //
  cout << "Read data from : " << Form("s035-fit-with-babar-belle%s.data",salephhcorr.data()) << endl;
  ifstream ifs(Form("s035-fit-with-babar-belle%s.data",salephhcorr.data())) ;
  if (!ifs.good()) {cout << "Cannot open input file : " << Form("s035-fit-with-babar-belle%s.data",salephhcorr.data()) << endl ; exit(1) ;}
  while(ifs.good()) {
    if (ifs.eof()) break;
    char firstch(' ') ; ifs.get(firstch) ;
    if (firstch=='#'||firstch=='\n') { // Skip such lines
      ifs.ignore(256,'\n') ;
    } else if (firstch=='*') {  // measurement line
      ifs >> measnum[nmeas] >> measgammaname[nmeas] >> measnodename[nmeas]
	  >> measvalue[nmeas] >> measerror[nmeas] 
	  >> expname[nmeas] >> author[nmeas] >> year[nmeas];
      //
      vector<string>::iterator it=find(nodename.begin(),nodename.end(),measnodename[nmeas]);      
      inode=it-nodename.begin();
      measnode[nmeas]=inode;
      //
      ifs.getline(buffer,256,'\n');
      string stringbuffer=string(buffer);
      meastitle[nmeas]="";
      bool first=false;
      for (string::iterator si=stringbuffer.begin();si!=stringbuffer.end();++si){
	if (13 == (int)*si) break; // ^M is 13
	if (' ' != *si) first=true;
	if (first) meastitle[nmeas]+=*si;
      }
      //      cout << measnum[nmeas] << " " << measgammaname[nmeas] << " " << measnodename[nmeas] << " " 
      //	   << measvalue[nmeas] << " " << measerror[nmeas] << " " 
      //	   << expname[nmeas] << " " << author[nmeas] << " " << year[nmeas] << " " << meastitle[nmeas] << endl;
      ++nmeas;
    } else if (firstch=='%') {  // correlation line
      ifs >> imeas1 >> imeas2 >> corrtemp;       ifs.ignore(256,'\n') ;
      if (corrmat[imeas1-1][imeas2-1] != 0) {cout << "WATCH OUT" << endl; exit(1);}
      corrmat[imeas1-1][imeas2-1] = corrmat[imeas2-1][imeas1-1] = corrtemp*1e-2;
      //      cout << "Correlation between imeas1 = " << measnum[imeas1-1] << " and imeas2 = " << measnum[imeas2-1] << " is = " << corrmat[imeas1-1][imeas2-1] << endl;
    }
  }
  //
  // Process Measurements
  //
  TMatrixD MeasErrorMatrix(nmeas,nmeas);
  TMatrixD InvMeasErrorMatrix(nmeas,nmeas);
  vector<int> icorrj[200]; // vector of measurements correlated with this measurement
  int ncorrij; // number of cycles of correlated measurements
  int ifirstj[200]; // first measurement in this cycle
  vector<int> veccorrij[200]; // vector of correlated measurements per cycle
  int newnode_all;
  int nodegroup_all[200]; // mapping of each measurement to groups of nodes
  int ncycle[200]; // mapping of each measurement to the cycle number
  int n_nodegroup[200]; // of newnode_all
  int ngroup[200]; // of measurements
  //
  process_measurements(nmeas, measnode, measerror, corrmat,
		       MeasErrorMatrix,  // output
		       InvMeasErrorMatrix, //output
		       icorrj, //output
		       ncorrij, //output
		       ifirstj, //output
		       veccorrij, //output
		       newnode_all, //output
		       nodegroup_all, //output
		       ncycle, // output
		       n_nodegroup, // output
		       ngroup //output
		       );
  //
  // PRINT INFORMATION ABOUT INPUT MEASUREMENTS
  //
  FILE *measinfofile=fopen(Form("s035-fit-with-babar-belle%s.info",salephhcorr.data()),"w");  
  print_measinfo(measinfofile,
		 nmeas, measnode, measvalue, measerror, corrmat,
		 expname, meastitle, measgammaname,
		 a_nodename, node_num_npar, node_den_npar, node_quan,
		 nbase, baseparm, basegamma, basenode, basetitle,
		 newnode_all, nodegroup_all, ncycle, ncorrij, ifirstj, icorrj, veccorrij, ngroup);
  fclose(measinfofile);
  //
  // Derivative of R [node] as a linearized sum of variables P_i [base quantities]
  //
  double node_num[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double node_den[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double node_val[nnode]; // value of each node
  vector<double> node_part[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- may change
		   node_num, // output
		   node_den, // output
		   node_val, // output
		   node_part // output
		   );
  //
  //  Prepare to obtain Fit Values
  //
  double basevalue_fit[nbase];
  double baseerror_fit[nbase];
  double** basecov_fit = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit[i] = new double[nbase];
  double** basecorr_fit = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit[i] = new double[nbase];
  double chisquared;
  int weak_none[200];  for (i=0;i<nmeas;++i) weak_none[i]=0; // all measurements are treated to be non-weak here
  //
  // COMBINE 
  //
  combine(uconstrain,
	  nmeas, weak_none, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix,
	  node_num_npar, node_den_npar, node_quan, node_parm, 
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit, // output
	  baseerror_fit, // output
	  basecov_fit,   // output
	  basecorr_fit,  // output
	  chisquared     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  // 
  double chisquared_temp = chisquared;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    // COMBINE 
    combine(uconstrain,
	    nmeas, weak_none, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix,
	    node_num_npar, node_den_npar, node_quan, node_parm, 
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit, // output
	    baseerror_fit, // output
	    basecov_fit,   // output
	    basecorr_fit,  // output
	    chisquared     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
		 chisquared, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared - chisquared_temp) < 1e-3) break;
    chisquared_temp = chisquared;
  }
  //
  cout << endl << "Results from original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE           SEED   FITVAL     FITERR   TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6e %10.6f %10.6f %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),baseseed[ibase],
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue[nnode];     // value of each node
  vector<double> NodeValue_part[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num, // output
		   NodeValue_den, // output
		   NodeValue,     // output
		   NodeValue_part // output
		   );
  //
  double NodeError[nnode];
  TMatrixD NodeErrorMatrix(nnode,nnode);
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part[jnode].at(jpar);
	  NodeErrorMatrix[inode][jnode] += 
	    ipartial*jpartial*baseerror_fit[iquan-1]*baseerror_fit[jquan-1]*basecorr_fit[iquan-1][jquan-1];
	}
      }
    }
    //
    NodeError[inode] = TMath::Sqrt(NodeErrorMatrix[inode][inode]);
  }
  //
  int   nchisquare_tot=0; 
  double chisquare_tot=0; 
  double chisquare_meas[200]; // chi2 per mesurement
  int   nchisquare_meas[200];  for (im=0;im<nmeas;++im) nchisquare_meas[im]=0;// nodes per mesurement     
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_none[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak_none[jm]==weakcompare) continue;
	  chisquare_meas[im]+= 
	    (measvalue[im]-NodeValue[measnode[im]])*InvMeasErrorMatrix[im][jm]*(measvalue[jm]-NodeValue[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak_none[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot+=1;
	chisquare_tot+=chisquare_node;
	nchisquare_meas[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror[im], 
		  NodeValue[measnode[im]], NodeError[measnode[im]],
		  chisquare_meas[im], chisquare_node, chisquare_tot,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d \n\n", chisquare_tot, nchisquare_tot);
  //
  cout << Form("%s = %10.6f +- %10.6f (%10.6f sigma); %s = %10.6f +- %10.6f (%10.6f sigma)\n\n",
	       nodetitle[nnode-2].data(),
	       NodeValue[nnode-2],
	       uconstrain?0:NodeError[nnode-2],
	       uconstrain?0:(1-NodeValue[nnode-2])/NodeError[nnode-2],
	       nodetitle[nnode-1].data(),
	       NodeValue[nnode-1],
	       NodeError[nnode-1],
	       NodeValue[nnode-1]/NodeError[nnode-1]);
  //
  // WRITE OUT MEASUREMENT FILE
  //
  FILE *measfile[2];
  for (int p=0;p<2;++p){
    if (p==0) measfile[0]=fopen(Form("combos_measurements_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) measfile[1]=fopen(Form("alucomb_measurements_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    FILE *thisfile = measfile[p];
    print_measfile(thisfile, p, uconstrain,
		   nmeas, weak_none, weakcompare, measgammaname, measnode, 
		   expname, author, year, meastitle,
		   measvalue, measerror, corrmat,
		   a_nodename, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef, 
		   nbase, baseparm, basegamma, basetitle, first_quan, baseseed);
    fclose(thisfile);
  }
  //
  // WRITE OUT AVERAGE.INPUT FILE FOR COMBOS/ALUCOMB
  //
  FILE *avefile[2];
  for (int p=0;p<2;++p){
    if (p==0) avefile[0]=fopen(Form("combos_average_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) avefile[1]=fopen(Form("alucomb_average_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==0) fprintf (avefile[p], "INCLUDE combos_measurements_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    if (p==1) fprintf (avefile[p], "INCLUDE alucomb_measurements_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    FILE *thisfile = avefile[p];
    print_avefile(thisfile, p, uconstrain,
		  nmeas, weak_none, weakcompare, measgammaname, measnode,
		  node_parm, node_quan, 
		  node_num_npar, node_num_parm, node_num_coef,
		  node_den_npar, node_den_parm, node_den_coef,
		  nbase, baseparm, basegamma, basetitle, first_quan, 
		  baseseed, node_num, node_den, node_part);
    fclose(thisfile);
  }
  //
  //  Prepare to obtain Fit Values for non-weak measurements
  //
  double basevalue_fit_noweak[nbase];
  double baseerror_fit_noweak[nbase];
  double** basecov_fit_noweak = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit_noweak[i] = new double[nbase];
  double** basecorr_fit_noweak = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit_noweak[i] = new double[nbase];
  double chisquared_noweak;
  //
  // Define weak nodes
  //
  double nsig[200];
  int weak[200];
  for (i=0;i<nmeas;++i){
    if (nchisquare_meas[i] > 0 ) {
      nsig[i] = measerror[i]/((sqrt(nchisquare_meas[i]))*NodeError[measnode[i]]);
    }else{
      nsig[i] = 999999;
    }
    weak[i] = (icorrj[i].size()==1) && (nsig[i] > 3.); // only for uncorrelated measurements
  }
  //
  // COMBINE 
  //
  combine(uconstrain,
	  nmeas, weak, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix,
	  node_num_npar, node_den_npar, node_quan, node_parm, 
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit_noweak, // output
	  baseerror_fit_noweak, // output
	  basecov_fit_noweak,   // output
	  basecorr_fit_noweak,  // output
	  chisquared_noweak     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_noweak, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  // 
  double chisquared_noweak_temp = chisquared_noweak;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit_noweak[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    //
    // COMBINE 
    //
    combine(uconstrain,
	    nmeas, weak, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix,
	    node_num_npar, node_den_npar, node_quan, node_parm, 
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit_noweak, // output
	    baseerror_fit_noweak, // output
	    basecov_fit_noweak,   // output
	    basecorr_fit_noweak,  // output
	    chisquared_noweak     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
		 chisquared_noweak, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared_noweak - chisquared_noweak_temp) < 1e-3) break;
    chisquared_noweak_temp = chisquared_noweak;
  }
  //
  cout << endl << "Comparison of Results from fit with [non-weak measurements only] w.r.t original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE  ORIG_FITVAL ORIG_FITERR SCAL_FITVAL SCAL_FITERR   SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6f  %10.6f  %10.6f  %10.6f %6.2f  %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basevalue_fit_noweak[ibase],
		 baseerror_fit_noweak[ibase],
		 baseerror_fit_noweak[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit_noweak[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit_noweak[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num_noweak[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den_noweak[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue_noweak[nnode];     // value of each node
  vector<double> NodeValue_part_noweak[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num_noweak, // output
		   NodeValue_den_noweak, // output
		   NodeValue_noweak,     // output
		   NodeValue_part_noweak // output
		   );
  //
  double NodeError_noweak[nnode];
  TMatrixD NodeErrorMatrix_noweak(nnode,nnode);
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix_noweak[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part_noweak[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part_noweak[jnode].at(jpar);
	  NodeErrorMatrix_noweak[inode][jnode] += 
	    ipartial*jpartial*baseerror_fit_noweak[iquan-1]*baseerror_fit_noweak[jquan-1]*basecorr_fit_noweak[iquan-1][jquan-1];
	}
      }
    }
    //
    NodeError_noweak[inode] = TMath::Sqrt(NodeErrorMatrix_noweak[inode][inode]);
  }
  //
  int   nchisquare_tot_noweak=0; 
  double chisquare_tot_noweak=0; 
  double chisquare_meas_noweak[200]; // chi2 per mesurement
  int   nchisquare_meas_noweak[200]; for (im=0;im<nmeas;++im) nchisquare_meas_noweak[im]=0;// nodes per mesurement
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas_noweak[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak[jm]==weakcompare) continue;
	  chisquare_meas_noweak[im]+= 
	    (measvalue[im]-NodeValue_noweak[measnode[im]])*InvMeasErrorMatrix[im][jm]*(measvalue[jm]-NodeValue_noweak[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas_noweak[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot_noweak+=1;
	chisquare_tot_noweak+=chisquare_node;
	nchisquare_meas_noweak[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror[im], 
		  NodeValue_noweak[measnode[im]], NodeError_noweak[measnode[im]],
		  chisquare_meas_noweak[im], chisquare_node, chisquare_tot_noweak,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d \n\n", chisquare_tot_noweak, nchisquare_tot_noweak);
  //
  cout << Form("%s = %10.6f +- %10.6f (%10.6f sigma); %s = %10.6f +- %10.6f (%10.6f sigma)\n\n",
	       nodetitle[nnode-2].data(),
	       NodeValue_noweak[nnode-2],
	       uconstrain?0:NodeError_noweak[nnode-2],
	       uconstrain?0:(1-NodeValue_noweak[nnode-2])/NodeError_noweak[nnode-2],
	       nodetitle[nnode-1].data(),
	       NodeValue_noweak[nnode-1],
	       NodeError_noweak[nnode-1],
	       NodeValue_noweak[nnode-1]/NodeError_noweak[nnode-1]);
  //
  cout << "Summary of fit with [non-weak measurements only]:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE           SEED    FITVAL     FITERR   RESCALED    SFAC TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6e %10.6f %10.6f %10.6f %6.2f %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),baseseed[ibase],
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 baseerror_fit_noweak[ibase],
		 baseerror_fit_noweak[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  // PDG-style Scale Factor Calculation
  //
  printf ("ncorrij = %d \n",ncorrij);
  for (i=0;i<ncorrij;++i){
    printf ("i = %d ifirstj[i]+1 = %3d ",i,ifirstj[i]+1);
    printf ("veccorrij[i].size() = %2d veccorrij[i][j]+1 = ",veccorrij[i].size());
    for (j=0;j<veccorrij[i].size();++j) printf ("%d ",veccorrij[i][j]+1);
    printf ("GammaName = "); 
    for (j=0;j<veccorrij[i].size();++j) printf ("%s ",measgammaname[veccorrij[i][j]].data());
    printf ("\n");
  }
  //
  vector<double> VectorOfPullMag;
  vector<TMatrixD> VectorOfMeasScaledErrorMatrix;
  vector<double> VectorOfPullMag_noweak;
  vector<TMatrixD> VectorOfMeasScaledErrorMatrix_noweak;
  for (int iweak=0;iweak<2;++iweak) {
    for (int ij=0;ij<ncorrij;++ij) {
      int nij=veccorrij[ij].size();
      TMatrixD ThisFitErrorMatrix(nij,nij);
      TMatrixD ThisMeasErrorMatrix(nij,nij);
      for (i=0;i<nij;++i){
	im=veccorrij[ij][i];
	inode=measnode[im];
	for (j=0;j<nij;++j){
	  jm=veccorrij[ij][j];
	  jnode=measnode[jm];
	  ThisMeasErrorMatrix[i][j] = (i==j) ? measerror[im]*measerror[im] : corrmat[im][jm] * measerror[im] * measerror[jm];
	  ThisFitErrorMatrix[i][j]  = (iweak==0) ? NodeErrorMatrix[inode][jnode] : NodeErrorMatrix_noweak[inode][jnode];
	  //	  if (iweak==1&&(ij==1||ij==7)) {
	  //	    cout << i << " " << j << " " << inode << " " << jnode << " " << nodegammaname[inode] << " " << nodegammaname[jnode] << " "
	  //		 << "NodeErrorMatrix: " << NodeErrorMatrix[inode][jnode] << " " << NodeErrorMatrix_noweak[inode][jnode] << " "
	  //		 << NodeErrorMatrix[inode][jnode]/NodeErrorMatrix_noweak[inode][jnode] << " "
	  //		 << endl;
	}
      }
      //
      // CMatrix = (MeasErrorMatrix - FitErrorMatrix)^{-1}
      //
      TMatrixD InvCMatrix(nij,nij);
      for (i=0;i<nij;++i) {
	for (j=0;j<nij;++j) {
	  InvCMatrix[i][j] = ThisMeasErrorMatrix[i][j] - ThisFitErrorMatrix[i][j];
	}
      }
      TMatrixD CMatrix = InvCMatrix;
      Double_t c_det;
      CMatrix.Invert(&c_det);
      //
      TVectorD EigenVal;
      TMatrixD EigenVecT = CMatrix.EigenVectors(EigenVal);
      TMatrixD EigenVec(nij,nij); EigenVec.Transpose(EigenVecT);
      cout << "EigenVal : "; for (i=0;i<nij;++i) cout << EigenVal[i] << " ";  cout << endl;
      bool EigenVal_AllPositive=true; for (i=0;i<nij;++i) if (EigenVal[i]<0) EigenVal_AllPositive=false;
      //
      double PullMag=0;
      double PullMag2=0;
      double PullVector[nij];
      for (i=0;i<nij;++i) {
	PullVector[i] = 0;
      }
      TMatrixD ThisMeasScaledErrorMatrix(nij,nij);
      //
      if (EigenVal_AllPositive) {
	TMatrixD DiagonalMatrix(nij,nij);
	TMatrixD DiagonalMatrix_sqrt(nij,nij);
	for (i=0;i<nij;++i) {
	  for (j=0;j<nij;++j) {
	    if (i==j) {
	      DiagonalMatrix(i,j) = EigenVal[i];
	      DiagonalMatrix_sqrt(i,j) = TMath::Sqrt(EigenVal[i]);
	    } else {
	      DiagonalMatrix(i,j) = 0;
	      DiagonalMatrix_sqrt(i,j) = 0;
	    }
	  }
	}
	//
	TMatrixD QMatrix = EigenVecT * DiagonalMatrix_sqrt * EigenVec; // Q = sqrt(C)
	//
	double FMinusM[nij];
	for (i=0;i<nij;++i) {
	  im=veccorrij[ij][i];
	  FMinusM[i] = (iweak==0) ? NodeValue[measnode[im]]-measvalue[im] : NodeValue_noweak[measnode[im]]-measvalue[im]; 
	}
	//
	for (i=0;i<nij;++i) {
	  for (j=0;j<nij;++j) {
	    PullVector[i] += QMatrix[i][j] * FMinusM[j];
	  }
	  PullMag2 += PullVector[i]*PullVector[i];
	}
	PullMag = TMath::Sqrt(PullMag2);
	//
	double UnitPullVector[nij];
	for (i=0;i<nij;++i) {
	  UnitPullVector[i] = (PullMag>0)? PullVector[i]/PullMag : 0;
	}
	TMatrixD ThisMeasCorrMatrix(nij,nij);
	for (i=0;i<nij;++i) {
	  for (j=0;j<nij;++j) {
	    if (i==j) {
	      ThisMeasCorrMatrix[i][j] = 1; // by construction diagonal elements of corrmat has not been filled
	    } else {
	      im=veccorrij[ij][i];
	      jm=veccorrij[ij][j];
	      ThisMeasCorrMatrix[i][j]=corrmat[im][jm];
	    }
	  }
	}
	//
	double CM=0;
	for (i=0;i<nij;++i) {
	  for (j=0;j<nij;++j) {
	    CM += UnitPullVector[i]*ThisMeasCorrMatrix[i][j]*UnitPullVector[j];
	  }
	}
	double SFactor = (PullMag>1) ? (PullMag2 - 1) * CM : 0; 
	//
	for (i=0;i<nij;++i) {
	  im=veccorrij[ij][i];
	  for (j=0;j<nij;++j) {
	    jm=veccorrij[ij][j];
	    ThisMeasScaledErrorMatrix[i][j] = ThisMeasErrorMatrix[i][j] + UnitPullVector[i] * measerror[im] * UnitPullVector[j] * measerror[jm] * SFactor;
	  }
	}
      }
      //
      if (iweak==0) {
	VectorOfPullMag.push_back(PullMag);
	VectorOfMeasScaledErrorMatrix.push_back(ThisMeasScaledErrorMatrix);
      }else{
	VectorOfPullMag_noweak.push_back(PullMag);
	VectorOfMeasScaledErrorMatrix_noweak.push_back(ThisMeasScaledErrorMatrix);
      }
      //
      cout << "iweak = " << iweak << " ij = " << ij << " PullMag = " << PullMag << " PullVector : "; for (i=0;i<nij;++i) cout << PullVector[i] << " "; cout << endl;
      //
    } //end of loop over ij
  } // end of loop over iweak
  //
  // With all measurements
  //
  double pullav_meas[200]; for (im=0;im<nmeas;++im) pullav_meas[im]=0; // pullav per measurement
  for (inode=0;inode<=newnode_all;++inode){
    int   npullsq_node=0;
    double pullsq_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_none[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	double pullsquare=0;
	if ( icorrj[im].size() > 1 ) { // first find which cycle this correlated measurement belongs to, then just fetch the result
	  pullsquare = TMath::Power(VectorOfPullMag[ncycle[im]],2);
	} else {
	  pullsquare = TMath::Power((measvalue[im]-NodeValue[measnode[im]]),2)/
	    (TMath::Power(measerror[im],2) - TMath::Power(NodeError[measnode[im]],2));
	}
	pullsq_node+=pullsquare;
	npullsq_node+=1;
      }
    }
    pullsq_node/=npullsq_node;
    //
    for (im=0;im<nmeas;++im){
      if (weak_none[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	pullav_meas[im]=TMath::Sqrt(pullsq_node);
	//cout << "inode = " << inode << " im+1 = " << im+1 << " npullsq_node = " << npullsq_node << " pullav_meas = " << pullav_meas[im] << endl;
      }
    }
  }
  //
  double measerror_scaled[200];
  for (i=0;i<nmeas;++i){
    if (icorrj[i].size() > 1 && VectorOfPullMag[ncycle[i]]>1) {
      double temperr2=measerror[i]*measerror[i]; // over-written below
      for (int itemp=0; itemp < veccorrij[ncycle[i]].size(); ++itemp) {
	if (i == veccorrij[ncycle[i]][itemp]) {
	  temperr2=((TMatrixD)(VectorOfMeasScaledErrorMatrix[ncycle[i]]))(itemp,itemp);
	  break;
	}
      }
      measerror_scaled[i] = TMath::Sqrt(temperr2);
    }else if (icorrj[i].size()==1 && pullav_meas[i]>1 ) { 
      measerror_scaled[i] = measerror[i]*pullav_meas[i];
    }else{
      measerror_scaled[i] = measerror[i];
    }
  }
  //
  double** corrmat_scaled = new double*[200]; for (i=0;i<200;++i) corrmat_scaled[i] = new double[200];
  for (int i=0;i<nmeas;++i){
    for (int j=0;j<nmeas;++j) {
      double tempcorrij=corrmat[i][j]; // over-written below
      if (icorrj[i].size() > 1 && VectorOfPullMag[ncycle[i]]>1) {
	for (int itemp=0; itemp < veccorrij[ncycle[i]].size(); ++itemp) {
	  if (i == veccorrij[ncycle[i]][itemp]) {
	    for (int jtemp=0; jtemp < veccorrij[ncycle[i]].size(); ++jtemp) {
	      if (j == veccorrij[ncycle[i]][jtemp]) {
		double temperr_i=((TMatrixD)(VectorOfMeasScaledErrorMatrix[ncycle[i]]))(itemp,itemp);
		double temperr_j=((TMatrixD)(VectorOfMeasScaledErrorMatrix[ncycle[i]]))(jtemp,jtemp);
		double temperr_ij=((TMatrixD)(VectorOfMeasScaledErrorMatrix[ncycle[i]]))(itemp,jtemp);
		tempcorrij = temperr_ij / TMath::Sqrt(temperr_i * temperr_j);
		break;
	      }
	    }
	    break;
	  }
	}
      }
      corrmat_scaled[i][j] = i==j ? 0 : tempcorrij;
      // cout << " i+1 = " << i+1 << " j+1 = " << j+1 
      //   << " icorrij = " << icorrj[i].size() << " pullav_meas = " << pullav_meas[i] << " PullMag = " << VectorOfPullMag[ncycle[i]] 
      //   << " corrmat = " << corrmat[i][j] << " corrmat_scaled = " << corrmat_scaled[i][j]
      //   << " measerror^2 = " << measerror[i] * measerror[i] 
      //   << " measerror_scaled^2 = " << measerror_scaled[i]*measerror_scaled[i]
      //   << endl;
    }
  }
  //
  TMatrixD MeasErrorMatrix_scaled(nmeas,nmeas);
  TMatrixD InvMeasErrorMatrix_scaled(nmeas,nmeas);
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) { // by construction diagonal elements of corrmat has not been filled
	MeasErrorMatrix_scaled(i,j)=measerror_scaled[i]*measerror_scaled[i];
      }else{
	MeasErrorMatrix_scaled(i,j)=corrmat_scaled[i][j]*measerror_scaled[i]*measerror_scaled[j]; 
      }
    }
  }
  InvMeasErrorMatrix_scaled = MeasErrorMatrix_scaled;
  Double_t determinant_scaled;
  InvMeasErrorMatrix_scaled.Invert(&determinant_scaled);  
  //  cout << Form("Determinant of MeasErrorMatrix_scaled = %10.5g\n", determinant_scaled);
  //
  // With non-weak measurements
  //
  double pullav_meas_noweak[200]; for (im=0;im<nmeas;++im) pullav_meas_noweak[im]=0; // pullav per measurement
  for (inode=0;inode<=newnode_all;++inode){
    int   npullsq_node_noweak=0;
    double pullsq_node_noweak=0;
    for (im=0;im<nmeas;++im){
      if (weak[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	double pullsquare=0;
	if ( icorrj[im].size() > 1 ) { // first find which cycle this correlated measurement belongs to, then just fetch the result
	  pullsquare = TMath::Power(VectorOfPullMag_noweak[ncycle[im]],2);
	} else {
	  pullsquare = TMath::Power((measvalue[im]-NodeValue_noweak[measnode[im]]),2)/
	    (TMath::Power(measerror[im],2) - TMath::Power(NodeError_noweak[measnode[im]],2));
	}
	pullsq_node_noweak+=pullsquare;
	npullsq_node_noweak+=1;
      }
    }
    pullsq_node_noweak/=npullsq_node_noweak;
    //
    for (im=0;im<nmeas;++im){
      if (weak[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	pullav_meas_noweak[im]=TMath::Sqrt(pullsq_node_noweak);
	//cout << "inode = " << inode << " im+1 = " << im+1 << " npullsq_node_noweak = " << npullsq_node_noweak << " pullav_meas_noweak = " << pullav_meas_noweak[im] << endl;
      }
    }
  }
  //
  double measerror_scaled_noweak[200];
  for (i=0;i<nmeas;++i){
    if (icorrj[i].size() > 1 && VectorOfPullMag_noweak[ncycle[i]]>1) {
      double temperr2=measerror[i]*measerror[i]; // over-written below
      for (int itemp=0; itemp < veccorrij[ncycle[i]].size(); ++itemp) {
	if (i == veccorrij[ncycle[i]][itemp]) {
	  temperr2=((TMatrixD)(VectorOfMeasScaledErrorMatrix_noweak[ncycle[i]]))(itemp,itemp);
	  break;
	}
      }
      measerror_scaled_noweak[i] = TMath::Sqrt(temperr2);
    }else if (icorrj[i].size()==1 && pullav_meas_noweak[i]>1 ) { 
      measerror_scaled_noweak[i] = measerror[i]*pullav_meas_noweak[i];
    }else{
      measerror_scaled_noweak[i] = measerror[i];
    }
  }
  //
  double** corrmat_scaled_noweak = new double*[200]; for (i=0;i<200;++i) corrmat_scaled_noweak[i] = new double[200];
  for (int i=0;i<nmeas;++i){
    for (int j=0;j<nmeas;++j) {
      double tempcorrij=corrmat[i][j]; // over-written below
      if (icorrj[i].size() > 1 && VectorOfPullMag_noweak[ncycle[i]]>1) {
	for (int itemp=0; itemp < veccorrij[ncycle[i]].size(); ++itemp) {
	  if (i == veccorrij[ncycle[i]][itemp]) {
	    for (int jtemp=0; jtemp < veccorrij[ncycle[i]].size(); ++jtemp) {
	      if (j == veccorrij[ncycle[i]][jtemp]) {
		double temperr_i=((TMatrixD)(VectorOfMeasScaledErrorMatrix_noweak[ncycle[i]]))(itemp,itemp);
		double temperr_j=((TMatrixD)(VectorOfMeasScaledErrorMatrix_noweak[ncycle[i]]))(jtemp,jtemp);
		double temperr_ij=((TMatrixD)(VectorOfMeasScaledErrorMatrix_noweak[ncycle[i]]))(itemp,jtemp);
		tempcorrij = temperr_ij / TMath::Sqrt(temperr_i * temperr_j);
		break;
	      }
	    }
	    break;
	  }
	}
      }
      corrmat_scaled_noweak[i][j] = i==j ? 0 : tempcorrij;
      // cout << " i+1 = " << i+1 << " j+1 = " << j+1 
      //   << " icorrij = " << icorrj[i].size() << " pullav_meas_noweak = " << pullav_meas_noweak[i] << " PullMag = " << VectorOfPullMag[ncycle[i]] 
      //   << " corrmat = " << corrmat[i][j] << " corrmat_scaled_noweak = " << corrmat_scaled_noweak[i][j]
      //   << " measerror^2 = " << measerror[i] * measerror[i] 
      //   << " measerror_scaled_noweak^2 = " << measerror_scaled_noweak[i]*measerror_scaled_noweak[i]
      //   << endl;
    }
  }
  //
  TMatrixD MeasErrorMatrix_scaled_noweak(nmeas,nmeas);
  TMatrixD InvMeasErrorMatrix_scaled_noweak(nmeas,nmeas);
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) { // by construction diagonal elements of corrmat has not been filled
	MeasErrorMatrix_scaled_noweak(i,j)=measerror_scaled_noweak[i]*measerror_scaled_noweak[i];
      }else{
	MeasErrorMatrix_scaled_noweak(i,j)=corrmat_scaled_noweak[i][j]*measerror_scaled_noweak[i]*measerror_scaled_noweak[j]; 
      }
    }
  }
  InvMeasErrorMatrix_scaled_noweak = MeasErrorMatrix_scaled_noweak;
  Double_t determinant_scaled_noweak;
  InvMeasErrorMatrix_scaled_noweak.Invert(&determinant_scaled_noweak);  
  //  cout << Form("Determinant of MeasErrorMatrix_scaled = %10.5g\n", determinant_scaled_noweak);
  //
  //
  //  Prepare to obtain Fit Values for non-weak measurements [scaled]
  //
  double basevalue_fit_noweak_scaled[nbase];
  double baseerror_fit_noweak_scaled[nbase];
  double** basecov_fit_noweak_scaled = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit_noweak_scaled[i] = new double[nbase];
  double** basecorr_fit_noweak_scaled = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit_noweak_scaled[i] = new double[nbase];
  double chisquared_noweak_scaled;
  //
  // Define weak nodes
  //
  double nsig_scaled[200];
  int weak_scaled[200];
  for (i=0;i<nmeas;++i){
    if (nchisquare_meas_noweak[i] > 0 ) {
      nsig_scaled[i] = measerror[i]/((sqrt(nchisquare_meas_noweak[i]))*NodeError_noweak[measnode[i]]);
    }else{
      nsig_scaled[i] = 999999;
    }
    weak_scaled[i] = (icorrj[i].size()==1) && (nsig_scaled[i] > 3.); // only for uncorrelated measurements
    //    weak_scaled[i] = weak[i];
    //
    //    cout << i 
    //	 << " nchisquare_meas = " << nchisquare_meas[i] << " nchisquare_meas_noweak = " << nchisquare_meas_noweak[i] 
    //	 << " nsig = " << nsig[i] << " nsig_scaled = " << nsig_scaled[i] 
    //	 << " weak = " << weak[i] << " weak_scaled = " << weak_scaled[i] 
    //	 << endl;
  }
  //
  // COMBINE
  //
  combine(uconstrain,
	  nmeas, weak_scaled, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix_scaled_noweak,
	  node_num_npar, node_den_npar, node_quan, node_parm, 
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit_noweak_scaled, // output
	  baseerror_fit_noweak_scaled, // output
	  basecov_fit_noweak_scaled,   // output
	  basecorr_fit_noweak_scaled,  // output
	  chisquared_noweak_scaled     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_noweak_scaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak_scaled,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  // 
  double chisquared_noweak_scaled_temp = chisquared_noweak_scaled;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit_noweak_scaled[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    //
    // COMBINE
    //
    combine(uconstrain,
	    nmeas, weak_scaled, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix_scaled_noweak,
	    node_num_npar, node_den_npar, node_quan, node_parm, 
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit_noweak_scaled, // output
	    baseerror_fit_noweak_scaled, // output
	    basecov_fit_noweak_scaled,   // output
	    basecorr_fit_noweak_scaled,  // output
	    chisquared_noweak_scaled     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_noweak_scaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak_scaled,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared_noweak_scaled - chisquared_noweak_scaled_temp) < 1e-3) break;
    chisquared_noweak_scaled_temp = chisquared_noweak_scaled;
  }
  //
  cout << endl << "Comparison of Results from fit with [non-weak measurements only] and [errors inflated with PDG-style scale factors] w.r.t original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE  ORIG_FITVAL ORIG_FITERR SCAL_FITVAL SCAL_FITERR   SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6f  %10.6f  %10.6f  %10.6f %6.2f  %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basevalue_fit_noweak_scaled[ibase],
		 baseerror_fit_noweak_scaled[ibase],
		 baseerror_fit_noweak_scaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit_noweak_scaled[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit_noweak_scaled[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num_noweak_scaled[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den_noweak_scaled[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue_noweak_scaled[nnode];     // value of each node
  vector<double> NodeValue_part_noweak_scaled[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num_noweak_scaled, // output
		   NodeValue_den_noweak_scaled, // output
		   NodeValue_noweak_scaled,     // output
		   NodeValue_part_noweak_scaled // output
		   );
  //
  double NodeError_noweak_scaled[nnode];
  TMatrixD NodeErrorMatrix_noweak_scaled(nnode,nnode);
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix_noweak_scaled[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part_noweak_scaled[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part_noweak_scaled[jnode].at(jpar);
	  NodeErrorMatrix_noweak_scaled[inode][jnode] += 
	    ipartial*jpartial*baseerror_fit_noweak_scaled[iquan-1]*baseerror_fit_noweak_scaled[jquan-1]*basecorr_fit_noweak_scaled[iquan-1][jquan-1];
	}
      }
    }
    //
    NodeError_noweak_scaled[inode] = TMath::Sqrt(NodeErrorMatrix_noweak_scaled[inode][inode]);
  }
  //
  int   nchisquare_tot_noweak_scaled=0; 
  double chisquare_tot_noweak_scaled=0; 
  double chisquare_meas_noweak_scaled[200]; // chi2 per mesurement
  int   nchisquare_meas_noweak_scaled[200]; for (im=0;im<nmeas;++im) nchisquare_meas_noweak_scaled[im]=0;// nodes per mesurement
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_scaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas_noweak_scaled[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak_scaled[jm]==weakcompare) continue;
	  chisquare_meas_noweak_scaled[im]+= 
	    (measvalue[im]-NodeValue_noweak_scaled[measnode[im]])*InvMeasErrorMatrix_scaled_noweak[im][jm]*(measvalue[jm]-NodeValue_noweak_scaled[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas_noweak_scaled[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak_scaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot_noweak_scaled+=1;
	chisquare_tot_noweak_scaled+=chisquare_node;
	nchisquare_meas_noweak_scaled[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror_scaled_noweak[im], 
		  NodeValue_noweak_scaled[measnode[im]], NodeError_noweak_scaled[measnode[im]],
		  chisquare_meas_noweak_scaled[im], chisquare_node, chisquare_tot_noweak_scaled,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d \n\n", chisquare_tot_noweak_scaled, nchisquare_tot_noweak_scaled);
  //
  cout << Form("%s = %10.6f +- %10.6f (%10.6f sigma); %s = %10.6f +- %10.6f (%10.6f sigma)\n\n",
	       nodetitle[nnode-2].data(),
	       NodeValue_noweak_scaled[nnode-2],
	       uconstrain?0:NodeError_noweak_scaled[nnode-2],
	       uconstrain?0:(1-NodeValue_noweak_scaled[nnode-2])/NodeError_noweak_scaled[nnode-2],
	       nodetitle[nnode-1].data(),
	       NodeValue_noweak_scaled[nnode-1],
	       NodeError_noweak_scaled[nnode-1],
	       NodeValue_noweak_scaled[nnode-1]/NodeError_noweak_scaled[nnode-1]);
  //
  cout << "Summary of fit with [non-weak measurements only] and [errors inflated with PDG-style scale factors]:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE           SEED    FITVAL     FITERR   RESCALED    SFAC TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6e %10.6f %10.6f %10.6f %6.2f %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),baseseed[ibase],
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 baseerror_fit_noweak_scaled[ibase],
		 baseerror_fit_noweak_scaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  // WRITE OUT MEASUREMENT FILE
  //
  FILE *measfile_scaled[2];
  for (int p=2;p<2;++p){
    if (p==0) measfile_scaled[0]=fopen(Form("combos_measurements_scaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) measfile_scaled[1]=fopen(Form("alucomb_measurements_scaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    FILE *thisfile = measfile_scaled[p];
    print_measfile(thisfile, p, uconstrain,
		   nmeas, weak_scaled, weakcompare, measgammaname, measnode, 
		   expname, author, year, meastitle,
		   measvalue, measerror_scaled_noweak, corrmat_scaled_noweak,
		   a_nodename, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef, 
		   nbase, baseparm, basegamma, basetitle, first_quan, baseseed);
    fclose(thisfile);
  }
  //
  // WRITE OUT AVERAGE.INPUT FILE FOR COMBOS/ALUCOMB
  //
  FILE *avefile_scaled[2];
  for (int p=2;p<2;++p){
    if (p==0) avefile_scaled[0]=fopen(Form("combos_average_scaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) avefile_scaled[1]=fopen(Form("alucomb_average_scaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==0) fprintf (avefile_scaled[p], "INCLUDE combos_measurements_scaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    if (p==1) fprintf (avefile_scaled[p], "INCLUDE alucomb_measurements_scaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    FILE *thisfile = avefile_scaled[p];
    print_avefile(thisfile, p, uconstrain,
		  nmeas, weak_scaled, weakcompare, measgammaname, measnode,
		  node_parm, node_quan, 
		  node_num_npar, node_num_parm, node_num_coef,
		  node_den_npar, node_den_parm, node_den_coef,
		  nbase, baseparm, basegamma, basetitle, first_quan, 
		  baseseed, node_num, node_den, node_part);
    fclose(thisfile);
  }
  //
  // Applying Ad-Hoc Scale Factors
  //
  double measerror_rescaled[200];
  for (i=0;i<nmeas;++i){
    if (measnode[i]==102) {// NAME = S035C9 GAMMA = 96 TITLE = G(K- K+ K- nu(tau))/G(total)
      double rescale = 5.435276;
      //if (measnode[i]==47) {// NAME = S035C21 GAMMA = 85 TITLE = G(K- pi+ pi- nu(tau) (ex.K0))/G(total)
      //  double rescale = 1;//1.59;
      measerror_rescaled[i] = measerror[i] * rescale;
      cout << i << " " << measnode[i]  << " " << measgammaname[i] << " " << measnodename[i] << " " << expname[i] << " " << meastitle[i] << " " 
	   << NodeValue[measnode[i]] << " " << NodeError[measnode[i]] << " " << measvalue[i] << " " << measerror[i]*rescale << " inflated by Scale Factor = " << rescale << endl;
    }else{
      measerror_rescaled[i] = measerror[i];
    }
  }
  //
  TMatrixD MeasErrorMatrix_rescaled(nmeas,nmeas);
  TMatrixD InvMeasErrorMatrix_rescaled(nmeas,nmeas);
  for (i=0;i<nmeas;++i) {
    for (j=0;j<nmeas;++j) {
      if (i==j) {
	MeasErrorMatrix_rescaled(i,j)=measerror_rescaled[i]*measerror_rescaled[i];
      }else{
	MeasErrorMatrix_rescaled(i,j)=corrmat[i][j]*measerror_rescaled[i]*measerror_rescaled[j]; // by construction diagonal elements of corrmat has not been filled
      }
    }
  }
  InvMeasErrorMatrix_rescaled = MeasErrorMatrix_rescaled;
  Double_t determinant_rescaled;
  InvMeasErrorMatrix_rescaled.Invert(&determinant_rescaled);  
  //  cout << Form("Determinant of MeasErrorMatrix_rescaled = %10.5g\n", determinant_rescaled);
  //
  //  Prepare to obtain Fit Values [rescaled]
  //
  double basevalue_fit_rescaled[nbase];
  double baseerror_fit_rescaled[nbase];
  double** basecov_fit_rescaled = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit_rescaled[i] = new double[nbase];
  double** basecorr_fit_rescaled = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit_rescaled[i] = new double[nbase];
  double chisquared_rescaled;
  int weak_none_rescaled[200];  for (i=0;i<nmeas;++i) weak_none_rescaled[i]=0; // all measurements are treated to be non-weak here
  //
  // Derivative of R [node] as a linearized sum of variables P_i [base quantities]
  //
  for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=baseseed_orig[ibase]; // <-- original baseseed re-used
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- may change
		   node_num, // output
		   node_den, // output
		   node_val, // output
		   node_part // output
		   );
  //
  // COMBINE
  //
  combine(uconstrain,
	  nmeas, weak_none_rescaled, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix_rescaled,
	  node_num_npar, node_den_npar, node_quan, node_parm, 
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit_rescaled, // output
	  baseerror_fit_rescaled, // output
	  basecov_fit_rescaled,   // output
	  basecorr_fit_rescaled,  // output
	  chisquared_rescaled     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_rescaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_rescaled,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  //
  double chisquared_rescaled_temp = chisquared_rescaled;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit_rescaled[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    //
    // COMBINE
    //
    combine(uconstrain,
	    nmeas, weak_none_rescaled, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix_rescaled,
	    node_num_npar, node_den_npar, node_quan, node_parm, 
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit_rescaled, // output
	    baseerror_fit_rescaled, // output
	    basecov_fit_rescaled,   // output
	    basecorr_fit_rescaled,  // output
	    chisquared_rescaled     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
		 chisquared_rescaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_rescaled,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared_rescaled - chisquared_rescaled_temp) < 1e-3) break;
    chisquared_rescaled_temp = chisquared_rescaled;
  }
  //
  cout << endl << "Comparison of Results from fit with [errors rescaled in a Ad-Hoc style for kkk only] w.r.t original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE  ORIG_FITVAL ORIG_FITERR SCAL_FITVAL SCAL_FITERR   SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6f  %10.6f  %10.6f  %10.6f %6.2f  %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basevalue_fit_rescaled[ibase],
		 baseerror_fit_rescaled[ibase],
		 baseerror_fit_rescaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit_rescaled[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit_rescaled[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num_rescaled[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den_rescaled[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue_rescaled[nnode];     // value of each node
  vector<double> NodeValue_part_rescaled[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num_rescaled, // output
		   NodeValue_den_rescaled, // output
		   NodeValue_rescaled,     // output
		   NodeValue_part_rescaled // output
		   );
  //
  double NodeError_rescaled[nnode];
  TMatrixD NodeErrorMatrix_rescaled(nnode,nnode);
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix_rescaled[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part_rescaled[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part_rescaled[jnode].at(jpar);
	  NodeErrorMatrix_rescaled[inode][jnode] += 
	    ipartial*jpartial*baseerror_fit_rescaled[iquan-1]*baseerror_fit_rescaled[jquan-1]*basecorr_fit_rescaled[iquan-1][jquan-1];
	}
      }
    }
    //
    NodeError_rescaled[inode] = TMath::Sqrt(NodeErrorMatrix_rescaled[inode][inode]);
  }
  //
  int   nchisquare_tot_rescaled=0; 
  double chisquare_tot_rescaled=0; 
  double chisquare_meas_rescaled[200]; // chi2 per mesurement
  int   nchisquare_meas_rescaled[200];  for (im=0;im<nmeas;++im) nchisquare_meas_rescaled[im]=0;// nodes per mesurement     
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_none_rescaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas_rescaled[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak_none_rescaled[jm]==weakcompare) continue;
	  chisquare_meas_rescaled[im]+= 
	    (measvalue[im]-NodeValue_rescaled[measnode[im]])*InvMeasErrorMatrix_rescaled[im][jm]*(measvalue[jm]-NodeValue_rescaled[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas_rescaled[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak_none_rescaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot_rescaled+=1;
	chisquare_tot_rescaled+=chisquare_node;
	nchisquare_meas_rescaled[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror_rescaled[im], 
		  NodeValue_rescaled[measnode[im]], NodeError_rescaled[measnode[im]],
		  chisquare_meas_rescaled[im], chisquare_node, chisquare_tot_rescaled,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d \n\n", chisquare_tot_rescaled, nchisquare_tot_rescaled);
  //
  cout << Form("%s = %10.6f +- %10.6f (%10.6f sigma); %s = %10.6f +- %10.6f (%10.6f sigma)\n\n",
	       nodetitle[nnode-2].data(),
	       NodeValue_rescaled[nnode-2],
	       uconstrain?0:NodeError_rescaled[nnode-2],
	       uconstrain?0:(1-NodeValue_rescaled[nnode-2])/NodeError_rescaled[nnode-2],
	       nodetitle[nnode-1].data(),
	       NodeValue_rescaled[nnode-1],
	       NodeError_rescaled[nnode-1],
	       NodeValue_rescaled[nnode-1]/NodeError_rescaled[nnode-1]);
  //
  cout << "Summary of fit with [errors rescaled in a Ad-Hoc style for kkk only]:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE           SEED    FITVAL     FITERR   RESCALED    SFAC TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6e %10.6f %10.6f %10.6f %6.2f %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),baseseed[ibase],
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 baseerror_fit_rescaled[ibase],
		 baseerror_fit_rescaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  // WRITE OUT MEASUREMENT FILE
  //
  FILE *measfile_rescaled[2];
  for (int p=0;p<2;++p){
    if (p==0) measfile_rescaled[0]=fopen(Form("combos_measurements_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) measfile_rescaled[1]=fopen(Form("alucomb_measurements_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    FILE *thisfile = measfile_rescaled[p];
    print_measfile(thisfile, p, uconstrain,
		   nmeas, weak_none_rescaled, weakcompare, measgammaname, measnode, 
		   expname, author, year, meastitle,
		   measvalue, measerror_rescaled, corrmat,
		   a_nodename, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef, 
		   nbase, baseparm, basegamma, basetitle, first_quan, baseseed);
    fclose(thisfile);
  }
  //
  // WRITE OUT AVERAGE.INPUT FILE FOR COMBOS/ALUCOMB
  //
  FILE *avefile_rescaled[2];
  for (int p=0;p<2;++p){
    if (p==0) avefile_rescaled[0]=fopen(Form("combos_average_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) avefile_rescaled[1]=fopen(Form("alucomb_average_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==0) fprintf (avefile_rescaled[p], "INCLUDE combos_measurements_rescaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    if (p==1) fprintf (avefile_rescaled[p], "INCLUDE alucomb_measurements_rescaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    FILE *thisfile = avefile_rescaled[p];
    print_avefile(thisfile, p, uconstrain,
		  nmeas, weak_none_rescaled, weakcompare, measgammaname, measnode,
		  node_parm, node_quan, 
		  node_num_npar, node_num_parm, node_num_coef,
		  node_den_npar, node_den_parm, node_den_coef,
		  nbase, baseparm, basegamma, basetitle, first_quan, 
		  baseseed, node_num, node_den, node_part);
    fclose(thisfile);
  }
  //
  //  Prepare to obtain Fit Values for non-weak measurements [rescaled]
  //
  double basevalue_fit_noweak_rescaled[nbase];
  double baseerror_fit_noweak_rescaled[nbase];
  double** basecov_fit_noweak_rescaled = new double*[nbase]; for (i=0;i<nbase;++i) basecov_fit_noweak_rescaled[i] = new double[nbase];
  double** basecorr_fit_noweak_rescaled = new double*[nbase]; for (i=0;i<nbase;++i) basecorr_fit_noweak_rescaled[i] = new double[nbase];
  double chisquared_noweak_rescaled;
  //
  // Define weak nodes
  //
  double nsig_rescaled[200];
  int weak_rescaled[200];
  for (i=0;i<nmeas;++i){
    if (nchisquare_meas_rescaled[i] > 0 ) {
      nsig_rescaled[i] = measerror[i]/((sqrt(nchisquare_meas_rescaled[i]))*NodeError_rescaled[measnode[i]]);
    }else{
      nsig_rescaled[i] = 999999;
    }
    weak_rescaled[i] = (icorrj[i].size()==1) && (nsig_rescaled[i] > 3.); // only for uncorrelated measurements
  }
  //
  // COMBINE
  //
  combine(uconstrain,
	  nmeas, weak_rescaled, weakcompare, nmeas_noweak,
	  measnode, measvalue, InvMeasErrorMatrix_rescaled,
	  node_num_npar, node_den_npar, node_quan, node_parm, 
	  nbase, basegamma, 
	  baseseed, node_num, node_den, node_part, // <-- input may change
	  basevalue_fit_noweak_rescaled, // output
	  baseerror_fit_noweak_rescaled, // output
	  basecov_fit_noweak_rescaled,   // output
	  basecorr_fit_noweak_rescaled,  // output
	  chisquared_noweak_rescaled     // output
	  );
  //
  cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
	       chisquared_noweak_rescaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak_rescaled,nmeas_noweak-nbase_u));
  //
  // Iterate ...
  // 
  double chisquared_noweak_rescaled_temp = chisquared_noweak_rescaled;
  while (1) {
    for (ibase=0;ibase<nbase;++ibase) baseseed[ibase]=basevalue_fit_noweak_rescaled[ibase]; // <-- baseseed updated
    get_num_den_part(nnode, node_parm,
		     node_num_npar, node_num_parm, node_num_coef,
		     node_den_npar, node_den_parm, node_den_coef,
		     baseparm, 
		     baseseed, // <-- input may change
		     node_num, // output
		     node_den, // output
		     node_val, // output
		     node_part // output
		     );
    //
    // COMBINE
    //
    combine(uconstrain,
	    nmeas, weak_rescaled, weakcompare, nmeas_noweak,
	    measnode, measvalue, InvMeasErrorMatrix_rescaled,
	    node_num_npar, node_den_npar, node_quan, node_parm, 
	    nbase, basegamma, 
	    baseseed, node_num, node_den, node_part, // <-- input may change
	    basevalue_fit_noweak_rescaled, // output
	    baseerror_fit_noweak_rescaled, // output
	    basecov_fit_noweak_rescaled,   // output
	    basecorr_fit_noweak_rescaled,  // output
	    chisquared_noweak_rescaled     // output
	    );
    //
    cout << Form("ChiSquared = %8.4f  NMeas = %3d Nbase_u = %2d NDOF = %3d CL = %8.4g\n",
		 chisquared_noweak_rescaled, nmeas_noweak, nbase_u, nmeas_noweak - nbase_u, TMath::Prob(chisquared_noweak_rescaled,nmeas_noweak-nbase_u));
    //
    if (TMath::Abs(chisquared_noweak_rescaled - chisquared_noweak_rescaled_temp) < 1e-3) break;
    chisquared_noweak_rescaled_temp = chisquared_noweak_rescaled;
  }
  //
  cout << endl << "Comparison of Results from fit with [non-weak measurements only] and [errors rescaled in a Ad-Hoc style for kkk only] w.r.t original fit:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE  ORIG_FITVAL ORIG_FITERR SCAL_FITVAL SCAL_FITERR   SFAC  TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6f  %10.6f  %10.6f  %10.6f %6.2f  %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 basevalue_fit_noweak_rescaled[ibase],
		 baseerror_fit_noweak_rescaled[ibase],
		 baseerror_fit_noweak_rescaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  //  for (ibase=0;ibase<nbase;++ibase) {
  //    cout << "basecorr["<<ibase<<"] = ";
  //    double checkzero=0;
  //    for (jbase=0;jbase<nbase;++jbase) {
  //      checkzero+=basecov_fit_noweak_rescaled[ibase][jbase];
  //      cout << Form("%8.4g ",basecorr_fit_noweak_rescaled[ibase][jbase]);
  //      if (jbase%16==15) cout << endl;
  //    }
  //    cout << " checkzero = " << checkzero << endl; 
  //  }
  //
  //
  double NodeValue_num_noweak_rescaled[nnode]; //   numerator sum [calculated using seed values of quantities] for each node
  double NodeValue_den_noweak_rescaled[nnode]; // denominator sum [calculated using seed values of quantities] for each node
  double NodeValue_noweak_rescaled[nnode];     // value of each node
  vector<double> NodeValue_part_noweak_rescaled[nnode]; // vector of partial derivatives w.r.t quantities for each node
  get_num_den_part(nnode, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef,
		   baseparm, 
		   baseseed, // <-- input may change [Using same input as in previous combine command]
		   NodeValue_num_noweak_rescaled, // output
		   NodeValue_den_noweak_rescaled, // output
		   NodeValue_noweak_rescaled,     // output
		   NodeValue_part_noweak_rescaled // output
		   );
  //
  double NodeError_noweak_rescaled[nnode];
  TMatrixD NodeErrorMatrix_noweak_rescaled(nnode,nnode);
  for (inode=0;inode<nnode;++inode) {
    for (jnode=0;jnode<nnode;++jnode) {
      NodeErrorMatrix_noweak_rescaled[inode][jnode] = 0; 
      for (ipar = 0; ipar < node_parm[inode].size(); ++ipar) {
	int iquan=node_quan[inode].at(ipar);
	double ipartial=NodeValue_part_noweak_rescaled[inode].at(ipar);
	for (jpar = 0; jpar < node_parm[jnode].size(); ++jpar) {
	  int jquan=node_quan[jnode].at(jpar);
	  double jpartial=NodeValue_part_noweak_rescaled[jnode].at(jpar);
	  NodeErrorMatrix_noweak_rescaled[inode][jnode] += 
	    ipartial*jpartial*baseerror_fit_noweak_rescaled[iquan-1]*baseerror_fit_noweak_rescaled[jquan-1]*basecorr_fit_noweak_rescaled[iquan-1][jquan-1];
	}
      }
    }
    //
    NodeError_noweak_rescaled[inode] = TMath::Sqrt(NodeErrorMatrix_noweak_rescaled[inode][inode]);
  }
  //
  int   nchisquare_tot_noweak_rescaled=0; 
  double chisquare_tot_noweak_rescaled=0; 
  double chisquare_meas_noweak_rescaled[200]; // chi2 per mesurement
  int   nchisquare_meas_noweak_rescaled[200];  for (im=0;im<nmeas;++im) nchisquare_meas_noweak_rescaled[im]=0;// nodes per mesurement     
  for (inode=0;inode<=newnode_all;++inode){
    int   nchisquare_node=0;
    double chisquare_node=0;
    for (im=0;im<nmeas;++im){
      if (weak_rescaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	chisquare_meas_noweak_rescaled[im]=0;
	for (jm=0;jm<nmeas;++jm) {
	  if (weak_rescaled[jm]==weakcompare) continue;
	  chisquare_meas_noweak_rescaled[im]+= 
	    (measvalue[im]-NodeValue_noweak_rescaled[measnode[im]])*InvMeasErrorMatrix_rescaled[im][jm]*(measvalue[jm]-NodeValue_noweak_rescaled[measnode[jm]]);
	}
	nchisquare_node+=1;
	chisquare_node+=chisquare_meas_noweak_rescaled[im];
      }
    }
    chisquare_node/=nchisquare_node;
    for (im=0;im<nmeas;++im){
      if (weak_rescaled[im]==weakcompare) continue;
      if (nodegroup_all[im]==inode) {
	nchisquare_tot_noweak_rescaled+=1;
	chisquare_tot_noweak_rescaled+=chisquare_node;
	nchisquare_meas_noweak_rescaled[im]=nchisquare_node;
	cout 
          << Form("i+1 = %3d group = %2d ngroup = %2d ncycle = %2d meas = %8.4e +- %8.4e fit = %8.4e +- %8.4e chi2 = %8.4f chi2node = %8.4f chi2tot = %8.4f %6s %s\n",
		  im+1,nodegroup_all[im],nchisquare_node,ncycle[im], 
		  measvalue[im], measerror_rescaled[im], 
		  NodeValue_noweak_rescaled[measnode[im]], NodeError_noweak_rescaled[measnode[im]],
		  chisquare_meas_noweak_rescaled[im], chisquare_node, chisquare_tot_noweak_rescaled,
		  expname[im].data(), measgammaname[im].data());
      }
    }
  }
  //
  cout << Form("chisquare_tot = %8.4f nchisquare_tot = %3d \n\n", chisquare_tot_noweak_rescaled, nchisquare_tot_noweak_rescaled);
  //
  cout << Form("%s = %10.6f +- %10.6f (%10.6f sigma); %s = %10.6f +- %10.6f (%10.6f sigma)\n\n",
	       nodetitle[nnode-2].data(),
	       NodeValue_noweak_rescaled[nnode-2],
	       uconstrain?0:NodeError_noweak_rescaled[nnode-2],
	       uconstrain?0:(1-NodeValue_noweak_rescaled[nnode-2])/NodeError_noweak_rescaled[nnode-2],
	       nodetitle[nnode-1].data(),
	       NodeValue_noweak_rescaled[nnode-1],
	       NodeError_noweak_rescaled[nnode-1],
	       NodeValue_noweak_rescaled[nnode-1]/NodeError_noweak_rescaled[nnode-1]);
  //
  cout << "Summary of fit with [non-weak measurements only] and [errors rescaled in a Ad-Hoc style for kkk only]:" << endl;
  cout << Form("# QUAN GAMMA  PARM   NODE           SEED    FITVAL     FITERR   RESCALED    SFAC TITLE\n");
  for (ibase=0;ibase<nbase;++ibase) {
    cout << Form("*%5d %5d %5d   %-5s  %10.6e %10.6f %10.6f %10.6f %6.2f %-48s\n",
		 ibase+1,basegamma[ibase],baseparm[ibase],basenode[ibase].data(),baseseed[ibase],
		 basevalue_fit[ibase],
		 baseerror_fit[ibase],
		 baseerror_fit_noweak_rescaled[ibase],
		 baseerror_fit_noweak_rescaled[ibase]/baseerror_fit[ibase],
		 basetitle[ibase].data());
  }
  cout << endl;
  //
  // WRITE OUT MEASUREMENT FILE
  //
  FILE *measfile_noweak_rescaled[2];
  for (int p=2;p<2;++p){
    if (p==0) measfile_noweak_rescaled[0]=fopen(Form("combos_measurements_noweak_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) measfile_noweak_rescaled[1]=fopen(Form("alucomb_measurements_noweak_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    FILE *thisfile = measfile_noweak_rescaled[p];
    print_measfile(thisfile, p, uconstrain,
		   nmeas, weak_rescaled, weakcompare, measgammaname, measnode, 
		   expname, author, year, meastitle,
		   measvalue, measerror_rescaled, corrmat,
		   a_nodename, node_parm,
		   node_num_npar, node_num_parm, node_num_coef,
		   node_den_npar, node_den_parm, node_den_coef, 
		   nbase, baseparm, basegamma, basetitle, first_quan, baseseed);
    fclose(thisfile);
  }
  //
  // WRITE OUT AVERAGE.INPUT FILE FOR COMBOS/ALUCOMB
  //
  FILE *avefile_noweak_rescaled[2];
  for (int p=2;p<2;++p){
    if (p==0) avefile_noweak_rescaled[0]=fopen(Form("combos_average_noweak_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==1) avefile_noweak_rescaled[1]=fopen(Form("alucomb_average_noweak_rescaled_%s%s.input",sconstrain.data(),salephhcorr.data()),"w");
    if (p==0) fprintf (avefile_noweak_rescaled[p], "INCLUDE combos_measurements_noweak_rescaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    if (p==1) fprintf (avefile_noweak_rescaled[p], "INCLUDE alucomb_measurements_noweak_rescaled_%s%s.input \n\n",sconstrain.data(),salephhcorr.data()); 
    FILE *thisfile = avefile_noweak_rescaled[p];
    print_avefile(thisfile, p, uconstrain,
		  nmeas, weak_rescaled, weakcompare, measgammaname, measnode,
		  node_parm, node_quan, 
		  node_num_npar, node_num_parm, node_num_coef,
		  node_den_npar, node_den_parm, node_den_coef,
		  nbase, baseparm, basegamma, basetitle, first_quan, 
		  baseseed, node_num, node_den, node_part);
    fclose(thisfile);
  }
  //
  delete [] corrmat_scaled_noweak;
  delete [] corrmat_scaled;
  delete [] basecorr_fit_noweak;
  delete [] basecov_fit_noweak;
  delete [] basecorr_fit;
  delete [] basecov_fit;
  delete [] corrmat;
  delete [] a_nodename;
}
