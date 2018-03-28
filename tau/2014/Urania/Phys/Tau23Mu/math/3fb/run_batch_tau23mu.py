#! /usr/bin/env python
__author="Marcin Chrzaszcz"

import os
import sys
from ROOT import *
import time
from XTuple import *
from MCLimit import *
NTOYS=2000000
#NTOYS=1000

print "Welcome!"
print "To use me: "
print "run <DECAY NUBER>"
#print "for DECAY NUMBERS see HFAG report"
print "Author: "+__author

def parse(filepath, number):
       parsed_file = open(filepath)
       lines = [line.split()[1:] for line in parsed_file if line.startswith(str(number))]         
       
       return lines
################################################################
def table(lines, number):
       size=len(lines)
       column=[]
       for i in range(0,size):
              if(number>0): column.append( float(lines[i][number] ))
              if(number==0): column.append( lines[i][number] )
       return column

################################################################
def get_BKG(n_bck, b_bck_err, names):
       size=len(n_bck)
       size2=len(b_bck_err)
       if(size!=size2):
              print "ERROR line 30!"
       bkg_hist =[]
       bkgerr_hist =[]
       numberBkg= []
       numberBkgErr=[]
       for i in range(0, size):
              nkg_tmp= TH1F(names[i]+"bkg", names[i]+"bkg",  1,0,10)
              nkg_tmp.SetBinContent(1, n_bck[i])
              nkg_err= TH1F(names[i]+"bkg_err", names[i]+"bkg_err", 1,0,10)
              nkg_err.SetBinContent(1,b_bck_err[i] )

              numberBkg.append(nkg_tmp.GetSum())
              nkg_tmp.Scale(1./nkg_tmp.GetSum())
              bkg_hist.append(nkg_tmp)

              numberBkgErr.append(nkg_err.GetSum())
              nkg_err.Scale(1./nkg_err.GetSum())
              bkgerr_hist.append(nkg_err)

       #bkg_up_hist.append(nkg_err_tmp_up)                           
       #bkg_down_hist.append(nkg_err_tmp_down)                       



       return bkg_hist, bkgerr_hist, numberBkg, numberBkgErr

################################################################
def get_SIG(names):
       size=len(names)
       sig_hist =[]
       for i in range(0, size):
              sig_tmp= TH1F(names[i]+"sig", names[i]+"sig",1,0,10)
              sig_tmp.SetBinContent(1,1.)
              sig_hist.append(sig_tmp)
       return sig_hist
              
################################################################
def get_obs(obs,names):
       size=len(names)
       obs_hist=[]
       for i in range(0, size):
              obs_tmp=  TH1F(names[i]+"obs", names[i]+"obs", 1,0,10)
              obs_tmp.SetBinContent(1,obs[i])
              print "Observed; ",obs[i]
              obs_hist.append(obs_tmp)
       return obs_hist

################################################################
def intefy(hist, norm):
  from scipy import random as rnd
  rethist = hist.Clone(hist.GetName()+"_random")
  rethist.Reset()
  for x in range(hist.GetNbinsX()):
    for y in range(hist.GetNbinsX()):
      bin = rethist.GetBin(x+1,y+1)
      rethist.SetBinContent(bin,float(int(norm * rnd.poisson(hist.GetBinContent(bin)))))
  return rethist


################################################################
# SWITCHES:
debug = True
use_BABAR= False
use_BELLE = False
use_LHCb = True

CHANNEL = "Tau -> 3 mu"

################################################################
# main declarations
#c1= TCanvas("c1", "c1", 800,600)
nul = csm_model()
nul_no_syst = csm_model()



################################################################
if( len(sys.argv)!=3):
       print "TO MANY PAREMTERS: " +str( len(sys.argv))
       sys.exit()
DECAY_NUMBER=sys.argv[1]
BR=float(sys.argv[2])
print "Running decay: "+str(DECAY_NUMBER)
params = parse('lfv_data.txt',DECAY_NUMBER)
names=table(params,0)
lumi=table(params,1)
x_section=table(params,2)
eff=table(params,3)
eff_err=table(params,4)
n_bck=table(params,5)
b_bck_err= table(params,6)
obs= table(params,7)
NUMBER_EXPERIMENTS=len(params)
#if(NUMBER_EXPERIMENTS!=2):
#       print "exit 1"
#       sys.exit()
#if(int(DECAY_NUMBER)!=183):
#       print "exit 2"
#       print "decay number: ",DECAY_NUMBER
#       sys.exit()

if(debug):
       print "names "+str( names )
       print "lumi "+str( lumi )
       print "xSection "+str( x_section )
       print "eff "+str( eff )
       print "eff_err "+str( eff_err )
       print "n_bck "+str( n_bck )
       print "n_bck_err "+str( b_bck_err )
       print "obs "+str( obs )
       
       

#print "Will combine"+str(NUMBER_EXPERIMENTS)+"Experiments: "+params[0][0]+"  "+params[1][0]

# CREATING NTAU SYSTEMATICS:
ntau = []
ntau_err = []
for i in range(0,NUMBER_EXPERIMENTS):
       #print params[i]
       ntau.append(lumi[i]*x_section[i]*2.*eff[i])
       ntau[i]=ntau[i]*(1e6);
       ntau_err.append(ntau[i]* (eff_err[i]/eff[i]))


# CREATING BACKGROUND SYSTEMATICS:
(bkg_hist, bkgerr_hist, numberBkg, numberBkgErr) =get_BKG(n_bck, b_bck_err, names)
print numberBkg
sig_hists=get_SIG(names)
print sig_hists
ob_hists=get_obs(obs,names)

### OK WE have got everything ready for limit computation

npb_plist_BABAR = []
bdevup_list_BABAR = []
npb_mlist_BABAR = []
bdevdown_list_BABAR = []
ffblist_BABAR = []

npb_plist_BELLE = []
bdevup_list_BELLE = []
npb_mlist_BELLE = []
bdevdown_list_BELLE = []
ffblist_BELLE = []



# loop over experiments
# experiment BaBar:
if use_BABAR:
       BABAR=-1
       BABAR_BKG_SYS=1
       if "BaBar" in names[0]:
              BABAR=0
       if "BaBar"  in names[1]:
              BABAR=1
       print "BaBar= ",BABAR
       bkg_BABAR = csm_template(BABAR_BKG_SYS)
       obs_BABAR= ob_hists[BABAR]
       sigpdf_BABAR=sig_hists[BABAR]
       bkgpdf_BABAR=bkg_hist[BABAR]  # This is scaled to 1                                            
       bkgERR_BABAR=numberBkgErr[BABAR]/numberBkg[BABAR]
       print "bkg err BaBar:",bkgERR_BABAR

       for ii in range(0,BABAR_BKG_SYS):
              bkg_BABAR.set_np("error pdf bin"+names[BABAR],bkgERR_BABAR,-bkgERR_BABAR, NULL,1,NULL,-1)
              
       bkg_BABAR.make(bkgpdf_BABAR,numberBkg[BABAR],0,0, names[BABAR])
       bkg_BABAR.add_to(nul)
       nul.set_interpolation_style(names[BABAR],1)
       bkg_no_syst_BABAR = csm_template(0)
       bkg_no_syst_BABAR.make(bkgpdf_BABAR,numberBkg[BABAR] ,0,0, names[BABAR])
       bkg_no_syst_BABAR.add_to(nul_no_syst)
       nul_no_syst.set_interpolation_style(names[BABAR],1)

       

if use_BELLE:
       BELLE=-1
       BELLE_BKG_SYS=1
       if "Belle" in names[0]:
              BELLE=0
       if "Belle"  in names[1]:
              BELLE=1
              #print BELLE
       bkg_BELLE = csm_template(BELLE_BKG_SYS)
       obs_BELLE= ob_hists[BELLE]
       sigpdf_BELLE=sig_hists[BELLE]
       bkgpdf_BELLE=bkg_hist[BELLE]  # This is scaled to 1                                             
       bkgERR_BELLE=numberBkgErr[BELLE]/numberBkg[BELLE]
       print "bkg err Belle:",bkgERR_BELLE

       for ii in range(0,BELLE_BKG_SYS):
              bkg_BELLE.set_np("error pdf bin"+names[BELLE],bkgERR_BELLE,-bkgERR_BELLE, NULL,1,NULL,-1)

       bkg_BELLE.make(bkgpdf_BELLE,numberBkg[BELLE] ,0,0, names[BELLE])
       bkg_BELLE.add_to(nul)
       nul.set_interpolation_style(names[BELLE],1)

       bkg_no_syst_BELLE = csm_template(0)
       bkg_no_syst_BELLE.make(bkgpdf_BELLE,numberBkg[BELLE] ,0,0, names[BELLE])
       bkg_no_syst_BELLE.add_to(nul_no_syst)
       nul_no_syst.set_interpolation_style(names[BELLE],1)
       
#if(BELLE==BABAR):
#       print "Belle "+str(BELLE)
#       print "BABAR "+str(BABAR)
#       print "You want to die ?"
#       sys.exit()
#print "afetr"
####################################################################################################
# ADDING LHCB!!!
####################################################################################################
polarity="down"
stupidbuffer = []

ontop    = [[1.,1.],[1.,1.],[1.,1.],[1.,1.],[3.9/4.05,4.15/4.05],[4.39/4.05,4.01/4.05],[4.23/4.05,3.8/4.05],[3.97/4.05,4.08/4.05]]

sigdofs  = [[20,21],[14,15],[16,17],[18,19],[ 1,5 ],[ 2,6 ],[ 3,7 ],[ 4,8 ]]
use_2011_tau23mu= True
use_2012_tau23mu= True
use_eta_PDF=False

use_s_alpha = True        
use_bck_sys = True        
use_sig_sys = True

bck_sys_2011= 0
bck_sys_2012= 0
alpha_2011 = 0.
s_alpha_2011 = 0.
alpha_2012 = 0.
s_alpha_2012 = 0.

if use_eta_PDF:
  alpha_2011 =  3.10e-9
  s_alpha_2011 =  0.37-9
else:
  alpha_2011 = 3.81e-9
  s_alpha_2011 = 0.46e-9

if use_eta_PDF:
  alpha_2012 = 1.39e-9
  s_alpha_2012 = 0.23e-9
else:
  alpha_2012 = 1.72e-9
  s_alpha_2012 = 0.19e-9

alpha_2011 = alpha_2011 * 9.607e-01 # improve by my ctau to avoid double accounting        
s_alpha_2011 = s_alpha_2011 * 9.607e-01 # improve by my ctau to avoid double accounting    
alpha_2012 = alpha_2012 * 9.607e-01 # improve by my ctau to avoid double accounting        
s_alpha_2012 = s_alpha_2012 * 9.607e-01 # improve by my ctau to avoid double accounting    

if use_bck_sys:
  bck_sys_2011 = bck_sys_2011+30
  bck_sys_2012 = bck_sys_2012+35

if alpha_2011==0:
    print "alpha shit!"
    sys.exit()
if alpha_2012==0:
    print "alpha shit!"
    sys.exit()

if s_alpha_2011==0:
      print "s_alpha shit!"
      sys.exit()
if s_alpha_2012==0:
      print "s_alpha shit!"
      sys.exit()

addSM_toBkg = 0
SYST = 0
ISTYLE = 1 # VERTICAL
npb_plist_2011 = []
bdevup_list_2011 = []
npb_mlist_2011 = []
bdevdown_list_2011 = []
ffblist_2011 = []

npb_plist_2012 = []
bdevup_list_2012 = []
npb_mlist_2012 = []
bdevdown_list_2012 = []
ffblist_2012 = []



'''
if use_2011_tau23mu:
  outname += "_2011"
if use_2012_tau23mu:
  outname += "_2012"
if external_BR is not None:
  outname += "_" + str(external_BR)
if use_eta_PDF:
  outname += "with_ETA"
else:
   outname += "with_ETA_VETO"

'''

if not (use_2011_tau23mu or use_2012_tau23mu):
  print "no year chosen"
  sys.exit()

#print outname

print "Using Tau23Mu 2011:", use_2011_tau23mu
print "Using Tau23Mu 2012:", use_2012_tau23mu



if use_LHCb:
       bkg_2011 = csm_template(bck_sys_2011)

       if use_eta_PDF:
              bkgfile_2011="../../params/3fb/2011_PDF_BCK_FINAL/bck_eta_notrash.root"
              f_bkg_2011=TFile(bkgfile_2011)
              bkgpdf_2011 = f_bkg_2011.Get("expected_exp").Clone("bkg_2011")
       else:
              bkgfile_2011="../../params/3fb/2011_PDF_BCK_FINAL/bck_eta_cut450_notrash.root"
              print "NOT TAKING THIS"
              f_bkg_2011=TFile(bkgfile_2011)
              bkgpdf_2011 = f_bkg_2011.Get("expected_exp").Clone("bkg_2011")

       numberBkg_2011 = bkgpdf_2011.GetSum()
       # bkgpdfsum_2011 = bkgpdf_2011.GetSum()                                             
       bkgpdf_2011.Scale(1./bkgpdf_2011.GetSum())
  ###############################################################                   
  ### BCK SYSTEMATICS                                                               
       if use_bck_sys:
       #nuissance                                                                        
              for b in range (0,bck_sys_2011):
                     fname = "../../params/3fb/2011_PDF_BCK_FINAL/"
                     if use_eta_PDF:
                            fname=fname+"Sys_ETA/Sys_"
                     else:
                            fname=fname+"Sys_ETAVETO/Sys_"
                     fname=fname+str(b)+".root"
                     f_now = TFile(fname)
                     name_plus="Sys_plus_"+str(b)
                     name_minus="Sys_minus_"+str(b)
                     pl = f_now.Get(name_plus)
                     mi = f_now.Get(name_minus)
                     ffblist_2011+=[f_now]
                     npb_plist_2011+=[pl]
                     npb_mlist_2011+=[mi]
                     bdevup_list_2011+=[(pl.GetSum()-numberBkg_2011)/numberBkg_2011/3.]
                     bdevdown_list_2011+=[(mi.GetSum()-numberBkg_2011)/numberBkg_2011/3.]
              for ii in range(0,bck_sys_2011):
                     npb_plist_2011[ii].Scale(1./npb_plist_2011[ii].GetSum())
                     npb_mlist_2011[ii].Scale(1./npb_mlist_2011[ii].GetSum())
                     bkg_2011.set_np("error pdf bin"+str(ii),bdevup_list_2011[ii],bdevdown_list_2011[ii],npb_plist_2011[ii],9.,npb_mlist_2011[ii],-9.)

       if use_eta_PDF:
              sigfile_2011 = "../../params/3fb/base_2011_calibV5_down_noveto_trash.root"
       else:
              sigfile_2011 = "../../params/3fb/base_2011_calibV4_down_trash.root"  # ORYGINAL FILE    
              sigfile_2011 = "../../params/3fb/base_2011_calibV4_up_sys1_trash.root"                  
              sigfile_2011 = "../../params/3fb/base_2011_calibV5_down_sys1_trash.root"

       sig_sys_pdfs_2011 = []
       sig_sys_corr_2011 = []
       
       f_sig_2011 = TFile(sigfile_2011)
       sigpdf_2011 = f_sig_2011.Get("central").Clone("sig_2011")
       sigpdfsum_2011 = sigpdf_2011.GetSum()
       print "trash efficiency: ", sigpdfsum_2011

       alpha_2011=alpha_2011/sigpdfsum_2011
       s_alpha_2011=s_alpha_2011/sigpdfsum_2011
       sigpdf_2011.Scale(1./sigpdf_2011.GetSum())

       #################   SIGNAL SYSTEMATICS                                                      
       if use_sig_sys:
              for ji in range(len(sigdofs)):  #FIXME what if trashing                                   
                     j = sigdofs[ji]
                     pdfbuffer = []
                     corrbuffer = []
                     for ii in range(len(j)):
                            jj = j[ii]
                            fname = "../../params/3fb/base_2011_calibV5_"+polarity+"_sys"+str(jj)+"_trash.root"
                            print "opening ",fname
                            f_now = TFile(fname)
                            print "ls: ",f_now.ls()
                            pdf = f_now.Get("central").Clone("central_"+str(jj))
                            stupidbuffer += [pdf,f_now]
                            print type(pdf)
                            pdfsum = pdf.GetSum()
                            pdf.Scale(1./pdfsum)
                            pdfbuffer += [pdf]
                            corrbuffer += [ ontop[ji][ii]*(pdfsum-sigpdfsum_2011 )/sigpdfsum_2011 ]
                     sig_sys_pdfs_2011 += [ pdfbuffer ]
                     sig_sys_corr_2011 += [ corrbuffer ]

       f_unbhist_2011=TFile('../../params/3fb/unblind/bck_eta_cut450_OBSERVED.root')
       unbhist_2011=f_unbhist_2011.Get('observed').Clone()

       print "normalized bkgpdf_2011.GetSum(): ",bkgpdf_2011.GetSum()
       #  bkg_2011 = csm_template(0)                                 
       bkg_2011.make(bkgpdf_2011,numberBkg_2011 ,0,0, "2011")
       bkg_2011.add_to(nul)
       nul.set_interpolation_style("2011",ISTYLE)
       
       bkg_no_syst_2011 = csm_template(0)
       bkg_no_syst_2011.make(bkgpdf_2011,numberBkg_2011 ,0,0, "2011")
       bkg_no_syst_2011.add_to(nul_no_syst)
       nul_no_syst.set_interpolation_style("2011",ISTYLE)
       
       ##########################################
       ### 2012 PDF:
       ##########################################
       bkg_2012 = csm_template(bck_sys_2012)
       if use_eta_PDF:
              bkgfile_2012="../../params/3fb/2012_PDF_BCK_FINAL/bck_eta_notrash.root"
              f_bkg_2012=TFile(bkgfile_2012)
              bkgpdf_2012 = f_bkg_2012.Get("expected_exp").Clone("bkg_2012")
       else:
              bkgfile_2012="../../params/3fb/2012_PDF_BCK_FINAL/bck_eta_cut450_notrash.root"
              f_bkg_2012=TFile(bkgfile_2012)
              bkgpdf_2012 = f_bkg_2012.Get("expected_exp").Clone("bkg_2012")

              
       numberBkg_2012=bkgpdf_2012.GetSum()
       bkgpdf_2012.Scale(1./bkgpdf_2012.GetSum())
       if use_bck_sys:

              for b in range (0,bck_sys_2012):
                     fname = "../../params/3fb/2012_PDF_BCK_FINAL/"
                     if use_eta_PDF:
                            fname=fname+"Sys_ETA/Sys_"
                     else:
                            fname=fname+"Sys_ETAVETO/Sys_"
                     fname=fname+str(b)+".root"
                     f_now = TFile(fname)
                     name_plus="Sys_plus_"+str(b)
                     name_minus="Sys_minus_"+str(b)
                     pl = f_now.Get(name_plus)
                     mi = f_now.Get(name_minus)
                     ffblist_2012+=[f_now]
                     npb_plist_2012+=[pl]
                     npb_mlist_2012+=[mi]
                     bdevup_list_2012+=[(pl.GetSum()-numberBkg_2012)/numberBkg_2012/3.]
                     bdevdown_list_2012+=[(mi.GetSum()-numberBkg_2012)/numberBkg_2012/3.]
              for ii in range(0,bck_sys_2012):
                     npb_plist_2012[ii].Scale(1./npb_plist_2012[ii].GetSum())
                     npb_mlist_2012[ii].Scale(1./npb_mlist_2012[ii].GetSum())
                     bkg_2012.set_np("error pdf bin"+str(ii),bdevup_list_2012[ii],bdevdown_list_2012[ii],npb_plist_2012[ii],9.,npb_mlist_2012[ii],-9.)

       if use_eta_PDF:
              sigfile_2012 = "../../params/3fb/base_2012_calibV5_down_noveto_trash.root"
       else:
              sigfile_2012 = "../../params/3fb/base_2012_calibV5_down_trash.root"

       sig_sys_pdfs_2012 = []
       sig_sys_corr_2012 = []


       f_sig_2012 = TFile(sigfile_2012)
       sigpdf_2012 = f_sig_2012.Get("central").Clone("sig_2012")
       sigpdfsum_2012 = sigpdf_2012.GetSum()
       print "trash efficiency: ",sigpdfsum_2012

       alpha_2012=alpha_2012/sigpdfsum_2012
       s_alpha_2012=s_alpha_2012/sigpdfsum_2012
       sigpdf_2012.Scale(1./sigpdf_2012.GetSum())
       
       if use_sig_sys:
            #################   SIGNAL SYSTEMATICS                                                      
              for ji in range(len(sigdofs)):  #FIXME what if trashing                                     
                     j = sigdofs[ji]
                     pdfbuffer = []
                     corrbuffer = []
                     for ii in range(len(j)):
                            jj = j[ii]
                            fname = "../../params/3fb/base_2012_calibV5_"+polarity+"_sys"+str(jj)+"_trash.root"
                            print "opening ",fname
                            f_now = TFile(fname)
                            print "ls: ",f_now.ls()
                            pdf = f_now.Get("central").Clone("central_"+str(jj))
                            stupidbuffer += [pdf,f_now]
                            print type(pdf)
                            pdfsum = pdf.GetSum()
                            pdf.Scale(1./pdfsum)
                            pdfbuffer += [pdf]
                            corrbuffer += [ ontop[ji][ii]*(pdfsum-sigpdfsum_2012 )/sigpdfsum_2012 ]
                     sig_sys_pdfs_2012 += [ pdfbuffer ]
                     sig_sys_corr_2012 += [ corrbuffer ]

       f_unbhist_2012=TFile('../../params/3fb/unblind/bck_eta_cut450_OBSERVED_2012.root')
       unbhist_2012=f_unbhist_2012.Get('observed').Clone()


       print "bkgpdf.GetSum(): ",bkgpdf_2012.GetSum()

       print "normalized bkgpdf_2012.GetSum(): ",bkgpdf_2012.GetSum()


       print "set bkg model"

       bkg_2012.make(bkgpdf_2012,numberBkg_2012 ,0,0, "2012")
       bkg_2012.add_to(nul)
       nul.set_interpolation_style("2012",ISTYLE)
       
       bkg_no_syst_2012 = csm_template(0)
       bkg_no_syst_2012.make(bkgpdf_2012,numberBkg_2012 ,0,0, "2012")
       bkg_no_syst_2012.add_to(nul_no_syst)
       nul_no_syst.set_interpolation_style("2012",ISTYLE)





              

#######################################       
def DoTestHyp(internal_BR):

           NSIGSYS=NUMBER_EXPERIMENTS

           test = csm_model()
           test_no_syst = csm_model()
           returnparameters = []
                      

           NSIGSYS_BABAR=1
           if use_BABAR:
                  bkg_BABAR.add_to(test)
                  bkg_no_syst_BABAR.add_to(test_no_syst)
                  dc_BABAR = csm_template(NSIGSYS_BABAR)
                  dc_no_syst_BABAR = csm_template(0)
                  dc_BABAR.set_np("nomlization of sig "+names[BABAR], ntau_err[BABAR]/ntau[BABAR], -ntau_err[BABAR]/ntau[BABAR] , NULL,1,NULL,-1)
                  dc_BABAR.make(sigpdf_BABAR,internal_BR*ntau[BABAR],0,1, names[BABAR])
                  dc_no_syst_BABAR.make(sigpdf_BABAR,internal_BR*ntau[BABAR],0,1, names[BABAR])
                  dc_BABAR.add_to(test)
                  dc_no_syst_BABAR.add_to(test_no_syst)
                  test.set_interpolation_style(names[BABAR],1)
                  test_no_syst.set_interpolation_style(names[BABAR],1)
                                    
                  returnparameters += [dc_BABAR,dc_no_syst_BABAR]
           
           NSIGSYS_BELLE=1
           if use_BELLE:
                  bkg_BELLE.add_to(test)
                  bkg_no_syst_BELLE.add_to(test_no_syst)
                  dc_BELLE = csm_template(NSIGSYS_BELLE)
                  dc_no_syst_BELLE = csm_template(0)
                  dc_BELLE.set_np("nomlization of sig "+names[BELLE], ntau_err[BELLE]/ntau[BELLE], -ntau_err[BELLE]/ntau[BELLE] , NULL,1,NULL,-1)
                  dc_BELLE.make(sigpdf_BELLE,internal_BR*ntau[BELLE],0,1, names[BELLE])
                  dc_no_syst_BELLE.make(sigpdf_BELLE,internal_BR*ntau[BELLE],0,1, names[BELLE])
                  dc_BELLE.add_to(test)
                  dc_no_syst_BELLE.add_to(test_no_syst)
                  test.set_interpolation_style(names[BELLE],1)
                  test_no_syst.set_interpolation_style(names[BELLE],1)
                  
                  returnparameters += [dc_BELLE,dc_no_syst_BELLE]
           #################################################################
           #### LHCb stuff
           if use_LHCb:
                  ############ 2011 #############################################
                  NSIGSYS_2011=0
                  if use_s_alpha:
                         NSIGSYS_2011 = NSIGSYS_2011 + 1
                  if use_sig_sys:
                         for j in sigdofs:
                                NSIGSYS_2011 += 1
                  bkg_2011.add_to(test)
                  bkg_no_syst_2011.add_to(test_no_syst)
                  dc_2011 = csm_template(NSIGSYS_2011)
                  dc_no_syst_2011 = csm_template(0)

                  if use_s_alpha:
                         dc_2011.set_np("nomlization of sig 2011", s_alpha_2011/alpha_2011, -s_alpha_2011/alpha_2011, NULL,1,NULL,-1)
                  if use_sig_sys:
                         for j in range(len(sig_sys_pdfs_2011)):
                                counter = 0
                                dc_2011.set_np("sigsys"+str(counter),sig_sys_corr_2011[j][0],sig_sys_corr_2011[j][1],sig_sys_pdfs_2011[j][0],1,sig_sys_pdfs_2011[j][1],-1)
                  

                  # END OF BCK SYSMATEDICS                                                                                                                         
                  dc_2011.make(sigpdf_2011,internal_BR/alpha_2011,0,1, "2011")
                  dc_no_syst_2011.make(sigpdf_2011,internal_BR/alpha_2011,0,1, "2011")
                  dc_2011.add_to(test)
                  dc_no_syst_2011.add_to(test_no_syst)

                  test.set_interpolation_style("2011",ISTYLE)
                  test_no_syst.set_interpolation_style("2011",ISTYLE)
                  returnparameters += [dc_2011,dc_no_syst_2011]
                  ########## 2012 #################################################
                  NSIGSYS_2012=0
                  if use_s_alpha:
                         NSIGSYS_2012 = NSIGSYS_2012 + 1
                  if use_sig_sys:
                         for j in sigdofs:
                                NSIGSYS_2012 += 1
                  

                  bkg_2012.add_to(test)
                  bkg_no_syst_2012.add_to(test_no_syst)

                  dc_2012 = csm_template(NSIGSYS_2012)
                  dc_no_syst_2012 = csm_template(0)
                  if use_s_alpha:
                         dc_2012.set_np("nomlization of sig 2012",  s_alpha_2012/alpha_2012, -s_alpha_2012/alpha_2012, NULL,1,NULL,-1)
                         
                  if use_sig_sys:
                         for j in range(len(sig_sys_pdfs_2012)):
                                counter = 0
                                dc_2012.set_np("sigsys"+str(counter),sig_sys_corr_2012[j][0],sig_sys_corr_2012[j][1],sig_sys_pdfs_2012[j][0],1,sig_sys_pdfs_2012[j][1],-1)
                                
                  dc_2012.make(sigpdf_2012,internal_BR/alpha_2012,0,1, "2012")
                  dc_no_syst_2012.make(sigpdf_2012,internal_BR/alpha_2012,0,1, "2012")
                  dc_2012.add_to(test)
                  dc_no_syst_2012.add_to(test_no_syst)
                  
                  test.set_interpolation_style("2012",ISTYLE)
                  test_no_syst.set_interpolation_style("2012",ISTYLE)
                  returnparameters += [dc_2012,dc_no_syst_2012]






           print "Before things get nasty"
           CL = mclimit_csm()
           CL.set_null_hypothesis(nul_no_syst)
           CL.set_test_hypothesis(test_no_syst)

           CL.set_null_hypothesis_pe(nul)
           CL.set_test_hypothesis_pe(test)

           if use_BABAR:
                  CL.set_datahist(obs_BABAR, names[BABAR])
           if use_BELLE:
                  CL.set_datahist(obs_BELLE, names[BELLE])
           if use_LHCb:
                  CL.set_datahist(unbhist_2011, "2011")
                  CL.set_datahist(unbhist_2012, "2012")


           CL.set_npe(NTOYS)
           print "Pseudo experiments: ",NTOYS
           CL.run_pseudoexperiments()
           print "After experiments"
           
           cls = CL.cls()
           clb = CL.clb()
           clsb = CL.clsb()

           print "Problem?"
           return CL, test, test_no_syst, returnparameters
                      

#################################################
## SCAN
##################################################
# THIS FUNCTION IS FOR SINGLE JOB!!!
##################################################
def do_scan1(filename):
       tup=XTuple(filename, labels = ["BR/F","ts/F", "chi2/F", "CLs/F","CLb/F", "CLsb/F","CLs_exp_b_med/F" , "CLs_exp_b_p1/F", "CLs_exp_b_p2/F","CLs_exp_b_m1/F" ,"CLs_exp_b_m2/F","CLb_exp_s_med/F" , "CLb_exp_s_p1/F", "CLb_exp_s_p2/F","CLb_exp_s_m1/F" ,"CLb_exp_s_m2/F"])
       print CHANNEL
       c90, c95 = 0.,0
       print filename
       buffer = []
       brloop = []
        #    if external_BR is None
        #  brloop = range(1,20)
        #else:
        #brloop = [external_BR]
        #for j in brloop:
       for j in range(1,2000):
               i =  0.7 + 0.05 * float(j)
               i=i*1.e-9
               print "ps for BR:", i
               if i>5.1e-8:
                        break
               CL,a,b,c = DoTestHyp(i) 
               buffer += [a,b,c]
               print "got CL"
                 #CL.setpxprintflag(1)
               cl = CL.cls()
               clb = CL.clb()
               clsb = CL.clsb()
               tup.fillItem("BR", i)
               tup.fillItem("ts", CL.ts())
               #tup.fillItem("chi2", CL.calc_chi2(testpes[i*1e-9],a2012.DataHist))
               tup.fillItem("CLs",cl)
               tup.fillItem("CLb",clb)
               tup.fillItem("CLsb",clsb)
               tup.fillItem("CLs_exp_b_med",CL.clsexpbmed())
               tup.fillItem("CLs_exp_b_p1",CL.clsexpbp1())
               tup.fillItem("CLs_exp_b_p2",CL.clsexpbp2())
               tup.fillItem("CLs_exp_b_m1",CL.clsexpbm1())
               tup.fillItem("CLs_exp_b_m2",CL.clsexpbm2())
               
               tup.fillItem("CLb_exp_s_med",CL.clbexpsmed())
               tup.fillItem("CLb_exp_s_p1",CL.clbexpsp1())
               tup.fillItem("CLb_exp_s_p2",CL.clbexpsp2())
               tup.fillItem("CLb_exp_s_m1",CL.clbexpsm1())
               tup.fillItem("CLb_exp_s_m2",CL.clbexpsm2())
                 
               print i,cl, clb, clsb, CL.clsexpbmed()
               tup.fill()
       tup.close()
       print "Done? ", filename
##################################################  
# THIS FUNCTION IS FOR PARRAREL JOBS
##################################################
def do_scan(filename, BR):
       tup=XTuple(filename, labels = ["BR/F","ts/F", "chi2/F", "CLs/F","CLb/F", "CLsb/F","CLs_exp_b_med/F" , "CLs_exp_b_p1/F", "CLs_exp_b_p2/F","CLs_exp_b_m1/F" ,"CLs_exp_b_m2/F","CLb_exp_s_med/F" , "CLb_exp_s_p1/F", "CLb_exp_s_p2/F","CLb_exp_s_m1/F" ,"CLb_exp_s_m2/F"])
       print CHANNEL
       c90, c95 = 0.,0
       print filename
       buffer = []
       brloop = []
       #    if external_BR is None  
       #  brloop = range(1,20) 
#       i=2e-9+((5e-8)-(2e-9))*(BR/480)
       i=4.6*1.e-8
       print "ps for BR:", i
       if i>0.9e-7:
              print "Here exiting"+ str(BR)
              sys.exit()
       CL,a,b,c = DoTestHyp(i)
       buffer += [a,b,c]
       print "got CL"
       cl = CL.cls()
       clb = CL.clb()
       clsb = CL.clsb()
       tup.fillItem("BR", i)
       tup.fillItem("ts", CL.ts())
       
       tup.fillItem("CLs",cl)
       tup.fillItem("CLb",clb)
       tup.fillItem("CLsb",clsb)
       tup.fillItem("CLs_exp_b_med",CL.clsexpbmed())
       tup.fillItem("CLs_exp_b_p1",CL.clsexpbp1())
       tup.fillItem("CLs_exp_b_p2",CL.clsexpbp2())
       tup.fillItem("CLs_exp_b_m1",CL.clsexpbm1())
       tup.fillItem("CLs_exp_b_m2",CL.clsexpbm2())
       tup.fillItem("CLb_exp_s_med",CL.clbexpsmed())
       tup.fillItem("CLb_exp_s_p1",CL.clbexpsp1())
       tup.fillItem("CLb_exp_s_p2",CL.clbexpsp2())
       tup.fillItem("CLb_exp_s_m1",CL.clbexpsm1())
       tup.fillItem("CLb_exp_s_m2",CL.clbexpsm2())

#       print i,cl, clb, clsb, CL.clsexpbmed()
       print "BR: ",i," Obs_Cls: ",cl, " Exp_Cls: ", CL.clsexpbmed()
       
       
       tup.fill()
       tup.close()
       print "Done Br= "+str(i)+  filename


NEXPERIMENT=NUMBER_EXPERIMENTS
if use_LHCb:
       NEXPERIMENT=NEXPERIMENT+1
       
outfile='results_'+str(NEXPERIMENT)+'/'+str(DECAY_NUMBER)+"/"+str(BR)
outfile1='mkdir -vp results_'+str(NEXPERIMENT)+'/'+str(DECAY_NUMBER)
import os
os.system(outfile1)
do_scan(outfile, BR)
