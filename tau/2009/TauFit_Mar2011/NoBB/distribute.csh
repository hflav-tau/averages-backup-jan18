#!/bin/sh

for dirname in \
   NoBB_nocorr \
   NoBB_lepcorr \
   NoBB_nocorr_addratio \
   NoBB_lepcorr_addratio \
   NoBB_nocorr_addratio_2 \
   NoBB_lepcorr_addratio_2 \
   BB_nocorr_addratio_3 \
   BB_lepcorr_addratio_3 \
   BB_lepcorr_addratio_4 \
   BB_bbrcorr_addratio_4 \
   BB_addpre-bbrcorr_addratio_4 \
   NoBB_addpre-bbrcorr_addratio_4 \
  ; do
 
  cd /afs/slac.stanford.edu/u/br/swaban/public/hfag/tau/2009/TauFit_Mar2011/NoBB/${dirname}
  pwd
#  cp -f /afs/slac.stanford.edu/u/br/swaban/public/hfag/tau/2009/TauFit_Mar2011/NoBB/readpdg.cc .
#  cp -f /afs/slac.stanford.edu/u/br/swaban/public/hfag/tau/2009/TauFit_Mar2011/NoBB/any.com .
#  cp -f /afs/slac.stanford.edu/u/br/swaban/public/hfag/tau/2009/TauFit_Mar2011/NoBB/run_readpdg .
  ./any.com readpdg
  source run_readpdg
  cp -p *.log log
  cp -p *.results log
#  cp -f /afs/slac.stanford.edu/u/br/swaban/public/hfag/tau/2009/TauFit_Mar2011/NoBB/summary0.tex .
#  cp -f /afs/slac.stanford.edu/u/br/swaban/public/hfag/tau/2009/TauFit_Mar2011/NoBB/summary0.csh .
#  cp -f /afs/slac.stanford.edu/u/br/swaban/public/hfag/tau/2009/TauFit_Mar2011/NoBB/maketable0.csh . 
  rm -f what_is_this.tex 
  if [ "${dirname}" == "NoBB_nocorr" ] ; then
    echo '\noindent Using inputs from: the non-B-Factory measurements of single quantities and excluding all correlations.' > what_is_this.tex
  fi
  if [ "${dirname}" == "NoBB_lepcorr" ] ; then
    echo '\noindent Using inputs from: the non-B-Factory measurements of single quantities and including leptonic correlations only.' > what_is_this.tex
  fi
  if [ "${dirname}" == "NoBB_nocorr_addratio" ] ; then
    echo '\noindent Using inputs from: the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$ from CLEO only) and excluding all correlations.' > what_is_this.tex
  fi
  if [ "${dirname}" == "NoBB_lepcorr_addratio" ] ; then
    echo '\noindent Using inputs from: the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$ from CLEO only) and including leptonic correlations only.' > what_is_this.tex
  fi
  if [ "${dirname}" == "NoBB_nocorr_addratio_2" ] ; then
    echo '\noindent Using inputs from: the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$ from CLEO and ARGUS) and excluding all correlation.' > what_is_this.tex
  fi
  if [ "${dirname}" == "NoBB_lepcorr_addratio_2" ] ; then
    echo '\noindent Using inputs from: the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$ from CLEO and ARGUS) and including leptonic correlations only.' > what_is_this.tex
  fi
  if [ "${dirname}" == "BB_nocorr_addratio_3" ] ; then
      echo '\noindent Using inputs from: $\mathrm{B_\mu/B_e}$ measurement from BaBar and the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$) and excluding all correlations.' > what_is_this.tex
  fi
  if [ "${dirname}" == "BB_lepcorr_addratio_3" ] ; then
    echo '\noindent Using inputs from: $\mathrm{B_\mu/B_e}$ measurement from BaBar and the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$) and including leptonic correlations only.' > what_is_this.tex
  fi
  if [ "${dirname}" == "BB_lepcorr_addratio_4" ] ; then
    echo '\noindent Using inputs from: $\mathrm{B_\mu/B_e}$,  $\mathrm{B_\pi/B_e}$ and $\mathrm{B_K/B_e}$ measurements from BaBar excluding all correlations; and the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$) and including leptonic correlations only.' > what_is_this.tex
  fi
  if [ "${dirname}" == "BB_bbrcorr_addratio_4" ] ; then
    echo '\noindent Using inputs from: $\mathrm{B_\mu/B_e}$,  $\mathrm{B_\pi/B_e}$ and $\mathrm{B_K/B_e}$ measurements from BaBar including their inter-correlations; and the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$) and including leptonic correlations only.' > what_is_this.tex
  fi
  if [ "${dirname}" == "BB_addpre-bbrcorr_addratio_4" ] ; then
    echo '\noindent Using inputs from: $\mathrm{B_\mu/B_e}$,  $\mathrm{B_\pi/B_e}$ and $\mathrm{B_K/B_e}$ measurements from BaBar including their inter-correlations; and the non-B-Factory measurements of single quantities and ratios (eg. $\mathrm{B_\mu/B_e}$) and including all inter-correlations.' > what_is_this.tex
  fi
  if [ "${dirname}" == "NoBB_addpre-bbrcorr_addratio_4" ] ; then
    echo '\noindent Using inputs from: all non-B-Factory measurements including all correlations.' > what_is_this.tex
  fi
  ./summary0.csh
  cp -p summary0.pdf log
  cd ..
done
