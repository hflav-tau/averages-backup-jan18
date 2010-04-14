#//////////////////////////////////////////////////////////////////////////////
#
# HFAG - Tau repository
#

#//////////////////////////////////////////////////////////////////////////////
# how to create your private working area

svn --username <google username> co https://hfag.googlecode.com/svn/trunk hfag
cd hfag

#//////////////////////////////////////////////////////////////////////////////
# organization of the working space

- the subversion repository contains code and cards
- log files and plots are stored in a data directory

a reference working space is kept at SLAC at:
/afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-googlecode/

a reference data directory is kept at SLAC at:
/afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data

the must be a "Data" directory at the root of the working space,
i.e. hfag/Data in you private working space
at SLAC you can set "Data" to a link to the reference data directory:

cd /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-googlecode/
ls -l Data
... Data -> /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data

#//////////////////////////////////////////////////////////////////////////////
# shell environment configuration

please do "source Common/config.csh" to configure your C shell
this script will add to your PATH
- the Root path in flora and some other SLAC nodes
- the R path in flora, iris, and other SLAC nodes

#//////////////////////////////////////////////////////////////////////////////
# examples

# --- To compile combos
cd combos ; gmake clean ; gmake ; cd -

# --- create a private data directory
mkdir Data

# --- link to the reference data directory
ln -s /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data Data

# --- to get information on Makefile targets
cd tau/2009; make; cd -

# --- to run combos for all tau 2009 HFAG averages
cd tau/2009; make combos; cd -

# --- To run combos in a single directory
cd tau/2009/TauToPimKzsNu ; make ; cd -

# --- To obtain PDG style averages
cd tau/2009/TauToPimKzsNu ; make pdg ; cd -

# --- it is useful to have a "log" link pointing to the Data directory
cd tau/2009/TauToPimKzsNu ; make link ; cd -

#//////////////////////////////////////////////////////////////////////////////
#
# Naming Convention adopted so that we can easily assimiliate files into latex: 
#

For daughters of tau decays:
First negative charge, then positive charge, then neutrals, then neutrino's; within each category : low to high mass.

Eg. tau- -> e- mu- pi- K- rho- K*-    e+ mu+ pi+ K+ rho+ K*+    gamma pi0 K0S K0L eta rho0 omega K*0    eta' phi nuebar numbar nut
    Tau  To Em Mm  Pim Km Rhom Kstarm Ep Mp  Pip Kp Rhop Kstarp Gamma Piz Kzs Kzl Eta Rhoz Omega Kstarz Etap Phi Nuebar Numbar Nu

When writing out decaymodes, no space needs to be inserted between particle names, but first letter is capitalized.

Avoid latex conflicts by replacing numbers [0-9] in names with something appropriate.

We follow a style similar to hepnames/pennames. 
http://mirror.ctan.org/macros/latex/contrib/hepnames/hepnames.pdf

#//////////////////////////////////////////////////////////////////////////////
#
# alucomb.r (by A.Lusiani)
#

The script alucomb.r is an almost complete replacement for Combos,
with some additional features.  From the tau/2009/* directories it can
be executed with:

shell> ../../../Common/bin/alucomb.r <combos .input file>

- can read Combos .input files

Additional features:
- uses the card SUMOFMEAS to indicate measurements that correspond to
  sum of quantities we want to average
- uses the card ERROR_CORR_WITH to indicate total (stat. + syst.)
  correlations between measurements
- computes S-factors to inflate the errors to account for a large chi-square

If available, it uses the "maxLik" R package, that can be installed by
the user as follows:
shell> R
> install.packages("maxLik")
> q()

#//////////////////////////////////////////////////////////////////////////////
#
# aluplot.r (by A.Lusiani)
#

Produces a plot.input file using information from the current directory and from
global files. It relies on the following files:

- ./average.input
- ./log/average.log
- ../../../tau/2009/Common/pdg_averages.input
- ../../../tau/2009/Common/results_asymm_errors.input

Invoke either as follows:

shell> ../../../Common/bin/aluplot.r

or as follows:

shell> make aluplot

To produce plots from the plot*.input files please do:

shell> make plot

#//////////////////////////////////////////////////////////////////////////////
#
# tip: remotely mount the slac reference directory
#

shell> sshfs -oworkaround=rename flora.slac.stanford.edu:/ ~/local/mounts/slac
shell> ln -s ~/local/mounts/slac/afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data ../../../Data

#//////////////////////////////////////////////////////////////////////////////
#
# end of 2009 paper CVS access
#

outside of SLAC AFS:
- env CVS_RSH=ssh cvs -d :ext:flora.slac.stanford.edu:/afs/slac.stanford.edu/www/xorg/hfag/cvs co Papers/EndOfYear09

in SLAC AFS:
- cvs -d /afs/slac.stanford.edu/www/xorg/hfag/cvs co Papers/EndOfYear09
