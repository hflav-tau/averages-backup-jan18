#############################################################################
## COMBOS Makefile - created 6 Feb 2003 - dkirkby@uci.edu
##
## Assumes a correctly-configured cernlib script in the path, or else
## set CERNLIBS on the command line.
##
## The choice of C and Fortran compilers and their flags is archictecture
## dependent. The following architectures (as reported by `uname -s`) are
## supported:
##
##  Architecture   Default CC   Default FC
## ----------------------------------------
##  Linux          gcc          g77
##  SunOS          gcc          g77
## ----------------------------------------
##
## Set FCOMP or CCOMP on the command line to select different compilers.
## Set FFLAGS or CFLAGS to override the default options for a compiler. Set
## FEXTRA or CEXTRA to add options to the defaults.
##
#############################################################################
##
## Scientific Linux 5 (A.Lusiani)
## - cernlib as proviced by RPM packages: cernlib, cernlin-devel, cernlib-utils
## - for g77 use RPM package compat-gcc-34-g77
#############################################################################

default : combos

# this command generates the linker options for CERNLIB
CERNLIBS = `cernlib genlib packlib mathlib kernlib`

# the objects we need to build COMBOS
OBJS = cmain.o master.o chi2_n_sym.o user_dms_sens.o adju_a_for_fbs.o avg_amplitude.o blue.o blue_asym.o chi2_asym_min.o chi2_sym.o chi2_sym_circ.o chi2_sym_circ_g.o chi2_sym_indep.o chi2_sym_min.o chi2_sym_rvk.o circularity.o circularity_g.o combos.o comby_lifetime.o final_fractions.o minuit_achi2.o minuit_chi2.o new_combination.o new_preparation.o new_termination.o reader.o set_dms_limit.o single.o writer.o

# get our architecture from uname
OSNAME:= $(shell uname -s)

# set default compilers and flags for recognized platforms

ifeq ($(OSNAME),Linux)
  ifndef CCOMP
    CCOMP = gcc -m32
# Option -m32 added on Dec 9, 2007 by O. Schneider to allow use on lxplus.cern.ch with Scientific Linux SLC 4
  endif
  ifndef FCOMP
    FCOMP = g77 -m32
# Option -m32 added on Dec 9, 2007 by O. Schneider to allow use on lxplus.cern.ch with Scientific Linux SLC 4
  endif
  # set default flags
  ifeq ($(FCOMP),g77 -m32)
    FFLAGS = -g -w -fno-automatic -fno-globals -Wno-globals -ffixed-line-length-none -funderscoring -fcase-lower -ffloat-store
  endif
# add default options for gfortran (A.Lusiani, March 2010)
  ifeq ($(shell expr '$(FCOMP)' : '\(gfortran\)'),gfortran)
    FFLAGS = -g -w -m32 -fno-automatic -ffixed-line-length-none -funderscoring -ffloat-store
  endif
endif

ifeq ($(OSNAME),SunOS)
  # set default compilers
  ifndef CCOMP
    CCOMP = gcc
  endif
  ifndef FCOMP
    FCOMP = f77
  endif
  # set default flags
  ifeq ($(FCOMP),g77)
    FFLAGS = -g -w -fno-automatic -fno-globals -Wno-globals -ffixed-line-length-none -funderscoring -fcase-lower -ffloat-store
  endif
  # the following has been tested to work on flora.slac.stanford.edu
  ifeq ($(FCOMP),f77)
    FFLAGS = -g -w -e -f77 -fnonstd
  endif
  # the following has been tested to work on flora.slac.stanford.edu
  LDFLAGS = -L/opt/SUNWspro/lib/ -lf77compat -lfsu -lsunmath
endif

# copy CCOMP,FCOMP to CC,FC (if set)
ifdef CCOMP
  CC = $(CCOMP)
endif
ifdef FCOMP
  FC = $(FCOMP)
endif

# add any extra compiler flags
ifdef CEXTRA
  CFLAGS+= $(CEXTRA)
endif
ifdef FEXTRA
  FFLAGS+= $(FEXTRA)
endif

combos : $(OBJS)
	$(FC) -o combos $(OBJS) $(CERNLIBS) $(LDFLAGS)

# --- regression testing ----------------------------------------------

#DIFF = gdiff -I"COMBOS version"
DIFF = diff

test :
	@echo "##################################################"
	@echo "Running the program with the demo input..."
	@cd demo ; ../combos demo.input > demo.output
	@echo "##################################################"
	@echo "Comparing new outputs to reference outputs..."
	-cd demo ; $(DIFF) fort.10 fort.10.ref
	-cd demo ; $(DIFF) fort.11 fort.11.ref
	-cd demo ; $(DIFF) fort.21 fort.21.ref
	-cd demo ; $(DIFF) demo.output demo.output.ref
	@echo "##################################################"

# --- housekeeping ---------------------------------------------------

clean :
	rm -f *~ *.o combos demo/fort.10 demo/fort.11 demo/fort.21 demo/demo.output

tarfile :
	tar -zcf combos.tar.gz GNUmakefile *.inc *.f *.c demo/*.input demo/*.ref kumac/*.kumac
