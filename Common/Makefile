.PHONY:	all combos pdg plot update_log update_pdglog update_plotlog update_all link alucomb aluplot

PLOTS = $(patsubst %.input,%.eps,$(wildcard plot*.input))

CURDIR = $(shell expr $$PWD : '.*/\([^/]*/[^/]*/[^/]*\)')
#
# we store data files such as plots and log files in a separate directory tree
# which is either the:
# - ../../../Data directory if it exists and is writable
# - SLAC HFAG-tau area /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data/
# - if none of the above is writable, a new ../../../Data is created and used
#
DATADIR =     $(shell \
	if [ -w ../../../Data/. ] ; then \
	  echo ../../../Data ; \
	elif [ -w /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data/. ] ; then \
	  echo /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data ; \
	else \
	  mkdir ../../../Data ; \
	  echo ../../../Data ; \
	fi ; \
	)
HFAGPLOTDIR = /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009

test:
	@echo $(CURDIR)
	@echo $(DATADIR)
	@if [ -w ../../../Data/. ] ; then \
	  _datadir=../../../Data/. ; \
	elif [ -w $(DATADIR)/. ] ; then \
	  _datadir=$(DATADIR) ; \
	fi

combos:
	@if [ -f average.input ] ; then \
	  /bin/rm -f average.log ; \
	  echo "../../../combos/combos average.input > average.log" ; \
	  ../../../combos/combos average.input > average.log ; \
	  echo "file $(CURDIR)/average.log created" ; \
	fi

alucomb:
	@if [ -f average.input ] ; then \
	  /bin/rm -f average_alucomb.log ; \
	  echo "../../../Common/bin/alucomb.r average.input > average_alucomb.log 2>&1" ; \
	  ../../../Common/bin/alucomb.r average.input > average_alucomb.log ; \
	  echo "file $(CURDIR)/average_alucomb.log created" ; \
	fi

aluplot:
	@if [ -f average.input ] ; then \
	  /bin/rm -f aluplot.log ; \
	  echo "../../../Common/bin/aluplot.r average.input > aluplot.log 2>&1" ; \
	  ../../../Common/bin/aluplot.r average.input > aluplot.log ; \
	  echo "file $(CURDIR)/aluplot.log created" ; \
	fi

pdg:
	@/bin/rm -f pdg_average.log
	@if [ -f pdg_average.input ] ; then \
	  /bin/rm -f pdg_average.log; \
	  echo "root -l -b -q ../../../Common/pdg_average.cc > pdg_average.log" ; \
	  root -l -b -q ../../../Common/pdg_average.cc > pdg_average.log; \
	  echo "file $(CURDIR)/pdg_average.log created" ; \
	fi

plot:	$(PLOTS)

$(PLOTS): %.eps: %.input
	@/bin/rm -f $(patsubst %.input,%.log, $<)
	root -l -b -q ../../../Common/plot.cc'("$<")' > $(patsubst %.input,%.log, $<)
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/$(patsubst %.input,%.log, $<) created:'
	@convert $@ $(patsubst %.input,%.gif, $<)
	@if [ -w /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/. ] ; then \
	  /bin/cp -f $(patsubst %.input,%.gif, $<) \
	  /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<)); \
	  echo "created" /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<)); \
	elif [ -w $(HOME)/local/mounts/slac/afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/. ] ; then \
	  /bin/cp -f $(patsubst %.input,%.gif, $<) \
	  $(HOME)/local/mounts/slac/afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<)); \
	  echo "created" /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<)) ; \
	else \
	  echo "#### scp $(patsubst %.input,%.gif, $<) flora.slac.stanford.edu:/afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<))"; \
	fi
#
# exploit remote access to SLAC
# sshfs -oworkaround=rename iris.slac.stanford.edu:/ ~/local/mounts/slac
#
plot2web: $(PLOTS)
	@if [ -w $(HOME)/local/mounts/slac//afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/. ] ; then \
	  echo $(patsubst %.eps,%.gif, $^) | tr ' ' "\n" | xargs -i \
	    echo "{}" | sed -e 's@plot\(.*\)@/bin/cp plot\1 $(HOME)/local/mounts/slac//afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/'`basename ${PWD}`'\1@' | sh -v ; \
	fi

update_log:
	@if [ -f average.log ] ; then \
	  /bin/mkdir -p $(DATADIR)/$(CURDIR) ; \
	  echo "/bin/cp -f average.log $(DATADIR)/$(CURDIR)" ; \
	  /bin/cp -f average.log $(DATADIR)/$(CURDIR)/. ; \
	fi

update_pdglog:
	@if [ -f pdg_average.log ] ; then \
	  /bin/mkdir -p $(DATADIR)/$(CURDIR) ; \
	  echo "/bin/cp -f pdg_average.log $(DATADIR)/$(CURDIR)" ; \
	  /bin/cp -f pdg_average.log $(DATADIR)/$(CURDIR)/. ; \
	fi

update_plotlog:
	@if [ -f plot.log ] ; then \
	  /bin/mkdir -p $(DATADIR)/$(CURDIR) ; \
	  echo "/bin/cp -f plot.log $(DATADIR)/$(CURDIR)" ; \
	  /bin/cp -f plot.log $(DATADIR)/$(CURDIR)/. ; \
	fi

update_alucomblog:
	@if [ -f average_alucomb.log ] ; then \
	  /bin/mkdir -p $(DATADIR)/$(CURDIR) ; \
	  echo "/bin/cp -f average_alucomb.log $(DATADIR)/$(CURDIR)" ; \
	  /bin/cp -f average_alucomb.log $(DATADIR)/$(CURDIR)/. ; \
	fi

update_all: update_log update_pdglog update_plotlog update_alucomblog

all: combos pdg plot update_all

link:
	@if [ -d "../../../Data" ] ; then \
	  /bin/mkdir -p "../../../Data/$(CURDIR)" ;\
	  ln -sf "../../../Data/$(CURDIR)" log ;\
	  echo "log is soft link to ../../../Data/$(CURDIR)" ;\
	elif [ -d /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data ] ; then \
	  /bin/mkdir -p "/afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data/$(CURDIR)" ;\
	  ln -fs /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data/$(CURDIR) log ;\
	else \
	  echo "no data directory found" ;\
	  echo "either create ../../../Data" ;\
	  echo "or make accessible /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data" ;\
	fi
