.PHONY:	all link combos alucomb pdg aluplot plot 
.PHONY:	update_combos update_alucomb update_pdg update_plot update_all

PLOTS = $(patsubst %.input,%.eps,$(wildcard plot*.input))

CURDIR = $(shell expr $$PWD : '.*/\([^/]*/[^/]*/[^/]*\)')
CURDIRBASE = $(shell expr $$PWD : '.*/\([^/]*\)')
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
	convert $@ $(patsubst %.input,%.gif, $<)
	@echo "file $(CURDIR)/$(patsubst %.input,%.log, $<) created"
	@if [ -w $(HFAGPLOTDIR)_disable_/. ] ; then \
	  /bin/cp -f $(patsubst %.input,%.gif, $<) \
	  $(HFAGPLOTDIR)/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<)); \
	  echo "created" $(HFAGPLOTDIR)/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<)); \
	elif [ -w $(HOME)/local/mounts/slac$(HFAGPLOTDIR)_disable_/. ] ; then \
	  /bin/cp -f $(patsubst %.input,%.gif, $<) \
	  $(HOME)/local/mounts/slac$(HFAGPLOTDIR)/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<)); \
	  echo "created" $(HFAGPLOTDIR)/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<)) ; \
	else \
	  : _disable_ echo "#### scp $(patsubst %.input,%.gif, $<) flora.slac.stanford.edu:$(HFAGPLOTDIR)/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<))"; \
	fi

update_plot: update_plot_all
	@/bin/mkdir -p $(DATADIR)/$(CURDIR)
	@if [ -f plot.log ] ; then \
	  echo "/bin/cp -f plot.log $(DATADIR)/$(CURDIR)" ;\
	  /bin/cp -f plot.log $(DATADIR)/$(CURDIR)/. ;\
	fi

.PHONY: update_plot_all
update_plot_all:  $(PLOTS)
	@/bin/mkdir -p $(DATADIR)/$(CURDIR)
	/bin/cp -f $^ $(DATADIR)/$(CURDIR)
	/bin/cp -f $(patsubst %.eps,%.gif, $^) $(DATADIR)/$(CURDIR)
	@if [ -w $(HFAGPLOTDIR)/. ] ; then \
	  for file in $(patsubst %.eps,%.gif, $^) ; do \
	    echo $${file} | /bin/sed -e 's@plot\(.*\)@/bin/cp plot\1 $(HFAGPLOTDIR)/$(CURDIRBASE)\1@' | sh -v ;\
	  done ;\
	  if [ `uname -s` = "Linux" ] ; then \
	  echo "updating plots gallery with albatross" ;\
	  (cd $(HFAGPLOTDIR); \
	   ../../../tau/bin/albatross --indexname gallery.html --clear; \
	   ../../../tau/bin/albatross --indexname gallery.html --light --noexif) ;\
	  else \
	  echo "you need to be on Linux to use albatross" ;\
	  echo "should do: cd $(HFAGPLOTDIR)" ;\
	  echo "should do: ../../../tau/bin/albatross --indexname gallery.html --clear" ;\
	  echo "should do: ../../../tau/bin/albatross --indexname gallery.html --light --noexif" ;\
	  fi ;\
	elif [ -w $(HOME)/local/mounts/slac$(HFAGPLOTDIR)/. ] ; then \
	  for file in $(patsubst %.eps,%.gif, $^) ; do \
	    echo $${file} | /bin/sed -e 's@plot\(.*\)@/bin/cp plot\1 $(HOME)/local/mounts/slac$$(HFAGPLOTDIR)/$(CURDIRBASE)\1@' | sh -v;\
	  done ;\
	  echo "should do: cd $(HFAGPLOTDIR)" ;\
	  echo "should do: ../../../tau/bin/albatross --indexname gallery.html --clear" ;\
	  echo "should do: ../../../tau/bin/albatross --indexname gallery.html --light --noexif" ;\
	else \
	  echo "cannot write to $(HFAGPLOTDIR)" ;\
	fi

update_combos:
	@if [ -f average.log ] ; then \
	  /bin/mkdir -p $(DATADIR)/$(CURDIR) ; \
	  echo "/bin/cp -f average.log $(DATADIR)/$(CURDIR)" ; \
	  /bin/cp -f average.log $(DATADIR)/$(CURDIR)/. ; \
	fi

update_pdg:
	@if [ -f pdg_average.log ] ; then \
	  /bin/mkdir -p $(DATADIR)/$(CURDIR) ; \
	  echo "/bin/cp -f pdg_average.log $(DATADIR)/$(CURDIR)" ; \
	  /bin/cp -f pdg_average.log $(DATADIR)/$(CURDIR)/. ; \
	fi

update_alucomb:
	@if [ -f average_alucomb.log ] ; then \
	  /bin/mkdir -p $(DATADIR)/$(CURDIR) ; \
	  echo "/bin/cp -f average_alucomb.log $(DATADIR)/$(CURDIR)" ; \
	  /bin/cp -f average_alucomb.log $(DATADIR)/$(CURDIR)/. ; \
	fi

update_all: update_combos update_pdg update_plot update_alucomb

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
