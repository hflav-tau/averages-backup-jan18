.PHONY:	all combos alucomb pdg aluplot plot 
.PHONY:	update_combos update_alucomb update_pdg update_plot update_all

PLOTS = $(patsubst %.input,%.eps,$(wildcard plot*.input))

CURDIR = $(shell expr $$PWD : '.*/\([^/]*/[^/]*/[^/]*\)')
CURDIRBASE = $(shell expr $$PWD : '.*/\([^/]*\)')
#
# we store data files such as plots and log files in a separate directory tree
# in ../../..Data which is:
# - either a link to local private directory like ../../../Data-local/
# - or a link to the SLAC data reference directory /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data/
#
DATADIR =     $(shell \
	if [ ! -w "../../../Data/." ] ; then \
	  /bin/mkdir -p "../../../Data-local" ;\
	  /bin/ln -sTf Data-local "../../../Data" ;\
	fi ;\
	echo "../../../Data" ;\
	)
HFAGPLOTDIR = /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009
HFAGDATADIR = /afs/slac.stanford.edu/www/xorg/hfag/tau/hfag-data
SSHFSDIR = $(HOME)/local/mounts/slac

# --- show what data directory you are using
.PHONY: test
test:
	@echo "you are in: $(CURDIR)"
	@ll="`/bin/ls -ld "$(DATADIR)"`" ;\
	if ll="`expr "$$ll" : '.* -> \(.*\)'`" ; then \
	  : ;\
	else \
	  echo else ;\
	  ll="$(DATADIR)" ;\
	fi ;\
	if expr "$$ll" : '/.*' >/dev/null ; then \
	  : ;\
	else \
	  ll="../../../$$ll" ;\
	fi ;\
	echo "data dir: $$ll"

# --- set data directory to local private area
.PHONY: local
local:
	@/bin/mkdir -p "../../../Data-local"
	@/bin/ln -sTf "Data-local" "../../../Data"
	@echo "Data reference directory now local in ../../../Data-local"

# --- set data directory to the SLAC HFAG refeerence data directory
.PHONY: slac
slac:
	@if [ -w $(HFAGDATADIR)/. ] ; then \
	  /bin/ln -sTf $(HFAGDATADIR) "../../../Data" ;\
	  echo "data directory set to SLAC" ;\
	elif [ -w $(SSHFSDIR)$(HFAGDATADIR) ] ; then \
	  /bin/ln -sTf $(SSHFSDIR)$(HFAGDATADIR) "../../../Data" ;\
	  echo "data directory set to SLAC via SSHFS" ;\
	else \
	  echo "no write access to $(HFAGDATADIR)" ;\
	fi

# --- mount SLAC AFS via sshfs
.PHONY: slac_mount
slac_mount:
	sshfs -oworkaround=rename flora.slac.stanford.edu:/ $(HOME)/local/mounts/slac

# --- un-mount SLAC AFS via sshfs
.PHONY: slac_umount
slac_umount:
	fusermount -u $(HOME)/local/mounts/slac

# --- sync all hfag-data tree into ../../../Data-local/.
.PHONY: data_download
data_download:
	@if [ -r $(HFAGDATADIR)/. ] ; then \
	  echo rsync -avz --force --delete $(HFAGDATADIR)/. ../../../Data-local/. ;\
	  rsync -avz --force --delete $(HFAGDATADIR)/. ../../../Data-local/. ;\
	else \
	  echo rsync -essh -avz --force --delete flora.slac.stanford.edu:$(HFAGDATADIR)/. ../../../Data-local/. ;\
	  rsync -essh -avz --force --delete flora.slac.stanford.edu:$(HFAGDATADIR)/. ../../../Data-local/. ;\
	fi

# --- (currently disabled) upload local ../../../Data-local/. tree into hfag-data/.
.PHONY: data_upload
data_upload:
	@if [ -r $(HFAGDATADIR)/. ] ; then \
	  echo rsync -avz $(HFAGDATADIR)/. ../../../Data-local/. ;\
	  read -a REPLY -n1 -p "  <y> to confirm, anything else to abort: " ;\
	  if [ $$REPLY = "y" -o $$REPLY = "Y" ] ; then \
	    echo " OK" ;\
	    rsync -n -avz ../../../Data-local/. $(HFAGDATADIR)/. ;\
	  fi ;\
	else \
	  echo rsync -essh -avz ../../../Data-local/. flora.slac.stanford.edu:$(HFAGDATADIR)/. ;\
	  read -a REPLY -n1 -p "  <y> to confirm, anything else to abort: " ;\
	  if [ $$REPLY = "y" -o $$REPLY = "Y" ] ; then \
	    echo " OK" ;\
	    rsync -n -essh -avz ../../../Data-local/. flora.slac.stanford.edu:$(HFAGDATADIR)/. ;\
	  fi ;\
	fi

# --- create "log" as link to the respective dir in the data directory
.PHONY: link
link:
	@if [ -d "../../../Data/." ] ; then \
	  /bin/mkdir -p "../../../Data/$(CURDIR)" ;\
	  /bin/ln -sTf "../../../Data/$(CURDIR)" log ;\
	  echo "log is soft link to ../../../Data/$(CURDIR)" ;\
	else \
	  echo "no data directory found" ;\
	  echo "either \"make local\" or \"make slac\"" ;\
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

alucomb_sum:
	@if [ -f average.input ] ; then \
	  /bin/rm -f average_alucomb.log ; \
	  echo "../../../Common/bin/alucomb.r average.input > average_alucomb.log 2>&1" ; \
	  ../../../Common/bin/alucomb.r average.input > average_alucomb.log ; \
	  echo "file $(CURDIR)/average_alucomb.log created" ; \
	  echo "" >> average_alucomb.log ; \
	  echo "sum of all measurements" >> average_alucomb.log ; \
	  echo "" >> average_alucomb.log ; \
	  ../../../Common/bin/aluelab-sum.r average_alucomb.log >> average_alucomb.log ; \
	  echo "file $(CURDIR)/average_alucomb.log: added sum of measurements" ; \
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
	@echo "do \"make update_plot\" to save plots in reference data dir"

$(PLOTS): %.eps: %.input ../../../Common/plot.cc
	@/bin/rm -f $(patsubst %.input,%.log, $<)
	root -l -b -q ../../../Common/plot.cc'("$<")' > $(patsubst %.input,%.log, $<)
	convert $@ $(patsubst %.input,%.gif, $<)
	@echo "file $(CURDIR)/$(patsubst %.input,%.log, $<) created"

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
	elif [ -w $(SSHFSDIR)$(HFAGPLOTDIR)/. ] ; then \
	  for file in $(patsubst %.eps,%.gif, $^) ; do \
	    echo $${file} | /bin/sed -e 's@plot\(.*\)@/bin/cp plot\1 $(SSHFSDIR)$(HFAGPLOTDIR)/$(CURDIRBASE)\1@' | sh -v ;\
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
