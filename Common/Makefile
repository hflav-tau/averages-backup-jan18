.PHONY:	combos pdg plot update_log update_pdglog update_plotlog update_all

PLOTS = $(patsubst %.input,%.eps,$(wildcard plot*.input))

combos:
	@/bin/rm -f average.log
	../../../combos/combos average.input > average.log
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/average.log created:'

pdg:
	@/bin/rm -f pdg_average.log
	@if [ -f pdg_average.input ] ; then \
	  /bin/rm -f pdg_average.log; \
	  root -l -b -q ../../../Common/pdg_average.cc > pdg_average.log; \
	  /bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/pdg_average.log created:'; \
	fi

plot:	$(PLOTS)

$(PLOTS): %.eps: %.input
	@/bin/rm -f $(patsubst %.input,%.log, $<)
	root -l -b -q ../../../Common/plot.cc'("$<")' > $(patsubst %.input,%.log, $<)
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/$(patsubst %.input,%.log, $<) created:'
	@/usr/bin/convert $@ $(patsubst %.input,%.gif, $<)
	@if [ -w /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009 ] ; then \
	  /bin/cp -f $(patsubst %.input,%.gif, $<) \
	    /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/`basename ${PWD}`$(patsubst %.input,%.gif, $<); \
	  echo "copied `basename ${PWD}`$(patsubst %.input,%.gif, $<) to SLAC HFAG tau/plots/2009/" ; \
	else \
	  echo "#### scp $(patsubst %.input,%.gif, $<) flora.slac.stanford.edu:/afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/`basename ${PWD}`$(subst plot,,$(patsubst %.input,%.gif, $<))"; \
	fi

update_log:
	@/bin/cp -f average.log log/
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/average.log copied to log/:'

update_pdglog:
	@/bin/cp -f pdg_average.log log/
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/pdg_average.log copied to log/:'

update_plotlog:
	@/bin/cp -f plot.log log/
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/plot.log copied to log/:'

update_all: update_log update_pdglog update_plotlog
