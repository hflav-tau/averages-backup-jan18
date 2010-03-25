.PHONY:	combos pdg update_log update_pdglog update_plotlog update_all

combos:
	@/bin/rm -f average.log
	../../../combos/combos average.input > average.log
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/average.log created:'

pdg:
	@/bin/rm -f pdg_average.log
	root -l -b -q ../../../Common/pdg_average.cc > pdg_average.log
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/pdg_average.log created:'

plot:
	@/bin/rm -f plot.log
	root -l -b -q ../../../Common/plot.cc > plot.log
	@/bin/pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:file \1/plot.log created:'
	@/usr/bin/convert plot.eps plot.gif 
	@/bin/cp -f plot.gif /afs/slac.stanford.edu/www/xorg/hfag/tau/plots/2009/`basename ${PWD}`.gif

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
