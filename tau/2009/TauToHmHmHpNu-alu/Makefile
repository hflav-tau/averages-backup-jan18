combos: average.input babar.input belle.input
	../../../combos/combos average.input > average.log
	@echo "file average.log has been produced"

pdg:
	bbrroot -l -b -q ../../../Common/pdg_average.cc

babar.input: babar.r
	@if which Rscript > /dev/null 2>&1; then ./babar.r > babar.input; echo converted babar.r to babar.input; fi

belle.input: belle.r
	@if which Rscript > /dev/null 2>&1; then ./belle.r > belle.input; echo converted belle.r to belle.input; fi
