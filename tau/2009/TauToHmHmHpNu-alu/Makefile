combos: average.input babar.input belle.input
	../../../combos/combos average.input > average.log
	@echo "file average.log has been produced"

pdg:
	bbrroot -l -b -q ../../../Common/pdg_average.cc

babar.input: babar.r
	./babar.r > babar.input

belle.input: belle.r
	./belle.r > belle.input
