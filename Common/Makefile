combos:
	../../../combos/combos average.input > average.log
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/average.log created"

pdg:
	root -l -b -q ../../../Common/pdg_average.cc

plot:
	root -l -b -q ../../../Common/plot.cc

update_log:
	@cp -f average.log log/
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/average.log copied to log/"
