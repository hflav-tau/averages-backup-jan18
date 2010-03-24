combos:
	rm -f average.log ; ../../../combos/combos average.input > average.log
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/average.log created"

pdg:
	rm -f pdg_average.log ; root -l -b -q ../../../Common/pdg_average.cc > pdg_average.log
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/pdg_average.log created"

plot:
	rm -f plot.log ; root -l -b -q ../../../Common/plot.cc > plot.log
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/plot.log created"

update_log:
	@cp -f average.log log/
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/average.log copied to log/"

update_pdglog:
	@cp -f pdg_average.log log/
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/pdg_average.log copied to log/"

update_plotlog:
	@cp -f plot.log log/
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/plot.log copied to log/"

update_all:
	@cp -f average.log log/
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/average.log copied to log/"
	@cp -f pdg_average.log log/
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/pdg_average.log copied to log/"
	@cp -f plot.log log/
	@pwd | sed -e 's:.*/\([^/]*/[^/]*/[^/]*\)$$:\1:' | xargs -i echo "file {}/plot.log copied to log/"
