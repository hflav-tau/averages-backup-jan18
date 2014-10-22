LATEX = latex
HEVEA = hevea
HEVEAOPTS = -fix -O
BIBTEX = bibtex
BIBHVA = bibhva
HACHA = hacha
DIR = tau-report-web
DOC = tau-report-web

PLOTS = \
  plot-taulife-hfag-summer2014.pdf \
  hfag-tau-vus-plot.pdf \
  hfag-tau-lfv-plot.pdf \
  TauLFV_UL_2013001_averaged.pdf

PLOTSPDF = $(patsubst %.pdf, ../figures/tau/%.pdf, $(PLOTS))
PLOTSPNG = $(patsubst %.pdf, ../figures/tau/%.png, $(PLOTS))

#
# latex -recorder hfag-tau-web-report.tex
# grep INPUT hfag-tau-web-report.fls | grep -v "INPUT /" | grep -v aux | grep -v bbl | sed -e "s/INPUT //g" | uniq
#

DOCINCL = \
  ../symbols.tex \
  tau-defs.tex \
  tau-br-fit-aleph-hcorr_13.tex \
  tau-br-fit-aleph-hcorr_13-vadirect-elab.tex \
  tau-sections.tex \
  tau-lfv.tex \
  tau-lfv-combinations.tex \
  tau-br-fit-meas-byref.tex \
  tau-lfv-limits-plot.tex \
  tau-lfv-combinations-plot.tex \
  $(PLOTSPDF)

BASE = $(DIR)/$(DOC)
INSTALLDIR = www-slac/summer-2014

.PHONY: html
html: $(DIR)/index.html $(BASE)-style.css

$(DIR)/index.html: $(BASE).html
	$(HACHA) -o $(DIR)/index.html $(BASE).html
	/bin/rm -f $(DIR)/{previous_motif,next_motif,contents_motif}.gif
	/bin/rm -f $(DOC).image.out

$(BASE).html: $(DOC).hva $(DOC).tex $(DOCINCL) $(PLOTSPNG)
	@/bin/mkdir -p $(DIR)
	$(HEVEA) $(HEVEAOPTS) png.hva $(DOC).hva -o $(BASE).html $(DOC).tex
	$(BIBHVA) $(BASE)
	$(HEVEA) $(HEVEAOPTS) png.hva $(DOC).hva -o $(BASE).html $(DOC).tex
	/bin/cp -f $(PLOTSPDF) $(PLOTSPNG) $(DIR)/

$(DIR)/%-style.css: %-style.css
	/bin/cp $< $@

.PHONY: text
text: $(BASE).txt

$(BASE).txt: $(DOC).hva $(DOC).tex $(DOCINCL)
	@/bin/mkdir -p $(DIR)
	$(LATEX) $(DOC)
	$(BIBTEX) $(DOC)
	/bin/mv $(DOC).bbl $(DIR)/
	$(HEVEA) $(HEVEAOPTS) $(DOC).hva -text -o $(BASE).txt $(DOC).tex
	/bin/rm -f $(DOC).{haux,htoc,bbl,aux,dvi,log,blg,out} $(DOC).image.out

.PHONY: pdf
pdf: $(DOC).pdf

$(DOC).pdf: $(DOC).tex $(DOCINCL)
	latexmk -pdf $(DOC).tex
	@/bin/mkdir -p $(DIR)
	@/bin/mv -f $(DOC).pdf $(DIR)/
	@/bin/rm -f $(DOC).pdf
	@/bin/ln $(DIR)/$(DOC).pdf $(DOC).pdf
	@echo file $(DOC).pdf produced

.PHONY: purge
purge:
	/bin/rm -f $(BASE).{html,haux,htoc,hind,blg,bbl,hbbl} $(BASE).image.* $(DOC).pdf

.PHONY: clean
clean: slac_umount
	-/bin/rm -f $(DIR)/* $(DOC).pdf

.PHONY: slac_mount
slac_mount:
	/bin/mkdir -p www-slac
	sshfs -oworkaround=rename iris.slac.stanford.edu:/afs/slac.stanford.edu/www/xorg/hfag/tau www-slac

.PHONY: slac_umount
slac_umount:
	-fusermount -u www-slac
	-/bin/rmdir www-slac

.PHONY: install
install: pdf html purge
	/bin/rm -f $(DIR)/$(DOC)-style.css
	rsync -rv --force --delete $(DIR)/. $(INSTALLDIR)/.
