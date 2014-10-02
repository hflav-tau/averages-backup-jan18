DOC = Summer14
FIGURES = $(wildcard figures/*/*.pdf)
BIBS = $(wildcard */*.bib)
TEXFILES = $(wildcard *.tex)

all:	$(DOC).pdf

$(DOC).bbl $(DOC).pdf : $(DOC).bib $(BIBS) $(FIGURES) $(DOC).tex
	@echo Remaking $@
	@pdflatex $(DOC)
	@pdflatex $(DOC)
	@bibtex $(DOC) 
	@./fix_bbl < $(DOC).bbl > $(DOC).bbl.fixed; mv $(DOC).bbl{.fixed,}
	@pdflatex $(DOC)
	@pdflatex $(DOC)

clean:
	$(RM) $(addprefix $(DOC), .aux .cb .dvi .lof .log .lot .ps .toc .bbl .blg .pdf)
	@$(RM) $(wildcard *~)
	$(RM) $(addprefix $(DOC), -pics.pdf .out .fls .fdb_latexmk -autopp.toc -autopp.out)
