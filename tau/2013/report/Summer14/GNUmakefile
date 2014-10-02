DOC = Summer14
FIGURES = $(wildcard figures/*/*.pdf)
BIBS = $(wildcard */*.bib)
TEXFILES = $(wildcard *.tex)

all:	$(DOC).pdf

$(DOC).bbl $(DOC).pdf : $(DOC).bib $(BIBS) $(FIGURES) $(DOC).tex
	@echo Remaking $@
	@pdflatex -shell-escape $(DOC)
	@pdflatex -shell-escape $(DOC)
	@bibtex $(DOC) 
	@./fix_bbl < $(DOC).bbl > $(DOC).bbl.fixed; mv $(DOC).bbl{.fixed,}
	@pdflatex -shell-escape $(DOC)
	@pdflatex -shell-escape $(DOC)

clean:
	$(RM) $(addprefix $(DOC), .aux .cb .dvi .lof .log .lot .ps .toc .bbl .blg .pdf)  $(wildcard *~)

