
all:
	pdflatex hfag-standalone.tex
	bibtex   hfag-standalone
	pdflatex hfag-standalone.tex
	pdflatex hfag-standalone.tex

clean:
	/bin/rm -f *.pdf
	/bin/rm -f *~  *.aux  *.bbl  *.blg *.dvi *.log *.spl .*~

