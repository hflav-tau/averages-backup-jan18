all:
	pdflatex tau-report-lfv.tex
	bibtex   tau-report-lfv
	pdflatex tau-report-lfv.tex
	pdflatex tau-report-lfv.tex

clean:
	/bin/rm -f tau-report-lfv.pdf
	/bin/rm -f *~  *.aux  *.bbl  *.blg *.dvi *.log *.spl *.out *.toc .*~
