build:
	R CMD Rscript -e "library(knitr); knit('Bluejay.Rnw')"
	latex Bluejay.tex
	bibtex Bluejay.aux
	latex Bluejay.tex
	latex Bluejay.tex
	dvipdf Bluejay.dvi

clean:
	rm -rf *.tex *.dvi *.pdf

