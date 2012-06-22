# $Id: $

.PHONY: all
all: UserManual.pdf

UserManual.pdf: UserManual.tex 
	pdflatex UserManual.tex
	pdflatex UserManual.tex

tree.pdf: tree.tex
	latex tree.tex
	latex tree.tex
	latex tree.tex
	dvips tree.dvi
	ps2pdf tree.ps

UserManual.aux: UserManual.tex
	pdflatex UserManual.tex

UserManual.bbl: UserManual.aux
	bibtex UserManual.aux

.PHONY: clean
clean:
	rm -f UserManual.ps *.dvi UserManual.bbl UserManual.blg UserManual.toc *.aux *.log *~ 

.PHONY: purge
purge: clean
	rm -f UserManual.pdf

