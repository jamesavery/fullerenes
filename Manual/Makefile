
.PHONY: all
all: UserManual.pdf

UserManual.pdf: UserManual.tex 
	pdflatex UserManual.tex
	pdflatex UserManual.tex

tree.pdf: tree.tex
	latex tree.tex
	latex tree.tex
	dvips tree.dvi
	ps2pdf tree.ps

.PHONY: clean
clean:
	rm -f *.log *.aux *.dvi *.ps *.bib

.PHONY: distclean
distclean: clean
	rm -f UserManual.pdf tree.pdf

