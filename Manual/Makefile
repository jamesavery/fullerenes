
.PHONY: all
all: UserManual.pdf

UserManual.pdf: UserManual.tex Makefile tree.pdf
	pdflatex UserManual.tex
	pdflatex UserManual.tex

tree.pdf: tree.tex
	latex tree.tex
	latex tree.tex
	dvips tree.dvi
	ps2pdf tree.ps tree_tmp.pdf
	pdfcrop tree_tmp.pdf tree.pdf

.PHONY: clean
clean:
	rm -f *.log *.aux *.dvi *.ps *.bib tree_tmp.pdf

.PHONY: distclean
distclean: clean
	rm -f UserManual.pdf tree.pdf 

