
.PHONY: all
all: UserManual.pdf

UserManual.pdf: UserManual.tex Makefile workflow.pdf tree.pdf gc.pdf
	pdflatex UserManual.tex
	pdflatex UserManual.tex

workflow.pdf: workflow.tex
	latex workflow.tex
	latex workflow.tex
	dvips workflow.dvi
	ps2pdf workflow.ps workflow_tmp.pdf
	pdfcrop workflow_tmp.pdf workflow.pdf

tree.pdf: tree.tex
	latex tree.tex
	latex tree.tex
	dvips tree.dvi
	ps2pdf tree.ps tree_tmp.pdf
	pdfcrop tree_tmp.pdf tree.pdf

gc.pdf: gc.tex
	latex gc.tex
	latex gc.tex
	dvips gc.dvi
	ps2pdf gc.ps gc_tmp.pdf
	pdfcrop gc_tmp.pdf gc.pdf

.PHONY: clean
clean:
	rm -f *.log *.aux *.dvi *.ps *.bib *.out *_tmp.pdf

.PHONY: distclean
distclean: clean
	rm -f UserManual.pdf tree.pdf gc.pdf

