APJDIR= ApJ_v1
ARXIVDIR= arXiv_v1

all: ms.pdf

# Automatically copies the figures, but needs to be edited for plotone
apj:
	cp ms.tex $(APJDIR)
	echo "#!/bin/sh" > tmp_cp
	grep \includegraphics ms.tex | awk -F"[{}]" '{print "cp " $$(2) " $(APJDIR)"}' >> tmp_cp
	chmod a+x tmp_cp
	./tmp_cp
	$(RM) tmp_cp

ms.pdf: ms.tex ms.bib
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $<
	pdflatex $<
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "


.PHONY: clean spotless apj

clean:
	$(RM) *.log *.aux *.dvi *.out *.bbl
	$(RM) *~

spotless: clean
	$(RM) ms.pdf
