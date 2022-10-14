.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:miniprot.pdf

miniprot.pdf:miniprot.tex miniprot.bib
		pdflatex miniprot; bibtex miniprot; pdflatex miniprot; pdflatex miniprot;

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ miniprot.pdf
