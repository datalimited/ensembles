TEXT = anderson-etal-ensembles
MS = ms
SOM = supp
VERSION = 0.1
REFS = fisheries-ensembles

all: $(TEXT).pdf

$(MS).tex: $(MS).md
	pandoc -S --no-wrap --bibliography=$(REFS).bib --natbib -fmarkdown-all_symbols_escapable $< -o $@
	perl -p -i -e "s/Figure /Figure~/g" $(MS).tex
	perl -p -i -e "s/Fig\. /Fig.~/g" $(MS).tex
	perl -p -i -e "s/vs\. /vs.~/g" $(MS).tex
	perl -p -i -e "s/Table /Table~/g" $(MS).tex
	perl -p -i -e "s/Figs\. /Figs.~/g" $(MS).tex

$(TEXT).pdf: $(MS).tex $(TEXT).tex $(REFS).bib
	latexmk $(TEXT)

quick:
	pdflatex $(TEXT)

dropbox: text rtf
	cp $(TEXT).pdf ~/Dropbox/Public/$(TEXT)-v$(VERSION).pdf
	mv $(TEXT).rtf ~/Dropbox/Public/$(TEXT)-v$(VERSION).rtf

rtf: text
	latex2rtf -E0 $(TEXT).tex

extractbib: text
	bibtool -x $(MS).aux -o $(MS).bib -- 'expand.macros = ON'
	bibtool -x $(SOM).aux -o $(SOM).bib -- 'expand.macros = ON'

clean:
	latexmk -c

continuous:
	latexmk -pvc $(TEXT)
