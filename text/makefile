TEXT = anderson-etal-ensembles
MS = ms
SOM = supp
FIGS = figs
VERSION = 0.1
REFS = fisheries-ensembles

all: $(TEXT).pdf

%.tex: %.md
	pandoc -S --no-wrap --bibliography=$(REFS).bib --natbib -fmarkdown-all_symbols_escapable $< -o $@
	perl -p -i -e "s/Figure /Figure~/g" $@
	perl -p -i -e "s/Fig\. /Fig.~/g" $@
	perl -p -i -e "s/vs\. /vs.~/g" $@
	perl -p -i -e "s/Table /Table~/g" $@
	perl -p -i -e "s/Figs\. /Figs.~/g" $@

$(TEXT).pdf: $(MS).tex $(SOM).tex $(FIGS).tex $(TEXT).tex $(REFS).bib
	latexmk $(TEXT)

dropbox: $(TEXT).pdf rtf
	cp $(TEXT).pdf ~/Dropbox/Public/$(TEXT)-v$(VERSION).pdf
	mv $(TEXT).rtf ~/Dropbox/Public/$(TEXT)-v$(VERSION).rtf

rtf: $(TEXT).pdf
	latex2rtf -E0 $(TEXT).tex

extractbib: text
	bibtool -x $(MS).aux -o $(MS).bib -- 'expand.macros = ON'
	bibtool -x $(SOM).aux -o $(SOM).bib -- 'expand.macros = ON'

clean:
	latexmk -c

continuous:
	latexmk -pvc $(TEXT)
