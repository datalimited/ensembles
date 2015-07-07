TEXT = anderson-etal-ensembles
MS = ms
SOM = supp
FIGS = figs
REFS = fisheries-ensembles

all: $(TEXT).pdf

%.md: %.Rmd
	Rscript -e "knitr::knit('$<')"

%.tex: %.md
	pandoc -S --no-wrap --bibliography=$(REFS).bib --natbib -fmarkdown-all_symbols_escapable $< -o $@
	#perl -p -i -e "s/\\\citep/\\\cite/g" $@
	#perl -p -i -e "s/\\\citet/\\\cite/g" $@
	perl -p -i -e "s/Figure /Figure~/g" $@
	perl -p -i -e "s/Fig\. /Fig.~/g" $@
	perl -p -i -e "s/vs\. /vs.~/g" $@
	perl -p -i -e "s/Table /Table~/g" $@
	perl -p -i -e "s/Figs\. /Figs.~/g" $@

$(TEXT).pdf: $(MS).tex $(SOM).tex $(FIGS).tex $(TEXT).tex $(REFS).bib
	latexmk $(TEXT)
	#pdflatex $(TEXT)
	#bibtex $(TEXT)
	#pdflatex $(TEXT)
	#pdflatex $(TEXT)

dropbox: $(TEXT).pdf rtf
	cp $(TEXT).pdf ~/Dropbox/Public/$(TEXT).pdf
	mv $(TEXT).rtf ~/Dropbox/Public/$(TEXT).rtf

rtf: $(TEXT).pdf
	latex2rtf -f0 -M6 -E0 $(TEXT).tex

$(TEXT).docx: ms.Rmd
	Rscript -e "knitr::knit('ms.Rmd')"
	pandoc -S -s --bibliography=$(REFS).bib --csl=prsb.csl ms.md -o $@

extractbib: text
	bibtool -x $(MS).aux -o $(MS).bib -- 'expand.macros = ON'
	bibtool -x $(SOM).aux -o $(SOM).bib -- 'expand.macros = ON'

clean:
	rm ms.md
	latexmk -c

continuous:
	latexmk -pvc $(TEXT)
