## Ensembles

To recreate the analysis, source the `.R` in the `analysis` files in order:

```R
source("0-read.R")
source("1-read-old-ram-fits.R")
# etc. ...
```

Or run `make` on the command line from the `analysis` folder:

```sh
cd analysis
make
```

Generated data ends up in `generated-data` and figures in `figs`. Intermediate versions of some raw data have been cached in the `generated-data` folder so that the scripts will run from a freshly cloned version of this repository. (Some of the raw data is too large for Git.)

The draft paper is in the `text` folder. See the Markdown file `ms.md`. To render the PDF version of the paper run `make` on the command line in the `text` folder. In addition to `make`, you'll need [Pandoc](http://pandoc.org), [LaTeX](http://latex-project.org), and [latexmk](https://www.ctan.org/pkg/latexmk/?lang=en) installed. You can view the rendered version of the PDF [here](https://dl.dropboxusercontent.com/u/254940/anderson-etal-ensembles.pdf).

To make edits: create a branch, name it something like your initials, make your edits, push them to GitHub, and click the big green button to submit a pull request. Or markup a copy of the PDF. Or edit the [RTF version](https://dl.dropboxusercontent.com/u/254940/anderson-etal-ensembles.rtf) of the paper with track changes, but ignore any formatting in that document — it's auto-generated.
