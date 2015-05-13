## Ensembles

To recreate the analysis, source the `.R` files in order:

```R
source("0-read.R")
source("1-read-old-ram-fits.R")
# etc. ...
```

Generated data ends up in `generated-data` and figures in `figs`. Intermediate versions of some raw data have been cached in the `generated-data` folder so that the script will run from a freshly cloned version of this repository. (Some of the raw data is too large for Git.)

The draft paper is in the `text` folder. See the Markdown file `ms.md`. To render the PDF version of the paper run `make` on the command line in the `text` folder. In addition to `make`, you'll need [Pandoc](http://pandoc.org), and [LaTeX](http://latex-project.org) installed. You can also just read the rendered version of `ms.md` [on GitHub](https://github.com/datalimited/ensembles/blob/master/text/ms.md). You can view the current version of the PDF [here](https://dl.dropboxusercontent.com/u/254940/anderson-etal-ensembles-v0.1.pdf).

To make edits: create a branch, make your edits, and submit a pull request.
