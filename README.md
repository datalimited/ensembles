## Ensembles

To recreate the analysis, source the `.R` files in order:

```R
source("0-read.R")
source("1-read-old-ram-fits.R")
# etc. ...
```

Generated data ends up in `generated-data/` and figures in `figures`.

The draft paper is in the `text` folder. See the Markdown file `ms.md`. To render the PDF version of the paper run `make` on the command line in the `text` folder. In addition to `make`, you'll need [Pandoc](http://pandoc.org), and [LaTeX](http://latex-project.org) installed. You can also just read the rendered version of `ms.md` [on GitHub](https://github.com/datalimited/ensembles/blob/master/text/ms.md).

To make edits: create a branch, make your edits, and submit a pull request.
