# Improving estimates of population status and trend with superensemble models

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

Generated data ends up in `generated-data` and figures in `figs`. Intermediate versions of some raw data have been cached in the `generated-data` folder so that the scripts will run from a freshly cloned version of this repository. (Some of the raw data files are too large for Git.)

The paper is in the `text` folder. To render the PDF version of the paper run `make` on the command line in the `text` folder. In addition to `make`, you'll need [Pandoc](http://pandoc.org), [LaTeX](http://latex-project.org), and [latexmk](https://www.ctan.org/pkg/latexmk/?lang=en) installed.
