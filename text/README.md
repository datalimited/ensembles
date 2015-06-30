The main text is the file `ms.Rmd`. The figure captions are in `figs.md`. The supplement is in the file `supp.md`. These are all Markdown files and can be edited as plain text files in any text editor. If you don't have a favourite editor, RStudio will render these files nicely.

References are in the file `fisheries-ensembles.bib` in BibTeX format.

The PDF of the paper gets generated via `makefile` and the wrapper LaTeX file `anderson-etal-ensembles.tex`.

To make any substantial edits to the paper using these source files:

* Create a Git branch for you edits
* Edit; try and avoid unnecessary changes of line breaks so it's easy to see where you've made changes
* You can add comments right in the text, say in italics (wrap them in `*`) or in capital letters. Or with HTML comments: `<!-- comment here -->`.
* Commit your changes and push them to GitHub.
* Make a pull request.
