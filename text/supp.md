# Supporting Material

<!--
Supporting Information (i.e., online appendices) should be cited in the text of
the paper. Every piece is cited as Supporting Information, not by specific
appendix number. Before Literature Cited, insert a paragraph in the exact
format shown below that provides a brief description of supporting information
elements.

Supporting Information

XXX (Appendix S1), XXX (Appendix S2), and a XXX translation of the article
(Appendix S3) are available online. The authors are solely responsible for the
content and functionality of these materials. Queries (other than absence of
the material) should be directed to the corresponding author.
-->

# Supplementary Tables

\renewcommand{\thetable}{S\arabic{table}}
\setcounter{table}{0}

Table S1: Covariates in the ensemble models. The first four variables are predictions of the mean or slope of $B/B_\mathrm{MSY}$ in the last five years from the models described in @rosenberg2014. The last two variables are additional covariates incorporated into the ensemble models. These two variables are derived from spectral analysis and represent spectral densities at long- and short-term frequencies.

\begin{longtable}{>{\RaggedRight}m{4.5cm}>{\RaggedRight}p{9.4cm}}
\toprule
Variable & Description \\
\midrule
CMSY & Median estimated $B/B_{\mathrm{MSY}}$ from
  CMSY method \citep{martell2013} \\

COM-SIR & Median estimated $B/B_{\mathrm{MSY}}$ from
  COM-SIR method \citep{vasconcellos2005} \\

SSCOM & Median estimated $B/B_{\mathrm{MSY}}$
  from SSCOM method \citep{thorson2013} \\

mPRM & Median estimated $B/B_{\mathrm{MSY}}$
  from modified panel regression method \citep{costello2012} \\

Spectral density 0.05 & Spectral density (of fraction of maximum catch) at 20 years\\

Spectral density 0.20 & Spectral density (of fraction of maximum catch) at 5 years\\

\bottomrule
\label{tab:predictors}
\end{longtable}

\clearpage

# Supplementary Figures

\renewcommand{\thefigure}{S\arabic{figure}}
\setcounter{figure}{0}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/hex-slope-sim.pdf}
\caption{Same as Fig.~\ref{hexagon} but with the slope (Theil-Sen median slope
\citep{theil1950}) of $B/B_\mathrm{MSY}$ in the last five years as the
response variable.}
\label{scatter-sim-slope}
\end{center}
\end{figure}

\clearpage

<!-- \begin{figure}[htbp] -->
<!-- \begin{center} -->
<!-- \includegraphics[width=\textwidth]{../figs/roc-sim.pdf} -->
<!-- \caption{Receiver-operating-characteristic (ROC) curves from repeated
three-fold cross-validation of the simulated data of known status. Shown are
ROC curves for (a) ensemble methods and (b) individual data-limited model
estimates based on estimates of $B/B_\mathrm{MSY}$ with the response variable
representing whether true $B/B_\mathrm{MSY}$ was above or below 0.5, the
threshold for declaring a stock overfished in the United States and Australia.
The diagonal dashed line represents performance that is no better than
flipping a coin. The area under the curve represents the probability that the
model would correctly rank two randomly chosen stocks in terms of their mean
$B/B_\mathrm{MSY}$ in the last five years. Sensitivity (y axis) represents the
true positive rate (correctly categorizing a stock as having a
$B/B_\mathrm{MSY}$ greater than 0.5). Specificity (x axis) refers to the true
negative rate (correctly categorizing a stock as having a $B/B_\mathrm{MSY}$
less than one). The sensitivity and specificity values are shown across all
possible decision thresholds (values of $B/B_\mathrm{MSY}$) at which one could
divide the stocks into these two categories.} -->
<!-- \label{roc-sim} -->
<!-- \end{center} -->
<!-- \end{figure} -->

<!-- \clearpage -->

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=5in]{../figs/performance-slope-sim-inkscape.pdf}
\caption{Same as Fig.~\ref{performance} but with the slope of $B/B_\mathrm{MSY}$ in
the last five years as the response variable. This is based on the data shown
in Fig.~\ref{scatter-sim-slope}.}
\label{performance-sim-slope}
\end{center}
\end{figure}

\clearpage


\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/partial-sim.pdf}
\caption{Partial dependence plots for GBM ensemble models fitted to the
simulation data. Lines represent the marginal non-linear effect of each
predictor on mean $B/B_\mathrm{MSY}$ after integrating out the other predictor
values. The marginal effect at a given predictor value is the average
predicted response across all data points holding the given predictor at a
specific value.}
\label{partial-sim}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.6\textwidth]{../figs/lm-coefs.pdf}
\caption{Standardized regression coefficients from a linear model ensemble
predicting log mean $B/B_\mathrm{MSY}$ in the last five years. Coefficients
are centered (mean of each predictor is subtracted) and scaled (divided by two
standard deviations). Thick and thin lines represent +/- one and two standard
errors from the mean.}
\label{lm-coefs}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/partial-sim-2d.pdf}
\caption{Two-dimensional partial dependence plots for GBM ensemble models
fitted to the simulated dataset of known status. Red shading indicates an
expected $B/B_\mathrm{MSY}$ above one and blue shading an expected value below
one. White shading represents an expected value of 1. For example, in panel e
$B/B_\mathrm{MSY}$ is estimated to be low if SSCOM estimates
$B/B_\mathrm{MSY}$ to be low regardless of the estimate from CMSY. On the
other hand, if CMSY estimates $B/B_\mathrm{MSY}$ to be between 1.0 and about
1.5, the ensemble estimates $B/B_\mathrm{MSY}$ to be high unless SSCOM
estimates $B/B_\mathrm{MSY}$ to be very low.}
\label{partial-2d-sim}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/hex-mean-sim-basic-cv.pdf}
\caption{Same as Fig.~\ref{hexagon} but with simulation ensembles that have no
additional covariates (spectral densities were not included in these models).}
\label{hexagon-sim-basic}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/hex-mean-ram-cv.pdf}
\caption{RAM stocks fit with individual data-limited assessment methods (a--d)
and ensemble models (e--h) that were trained on the simulated dataset. These
are based on 3-fold cross-validation of the mPRM model. Panels e--h duplicate
the lower row in Fig.~\ref{hexagon}.}
\label{hexagon-ram}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/partial-sim-slope.pdf}
\caption{Partial dependence plots for GBM ensemble models fitted to the
simulation data. Lines represent the marginal non-linear effect of each
predictor on the slope of $B/B_\mathrm{MSY}$ after integrating out the other predictor values.}
\label{partial-sim-slope}
\end{center}
\end{figure}

<!-- vim: set formatoptions=nroq tw=78: -->
