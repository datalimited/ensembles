# Supporting Material

<!--
# Simulation model

Brief description of the simulation model based on the FAO report

[@rosenberg2014]

# Assessment methods

## CMSY --- catch maximum sustainable yield

## COMSIR --- catch only model with sample-importance-resampling

\begin{equation}
\hat{C}_{t+1} = P_{t+1} \left(B_t + r B_t \left(1 - \frac{B_t}{K}\right) - \hat{C}_t \right)
\end{equation}

\begin{equation}
P_{t+1} = P \left(1 + x\left(\frac{B_t}{aK}-1\right)\right)
\end{equation}

\begin{equation}
P_0 = \frac{C_0}{B_0}
\end{equation}

\begin{equation}
B_\mathrm{MSY} = \frac{K}{2}
\end{equation}

Priors: 

$a \sim \mathrm{uniform}(0, 1)$, $x \sim \mathrm{uniform}(0.000001, 1)$, $\log K \sim \mathrm{uniform}(\mathrm{max catch}, \log (100 \cdot \mathrm{max catch}))$

## mPRM --- modified panel regression model

## SSCOM --- state space catch-only-model

Harvest dynamics model:

\begin{equation}
\hat{E}_{t+1} = E_t \left( \frac{B_t}{a \cdot B_0 / 2} \right) ^ 2
\end{equation}

\begin{equation}
\hat{C}_t = E_t \cdot B_t
\end{equation}

Table summarizing these based on FAO report

[@rosenberg2014]

# Ensemble model descriptions

Mean ensemble:

Linear model ensemble:

Random forest ensemble:

Boosted regression ensemble:
-->

\clearpage

# Supplementary Tables

Table of the predictors

\begin{longtable}{>{\RaggedRight}m{4.5cm}>{\RaggedRight}p{8.5cm}}
\toprule
Variable & Description \\ 
\midrule
CMSY & Median estimated $B/B_{\mathrm{MSY}}$ from 
  CMSY method \citep{martell2013} \\ 

COM-SIR & Median estimated $B/B_{\mathrm{MSY}}$ from 
  CMSY method \citep{vasconcellos2005} \\ 

SSCOM & Median estimated $B/B_{\mathrm{MSY}}$ 
  from SSCOM method (REF) \\ 

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
\includegraphics[width=\textwidth]{../figs/partial-sim.pdf}
\caption{Partial dependence plots for GBM ensemble models fitted to the
simulation data. Lines represent the marginal non-linear effect of each
predictor on \bbmsy\\ after integrating out the other predictor values.}
\label{fig:partial-sim}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/gbm-partial-residuals.pdf}
\caption{Partial residuals plot from the GBM model. The dots represent the residuals when predicting from the model with the predictor set to its mean value. The lines represent the prediction when all other predictors are set to their mean value. Note that these are on a scale of $\log$ \bbmsy\\ residuals.}
\label{fig:partial-residuals-sim-gbm}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/partial-sim-2d.pdf}
\caption{Two-dimensional partial dependence plots for GBM ensemble models fitted to the simulated dataset of known status. Red shading indicates an expected \bbmsy\\ above one and blue shading an expected value below one. White shading represents an expected value of 1.}
\label{fig:partial-2d-sim}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/roc-sim.pdf}
\caption{Receiver-operating-characteristic (ROC) curves from repeated three-fold cross-validation of the simulated data of known status. Shown are ROC curves for (a) ensemble models and (b) individual data-limited model estimates based on estimates of \bbmsy\\ with the response variable representing whether true \bbmsy\\ was above or below one. The diagonal dashed line represents performance that is no better than flipping a coin. The area under the curve represents the probability that the model would correctly rank two randomly chosen stocks in terms of their mean \bbmsy\\ in the last five years. Sensitivity represents the true positive rate (correctly categorizing a stock as having a \bbmsy\\ greater than one). Specificity refers to the true negative rate (correctly categorizing a stock as having a \bbmsy\\ less than one). The sensitivity and sensitivity values are shown across all possible splits (values of \bbmsy\\) at which one could divide the stocks into these two categories.}
\label{fig:roc-sim}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/hex-mean-ram-cv.pdf}
\caption{RAM stocks fit with individual data-limited assessment methods (a--d) and ensemble models (e--h) that were trained on the simulated dataset. These are based on 3-fold cross-validation of the mPRM model. Panels e--h duplicate the lower row in Fig.~\ref{fig:hexagon}.}
\label{fig:hexagon-ram}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/hex-mean-sim-basic-cv.pdf}
\caption{Same as Fig.~\ref{fig:hexagon} but with simulation ensembles that have no additional covariates (spectral densities were not included in these models).}
\label{fig:hexagon-sim-basic}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/hex-slope-sim.pdf}
\caption{Same as Fig.~\ref{fig:hexagon} but with the slope (Theil-Sen median slope) of \bbmsy\\ in the last five years as the response variable.}
\label{fig:scatter-sim-slope}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=5in]{../figs/performance-slope-sim.pdf}
\caption{Same as Fig.~\ref{fig:performance} but with the slope of \bbmsy\\ in
the last five years as the response variable. This is based on the data shown in Fig.~\ref{fig:scatter-sim-slope}.}
\label{fig:performance-sim-slope}
\end{center}
\end{figure}

