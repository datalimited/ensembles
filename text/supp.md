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

# Supplementary Tables

\renewcommand{\thetable}{S\arabic{table}}
\setcounter{table}{0}

Table S1: Covariates in the ensemble models. The first four variables are predictions of the mean or slope of $B/B_\mathrm{MSY}$ in the last five years from the models described in Rosenberg et al. [@rosenberg2014]. The last two variables are additional covariates incorporated into the ensemble models. These two variables are derived from spectral analysis and represent spectral densities at long- and short-term frequencies.

\begin{longtable}{>{\RaggedRight}m{4.5cm}>{\RaggedRight}p{9.4cm}}
\toprule
Variable & Description \\ 
\midrule
CMSY & Median estimated $B/B_{\mathrm{MSY}}$ from 
  CMSY method \citep{martell2013} \\ 

COM-SIR & Median estimated $B/B_{\mathrm{MSY}}$ from 
  CMSY method \citep{vasconcellos2005} \\ 

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
\caption{Same as Fig.~\ref{hexagon} but with the slope (Theil-Sen median slope \citep{theil1950}) of $B/B_\mathrm{MSY}$ in the last five years as the response variable.}
\label{scatter-sim-slope}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/roc-sim.pdf}
\caption{Receiver-operating-characteristic (ROC) curves from repeated three-fold cross-validation of the simulated data of known status. Shown are ROC curves for (a) ensemble models and (b) individual data-limited model estimates based on estimates of $B/B_\mathrm{MSY}$ with the response variable representing whether true $B/B_\mathrm{MSY}$ was above or below one. The diagonal dashed line represents performance that is no better than flipping a coin. The area under the curve represents the probability that the model would correctly rank two randomly chosen stocks in terms of their mean $B/B_\mathrm{MSY}$ in the last five years. Sensitivity (y axis) represents the true positive rate (correctly categorizing a stock as having a $B/B_\mathrm{MSY}$ greater than one). Specificity (x axis) refers to the true negative rate (correctly categorizing a stock as having a $B/B_\mathrm{MSY}$ less than one). The sensitivity and sensitivity values are shown across all possible splits (values of $B/B_\mathrm{MSY}$) at which one could divide the stocks into these two categories.}
\label{roc-sim}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=5in]{../figs/performance-slope-sim.pdf}
\caption{Same as Fig.~\ref{performance} but with the slope of $B/B_\mathrm{MSY}$ in
the last five years as the response variable. This is based on the data shown in Fig.~\ref{scatter-sim-slope}.}
\label{performance-sim-slope}
\end{center}
\end{figure}

\clearpage


\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/partial-sim.pdf}
\caption{Partial dependence plots for GBM ensemble models fitted to the
simulation data. Lines represent the marginal non-linear effect of each
predictor on $B/B_\mathrm{MSY}$ after integrating out the other predictor values. The marginal effect at a given predictor value is the average predicted response across all data points holding the given predictor at a specific value.}
\label{partial-sim}
\end{center}
\end{figure}

<!--
\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/gbm-partial-residuals.pdf}
\caption{Partial residuals plot from the GBM model. The dots represent the residuals when predicting from the model with the predictor set to its mean value. The lines represent the prediction when all other predictors are set to their mean value. Note that these are on a scale of $\log$ $B/B_\mathrm{MSY}$ residuals.}
\label{partial-residuals-sim-gbm}
\end{center}
\end{figure}
-->
\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.6\textwidth]{../figs/lm-coefs.pdf}
\caption{Standardized regression coefficients from a linear model ensemble predicting log mean $B/B_\mathrm{MSY}$ in the last five years. Coefficients are centered (mean of each predictor is subtracted) and scaled (divided by two standard deviations). Thick and thin lines represent +/- one and two standard errors from the mean.}
\label{lm-coefs}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/partial-sim-2d.pdf}
\caption{Two-dimensional partial dependence plots for GBM ensemble models fitted to the simulated dataset of known status. Red shading indicates an expected $B/B_\mathrm{MSY}$ above one and blue shading an expected value below one. White shading represents an expected value of 1.}
\label{partial-2d-sim}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/hex-mean-sim-basic-cv.pdf}
\caption{Same as Fig.~\ref{hexagon} but with simulation ensembles that have no additional covariates (spectral densities were not included in these models).}
\label{hexagon-sim-basic}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6.5in]{../figs/hex-mean-ram-cv.pdf}
\caption{RAM stocks fit with individual data-limited assessment methods (a--d) and ensemble models (e--h) that were trained on the simulated dataset. These are based on 3-fold cross-validation of the mPRM model. Panels e--h duplicate the lower row in Fig.~\ref{hexagon}.}
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


