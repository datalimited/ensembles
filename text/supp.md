# Supporting Material

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


\clearpage

# Supplementary Tables

Table of the predictors

\begin{longtable}{>{\RaggedRight}m{3.2cm}>{\RaggedRight}p{8.5cm}}
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

Maximum catch & Maximum observed catch\\

Spectral density 0.05 & Spectral density at 20 years\\

Spectral density 0.20 & Spectral density at 5 years\\

Lifehistory & A categorical description of lifehistory. 
  One of three values in the simulated dataset.\\

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
predictor on mean \bbmsy\\. TODO Add random forest model to this too if we keep
both around.}
\label{fig:partial-sim}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/partial-sim-2d.pdf}
\caption{Partial dependence two-dimensional plots for simulation GBM ensembles. Red is above one and blue is below one.}
\label{fig:partial-2d-sim}
\end{center}
\end{figure}

<!--\begin{figure}[htbp]-->
<!--\begin{center}-->
<!--\includegraphics[width=\textwidth]{../figs/cv-sim-mean-scatter.png}-->
<!--\caption{Cross-validated simulation ensemble scatter plot --- mean B/Bmsy}-->
<!--\label{fig:scatter-sim-mean}-->
<!--\end{center}-->
<!--\end{figure}-->

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/hex-slope-sim.pdf}
\caption{Same as Fig.~\ref{fig:sim-hexagon} but with the slope of \bbmsy\\ in
the last three years as the response variable.}
\label{fig:scatter-sim-slope}
\end{center}
\end{figure}

<!--
\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/performance-sim-scatter.pdf}
\caption{Performance for simulation ensembles.}
\label{fig:performance-sim}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.6\textwidth]{../figs/roc-sim.pdf}
\caption{ROC for simulation data ensembles.}
\label{fig:roc-sim}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/partial-ram.pdf}
\caption{Partial dependence plots for RAM ensembles.}
\label{fig:partial-ram}
\end{center}
\end{figure}
-->

<!--
\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/cv-ram-mean-scatter.png}
\caption{Cross-validated RAM ensemble scatter plot --- mean B/Bmsy}
\label{fig:scatter-ram-mean}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/cv-ram-slope-scatter.png}
\caption{Cross-validated RAM ensemble scatter plot --- slope B/Bmsy}
\label{fig:scatter-ram-slope}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/performance-ram-scatter.pdf}
\caption{Performance for RAM ensembles.}
\label{fig:performance-ram}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.6\textwidth]{../figs/roc-ram.pdf}
\caption{ROC for RAM data ensembles.}
\label{fig:roc-ram}
\end{center}
\end{figure}
-->
