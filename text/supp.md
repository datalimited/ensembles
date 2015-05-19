# Supporting Material

# Simulation model

Brief description of the simulation model based on the FAO report

[@rosenberg2014]

# Assessment methods

CMSY is...

COMSIR is...

mPRM is...

SSCOM is...

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
\caption{Partial dependence plots for simulation ensembles.}
\label{fig:partial-sim}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/partial-sim-2d.pdf}
\caption{Partial dependence 2d plots for simulation ensembles.}
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
\caption{Cross-validated simulation ensemble scatter plot --- slope B/Bmsy}
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
