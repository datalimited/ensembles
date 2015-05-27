# Figures

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=4in]{../figs/motivate.pdf}
\caption{Assessment methods can suggest conflicting population statuses. Example
trajectories of \bbmsy\\ estimated by four data-limited assessment methods
(colours) and a TODO data-rich stock assessment (black). Lines indicate median
fits and shaded regions interquartile ranges.}
\label{fig:motivate}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=0.9\textwidth]{../figs/didactic.pdf}
\caption{Using an ensemble model to predict population status.
Individual models are fit to populations of known or assumed status. Estimates
from these individual models, along with potential additional covariates, are
then used as coviarates in an ensemble model fitted the known or assumed
population status. Finally, the same individual models are fit to a population
of interest and combined with the fitted ensemble model to derive the ensemble
prediction of population status.}
\label{fig:didactic}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/fig2.pdf}
\caption{True vs. predicted \bbmsy\\ for a simulated dataset of known
population status. Upper panels represent four individual data-limited methods.
The output from these methods is combined, along with additional covariates, to
form the ensemble model estimates in the lower panels. These scatterplots
represent repeated three-fold cross-validation tests where the ensemble models
are built on two-thirds of the data and tested on the third. The data were
binned into hexagons for visual presentation. Darker areas indicate areas with
greater density of data.}
\label{fig:sim-hexagon}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=3.7in]{../figs/fig3.pdf}
\caption{Performance of individual and ensemble models fitted to a 
dataset with known population statuses. 
The x-axis represents within-population innaccuracy: median absolute relative
error. 
The y-axis represents across-population Spearman rank-order correlation. The
top-left corner contains methods with the best performance across the two
metrics. 
The colour shading represents bias in relative error: white points are
unbiased, blue points represent methods that predict \bbmsy\\ values that are
too high, red points represent methods that predict \bbmsy\\ values that are
too low. 
These performance metrics are derived from the data in
Fig.~\ref{fig:sim-hexagon} and based on repeated three-fold cross-validation
testing.}
\label{fig:performance-sim}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=5in]{../figs/hex-mean-ram-cv.pdf}
\includegraphics[width=3in]{../figs/ram-ensemble-performance-cv.pdf}
\caption{RAM stocks fit with ensemble models built from the simulated dataset. Early draft figures to be revised. Need to look into why linear model fails so badly. These are based on 3-fold cross-validation of the mPRM model.}
\label{fig:performance-ram}
\end{center}
\end{figure}

<!--
\begin{figure}[htbp]
\begin{center}
\includegraphics[width=3.7in]{../figs/performance-beanplots-sim.pdf}
\caption{Performance with simulation ensembles (alternate of previous figure).}
\label{fig:performance-sim}
\end{center}
\end{figure}
-->
