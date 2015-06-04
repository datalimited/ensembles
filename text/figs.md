# Figures

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=3.7in]{../figs/motivate.pdf}
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
\includegraphics[width=0.7\textwidth]{../figs/didactic2.pdf}
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
\includegraphics[width=\textwidth]{../figs/hex-mean-sim-ram-cv.pdf}
\caption{True vs. predicted \bbmsy\\ for a simulated dataset of known
population status (top two rows) and the RAM Legacy stock assessment database
(third row). Upper panels represent four individual data-limited methods. The
output from these methods is combined, along with additional covariates, to
form the ensemble model estimates in the lower panels. These scatterplots
represent repeated three-fold cross-validation tests where the ensemble models
are built on two-thirds of the data and tested on the third. The data were
binned into hexagons for visual presentation. Darker areas indicate areas with
greater density of data.}
\label{fig:hexagon}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6in]{../figs/performance-gg.pdf}
\caption{Performance of individual and ensemble models fitted to a dataset with
known population statuses (left panel) and the RAM Legacy stock assessment
database (right panel). The x-axis represents within-population innaccuracy:
median absolute relative error. The y-axis represents across-population
Spearman rank-order correlation. The top-left corner contains methods with the
best performance across the two metrics. The colour shading represents bias in
relative error: white points are unbiased, blue points represent methods that
predict \bbmsy\\ values that are too high, red points represent methods that
predict \bbmsy\\ values that are too low. These performance metrics are derived
from the data in Fig.~\ref{fig:sim-hexagon} and based on repeated three-fold
cross-validation testing.
}
\label{fig:performance}
\end{center}
\end{figure}
