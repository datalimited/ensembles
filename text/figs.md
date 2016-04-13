# Figures

<!--The legends for all figures should be grouped on a page that precedes the
figures. Do not place a figure and its legend on the same page.-->

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/didactic5.pdf}
\caption{Using a superensemble model to predict population status.
The process is illustrated didactically on the left and with R pseudocode on the right. 
(a) Individual models (red and blue lines) are fit to data (dots) 
from populations of known or assumed status (known status shown by black line).
The shaded gray boxes indicate the recent time period that we are ultimately interested in in this paper. 
Estimates from these individual models, potentially combined with additional
covariates, are then used as covariates in a statistical model fitted to the
known or assumed population status as the response (here represented as a simple linear model). 
The symbols $\beta$ and $X$ represent estimated parameters in the linear model and the status estimates from models 1 and 2, respectively.
The symbol $\epsilon$ represents error in this linear model. 
The $i$ subscripts represent individual fish stocks from $1$ to $n$, and $y$ represents the known status.
(b) The superensemble can then be used to predict on new stocks of interest.
The same individual models are fit to populations of interest and then combined using the previously fitted superensemble mode. 
}
\label{didactic}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=3.7in]{../figs/motivate.pdf}
\caption{Different models can suggest conflicting population statuses and
trends. Shown here are trajectories of estimated \bbmsy\\ from four
data-limited assessment methods (colours) and a data-rich stock assessment
(black). Lines indicate median fits and shaded regions interquartile ranges.
Dashed horizontal line indicates $B/B_\mathrm{MSY} = 1$.}
\label{motivate}
\end{center}
\end{figure}

\clearpage


\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/hex-mean-sim-ram-cv.pdf}
\caption{True (or assessed) population status (x axis) vs. predicted population status from individual models and ensemble methods with cross-validation (y axis).
These scatterplots represent the aggregate results of repeated three-fold cross-validation tests where the ensemble models are built on two-thirds of the data and tested on the remaining third.
(a--d) Individual data-limited model estimates of mean $B/B_\mathrm{MSY}$ (biomass divided by biomass at maximum sustainable yield) in the last five years for a simulated dataset of known population status.
(e--h) Ensemble estimates for the same populations. Shown are a mean, a linear model with two-way interactions (LM), a random forest ensemble (RF), and a generalised boosted regression model (GBM).
(i--l) The same ensemble models, which were trained on the simulated dataset, applied to the RAM Legacy stock assessment database and compared to data-rich stock-assessed status.
In the case of the RAM Legacy stock assessment data, the modified panel regression model (mPRM) was refit on each cross-validation split.
The data were binned into hexagons for visual presentation. Darker areas indicate areas with greater density of data. Yellow-red shading and yellow-blue shading distinguishes individual models from ensemble methods.}
\label{hexagon}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=6in]{../figs/performance-inkscape.pdf}
\caption{
Performance metrics of individual and ensemble models predicting
$B/B_\mathrm{MSY}$ (mean biomass divided by biomass at maximum sustainable
yield) in the last five years fitted to a dataset with (a) known population
statuses and (b) the RAM Legacy stock assessment database. The x-axis
represents within-population inaccuracy: median absolute proportional error
(MAPE). The y-axis represents across-population Spearman rank-order
correlation. The top-left corner contains methods with the best performance
across the two metrics. The colour shading represents bias (median
proportional error; MPE): white points are unbiased, blue points represent
methods that predict \bbmsy\\ values that are too high, red points represent
methods that predict \bbmsy\\ values that are too low. These performance
metrics are derived from the data in Fig.~\ref{hexagon} and based on repeated
three-fold cross-validation testing.}
\label{performance}
\end{center}
\end{figure}

<!-- vim: set formatoptions=nroq tw=78: -->
