# Figures

<!--The legends for all figures should be grouped on a page that precedes the
figures. Do not place a figure and its legend on the same page.-->

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=3.7in]{../figs/motivate.pdf}
\caption{Different models can suggest conflicting population statuses and trends. Shown here are trajectories of estimated \bbmsy\\ from four data-limited assessment methods (colours) and a data-rich stock assessment (black). Lines indicate median fits and shaded regions interquartile ranges. (\textit{Could add a panel with a simulated stock, change the stock that is shown, and/or show a number of small panels with different conflicting patterns but without the confidence intervals.})}
\label{motivate}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=4.2in]{../figs/didactic4.pdf}
\caption{Using a superensemble model to predict population status.
(a) Individual models are fit to populations of known or assumed status. (b) Estimates from these individual models, potentially combined with additional covariates, are then used as coviarates in an yet another statistical model fitted to the known or assumed population status as the response. (c) Finally, the same individual models are fit to a population of interest and (d) combined using the previously fitted ensemble model to derive the superensemble prediction of population status. \textit{There has been discussion about adding a figure that outlines the cross-validation procedure and/or a figure/table that describes the various superensemble methods. There isn't a lot of figure room left, but I may be able to squeeze another panel here.}}
\label{didactic}
\end{center}
\end{figure}

\clearpage

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/hex-mean-sim-ram-cv.pdf}
\caption{True (or assessed) population status (x axis) vs. predicted population status from individual models and ensemble methods with cross-validation (y axis).
These scatterplots represent the aggregate results of repeated three-fold cross-validation tests where the ensemble models are built on two-thirds of the data and tested on the third.
(a--d) Individual data-limited model estimates of $B/B_\mathrm{MSY}$ (biomass divided by biomass at maximum sustainable yield) in the last five years for a simulated dataset of known population status.
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
\includegraphics[width=6in]{../figs/performance.pdf}
\caption{
Performance metrics of individual and ensemble models predicting $B/B_\mathrm{MSY}$ (mean biomass divided by biomass at maximum sustainable yield) in the last five years fitted to a dataset with (a) known population statuses and (b) the RAM Legacy stock assessment database. 
The x-axis represents within-population innaccuracy: median absolute proportional error (MAPE). 
The y-axis represents across-population Spearman rank-order correlation. 
The top-left corner contains methods with the best performance across the two metrics. 
The colour shading represents bias (median proportional error; MPE): white points are unbiased, blue points represent methods that predict \bbmsy\\ values that are too high, red points represent methods that predict \bbmsy\\ values that are too low. 
These performance metrics are derived from the data in Fig.~\ref{hexagon} and based on repeated three-fold cross-validation testing. \textit{(TODO: I'll tweak the labels to avoid overlap near the end. Otherwise, they move slightly with every little thing we tweak in the analysis.)}.
}
\label{performance}
\end{center}
\end{figure}
