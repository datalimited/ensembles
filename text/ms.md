# Robust estimates of population status from ensemble models

# Combining data-limited fisheries models to derive robust estimates of stock status

# Improving estimates of ecological population status in the face of disparate models

# Improving ecological decision making in the face of conflicting output with ensemble models

Sean C. Anderson^1^, 
Jamie Afflerbach^n^, 
Andrew B. Cooper^n^, 
Mark Dickey-Collas^n^, 
Olaf P. Jensen^n^, 
Kristin M. Kleisner^n^, 
Catherine Longo^n^, 
C\'{o}il\'{i}n Minto^n^, 
Giacomo Chato Osio^n^, 
Dan Ovando^n^, 
Andrew A. Rosenberg^n^, 
Elizabeth R. Selig^n^, 
James T. Thorson^n^
(order to be determined and others may be added)

\noindent
^1^School of Resource and Environmental Management,
Simon Fraser University, Burnaby, BC, V5A 1S6, Canada

\noindent
^\*^Corresponding author: Sean C. Anderson,
School of Resource and Environmental Management,
Simon Fraser University,
Burnaby BC, V5A 1S6;
E-mail: sean_anderson@sfu.ca

## Figures

Figure 1: Illustration of ensemble methods with cartoon results

Figure 2: Cross-validation scatter plots

Figure 3: Cross-validation performance metrics (correlation across stocks, MARE
within stocks, relative error distributions to illustrate bias and precision)

Figure 4: Performance metrics or applying simulation to RAM and/or RAM to
simulation (scatter panel and performance panels)

## Supplementary Figures

- partial dependence plots
- 2D partial dependence plots
- scatterplots and performance metrics for slopes
- scatterplots and performance metrics for cross-validated RAM stocks
- ROC plots

We often have multiple models of ecological systems. For example, models might incorporate different data types, assume alternate population dynamics, or make contrasting assumptions about the starting state of a system. These models may suggest different conclusions about population status and often times we do not know which model is best. How can we reconcile multiple models to make robust management decisions about ecological resources?
## Main messages

When models give conflicting estimates of a population's status, there are a number of solutions. (1) We can pick a single model. However, it isn't always obvious which model is best, and even objectively second-best models may contain useful additional information. (2) We can maintain multiple models throughout the decision-making process. However, this can be complicated (REFs)... (3) We can combine the model outputs. This could be as simple as taking the average, potentially weighting the models by some performance metric. This could be as complicated as forming an additional 'hyper' model that draws on interactions between model outputs and potentially incorporates additional information to derive a best estimate.
Ensembles have great potential, but your training dataset needs to be diverse
and representative of the data you're apply to or they can go very wrong.

Ensemble models are widely used in the climate sciences, where thousands of models are combined across different structural ...
The idea behind ensemble models forms the backbone of many machine learning methods. For example...
Will never do better than out-of-bag measures.

In fisheries science, a common problem is estimating the status of an exploited fish population. For the majority of fish stocks, we have limited information to go on and stock status isn't known. In recent years, a number of methods have been proposed to derive population status based on limited information and a set of assumptions. However, we know that these models can give conflicting output and no one model is that accurate.
Ensembles not formulaic --- need to keep digging.

Here, we develop ensemble models for data-limited exploited fish populations. We explore a variety of ensemble approaches applied to both simulated and real-world fish stocks and compare performance...

# Methods

# Results

# Discussion

# Acknowledgements

Funding...

\bibliographystyle{apalike}
\bibliography{fisheries-ensembles}

# Figures

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

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/cv-sim-mean-scatter.png}
\caption{Cross-validated simulation ensemble scatter plot --- mean B/Bmsy}
\label{fig:scatter-sim-mean}
\end{center}
\end{figure}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=\textwidth]{../figs/cv-sim-slope-scatter.png}
\caption{Cross-validated simulation ensemble scatter plot --- slope B/Bmsy}
\label{fig:scatter-sim-slope}
\end{center}
\end{figure}

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


