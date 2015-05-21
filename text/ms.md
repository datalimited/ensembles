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

\clearpage

## Main messages

Ensembles have great potential, but your training dataset needs to be diverse
and representative of the data you're apply to or they can go very wrong.

Will never do better than out-of-bag measures.

Ensembles not formulaic --- need to keep digging.

# Abstract

# Introduction

We often have multiple models of ecological systems. For example, models might
incorporate different data types, assume alternate population dynamics, or make
contrasting assumptions about the starting state of a system. These models may
suggest different conclusions about population status and often times we do not
know which model is best. How can we reconcile multiple models to make robust
management decisions about ecological resources?

When models give conflicting estimates of a population's status, there are a
number of solutions. (1) We can pick a single model. However, it isn't always
obvious which model is best, and even objectively second-best models may contain
useful additional information. (2) We can maintain multiple models throughout
the decision-making process. However, this can be complicated (REFs)... (3) We
can combine the model outputs. This could be as simple as taking the average,
potentially weighting the models by some performance metric. This could be as
complicated as forming an additional 'hyper' model that draws on interactions
between model outputs and potentially incorporates additional information to
derive a best estimate.

Ensemble models are widely used in the climate sciences, where thousands of
models with varying structural combined across different structural ... The idea behind ensemble
models forms the backbone of many machine learning methods. For example...

In fisheries science, a common problem is estimating the status of an exploited
fish population. For the majority of fish stocks, we have limited information to
go on and stock status isn't known. In recent years, a number of methods have
been proposed to derive population status based on limited information and a set
of assumptions. However, we know that these models can give conflicting output
and no one model is that accurate.

Here, we develop ensemble models for data-limited exploited fish populations. We
explore a variety of ensemble approaches applied to both simulated and
real-world fish stocks and compare performance...

# Methods

We tested the ability of ensemble models to improve estimates of population status. Specifically, we applied ensemble models to a large-scale fully factorial dataset of simulated fisheries. We combined 

Repeated three-fold cross-validation [@hastie2009] to test out-of-sample
prediction error. Split into three, build models on two of the three, test on
the third, repeat for all splits, repeat entire procedure N times.

## Datasets

<!--TODO: github flr repo not available-->
We developed and tested ensemble models with both a fully factorial simulated
dataset [@rosenberg2014] and the RAM Legacy Stock Assessment Database
[@ricard2012]. @rosenberg2014 provide a full description of the simulation
model and the code to generate the simulations is available at
<https://github.com/flr/StockSims>. To summarize, the model included all
combinations of the following factors:

- Three life histories: small pelagic, demersal, large pelagic.
- Two time series lengths: 20 and 60 years.
- Three levels of initial depletion: 1, 0.7, and 0.4.
- Four harvest dynamics: (1) a constant harvest rate, (2) a harvest rate that
  is coupled with biomass to mimic an open-access single-species fishery, (3) a
  one-way trip where harvest rate increases five percent per year to 80% of the
  level at which the stock would crash, and (4) roller-coaster shaped where the
  harvest rate increases, maintains a steady level at 80% of the level at which
  the stock would crash, and then decreases to FMSY. 
- Two levels of multiplicative standard deviation around catch (0 and 0.2)
- Two levels of standard deviation of recruitment variability (0.2 and 0.6)
- Two levels of autoregressive correlation on recruitment variability (0 and 0.6)

\noindent
Ten stochastic draws of recruitment and catch-recording variability were run
for each combination for a total of 5760 stocks.

Our analysis with the RAM Legacy Stock Assessment Database was based on a
version downloaded on XX. This version includes XX stocks from XX countries
across XX taxonomic orders.

## Individual models of population status

COM-SIR
CMSY
SSCOM
mPRM

## Additional covariates

Ensemble methods allow us to incorporate additional covariates into our models
and potentially leverage interactions of how these covariates combine with the
main individual model estimates to estimate population status.
Examples might be life-history characteristics, additional information on
exploitation patterns, or additional statistical properties of the data going
into the individual models.
For simplicity, and to allow us to apply models developed with the simulated
dataset to the real-world dataset, we added only one set of additional
covariates: spectral properties of the catch time series itself. We fit
spectral models to the catch time series and recorded a representative
short-term and long-term spectral density: 0.05 (corresponding to a 20-year
cycle) and 0.2 (corresponding to a 5-year cycle).

## Ensemble models

Mean
Linear model (interactions chosen through cross-validation, or AIC)
Random forest
GBM

## Testing model performance

A critical component to any predictive modelling exercise is to evaluate the
predictive performance of a model on new data [@hastie2009]. Typically, data
are limited in availability, and so a common and effective tool is
cross-validation [@hastie2009]. We used repeated three-fold cross validation to
test predictive performance: we randomly divided the dataset into three sets,
build the model on two-thirds of the data, and evaluate predictive performance
on the remaining third of the data. We repeated this across each of the thirds
of the data and then repeated the whole procedure XX times to account for bias
that may result from any one set of validation splits. 

## Performance metrics

# Results

# Discussion

# Acknowledgements

Funding...

# Citation notes

<!--
- @breiner2015 ensemble models of species distribution models for rare species

- @jones2015 ensemble models of species distribution models - globally for
  marine biodiversity

- [@greene2006] example similar to ours but with climate models (Bayesian
  multilevel ensemble)

- [@kell2007] FLR

- multiple learners book [@alpaydin2010]

- [@caruana2004] nice paper on ensemble models; prob of overfitting increases
  with more models, bagging important

- [@tebaldi2007] key paper: review of 'multi-model ensembles for climate projections'

- famous ensemble is DEMETER: Development of a European Multi-model Ensemble
  System for Seasonal to Interannual Prediction @hagedorn2005 is a good
  reference

- Examples of where ensemble models are shown to be better than any one:
  @thomson2006 (public health - malaria), @cantelaube2005 (agriculture - crop
  yield) (citations taken from @tebaldi2007)

- simple averages are used very commonly in climate science - e.g. IPCC 2001

- @tebaldi2007: weighting obviously makes sense, but how do we define a
  performance metric that is based on past observations that is relevant to the
  future? 

- @tebaldi2007: model independence important

- for climate, weighted averages perform better than simple averages [@min2006],
  but are those same weights applicable to the future (or in our case other
  fisheries)?

- error cancellation is one but not the only reason for superiority of ensemble
  models [@hagedorn2005]

- ensemble model may be only marginally better than the best single model in
  any given case, but we don't usually know which is the best single model
  [@hagedorn2005]

- ensemble models useful for regional climate models too; as an example,
  @pierce2009 use 42 metrics to characterize model performance in regional
  climate model ensembles; found that ensemble models were superior to any one
  model --- especially when considering multiple metrics

- @knutti2009: key review paper on motivation and challenges of combining
  climate projection models; performance on testing / current data may only
  weakly relate to future / other datasets

- @murphy2004: Nature paper on ensembles of climate model simulations; weighted
  average better performance than unweighted average

- @dietterich2000: highly cited book chapter "Ensemble Methods in Machine
  Learning"; ensembles are often much more accurate than the "individual
  classifiers that make them up"; but components must be diverse and accurate

- @dietterich2000: a main justification for ensembles: they are
  representational --- no one model usually can contain all hypothetical
  functional forms, but many separate models can cover more hypotheses

- strong paper showing that ensemble models are most accurate when the various
  individual models make errors in uncorrelated ways [@ali1996]

- without substantial training-testing data, appropriate model weights can be
  very hard to deduce and can cause more harm than good (compared to equal
  weighting) [@weigel2010] (they give the example of seasonal forecasting where
  a ton of hind cast testing can be done, vs. long-term climate) (asymmetrical
  loss function)

- [@weigel2010] optimal weights are always as good or better than equal
  weights, but if you get the weights wrong, you can be better off just using
  equal weights; but averaging was almost always better than any one model

- @rykiel1996 "Testing ecological models: the meaning of validation" (see for
  performance criteria, theory on model testing and assessment)
-->
