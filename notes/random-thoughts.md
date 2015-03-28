---
title:  Combining Data-Limited Fisheries Models to Derive Robust Estimates of Stock Status
author: Sean Anderson
output:
  html_document:
    toc: true
    theme: united
date: 2015-03-19
bibliography: combining-models.bib
csl: cjfas.csl
---

Coilin: 2015-03-27:

- spectral: qualitative effort dynamics
- chris francis - averaging models in assessments
- no effort dynamics in CMSY - low catches always mean mean reduced effort

- look at FLR - do they have these methods...
- something to do with life-history... Demersal, large pelagic, small demersal, "missing"
(- time series length?)
- spectral predictor?
# - propensity score stuff -- max catch in the time series...
- the propensity score itself... RAM'iness
- log-total catch (summed across years)
- log(cmax)
- yrs since development...
- temperature...
- slope over 10 years with robust regression? - iucn categories

# Things to discuss
- travel dates?
- my travelling in May?
- how much I should be pursuing ETPS data?
- R package
- *train on simulated dataset or on training/testing RAM stocks?*
- if training/testing RAM: how many splits, cross-validation structure?
- costello method circularity is a major issue... only way around it I can see is by splitting the RAM stocks
- possible approaches: straight averaging, averaging with some form externally derived weights, combining via machine learning model, combining via parametric (perhaps Bayesian) model, either of the last two with additional covariates describing the 'scenario'
- a question will the assymetry of risk/loss: interesting climate paper found that fancy weighted models could be a bit better if weights were right, but also much more wrong
- 2 levels of assymetry: above below B/B_MSY and min-max style ensemble method level
- at least 2 levels of precision: precision of estimates coming out of each model, precision of estimates across models (a Bayesian-flavoured model would let us propogate the first through to the second)
- practical and keep it simple vs. theoretical best performance?
- a common theme seems to be that averaging (or weighted averaging etc.) may only perform best when evaluated across many performance metrics... it may lose out to any one model or any one performance metric
- maybe we don't *just* want to be comparing B/B_MSY in the last year or average of last 5 years
- next most important one, in my mind, is probably the trend over the last N years

# What attributes of the time series can we use to devise intelligent weights?
- if the time series looks like X and these conditions are Y, these are the weights...
- properties of the catch time series itself: 
  - absolute size (min, max, median...)
  - length
  - spectral properties
  - variability
  - general shape/trend
- properties of the species
  - need to be careful with this, since we could easily end up giving the Costello method everything it was trained with
  - pelagic, ...
  - ballpark intrinsic growth rate
  - ballpark K / B0?

# If we feed this into a machine learning model, what does this look like?
- each row of training data represents one simulated stock?
- columns represent measurable attributes of that time-series/stock?
- response column is a performance metric, or (in a multivariate analysis), a set of columns represent a suite of performance metrics?
- the output you want is a set of model weights
- or you condense all these metrics down to an average performance score and use these to derive weights
- alternatively, you could use each performance metric individually and derive a set of different possible weights... then explore the range of outputs across different possible weights
- or these metrics could represent different metrics (e.g. precision, bias) and you could then present the results in 2 or 3 dimensional space
- you want to reduce bias? look towards this corner
- if you want to balance bias and precision, look over here...
- this would be multidimensional weighting space
- then we could pick example weights from this space

- or does the machine learning model look like this:
  - you have B/Bmsy from each model as predictors and you have true B/Bmsy as the response
  - how do you incorporate other predictors (like spectral properties)?
  - surely this has been worked out — TODO look into

- What part of the time series do we care about? All? Most recent? Last 5 years? Weighted to favour recent years but also care about the past?

# Issues and caveats
- Costello method is trained on the RAM dataset. Solutions? One might be to split the RAM dataset into training and testing sets. If we want to get fancy, we could do this a number of times, say 10-fold cross-validation.
- 3 of the 4 models are very similar: probably are not "diverse learners" 

# What performance metrics could we use and what are the reasons for each?
- metrics of bias, precision, accuracy
- metrics of relative trend performance
- any of the above with equal penalty for over/under estimation
- any of the above with asymmetric penalty for over/under estimation

# What different ways can we weight the models?
- equal weighting
- based on bias/precision/accuracy metrics with respect to training data
- based on agreement/disagreement (convergence)
- based on estimated model precision (certainty)
- basing weights on linear combination of attributes of each stock
- basing weights on attributes of each stock that can interact (possibly in highly non-linear ways)

# Methods of combining
- ensemble all at once (equal vote counting)
- boosting (not feasible in reality where we don't know 'truth', but could be used in our training model weight development)

# Training datasets
- simulation dataset
- could use RAM (via bagging and possibly boosting) to train for other truly data-limited stocks
- worth comparing performance of these two on the RAM dataset?

# Conveying uncertainty
- there is uncertainty in each model estimate and across estimates
- also uncertainty in terms of difference weight choices

# Possible paper story lines
- Start with a review of the issues presented here, apply a few of them to the RAM dataset trained on the simulated data
- ... 

# How will we guage performance on the RAM dataset?
- same measures as used to calculate weights from simulated/training dataset?

# To do
- lay out with diagrams how various ensemble schemes would look
- look into weighting ensembles based on multiple performance measures
- chat with Coilin
- discuss with Ian S.?
- look into time series classification approaches
- compile table of all the various performance metrics we might use, examples of where they have been used, what they purport to represent, advantages, disadvantages
- outline potential paper story lines
- compile table of weighting options that have been used, where they were used, advantages, and disadvantages
