---
bibliography: fisheries-ensembles.bib
output: word_document
---

# Ensembles discussion

* 1st paragraph: summarize the findings
    * ensembles might not always be the best, but they're consistently pretty
      good (or at worst you're not far off) across many performance dimensions,
      response variables, datasets, and ensemble types
    * probably pick something non-linear and non-parametric, but even an
      average may be useful
    * ensembles are powerful: we should consider them more often in ecology

* What makes for a good ensemble?
    * diverse 'learners' that are good predictors themselves 
      [e.g. @ali1996; @tebaldi2007]
    * not overfit: benefits of getting fancy ensembles right vs the consequences 
      of getting them wrong [e.g. @weigel2010] 
        * probability of overfitting increases with more models [e.g. @caruana2004]
        * 'naive' mean performed fairly well here... never too wrong... especially
          on entirely new (RAM) data
        * this is a benefit to methods like random forest that internally 
          cross-validate to reduce overfitting risks
        * importance of representative training data
    
* Similarities and differences to the ubiquitous coefficient averaging in 
  ecology (i.e. compare to Burnham and Anderson)
    * based on similar ethos / many similar ideas
    * one is based on combining coefficients, the other on the predictions
    * ensembles a more general purpose tool: doesn't require information 
      theoretics, can combine completely different types of models (as done
      here) (e.g. combine parametric and non-parametric models, combine
      frequentist and Bayesian predictions)
    * ensembles can exploit non-linear interactions between individual models and
      additional covariates

* The importance of tailoring the response variable to the question of interest
    * example of mean and slope here, another example is above or below some 
      threshold
    * no reason you can't have multiple ensembles that use the same individual 
      models to answer different questions
    * individual models optimally combine in different ways depending on the 
      response

* How ensembles can help understand mechanistic underpinnings of individual 
  models (Coilin)
    * discover surprises: counterintuitive relationships and interactions
    * can dig into surprising conditions and learn about why certain models
    * fail; potentially learn how to improve models
    * can learn under what conditions certain models perform well
    * pick out an example

* Other possible applications of ensemble models outside fisheries status and
  trajectory
    * other taxa examples
    * other kinds of response variables
	
* Other possible topics:
    * how much better are the ensembles really? translate the results into
      some concrete examples in an absolute not relative sense
    * other types of ensembles: assigning weights, incorporating uncertainty
      around estimates, GAMs, ...
    * when wouldn't you want to use ensembles? perhaps if models represent 
      dichotomous assumptions where either one or the other is right and imply
      different management actions? ...
    * your ideas here...
	
# References
