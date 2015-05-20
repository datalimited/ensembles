Book: Introduction to Machine Learning, second edition Ch. 17 

@alpaydin2010

Awesome chapter on this topic

Important questions are:
1. How do we generate models that complement each other?
2. How do we combine their output?

Combining models doesn't necessarily improve accuracy, but it always requires more computational time and complexity

Important to generate "diverse learners"

Diverse learners can come from fundamentally different models (say in model structure or different underlying assumptions...)

Diverse learners can also come from different training datasets (or even different subsets of the same dataset, i.e. bagging)

Another way of thinking about diversity, is they should ideally target different subdomains of the problem

How much the various models can specialize also depends on how the models will be combined
e.g. if it's a vote counting scheme then they all vote on all data, so they can't be too specialized
but if you can figure out where some are good and others are bad, you can exploit that in how you combine them

Types of model combination schemes
- voting or stacking: known as 'learner fusion' in a 'global approach'
- gating: input is looked at and one (or a few) models are chosen
- multistage: use one model, then use other models only on instances that aren't fit well by previous models ('cascading' is an example, this sounds like 'boosting' to me)

If all models are given the same weight and averaged, this is a form of 'voting'
also called a 'linear opinion pool'

Outputs may need to be normalized first to the same scale

Could be combined by sum (average) weighted average, median, min, max, product...

In product combination, each model has veto power!

Accuracy weights can be derived from a separate testing dataset

'Bayesian model combination': weights become priors, model decisions are like likelihoods
- in this context, simple averaging is like using uniform priors

'Bagging' is short for bootstrap aggregating

Suggest to use median when combining bagging with continuous data

Original boosting algorithm reference: Schapire 1990
- a sequence of splitting your data, fitting one model and testing with another... see the chapter
- but this requires a lot of data

AdaBoost: adaptive boosting, see Freund and Schapire 1996
- uses same training set repeatedly
- learner must be sufficiently weak
- use misfit to increase probability of drawing that row of data
- at end, models are combined weighted by the models accuracy on each training set
- benefit comes from increasing the 'margin'
- basically this is bagging, but which rows you pick each time aren't (entirely) left to chance

'Mixture of experts' method
- different models are given more or less weight on different parts of the data

'Stacked generalization'
- see Wolpert 1992
- combine models in a not-necessarily linear way
- the combiner model should be trained on different unused data to the base-learning models

(Now at page 439)

'Cascading'
- sequence of models
- ordered in terms of cost/complexity (say time to run)
- and we use the later models only if the previous models are not confident
- this differs from, say, AdaBoost: not only are we prioritizing the estimates we got wrong, *we're also prioritizing the estimates we weren't certain about*
- in a sense: you're starting simple and then adding exceptions
- first ones might be linear, latter ones may be highly non-linear

What we're doing (and what is discussed above) is called 'late integration' where models are made separately and combined; this is probably superior to 'early integration'
