# causalGLM

## Semiparametric generalized linear models for causal inference using targeted machine-learning 

Remarkably, it is possible to get robust and efficient inference for causal quantities using machine-learning. In the search for causal answers, assuming parametric models can be dangerous. With even a little bit of confounding they can give remarkably incorrect answers. Rather than assuming a fully parametric model, instead assume a parametric model for only the feature of the data-generating distribution that you care about. That is, assume a semiparametric model! Let the data speak for itself and use machine-learning to model the nuisance features of the data that are not directly related to your causal question. Why worry about things that don't matter for your question (and can only hurt you)? It is not worth the risk of being wrong.

In this package, we utilize targeted machine-learning (TMLE) to generalize the parametric generalized linear models commonly used for treatment effect estimation (e.g. the R package glm) to the world of semi and nonparametric models. There is virtually no loss in precision/p-values/confidence-interval-widths with these methods relative to parametric generalized linear models, but the bias reduction from these methods can be substantial! These methods even work well with small sample sizes (robust inference even when sample sizes are as small as 50-100). We employ auto-machine-learning that adapts the aggressiveness of the ML algorithms with sample size, thereby allowing for robust and correct inference in all types of settings. 


So far, this package supports:

1. Conditional average treatment effect estimation with "spCATE". (Causal semiparametric linear regression with general link functions)
2. Conditional odds ratio estimation between two binary variables with "spOR". (Causal semiparametric logistic regression)
3. Conditional relative risk regression for nonnegative outcomes and a binary treatment with "spRR". (Causal semiparametric log-linear relative-risk regression with general link functions)

The functions are easy to use. The only required input from users is the data (of course) and a R-formula object (just like glm) for the CATE/OR/RR.

Outputs include:
1. Coefficient estimates
2. Z-scores and p-values for coefficients
3. 95% confidence intervals for coefficients
4. 95% prediction/confidence intervals for evaluations of the CATE/RR/OR

 

## Targeted learning for robust efficient inference using machine-learning

Targeted learning is a general framework for using machine-learning in real-world settings to estimate causal parameters and obtain efficient inference. Targeted learning works by first estimating the data-generating distribution and conditional mean functions using parametric or nonparametric black-box machine-learning, and then performing a targeted bias correction to obtain correct inference (targeted maximum likelihood estimation).

Targeted maximum likelihood estimation (TMLE) is a generalization of the well-known maximum likelihood estimation framework but it allows for inference even when using machine-learning and variable selection procedures.  

TMLE dates back to 2000-2006 and is a state-of-the-art method for efficient semiparametric estimation and inference using machine-learning tools.


### Data-structure
All functions utilize the data-structure (W,A,Y) where
1. "W = (W_1,W_2,W_3,...)" represents a vector of baseline variables/covariates/confounders for which to adjust (passed as a named matrix or data.frame)
2. "A" represents a binary treatment or exposure assignment whose effect on the outcome we are interested in.
3. "Y" is an arbitrary outcome.

### Conditional average treatment effects (CATE) with "spCATE"
"spCATE" implements causal linear regression for additive treatment effects.  
The default model is the so-called "partially linear regression model" defined as
E[Y|A,W] = A CATE(W) + E[Y|A=0,W] where CATE(W) = E[Y|A=1,W] - E[Y|A=0,W] is user-specified and E[Y|A=0,W] is learned nonparametrically using machine-learning.


Using the argument "formula_CATE", one can specify a linear model for the CATE of the form CATE(W) = a0 + a W_1 + b W_2 + c W_3 (or whatever you want).

This function only assumes a parametric model for CATE(W) and does not assume anything about E[Y|A=0,W]. We use robust machine-learning to learn E[Y|A=0,W] and use targeted learning for valid, robust, and efficient inference.


Useful models include:

1. Constant CATE: formula_CATE = ~ 1

This specifies a constant CATE and the function will return estimates and inference for a single coefficient which can be directly interpreted as the treatment effect.

2. Effect modification and subgroup effects: formula_CATE = ~ 1 + W_1

This specifies a CATE model with effect modification by the baseline variable W_1. Estimates and inference are returned for both the intercept and coefficient in front of W_1. This can be interpreted just like linear regression-based estimates are interpreted.

3. Crazy 8 and Crazy CATE: formula_CATE = ~ 1 + W_1 + W_1*W_2 + poly(W_1, degree = 2, raw = T)

Be crazy and parametric model the CATE to your hearts content. The above formula is equivalent to CATE(W) = a + b W_1 + c W_1*W_2 + d W_1^2

No model is (or can or should) be assumed for E[Y|A=0,W]. This is estimated using default or user-specified machine-learning via sl3.
By default, a theoretically understood, flexible, robust, sparsity and smoothness adapting smoothing spline is used to estimate this nonparametric component. Specifically, we employ the R package hal9001 which implements the highly adaptive lasso estimator (HAL). This allows you to focus all your energy and attention on making a good model for CATE. No need to worry about parts of the data distribution that don't matter for what you care about!

If one wants to use a different link function for the CATE, you can pass a family object using the argument "family_CATE" (by default identity link). For example, if you want "formula_CATE = ~ 1 + W_1" to imply the exponential CATE model CATE(W) = exp(a + b * W_1) then use "family_CATE = poisson()".  


### Conditional odds ratio (OR) with "spOR"
When Y is binary, the adjusted causal odds ratio between A and Y may be of interest. "spOR" implements causal logistic regression for odds ratio estimation.
The model used is the so-called "partially-linear logistic regression model" which *only* assumes

logOR(W) := log[ {P(Y=1|A=1,W)/P(Y=0|A=1,W)} / {P(Y=1|A=0,W)/P(Y=0|A=0,W)} ] ~ user-specified parametric model.
That is, the user specifies a parametric model for the log odds between A and Y and nothing else is assumed known.

This is equivalent to assuming the logistic regression model

P(Y=1|A,W) = expit{A*logOR(W) + logit(P(Y=1|A=0,W))}

where P(Y=1|A=0,W) is unspecified and learned using machine-learning.

Using the argument "formula_logOR", one can specify a linear model for the log conditional odds ratio between A and Y of the form logOR(W) = a0 + a W_1 + b W_2 + c W_3 (or whatever you want).

Useful models include:

1. Constant odds ratio: formula_logOR ~ 1
2. Odds ratio modification and subgroup effects: formula_logOR ~ 1 + W_1 + W_2

Machine-learning estimation is done just like "spCATE".

### Conditional Relative Risk regression (RR) with "spRR"
When Y is nonnegative (e.g. binary or a count), the causal relative risk or causal relative treatment effect may be of interest. "spRR" implements causal relative risk regression (using generalized causal log-linear/poisson regression).

The model used is the so-called "partially-linear relative risk/poisson regression model" which *only* assumes

log RR(W) := log{E[Y|A=1,W] / E[Y|A=0,W]} ~ user-specified parametric model.

That is, we only assume the user specified parametric model (at the exponential scale) for the relative risk of Y with respect to A.

This is equivalent to assuming the poisson-type regression model
E[Y|A,W] = E[Y|A=0,W] exp(log RR(W)) = E[Y|A=0,W] RR(W),
where log RR(W) is parametric and E[Y|A=0,W] is the background/placebo outcome model which is unspecified and learned using machine-learning.

Using the argument "formula_logRR", one can specify a linear model for the log relative risk of the form log RR(W) = a0 + a W_1 + b W_2 + c W_3 (or whatever you want).

Useful models include:

1. Constant odds ratio: formula_logRR ~ 1
2. Odds ratio modification and subgroup effects: formula_logRR ~ 1 + W_1 + W_2

Machine-learning estimation is done just like "spCATE".

This function also supports general link functions using the "family_RR" argument in the same way as used in "spCATE".

### User friendly interface
A minimalistic yet still very flexible front-end function for all routines is provided through the "causalGLM" function. Check out the vignette to see how to use it! The only necessary arguments are: 
1. A formula object for the CATE, OR, or RR
2. The data: W, A, Y
3. Choice of estimand: "CATE", "OR", "RR"
4. That's it! Feel free to customize the machine-learning routines available using the "learning_method" argument. Built in options are: auto-HAL, glm, glmnet, gam, earth (MARS), CV-autotuned-xgboost. Cross-fitting is performed automatically. If you want to make your own learner, use the sl3_Learner argument and the tlverse/sl3 package.

## Is this like double-machine-learning?
Yes, but better! TMLE, unlike double-machine-learning (a special case of the estimating equation methodology developed by James Robins (van der Laan, Robins, 2003)), is a substitution estimator and therefore respects all constraints of the statistical model. This leads to substantially improved finite-sample performance especially in real-world settings with model misspecification and positivity/imbalance issues.

We support sample-splitting/cross-fitting through the tlverse/sl3 machine-learning pipeline which can be passed into all the implemented methods to specify machine-learning algorithms. (By default robust machine-learning is performed so user specification is not necessary). 

Example code:
devtools::install_github("tlverse/sl3", ref="devel")
library(sl3)
lrnr <- Lrnr_glm$new()
lrnr <- Lrnr_xgboost$new()
lrnr <- Lrnr_gam$new()
lrnr_cross_fit <- make_learner(Pipeline, Lrnr_cv$new(), lrnr)

Relevant reads:
https://vanderlaan-lab.org/2019/12/24/cv-tmle-and-double-machine-learning/
https://pubmed.ncbi.nlm.nih.gov/31742333/
https://digitalassets.lib.berkeley.edu/etd/ucb/text/Porter_berkeley_0028E_11248.pdf

## Need a new or specialized method?

Any confusion? Questions? Don't know which method to use? None of the methods handle your problem? Need a custom/specialized method?

Just send me a message. I would be happy to develop and implement a new method to add to this package.

## References:
Most methods are based on theory and pseudo-code provided in the working paper van der Laan (2009), some of which is also published in journals: https://core.ac.uk/download/pdf/61320177.pdf


The relative risk method is treated in Targeted Maximum Likelihood Estimation of Conditional Relative Risk in a Semi-parametric Regression Model, Tuglus et al. (2011): https://biostats.bepress.com/ucbbiostat/paper283/.
The CATE method is treated in Statistical Inference for Variable Importance, van der Laan (2006): https://biostats.bepress.com/ucbbiostat/paper188/
For machine-learnng, the package tlverse/hal9001 and tlverse/sl3 are used: https://github.com/tlverse/hal9001

See also:

Estimation of a non-parametric variable importance measure of a continuous exposure, Chambaz et al. (2012): https://projecteuclid.org/journals/electronic-journal-Nonparametricof-statistics/volume-6/issue-none/Estimation-of-a-non-parametric-variable-importance-measure-of-a/10.1214/12-EJS703.full

Causal effects based on marginal structural models, Neugebauer, van der Laan (2007): 
https://www.researchgate.net/publication/222318646_Nonparametric_causal_effects_based_on_marginal_structural_models  

Related R packages: 

https://github.com/ck37/varimpact/tree/master/R
https://academic.oup.com/bioinformatics/article/31/18/3054/241218
https://cran.case.edu/web/packages/tmle.npvi/tmle.npvi.pdf

For fully nonparametric ATE-type methods, see the tmle package: https://cran.r-project.org/web/packages/tmle/index.html
Or tlverse/tmle3: https://tlverse.org


 
 


