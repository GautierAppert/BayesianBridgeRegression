# Bayesian Bridge Regression

* This repository contains R code to estimate the **Bridge regression in a Bayesian framework**,  using Monte Carlo methods. 
Here, the Bayesian formulation can be viewed as a trick to overcome the computational difficulties 
that arise from the constrained optimization problem defined by the Bridge model.
Our proposal is to use **Importance Sampling** to estimate the parameters instead of using the Gibbs sampler proposed in [The Bayesian Bridge , 2014].
The idea is to see if we can achieve similar results as in [The Bayesian Bridge , 2014]  with a reasonable CPU-Time
using this much more direct approach.

* The **importance sampling** method produces similar results as the **Gibbs sampler** on the same data set used in the article [The Bayesian Bridge , 2014].
Nonparametric estimation of marginal posteriors coming from the importance sampling method coincides very strongly with the ones coming from the Gibbs sampler.

* We also apply the **Sequential Monte Carlo method (SMC)** on the bayesian bridge model in the simplest case where all the hyperparameters are fixed
and the number of explanatory variables is relatively small (less that 10). We also obtain identical results as the ones given by the Gibbs sampler in
[The Bayesian Bridge , 2014].



