# BayesianBridgeRegression

* This repository contains R code to estimate the **Bridge regression in a Bayesian framework**,  using Monte Carlo methods to estimate the parameters. 
Here, the Bayesian formulation can be viewed as a trick to overcome the computational difficulties 
that arise from the constrained optimization problem defined by the Bridge model.
Our proposal is to use **Importance Sampling** to estimate the parameters instead of using the Gibbs sampler proposed in [The Bayesian Bridge , 2014].
The idea is to see if we can achieve similar results as the Gibbs sampler with a reasonable calculation time
using this much more direct method.

* The **importance sampling** method produces similar results as the **Gibbs sampler** on the same data set used in the article [The Bayesian Bridge , 2014].
Nonparametric estimation of marginal posterior coming from the importance sampling method coincides very strongly with the Gibbs sampler.
We also apply the Sequential Monte Carlo method on the bayesian bridge model in the simplest case where all the hyperparameters are fixed
and number of explanatory variables is relatively small (less that 10). We also obtain identical results with respect to the Gibbs sampler.
