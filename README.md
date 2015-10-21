# BayesianBridgeRegression
Bridge regression in a Bayesian framework,  using Monte Carlo methods in order to estimate the parameters. 
Estimation of the Bayesian Bridge regression using importance sampling method and Gibbs sampler (article The Bayesian Bridge (2014)).
Bayesian formulation can be an alternative against the constrained optimization problem
defined by the bridge model.
The importance sampling method is based on the simulation of the parameters via an instrumental law that is easy to simulate using a fast and reliable simulator.
The underlying idea is to see if we can get similar results as the Gibbs sampler with a reasonable calculation time
using this much more direct method.
The importance sampling method produces similar results in comparison with the method of Gibbs sampler
on the article data sets.
Specifically nonparametric estimation of marginal posterior densities of parameters according to the importance sampling method and Gibbs sampler coincides very strongly.
We apply the Sequential Monte Carlo method on the bayesian bridge model in the simplest case where all the hyperparameters are fixed and number of explanatory variables is relatively small ($\leq 10$).
We observe that we obtain identical results of the estimates of the Gibbs sampler method using the Monte Carlo sequential algorithm.
