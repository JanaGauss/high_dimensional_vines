# Properties of stepwise parameter estimation in high-dimensional vine copulas

Code and data to reproduce the results of the paper *Gauss and Nagler. Properties of stepwise parameter estimation in high-dimensional vine copulas*.

## Code and data:

The `code` folder contains all files to run the simulations. 
Aggregated results (generated using the `create_data.R` file) are provided in the `results` folder.
All figures can be reproduced using the `estimate_parameter_results.R` and `assumptions_results.R` files.

## Abstract

The increasing use of vine copulas in high-dimensional settings, where the number of parameters is often of the same order as the sample size, calls for asymptotic theory beyond the traditional fixed-$p$, large-$n$ framework. We establish consistency and asymptotic normality of the stepwise maximum likelihood estimator for vine copulas when the number of parameters diverges as $n \to \infty$. Our theoretical results cover both parametric and nonparametric estimation of the marginal distributions, as well as truncated vines, and are also applicable to general estimation problems, particularly other sequential procedures. Numerical experiments suggest that the derived assumptions are satisfied if the pair copulas in higher trees converge to independence copulas sufficiently fast. A simulation study substantiates these findings and identifies settings in which estimation becomes challenging. In particular, the vine structure strongly affects estimation accuracy, with D-vines being more difficult to estimate than C-vines, and estimates in Gumbel vines exhibit substantially larger biases than those in Gaussian vines.

