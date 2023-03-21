# BLOUCH
Bayesian Linear Ornstein-Uhlenbeck models for Comparative Hypotheses (BLOUCH) fits adaptive models of continuous trait evolution in a Bayesian framework based on categorical or continuous predictors, and incorporates measurement error following the approach of Hansen et al. (2008). Blouch can also make phylogenetically informed predictions of known or unknown traits from any clade, given a dataset of comparative measurements and a phylogeny including the taxa of interest.

Blouch is a Bayesian version of its frequentist brother, Slouch (Kopperud et al. 2020), which is available <a href="https://github.com/kopperud/slouch" title="here.">here</a>. While the front-end component of Blouch is written in R (R Core Team, 2015), the nuts and bolts are written in the language Stan (Carpenter et al., 2017), which allows estimation of Bayesian models using Markov chain Monte Carlo (MCMC) methods based on the Hamilton Monte Carlo sampler.

## Getting Started
If you are just getting started with blouch I recommend starting with the tutorial vignettes available on the package website. Blouch is based on an article currently in review:

+ Grabowski, M (in review). Bayesian Linear Ornstein-Uhlenbeck models for Comparative Hypotheses (BLOUCH).

## Instalation Instructions
To install the R and Stan functions associated with Blouch from github, first install the package devtools:
```{r}
install.packages("devtools")
library(devtools)
```
Then install blouch
```{r}
devtools::install_github("mark-grabowski/blouch")
library(blouch)
```

## Documentation
Please visit the package website <a href="https://mark-grabowski.github.io/blouch/" title="here.">here</a>.


## References
Carpenter, B., A. Gelman, M. D. Hoffman, D. Lee, B. Goodrich, M. Betancourt, M. Brubaker, J. Guo, P. Li, and A. Riddell. 2017. Stan: A Probabilistic Programming Language. Journal of Statistical Software 76:1–32.

Hansen, T. F., J. Pienaar, and S. H. Orzack. 2008. A comparative method for studying adaptation to a randomly evolving environment. Evolution 62:1965–1977.

Kopperud, B. T., J. Pienaar, K. L. Voje, S. H. Orzack, and T. F. Hansen. 2020. Slouch: Stochastic Linear Ornstein-Uhlenbeck Comparative Hypotheses. R package version 2.1.4.

