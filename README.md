# BLOUCH
Bayesian Linear Ornstein-Uhlenbeck models for Comparative Hypotheses (BLOUCH) fits adaptive models of continuous trait evolution in a Bayesian framework based on categorical or continuous predictors, and incorporates measurement error. Blouch can also make phylogenetically informed predictions of known or unknown traits from any clade, given a dataset of comparative measurements and a phylogeny including the taxa of interest.

# Getting Started
If you are just getting started with blouch I recommend starting with the tutorial vignettes available on the package website. Blouch is based on an article currently in review:

+ Grabowski, M (in review). Bayesian Linear Ornstein-Uhlenbeck models for Comparative Hypotheses (BLOUCH). (arXiv preprint, code on GitHub)

# Install Instructions
To install the R and Stan functions associated with Blouch from github, first install the package devtools:
```{r}
install.packages("devtools")
library(devtools)
```
Then install blouch
```{r}
devtools::install_github("Mark-Grabowski/blouch")
library(blouch)
```

# Documentation
Please visit the package website here.

https://mark-grabowski.github.io/blouch/

# References
Hansen et al. (2008)

Stan Development Team. 2021. Stan Modeling Language Users Guide and Reference Manual, VERSION. https://mc-stan.org

Stan Development Team (2020). “RStan: the R interface to Stan.” R package version 2.21.2, http://mc-stan.org/.

