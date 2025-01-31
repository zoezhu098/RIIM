# RIIM: Randomization-Based Inference Under Inexact Matching

## Author
Jianan Zhu and Siyu Heng

## Maintainer
Jianan Zhu (Email: jz4698@nyu.edu)

## Description
RIIM is an R package for randomization-based inference for average treatment effects under inexact matching introduced in Zhu and Heng (2024).

Before installing this R package, please ensure that you have installed the following R packages: xgboost, MASS, mvtnorm, VGAM, and optmatch. To install package REIM in R from GitHub, please run the following commands:
```{r}
install.packages("devtools") 
library(devtools) 
install_github("zoezhu098/RIIM")
```

## Reference
Zhu, J., Zhang, J., Guo, Z., & Heng, S. (2024). Randomization-Based Inference for Average Treatment Effect in Inexactly Matched Observational Studies. arXiv:2308.02005.

Hansen, B.B. and Klopfer, S.O. (2006), ‘ Optimal full matching and related designs via network flows’, Journal of Computational and Graphical Statistics, 15, 609–627.

Hansen, B.B. (2004), ‘Full Matching in an Observational Study of Coaching for the SAT’, Journal of the American Statistical Association, 99, 609–618.

Rosenbaum, P. (1991), ‘A Characterization of Optimal Designs for Observational Studies’, Journal of the Royal Statistical Society, Series B, 53, 597–610.

Colin B.F. (2018), 'On Mitigating the Analytical Limitations of Finely Stratified Experiments', Journal of the Royal Statistical Society Series B: Statistical Methodology, 80 (5), 1035–1056.

Kang H., Kreuels B., May J., & Small S.D. (2016), 'Full matching approach to instrumental variables estimation with application to the effect of malaria on stunting', The Annals of Applied Statistics, 10 (1), 335-364.
