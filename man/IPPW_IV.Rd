\name{IPPW_IV}
\alias{IPPW_IV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The bias-corrected Wald estimator for the complier average treatment effect}
\description{
The bias-corrected Wald estimator for randomization-based estimation and inference for the average treatment effect under inexact matching, proposed in Zhu, Zhang, Guo and Heng (2024). Optimal full matching is conducted before the estimation and inferencce. The propensity scores are estimated by XGBoost in function.
}
\usage{
IPPW_IV = function(Y, Z, X, D, min.controls = 0.001, max.controls = 10000, caliper = TRUE, calipersd = 0.2, classical = FALSE, gamma = 0.1, lower.bound, upper.bound, by,alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{The observed outcome vector.}
  \item{Z}{The instrument vector.}
  \item{X}{The covariates matrix.}
  \item{D}{The treatment indicator vector.}
  \item{min.control}{The minimum ratio of controls to treatments permitted within a matched set. The default is 0.001.}
  \item{max.control}{The maximum ratio of controls to treatments permitted within a matched set. The default is 10000.}
  \item{caliper}{Whether adding caliper or not. The default is TRUE.}
  \item{calipersd}{The standard deviation of the logit propensity score for caliper. The default is 0.2.}
  \item{classical}{Whether using the classical estimator (proposed by Kang H.(2016)) to estimate the complier average treatment effect. The default is FALSE.}
  \item{gamma}{The regularization threshold. The default is 0.1.}
  \item{lower.bound}{The starting value of the search region for the point estimate.}
  \item{upper.bound}{The end value of the search region for the point estimate.}
  \item{by}{The increment of the search region for the point estimate.}
  \item{alpha}{The prespecified level alpha for the CI.}

}

\value{
  \item{estimate}{The estimation for the average treatment effect on the IPPW estimator.}
  \item{CI}{The confidence interval for the average treatment effect on the IPPW estimator.}
  \item{value}{The corresponding z-scores for all the lambda provided.}
  \item{balance}{The pre- and post-matching covariate balance table.}



}
\references{
Zhu, J., Zhang, J., Guo, Z., & Heng, S. (2024). Randomization-Based Inference for Average Treatment Effect in Inexactly Matched Observational Studies. arXiv preprint, arXiv:2308.02005.

Hansen, B. B., & Klopfer, S. O. (2006). Optimal full matching and related designs via network flows. Journal of computational and Graphical Statistics, 15(3), 609-627.

Hansen, B. B. (2004). Full matching in an observational study of coaching for the SAT. Journal of the American Statistical Association, 99(467), 609-618.

Rosenbaum, P. R. (1991). A characterization of optimal designs for observational studies. Journal of the Royal Statistical Society: Series B (Methodological), 53(3), 597-610.

Fogarty, C. B. (2018). On mitigating the analytical limitations of finely stratified experiments. Journal of the Royal Statistical Society Series B: Statistical Methodology, 80(5), 1035-1056.

Kang, H., Kreuels, B., May, J., & Small, D. S. (2016). Full matching approach to instrumental variables estimation with application to the effect of malaria on stunting. The Annals of Applied Statistics, 10(1), 335-364.

}
\author{
Jianan Zhu (maintainer, \url{jz4698@nyu.edu}), Jeffrey Zhang, Zijian Guo and Siyu Heng.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(MASS)
library(VGAM)
library(xgboost)
library(mvtnorm)
library(optmatch)

# Generate data
set.seed(1)
d = 5
n = 400
sigma = diag(d)

# generate X
X_d = rmvnorm(n, mean = rep(0,d), sigma = sigma)
X_d[,4] = rlaplace(n, location = 0, scale = sqrt(2)/2)
X_d[,5] = rlaplace(n, location = 0, scale = sqrt(2)/2)

# generate Z
C = -2.5 
fx = 0.1*(X_d[,1])^3 + 0.3*(X_d[,2]) + 0.2*log((X_d[,3])^2) + 0.1*(X_d[,4]) + 0.2*X_d[,5] + abs(X_d[,1]*X_d[,2]) + (X_d[,3]*X_d[,4])^2 + 0.5*(X_d[,2]*X_d[,4])^2 + rnorm(n,0,1) + C
p = exp(fx)/(1+exp(fx)) # the probability of receiving the treatment
Z = rep(0,length(p))
for(i in seq_along(p)){
  Z[i] = rbinom(1,1,p[i])
}

#joint distribution
matrix = matrix(c(1,0.8,0.8,1),2,2)
sigma = mvrnorm(n,c(0,0),matrix)
  
# generate the treatment effect D:
fx_D0 = 0.7*X_d[,1] + 0.4*sin(X_d[,2]) + 0.4*abs(X_d[,3]) + 0.6*X_d[,4] 
+ 0.1*X_d[,5] + 0.3*X_d[,3]*X_d[,4] - 1 + sigma[,1]
ibu <- rnorm(n,0,1)
D_0 = ifelse(fx_D0>ibu,1,0)
  
fx_D1 = fx_D0 + 2 + 0.8*X_d[,2]^2
D_1 = ifelse(fx_D1>ibu,1,0)
D = (1-Z)*D_0 + Z*D_1
  
# generate continuous outcome Y:
Y_0 = 0.4*(X_d[,1])^2 + 0.1*abs(X_d[,2]) + 0.1*X_d[,3]^2 + 0.2*cos(X_d[,4]) 
+ 0.5*sin(X_d[,5]) + sigma[,2]
Y_1 = Y_0 + 1 + 0.1*X_d[,1] + 0.3*X_d[,3]^2
Y = (1-D)*Y_0 + D*Y_1

# The estimation
est = IPPW_IV(Y,Z,X,D,min.controls = 0.001, max.controls = 10000,caliper=TRUE, calipersd = 0.2,classical = FALSE,gamma = 0.1,lower.bound=0,upper.bound=3,by=0.01,alpha=0.05)$estimate

# The Confidene Intervel
CI = IPPW_IV(Y,Z,X,D,min.controls = 0.001, max.controls = 10000,caliper=TRUE, calipersd = 0.2,classical = FALSE,gamma = 0.1,lower.bound=0,upper.bound=3,by=0.01,alpha=0.05)$CI

# The balance table
balance = IPPW_IV(Y,Z,X,D,min.controls = 0.001, max.controls = 10000,caliper=TRUE, calipersd = 0.2,classical = FALSE,gamma = 0.1,lower.bound=0,upper.bound=3,by=0.01,alpha=0.05)$balance


