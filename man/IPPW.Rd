\name{IPPW}
\alias{IPPW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The inverse post-matching probability weighting (IPPW) estimator for the average treatment effect}
\description{
The inverse post-matching probability weighting (IPPW) estimator for randomization-based estimation and inference for the average treatment effect under inexact matching, proposed in Zhu and Heng (2024). Optimal full matching is conducted before the estimation and inferencce.
}

\usage{
IPPW = function(Y, Z, X, prob, caliper = TRUE, dim = FALSE, gamma = 0.1, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{The observed outcome vector.}
  \item{Z}{The treatment indicator vector.}
  \item{X}{The covariates matrix.}
  \item{prob}{The propensity score vector.}
  \item{min.control}{The minimum ratio of controls to treatments permitted within a matched set. The default is 0.001.}
  \item{max.control}{The maximum ratio of controls to treatments permitted within a matched set. The default is 10000.}
  \item{caliper}{Whether adding caliper or not. The default is TRUE.}
  \item{calipersd}{The standard deviation of the logit propensity score for caliper. The default is 0.2.}
  \item{dim}{Whether using difference-in-means estimator to estimate the average treatment effect. The default is FALSE.}
  \item{gamma}{The regularization threshold. The default is 0.1.}
  \item{alpha}{The prespecified level alpha for the CI.}

}

\value{
  \item{estimate}{The estimation for the average treatment effect on the IPPW estimator.}
  \item{var}{The variance for the average treatment effect on the IPPW estimator.}
  \item{low}{The lower bound for the average treatment effect on the IPPW estimator.}
  \item{up}{The upper bound for the average treatment effect on the IPPW estimator.}
  \item{CI}{The confidence interval for the average treatment effect on the IPPW estimator.}
  \item{balance}{The pre- and post-matching covariate balance table.}

}
\references{
Zhu, J., Zhang, J., Guo, Z., & Heng, S. (2024). Randomization-Based Inference for Average Treatment Effect in Inexactly Matched Observational Studies. arXiv:2308.02005.

Rosenbaum, P. (1991), ‘A Characterization of Optimal Designs for Observational Studies’, Journal of the Royal Statistical Society, Series B, 53, 597–610.
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

# generate Y 
Y_0 = 0.2*(X_d[,1])^3 + 0.2*abs(X_d[,2]) + 0.2*X_d[,3]^3 + 0.5*abs(X_d[,4]) + 0.3*X_d[,5] + rnorm(n,0,1)
Y_1 = Y_0 + 1 + 0.3*X_d[,1] + 0.2*X_d[,3]^3
Y = (1-Z)*Y_0 + Z*Y_1

# Use XGBoost to estimate propensity score
X = X_d
df = data.frame(Z,X)
index_model1 = sample(400,200)
df1 = df[index_model1,]
df2 = df[-index_model1,]
prob_xgb = rep(0,400)
xgb.model1 = xgboost(data = as.matrix(df1[2:6]), label = df1$Z, nrounds = 2, objective = "binary:logistic",verbose = 0)
prob_xgb[-index_model1] = predict(xgb.model1, as.matrix(df2[2:6]))
xgb.model2 = xgboost(data = as.matrix(df2[2:6]), label = df2$Z, nrounds = 2, objective = "binary:logistic",verbose = 0)
prob_xgb[index_model1] = predict(xgb.model2, as.matrix(df1[2:6]))

# The estimation
est = IPPW(Y,Z,X,prob_xgb,caliper=TRUE,dim=FALSE,gamma=0.1,0.05)$estimate

# The Confidene Intervel
CI = IPPW(Y,Z,X,prob_xgb,caliper=TRUE,dim=FALSE,gamma=0.1,0.05)$CI

# The balance table
balance = IPPW(Y,Z,X,prob_xgb,caliper=TRUE,dim=FALSE,gamma=0.1,0.05)$balance

