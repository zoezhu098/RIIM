\name{IV_CI}
\alias{IV_CI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The bias-corrected Wald estimator for for the complier average treatment effect}
\description{
The bias-corrected Wald estimator for randomization-based estimation and inference for the average treatment effect under inexact matching, proposed in Zhu and Heng (2024). Optimal full matching is conducted before the estimation and inferencce.
}
\usage{
IV_CI = function(Y, Z, X, D, prob, caliper = TRUE, gamma = 0.1, lambda, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{The observed outcome vector.}
  \item{Z}{The instrument vector.}
  \item{X}{The covariates matrix.}
  \item{D}{The treatment indicator vector.}
  \item{prob}{The propensity score vector.}
  \item{caliper}{Whether adding caliper or not. The default is TRUE. The caliper is 0.2 standard deviation of the logit propensity score.}
  \item{gamma}{Prespecified small number. The default is 0.1.}
  \item{lambda}{The candidate vector for the complier average treatment effect.}
  \item{alpha}{The prespecified level alpha for the CI.}

}

\value{
  \item{estimate}{The estimation for the average treatment effect on the IPPW estimator.}
  \item{CI}{The confidence interval for the average treatment effect on the IPPW estimator.}
  \item{value}{The corresponding z-scores for all the lambda provided.}
  \item{balance}{The pre- and post-matching covariate balance table.}



}
\references{
Zhu, J., Zhang, J., Guo, Z., & Heng, S. (2024). Randomization-Based Inference for Average Treatment Effect in Inexactly Matched Observational Studies. arXiv:2308.02005.

Rosenbaum, P. (1991), ‘A Characterization of Optimal Designs for Observational Studies’, Journal of the Royal Statistical Society, Series B, 53, 597–610.
}
\author{
Jianan Zhu (maintainer, \url{jz4698@nyu.edu}) and Siyu Heng.
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

#joint distribution
matrix = matrix(c(1,0.8,0.8,1),2,2)
sigma = mvrnorm(n,c(0,0),matrix)
  
# generate the treatment effect D:
fx_D0 = 0.7*X_d[,1] + 0.4*sin(X_d[,2]) + 0.4*abs(X_d[,3]) + 0.6*X_d[,4] + 0.1*X_d[,5] + 0.3*X_d[,3]*X_d[,4] - 1 + sigma[,1]
ibu <- rnorm(n,0,1)
D_0 = ifelse(fx_D0>ibu,1,0)
  
fx_D1 = fx_D0 + 2 + 0.8*X_d[,2]^2
D_1 = ifelse(fx_D1>ibu,1,0)
D = (1-Z)*D_0 + Z*D_1
  
# generate continuous outcome Y:
Y_0 = 0.4*(X_d[,1])^2 + 0.1*abs(X_d[,2]) + 0.1*X_d[,3]^2 + 0.2*cos(X_d[,4]) + 0.5*sin(X_d[,5]) + sigma[,2]
Y_1 = Y_0 + 1 + 0.1*X_d[,1] + 0.3*X_d[,3]^2
Y = (1-D)*Y_0 + D*Y_1

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
est = IV_CI(Y,Z,X,D,prob_xgb,caliper=TRUE,gamma=0.1,lambda = seq(0,3,0.01),0.05)$estimate

# The Confidene Intervel
CI = IV_CI(Y,Z,X,D,prob_xgb,caliper=TRUE,gamma=0.1,lambda = seq(0,3,0.01),0.05)$CI

# The balance table
balance = IV_CI(Y,Z,X,D,prob_xgb,caliper=TRUE,gamma=0.1,lambda = seq(0,3,0.01),0.05)$balance



























