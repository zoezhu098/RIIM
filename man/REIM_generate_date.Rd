\name{generate_data}
\alias{generate_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate data from two models in paper}
\description{
It generates data from two models described in Zhu and Heng (2023).
}
\usage{
generate_data(d = 5, n = 1000, C = -0.5,
              model.type = c("model1","model2"),
              hypothesis.type = c("sharp","weak"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{The number of covariates. The default is 5.}
  \item{n}{Sample size.}
  \item{C}{The intercept in the model. The default is -0.5.}
  \item{model.type}{"model1" is the nonlinear logistic model. "model2" is the nonlinear selection model.}
  \item{effect.type}{"constant" is the model with constant effect. "ATE" is the model for the average treatment effect.}

}

\value{
  \item{X}{All the covariates information in the dataset.}
  \item{Z}{The treatment assignment indicators vector in the dataset.}
  \item{R}{The observed outcomes vector in the dataset.}
  \item{P}{The true propensity scores vector in the dataset.}
  \item{R_t}{The potential outcomes vector under treatment in the dataset.}
  \item{R_c}{The potential outcomes vector under control in the dataset.}



}
\references{
Zhu, J., & Heng, S. (2023). Bias Correction for Randomization-Based Estimation in Inexactly Matched Observational Studies. arXiv preprint arXiv:2308.02005.
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
# Model 1 with constant effect
data1 <- generate_data(d = 5, n = 2000, C = -2.5, model.type = "model1", effect.type = "constant")

# Model 2 with average treatment effect
data2 <- generate_data(d = 5, n = 2000, C = -2, model.type = "model2", effect.type = "ATE")



