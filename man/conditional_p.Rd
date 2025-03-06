\name{conditional_p}
\alias{conditional_p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The function for calculating the post-matching treatment assignment probabilities}
\description{
This function calculates the regularized post-matching treatment assignment probabilities, as described in Zhu, Zhang, Guo, and Heng (2024).
}
\usage{
conditional_p = function(treated.subject.index,matched.control.subject.index,p,gamma=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treated.subject.index}{The index list for treated subjects after matching.}
  \item{matched.control.subject.index}{The index list for control subjects after matching.}
  \item{p}{The estimated propensity score vector.}
  \item{gamma}{Prespecified small number as the regularization threshold. The default is 0.1.}

}

\value{
  \item{prob}{The regularized post-matching treatment assignment probabilities vector.}


}
\references{
Zhu, J., Zhang, J., Guo, Z., and Heng, S. (2024). Randomization-Based Inference for Average Treatment Effect in Inexactly Matched Observational Studies. arXiv:2308.02005.
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

