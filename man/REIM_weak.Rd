\name{REIM_ATE}
\alias{REIM_ATE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The bias-corrected difference-in-means estimator for the average treatment effect}
\description{
The bias-corrected difference-in-means estimator for randomization-based estimation and inference for the average treatment effect under inexact matching, proposed in Zhu and Heng (2023).
}

\usage{
REIM_ATE(Y, Z, e, gamma = 0.1, group.index, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{The observed outcome vector.}
  \item{Z}{The treatment indicator vector.}
  \item{e}{The estimated propensity score vector.}
  \item{gamma}{Prespecified small number. The default is 0.1.}
  \item{group.index}{The matched set index vector.}
  \item{alpha}{The prespecified level alpha for the CI.}

}

\value{
  \item{Estimate}{The bias-corrected estimation of the average treatment effect.}
  \item{CI}{The bias-corrected confidence interval for the average treatment effect.}

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

# Model 2 with average treatment effect
set.seed(1)
data1 <- generate_data(d = 5, n = 2000, C = -2, model.type = "model2", effect.type = "ATE")

#function smahal
smahal=
  function(z,X){
  X<-as.matrix(X)
  n<-dim(X)[1]
  rownames(X)<-1:n
  k<-dim(X)[2]
  m<-sum(z)
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,m,n-m)
  Xc<-X[z==0,]
  Xt<-X[z==1,]
  rownames(out)<-rownames(X)[z==1]
  colnames(out)<-rownames(X)[z==0]
  library(MASS)
  icov<-ginv(cv)
  for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
  out
}

# Fit a propensity score using logistic regression
Z = data1$Z
X = data1$X
Y = data1$R
e = data1$P
propscore.model = glm(Z ~ X, family = 'binomial',x=TRUE,y=TRUE)
treated = propscore.model$y
Xmat=propscore.model$x[,-1]
distmat=smahal(treated,Xmat)
logit.propscore=predict(propscore.model)
subject.index=seq(1,length(treated),1)
rownames(distmat)=subject.index[treated==1]
colnames(distmat)=subject.index[treated==0]

# Full matching
library(optmatch)
matchvec=fullmatch(distmat,min.controls=0.001,max.controls=10000)
treated.index = which(data1$Z == 1)
treated.subject.index=vector("list",length(treated.index))
matched.control.subject.index=vector("list",length(treated.index))
matchedset.index=substr(matchvec,start=3,stop=10)
matchedset.index.numeric=as.numeric(matchedset.index)
subjects.match.order=as.numeric(names(matchvec))
matchedset_index = length(unique(matchedset.index.numeric))

# Total number in each set
l <- rep(0,length(treated.subject.index))
for(i in 1:length(treated.subject.index)){
  matched.set.temp=which(matchedset.index.numeric==i)
  matched.set.temp.indices=subjects.match.order[matched.set.temp]
  l[i] <- length(matched.set.temp.indices)
}

# The order of matchvec
for(i in 1:length(treated.index)){
  matched.set.temp=which(matchedset.index.numeric==i)
  matched.set.temp.indices=subjects.match.order[matched.set.temp]
  treated.temp.index=which(matched.set.temp.indices %in% treated.index)
  if(length(treated.temp.index) != 0){
    treated.subject.index[[i]]=matched.set.temp.indices[treated.temp.index]
    matched.control.subject.index[[i]]=matched.set.temp.indices[-treated.temp.index]
  }
}

#remove null
if(sum(sapply(treated.subject.index, is.null)) != 0){
  treated.subject.index<- treated.subject.index[-which(sapply(treated.subject.index, is.null))]
  matched.control.subject.index<-matched.control.subject.index[-which(sapply(matched.control.subject.index, is.null))]
}

# Get the group.index
set = NULL
for (i in 1:length(treated.subject.index)) {
  set[[i]] = c(unlist(treated.subject.index[[i]]),unlist(matched.control.subject.index[[i]]))
}

Ni = max(sapply(set,length))

group.index <- rep(0,length(Y))
for(j in 1:length(Y)){
  if(Z[j] == 0){
    group.index[j] <- which(sapply(matched.control.subject.index,function(x){
      is.element(j,x)
    }))
} else if(Z[j] == 1){
    group.index[j] <- which(sapply(treated.subject.index,function(x){
      is.element(j,x)
    }))
  }
}

# The estimator
est = REIM_ATE(Y, Z, e, gamma = 0.1, group.index, alpha = 0.05)[[2]]

# The Confidene Intervel
CI = REIM_ATE(Y, Z, e, gamma = 0.1, group.index, alpha = 0.05)[[1]]
