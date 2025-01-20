
# bias corrected estimator for the effect ratio
IPPW_IV <- function(Y, Z, X, D, prob, caliper = TRUE, gamma = 0.1, lambda, alpha){
  
  # Load optmatch
  library(optmatch)
  
  # Smahal function
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
  
  # Add caliper function
  addcaliper=function(dmat,z,logitp,calipersd=.2,penalty=1000){
    sd.logitp=sd(logitp)
    adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
    adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
    dmat=dmat+adif*penalty
    dmat
  }
  
  # Matching
  treated.index = which(Z == 1)
  propscore.model = glm(Z ~ X, family = 'binomial',x=TRUE,y=TRUE)
  treated = propscore.model$y
  Xmat=propscore.model$x[,-1]
  distmat=smahal(treated,Xmat)
  logit.propscore=predict(propscore.model)
  
  # adding caliper
  if(caliper == TRUE) {
    distmat=addcaliper(distmat,treated,logit.propscore,calipersd=.2)
    subject.index=seq(1,length(treated),1)
    rownames(distmat)=subject.index[treated==1]
    colnames(distmat)=subject.index[treated==0]
    
    matchvec=fullmatch(distmat,min.controls=0.001,max.controls=10000)
  } else {
    subject.index=seq(1,length(treated),1)
    rownames(distmat)=subject.index[treated==1]
    colnames(distmat)=subject.index[treated==0]
    
    matchvec=fullmatch(distmat,min.controls=0.001,max.controls=10000)
  }
  
  treated.subject.index=vector("list",length(treated.index))
  matched.control.subject.index=vector("list",length(treated.index))
  matchedset.index=substr(matchvec,start=3,stop=10)
  matchedset.index.numeric=as.numeric(matchedset.index)
  subjects.match.order=as.numeric(names(matchvec))
  matchedset_index = length(unique(matchedset.index.numeric))
  
  # total number in each set
  l <- rep(0,length(treated.subject.index))
  for(i in 1:length(treated.subject.index)){
    matched.set.temp=which(matchedset.index.numeric==i)
    matched.set.temp.indices=subjects.match.order[matched.set.temp]
    l[i] <- length(matched.set.temp.indices)
  }
  
  # the order of matchvec
  for(i in 1:length(treated.index)){
    matched.set.temp=which(matchedset.index.numeric==i)
    matched.set.temp.indices=subjects.match.order[matched.set.temp]
    treated.temp.index=which(matched.set.temp.indices %in% treated.index)
    if(length(treated.temp.index) != 0){
      treated.subject.index[[i]]=matched.set.temp.indices[treated.temp.index]
      matched.control.subject.index[[i]]=matched.set.temp.indices[-treated.temp.index]
    }
  }
  
  # remove null
  if(sum(sapply(treated.subject.index, is.null)) != 0){
    treated.subject.index<- treated.subject.index[-which(sapply(treated.subject.index, is.null))]
    matched.control.subject.index<-matched.control.subject.index[-which(sapply(matched.control.subject.index, is.null))]
  }
  
  # Calculate standardized differences
  treatedmat = X[Z == 1,]
  control.b = X[Z == 0,]
  controlmean.b = apply(control.b,2,mean)
  treatmean = apply(treatedmat,2,mean)
  treatvar = apply(treatedmat,2,var)
  controlvar = apply(control.b, 2, var)
  stand.diff.before=(treatmean-controlmean.b)/sqrt((treatvar+controlvar)/2)
  treatedmat.after=matrix(0,nrow=length(matched.control.subject.index),ncol=5)
  for (i in 1:length(matched.control.subject.index)) {
    if(length(treated.subject.index[[i]])>1){
      treatedmat.after[i,]=apply(X[treated.subject.index[[i]],],2,mean)
    } else {
      treatedmat.after[i,]=X[treated.subject.index[[i]],]
    }
  }
  controlmat.after=matrix(0,nrow=length(matched.control.subject.index),ncol=5)
  for (i in 1:length(matched.control.subject.index)) {
    if(length(matched.control.subject.index[[i]])>1){
      controlmat.after[i,]=apply(X[matched.control.subject.index[[i]],],2,mean)
    } else {
      controlmat.after[i,]=X[matched.control.subject.index[[i]],]
    }
  }
  controlmean.after=apply(controlmat.after,2,mean)
  treatedmean.after=apply(treatedmat.after,2,mean)
  stand.diff.after=(treatedmean.after-controlmean.after)/sqrt((treatvar+controlvar)/2)
  balance = cbind(stand.diff.before,stand.diff.after)
  
  p = conditional_p(treated.subject.index,matched.control.subject.index,prob,gamma)
  
  # Create set
  set = NULL
  for (i in 1:length(treated.subject.index)) {
    set[[i]] = c(unlist(treated.subject.index[[i]]),unlist(matched.control.subject.index[[i]]))
  }
  
  # Estimation
  nume = (Y*(Z-p))/(p*(1-p))
  deno = (D*(Z-p))/(p*(1-p))
  lambda_est_bc = sum(nume)/sum(deno)
  
  # Confidence interval
  value = rep(0,length(lambda))
  for (m in 1:length(lambda)) {
    V_i = rep(0,length(set))
    for (i in 1:length(set)) {
      index = set[[i]]
      n_i = length(set[[i]])
      V_i[i] = sum(Z[index]*(Y[index]-lambda[m]*D[index])/p[index])-sum((1-Z[index])*(Y[index]-lambda[m]*D[index])/(1-p[index])) 
    }
    
    T_lambda = sum(V_i)/length(set)
    variance_lambda = sum((V_i-T_lambda)^2)/(length(set)*(length(set)-1)) 
    value[m] = T_lambda/sqrt(variance_lambda)
  }
  
  make_interval = function(x,y){
    paste0("[",x,",",y,"]")
  }
  thre = qnorm(1-alpha/2)
  lower = lambda[which(value >= -thre & value <= thre)][1]
  upper = tail(lambda[which(value >= -thre & value <= thre)],1)
  CI = make_interval(format(lower,digits=3),format(upper,digits=4))
  
  return(list(estimate=lambda_est_bc,CI=CI,value=value,balance=balance))
}

  
  
  
  
  
  
  
  
  
  
  
  
