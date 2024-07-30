
# FULL MATCHING
ATE_CI = function(Y, Z, p1, p2, set.index, treated.index, beta){

  # Estimation
  tae_weight = rep(0,length(set.index))
  for (i in 1:length(set.index)) {
    index = set.index[[i]]
    n = length(set.index[[i]])
    tae = sum(Z[index]*Y[index]/(n*p1[index]))-sum((1-Z[index])*Y[index]/(n*(1-p1[index])))
    tae_weight[i] = n*tae/length(Z)
  }
  tae_all = sum(tae_weight)

  # CI
  N = length(Y)
  B = length(set.index)
  n_vector = sapply(set.index,length)
  w_vector = B*n_vector/N
  W = diag(w_vector,B,B)
  Q = matrix(0,B,1)
  Q[,1] = 1
  H_Q = Q%*%ginv((t(Q)%*%Q))%*%t(Q)
  y = rep(0,length(B))
  for(i in seq_along(1:B)){
    index = set.index[[i]]
    ni = length(index)
    tae = sum((Z[index]*Y[index])/(ni*p2[index]))-sum(((1-Z[index])*Y[index])/(ni*(1-p2[index])))
    y[i] = tae/sqrt(1-H_Q[i,i])
  }
  I = diag(1,B,B)
  variance = (1/B^2)*t(y)%*%W%*%(I-H_Q)%*%W%*%y

  make_interval = function(x,y){
    paste0("[",x,",",y,"]")
  }
  thre = qnorm(1-beta/2)
  low = tae_all-thre*sqrt(variance[1,1])
  up = tae_all+thre*sqrt(variance[1,1])
  CI = make_interval(format(low,digits=3),format(up,digits=4))
  variance1 = variance[1,1]

  return(list(CI=CI,Estimate=tae_all,var=variance1,low=low,up=up))

}






