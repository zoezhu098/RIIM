
# FULL MATCHING
REIM_ATE = function(Y, Z, p, set.index, treated.index, beta){

  # Estimation
  tae_weight = rep(0,length(set.index))
  for (i in 1:length(set.index)) {
    index = set.index[[i]]
    n = length(set.index[[i]])
    tae = sum(Z[index]*Y[index]/(n*p[index]))-sum((1-Z[index])*Y[index]/(n*(1-p[index])))
    tae_weight[i] = n*tae/length(Z)
  }
  tae_all = sum(tae_weight)

  # CI
  N = length(Y)
  B = length(set.index)
  n_vector = sapply(set.index,length)
  w_vector = B*n_vector/N
  W = diag(w_vector,B,B)
  Q = matrix(0,B,2)
  Q[,1] = 1
  Q[,2] = w_vector-1
  H_Q = Q%*%ginv((t(Q)%*%Q))%*%t(Q)
  y = rep(0,length(B))
  for(i in seq_along(1:B)){
    n1 = length(treated.index[[i]])
    n0 = length(set.index[[i]]) - n1
    t_hat = sum(Z[set.index[[i]]]*Y[set.index[[i]]]/n1)-sum((1-Z[set.index[[i]]])*Y[set.index[[i]]]/n0)
    y[i] = t_hat/sqrt(1-H_Q[i,i])
  }
  I = diag(1,B,B)
  variance = sqrt((1/B^2)*t(y)%*% W %*%(I-H_Q)%*%W%*%y)

  make_interval = function(x,y){
    paste0("[",x,",",y,"]")
  }
  thre = qnorm(1-beta/2)
  low = tae_all-thre*variance[1,1]
  up = tae_all+thre*variance[1,1]
  CI = make_interval(format(low,digits=3),format(up,digits=4))

  return(list(CI=CI,Estimate=tae_all))

}






