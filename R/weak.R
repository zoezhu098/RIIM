
#FULL MATCHING
REIM_ATE = function(Y, Z, e, gamma = 0.1, group.index, alpha){

  #create set.index and treated.index
  all_groups = unique(group.index)
  set.index <- treated.index <- NULL
  for (i in 1:length(all_groups)) {
    set.index[[i]] <- which(group.index == i)
    treated.index[[i]] <- set.index[[i]][Z[set.index[[i]]] == 1]
  }

  #calculate probability based on propensity score
  prob = rep(0,length(e))
  for(z in 1:length(treated.index)){
    if(length(treated.index[[z]]) == 1){
      index = set.index[[z]]
      n = length(index)
      gi_all = rep(0,n)
      for(j in 1:n){
        p_i1 = e[index[j]]
        p_i2_ni = e[index[-j]]
        g_i1 = p_i1*prod(1-p_i2_ni)
        gi_all[j] = g_i1
      }
      for(m in 1:n){
        prob[index[m]] =  gi_all[m]/sum(gi_all)
      }

      if(sum(prob[index]<gamma | prob[index]>(1-gamma))>0){
        prob[index] <- 1/length(index)
      }

    } else {
      index = set.index[[z]]
      n = length(index)
      hi_all = rep(0,n)
      for(j in 1:n){
        p_i1 = e[index[j]]
        p_i2_ni = e[index[-j]]
        h_i1 = (1-p_i1)*prod(p_i2_ni)
        hi_all[j] = h_i1
      }
      for(m in 1:n){
        prob[index[m]] =  hi_all[m]/sum(hi_all)
      }

      if(sum(prob[index]<gamma | prob[index]>(1-gamma))>0){
        prob[index] <- 1/length(index)
      }
    }
  }

  N = length(Y)
  est = rep(0,length(treated.index))
  for(i in 1:length(treated.index)){
    if(length(treated.index[[i]]) == 1){
      index = set.index[[i]]
      n_i = length(index)
      beta1 = 1/(n_i*prob[index])
      beta2 = (n_i-1)/(n_i*(1-prob[index]))
      est[i] = n_i*(mean(sum(Z[index]*Y[index]*beta1) - sum((1-Z[index])*Y[index])*beta2/(n_i-1)))/N

    } else {
      index = set.index[[i]]
      n_i = length(index)
      beta1 = (n_i-1)/(n_i*(1-prob[index]))
      beta2 = 1/(n_i*prob[index])
      est[i] = n_i*(mean(sum(Z[index]*Y[index]*beta1/(n_i-1)) - sum((1-Z[index])*Y[index])*beta2))/N

    }
  }
  weak_est = sum(est)


  #CI
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
  thre = qnorm(1-alpha/2)
  low = weak_est-thre*variance[1,1]
  up = weak_est+thre*variance[1,1]
  CI = make_interval(format(low,digits=3),format(up,digits=4))

  return(list(CI=CI,Estimate=weak_est))

}






