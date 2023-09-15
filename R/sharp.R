
#FULL MATCHING
#exact p-value/CLT
REIM_constant <- function(Y, Z, e, beta_low, beta_up, group.index, type = c("exact","CLT"), K=10000, test = c("t-test","rank-sum"), alpha, q = NULL){

  if(type == "exact"){

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
      }
    }

    #different p-value based on different beta
    beta_all = seq(beta_low,beta_up,0.01)
    exact_p <- function(Y, Z, prob, beta, test = c("t-test","rank-sum"), set.index, treated.index, K=10000){
      #adjusted response
      Y_c = Y - beta*Z

      if(test == "t-test"){
        q = Y_c
      } else if(test == "rank-sum"){
        #rank
        q = rank(Y_c,ties.method = "average")
      }

      t_obs = sum(Z*q)

      #calculate the expectation of t under null
      prob2 = prob

      for(i in 1:length(set.index)){
        if(length(treated.index[[i]]) == 1){
          index = set.index[[i]]
          prob2[index] = prob[index]
        } else {
          index = set.index[[i]]
          prob2[index] = 1-prob[index]
        }
      }

      t_beta = sum(prob2*q)

      t_k = rep(0,K)
      for (k in 1:K) {
        #cat("k = ", k, "\n")
        set.seed(k)
        Z_k = rep(0,length(Z))

        #treatment assignment
        for(i in 1:length(set.index)){
          if(length(treated.index[[i]]) == 1){
            t_index = sample(set.index[[i]],1,prob = prob[set.index[[i]]])
            Z_k[t_index] = 1
            Z_k[setdiff(set.index[[i]],t_index)] = 0
          } else {
            t_index = sample(set.index[[i]],1,prob = prob[set.index[[i]]])
            Z_k[t_index] = 0
            Z_k[setdiff(set.index[[i]],t_index)] = 1
          }
        }
        t_k[k] = sum(Z_k*q)
      }

      p_value <-  mean(ifelse(abs(t_k-t_beta) >= abs(t_obs-t_beta),1,0))
      return(p_value)
    }
    p_vall = rep(0,length(beta_all))
    for (j in seq_along(beta_all)) {
      if(test == "t-test"){
        p_vall[j] <- exact_p(Y, Z, prob, beta[j], test = "t-test", set.index, treated.index, K)
      } else {
        p_vall[j] <- exact_p(Y, Z, prob, beta[j], test = "rank-sum", set.index, treated.index, K)
      }
    }

    #calculate CI and point estimation
    CI_full = function(p_value,beta,alpha = 0.05){
      make_interval = function(x,y){
        paste0("[",x,",",y,"]")
      }
      #keep beta with the corresponding p_value > 0.05
      min_beta = min(beta[p_value >= alpha])
      max_beta = max(beta[p_value >= alpha])
      CI = make_interval(min_beta,max_beta)
      Estimate = median(beta[which.max(p_value)])
      list(CI,Estimate)
    }

    return(CI_full(p_vall,beta_all))
  }

  if(type == "CLT"){

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
      }
    }

    #different p-value based on different beta
    beta_all = seq(beta_low,beta_up,0.01)
    exact_p <- function(Y, Z, prob, beta, test = c("t-test","rank-sum"), set.index, treated.index){
      #adjusted response
      Y_c = Y - beta*Z

      if(test == "t-test"){
        q = Y_c
      } else if(test == "rank-sum"){
        #rank
        q = rank(Y_c,ties.method = "average")
      }

      #t observe
      t_obs = sum(Z*q)

      #ET
      prob = prob #for mi=1, prob is p(z=1);for ni-mi=1, prob is p(z=0)

      ET = 0
      for(i in 1:length(set.index)){
        if(length(treated.index[[i]]) == 1){
          index = set.index[[i]]
          ET = ET + sum(prob[index]*q[index])
        } else {
          index = set.index[[i]]
          ET = ET + sum((1-prob[index])*q[index])
        }
      }

      #VarT
      VarT = 0
      for(k in 1:length(set.index)){
          index = set.index[[k]]
          VarT = VarT + sum(prob[index]*(q[index]^2)) - (sum(prob[index]*q[index]))^2
      }

      normal_value = (t_obs-ET)/sqrt(VarT)
      p_value = 2*(1-pnorm(abs(normal_value)))
      return(p_value)
    }

    p_vall = rep(0,length(beta_all))
    for (l in seq_along(beta_all)) {
      if(test == "t-test"){
        p_vall[l] <- exact_p(Y, Z, prob, beta_all[l], test = "t-test", set.index, treated.index)
      } else {
        p_vall[l] <- exact_p(Y, Z, prob, beta_all[l], test = "rank-sum", set.index, treated.index)
      }
    }

  }

  #calculate CI and point estimation
  CI_full = function(p_value,beta,alpha = 0.05){
    make_interval = function(x,y){
      paste0("[",x,",",y,"]")
    }
    #keep beta with the corresponding p_value > 0.05
    min_beta = min(beta[p_value >= alpha])
    max_beta = max(beta[p_value >= alpha])
    CI = make_interval(min_beta,max_beta)
    Estimate = median(beta[which.max(p_value)])
    list(CI=CI,Estimate=Estimate)
  }

  return(CI_full(p_vall,beta_all,alpha))
}








