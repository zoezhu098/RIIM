#' Function for generating data from two models in paper
#' @param d The number of covariates. The default is 5.
#' @param n Sample size.
#' @param C The intercept in the model. The default is -0.5.
#' @param model.type "model1" is the nonlinear logistic model. "model2" is the nonlinear selection model.
#' @param effect.type "constant" is the model with constant effect. "ATE" is the model for the average treatment effect.
#' @return A list that contains: "X" - All the covariates information in the dataset.
#'                               "Z" - The treatment assignment indicators vector in the dataset.
#'                               "R" - The observed outcomes vector in the dataset.
#'                               "P" - The true propensity scores vector in the dataset.
#'                               "R_t" - The potential outcomes vector under treatment in the dataset.
#'                               "R_c" - The potential outcomes vector under control in the dataset.
#' @import mvtnorm, VGAM
#' @export



generate_data <- function(d = 5, n = 1000, C = -0.5, model.type = c("model1","model2"), effect.type = c("constant","ATE")){

  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package \"mvtnorm\" is required for this function.", call. = FALSE)
  }

  if (!requireNamespace("VGAM", quietly = TRUE)){
    stop("Package \"VGAM\" is required for this function.", call. = FALSE)
  }

  # X_1 to X_d are normal(0, 1)
  sigma = diag(d)
  X_d = rmvnorm(n, mean = rep(0,d), sigma = sigma)
  X_d[,4] = rlaplace(n, location = 0, scale = sqrt(2)/2)
  X_d[,5] = rlaplace(n, location = 0, scale = sqrt(2)/2)


  # Generate Z
  fx = 0.3*X_d[,1] + 0.3*X_d[,2] + 0.6*cos(X_d[,3]) + 0.4*abs(X_d[,4]) + 0.4*X_d[,5] + abs(X_d[,1]*X_d[,2]) + X_d[,4]*X_d[,5] + X_d[,1]*X_d[,3] + C

  if(model.type == "model1"){
    p = exp(fx)/(1+exp(fx)) # the probability of receiving the treatment
    Z = rep(0,length(p))
    for(i in seq_along(p)){
      Z[i] = rbinom(1,1,p[i])
    }
  } else if(model.type == "model2"){
    ibu <- rnorm(n,0,1)
    Z = ifelse(fx>ibu,1,0)
    p = pnorm(fx) # the probability of receiving the treatment
  }

  # Generate two potential outcomes for each unit
  R_c = X_d[,1] + abs(X_d[,2]) + 0.2*X_d[,3]^3 +  0.2*X_d[,5] + rnorm(n,0,1)

  if(effect.type == "constant"){
    R_t = R_c + 1
  } else if(effect.type == "ATE"){
    R_t = R_c + 1 + 0.3*X_d[,1] + 0.2*X_d[,3]^3
  }

  R = (Z == 0)*R_c + (Z == 1)*R_t + 0
  # return a list of two objects: X and Z
  return(list(X = X_d, Z = Z, R = R, P = p, R_t = R_t, R_c = R_c))
}





