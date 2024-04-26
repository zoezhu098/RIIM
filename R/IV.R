
# bias corrected estimator for the effect ratio
bias_corr_er <- function(Z,Y,D,p){
  nume = (Y*(Z-p))/(p*(1-p))
  deno = (D*(Z-p))/(p*(1-p))
  lambda_est_bc = sum(nume)/sum(deno)
  return(lambda_est_bc)
}

