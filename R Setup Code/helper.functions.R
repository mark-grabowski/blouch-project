calc_adaptive_dmX<-function(a,T_term,X){
  N<-length(T_term);
  if(is.null(dim(X))==FALSE){Z<-dim(X)[2]}else{Z<-1}
  rho<-(1 - (1 - exp(-a * T_term))/(a * T_term))
  rhos<-matrix(rep(rho,Z),nrow=N,ncol=Z)
  dmX<-X * rhos
  return(dmX)
}

calc_mixed_dmX<-function(a,T_term,X,Z_direct,Z_adaptive){
  N<-length(T_term);
  rho<-(1 - (1 - exp(-a * T_term))/(a * T_term))
  rhos<-matrix(rep(rho,Z_adaptive),nrow=N,ncol=Z_adaptive)
  dmX<-cbind(X[,1:Z_direct],X[,(Z_direct+1):(Z_adaptive+Z_direct)]*rhos)
  return(dmX)
}
calc_direct_V<-function(phy, sigma2_y, a){
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  Vt<-sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) * exp(-a * tij)) #ta - time from root to tips, tij  - total time separating spcies
  return(Vt)
}

calc_adaptive_V<-function(a, sigma2_y, ta,  tij,  tja,  T_term,  beta,  sigma2_x, Z_adaptive, n_reg){
  N<-dim(ta)[1];
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  T_term<-ts[[3]]
  tja<-ts[[4]]
  ones<-matrix(rep(1,Z_adaptive),nrow=Z_adaptive,ncol=1)
  #beta<-matrix(beta,nrow=n_reg,ncol=1)
  ti<-matrix(T_term,length(T_term),N);
  #if(Z_adaptive>1){
  #  var_opt<-t(beta)%*%sigma2_x%*%ones
  #}
  #else{
  var_opt<-sum(as.matrix(beta^2)%*%sigma2_x)#}
  term0<-(var_opt + sigma2_y) / (2 * a) * (1 - exp( -2 * a * ta)) * exp(-a * tij)
  term1<-(1 - exp(-a * ti)) / (a * ti)
  term2<-exp(-a * tja) * (1 - exp(-a * ti)) / (a * ti)
  Vt<-term0 + var_opt * (ta * term1 * t(term1) - ((1 - exp(-a * ta)) / a) * (term2 + t(term2)))
  return(Vt)
}

calc_multiadaptive_cov_plot<-function(a,sigma2_y,beta,x,Z_adaptive,n_reg){
  #Assuming n_reg>1
  ti<-1
  var_opt = sum(beta^2)
  term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * (1-x/2))) * exp(-a * x)
  term1 = (1 - exp(-a * ti)) / (a * ti)
  term2 = exp(-a * x) * (1 - exp(-a * ti)) / (a * ti)
  Vt<-term0 + var_opt * ((1-x/2) * term1 * term1 - ((1 - exp(-a * (1-x/2))) / a) * (term2 + term2))
  return(Vt)
}
