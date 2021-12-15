
# Save this file as `R/blouchOU.setup.v1.R`

#' This is the R setup file for blouchOU_v1.stan
#' 
#'
#' @param trdata Data formatted by make.treedata function from the R Pacakge treeplyr
#' @param names.traits Vector of trait names
#' @return An object of class list for use in blouchOU_v1.stan
#' @export

blouchOU.setup.v1<-function(trdata,names.traits){
  #blouchOU.setup.v1<-function(trdata,names.traits,dist.values){
  #Trait names are in the order trait 1 is adapting to trait 2, data should already be logged
  #Names.traits = trait 1, trait 2, mv trait 1, mv trait 2
  #Define tree parameters
  phy<-trdata$phy
  #max(branching.times(phy))
  phy$edge.length<-phy$edge.length/max(branching.times(phy)) ## rescale tree to height 1

  n<-length(phy$tip.label)
  mrca1 <- ape::mrca(phy)
  times <- ape::node.depth.edgelength(phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(phy$tip.label, phy$tip.label))
  #ta[ta==0]<-.Machine$double.eps #Added because dividing by ta in newest Blouch
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  if(is.na(names.traits[3])!=TRUE){
    mv.response<-as.vector(pull(trdata$dat,names.traits[3]))}
  else{
    mv.response<-as.vector(rep(0,n))}

  if(is.na(names.traits[4])!=TRUE){
    mv.pred<-matrix(pull(trdata$dat,names.traits[4]),nrow=n,ncol=1)}
  else{
    mv.pred<-matrix(0,nrow=n,ncol=1)}

  test<-sigma.X.estimate(phy,ta, predictor = (pull(trdata$dat,names.traits[2])), mv.predictor = mv.pred)

  brownian_mean<-test[1]
  sigma_squared_x<-test[2]

  coef<-lm(Y~X,data=data.frame(Y=(pull(trdata$dat,names.traits[1])),X=(pull(trdata$dat,names.traits[2]))))$coef
  ols.intercept<-coef[1]
  ols.slope<-coef[2]


  stan_constraint_data<-list(N = n,
                             Z = 1,
                             Y = as.vector((pull(trdata$dat,names.traits[1]))),
                             mv_response = mv.response,
                             random_cov = matrix(0,nrow = n,ncol=1),
                             mv_random_cov = matrix(0,nrow=n,ncol = 1),
                             direct_cov = as.matrix((pull(trdata$dat,names.traits[2]))),
                             mv_direct_cov = mv.pred,
                             ta = as.matrix(ta),
                             T_term = as.vector(T.term),
                             tia = as.matrix(tia),
                             tja = as.matrix(tja),
                             tij = as.matrix(tij),
                             brownian_mean = brownian_mean,
                             sigma_squared_x = sigma_squared_x,
                             ols_intercept = ols.intercept,
                             ols_slope = ols.slope
                             #mean_log = dist.values$mean_log,
                             #sd_log = dist.values$sd_log,
                             #intercept_sd = dist.values$intercept_sd,
                             #slope_sd = dist.values$slope_sd,
                             #sigma2_y_scale = dist.values$sigma2_y_scale
                             )

  #############################################################################################
  #Setting up data - random covariates
  stan_adaptive_data<-list(N = n,
                           Z = 1,
                           Y = as.vector((pull(trdata$dat,names.traits[1]))),
                           mv_response = mv.response,
                           random_cov = as.matrix((pull(trdata$dat,names.traits[2]))),
                           mv_random_cov = mv.pred,
                           direct_cov =  matrix(0,nrow = n,ncol=1),
                           mv_direct_cov = matrix(0,nrow = n,ncol=1),
                           ta = as.matrix(ta),
                           T_term = as.vector(T.term),
                           tia = as.matrix(tia),
                           tja = as.matrix(tja),
                           tij = as.matrix(tij),
                           brownian_mean = brownian_mean,
                           sigma_squared_x = sigma_squared_x,
                           ols_intercept = ols.intercept,
                           ols_slope = ols.slope
                           #mean_log = dist.values$mean_log,
                           #sd_log = dist.values$sd_log,
                           #intercept_sd = dist.values$intercept_sd,
                           #slope_sd = dist.values$slope_sd,
                           #sigma2_y_scale = dist.values$sigma2_y_scale
                           )

  stan_OU1_data<-list(N = n,
                          Z = 1,
                          Y = as.vector((pull(trdata$dat,names.traits[1]))),
                          mv_response = mv.response,
                          ta = as.matrix(ta),
                          T_term = as.vector(T.term),
                          tia = as.matrix(tia),
                          tja = as.matrix(tja),
                          tij = as.matrix(tij),
                          Y_mean = mean(as.vector((pull(trdata$dat,names.traits[1]))))
                      )
  
  return(list(stan_constraint_data,stan_adaptive_data,stan_OU1_data))
}


sigma.X.estimate <-
  function (phy, ta, predictor, mv.predictor) {
    predictor <- matrix(predictor, nrow = length(phy$tip.label))
    
    N <- length(phy$tip.label)
    v <- ta # Time from root to most recent ancestor
    w <- matrix(data = 1, nrow = N, ncol = 1)
    me <- diag(mv.predictor)
    dat <- predictor
    beta1 <- solve(t(w) %*% solve(v) %*% w) %*% (t(w) %*% solve(v) %*% dat)
    e <- dat - c(beta1)
    sigma_squared <- as.numeric((t(e) %*% solve(v) %*% e) / (N-1))
    repeat{
      beta1 <- solve(t(w) %*% solve(v + me/sigma_squared) %*% w) %*% (t(w) %*% solve(v + me/sigma_squared) %*% dat)
      e <- dat - c(beta1)
      sigma_squared1 <- (t(e) %*% solve(v + me/sigma_squared) %*% e) / (N-1)
      if (abs(as.numeric(sigma_squared1) - sigma_squared) <= 0.0000001 * sigma_squared){
        break
      }
      sigma_squared <- as.numeric(sigma_squared1)
    }
    return(list(mean = as.numeric(beta1),
                sigma_squared = as.numeric(sigma_squared)))
  }


