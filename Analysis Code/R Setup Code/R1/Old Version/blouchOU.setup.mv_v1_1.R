
# Save this file as `R/blouchOU.setup.v1.R`

#' This is the R setup file for blouchOU_v1.stan
#' 
#' @importFrom stats lm
#' @importFrom utils head tail
#' @importFrom dplyr pull
#' @importFrom ape branching.times
#' @param trdata Data formatted by make.treedata function from the R Pacakge treeplyr
#' @param names.traits Vector of trait names
#' @return An object of class list for use in blouchOU_v1.stan
#' @export

blouchOU.setup.mv<-function(trdata,names.traits,names.me.traits=NULL){
  #Trait names are in the order: trait 1 is adapting to trait 2, 3 ...
  #y ~ b1X2+b2X2 ...
  #Data should already be logged
  #names.traits = trait 1, trait 2, trait 3 ...
  #names.me.traits = mv trait 1, mv trait 2, mv trait 3 ...
  n<-length(trdata$phy$tip.label)
  y<-as.data.frame(trdata$dat %>%
    select(any_of(names.traits[1])))
  Xs<-as.data.frame(trdata$dat %>%
    select(any_of(names.traits[-1])))
  if(!is.null(names.me.traits)){
    mv.response<-as.data.frame(trdata$dat %>%select(any_of(names.me.traits[1])))
    mv.predictors<-as.data.frame(trdata$dat %>%select(any_of(names.me.traits[-1])))
  }
  else{
    mv.response<-data.frame(rep(0,n))
    mv.predictors<-matrix(0,nrow=n,ncol=dim(Xs)[2])
  }

  ######################
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
  ########################
  
  #mv.pred<-matrix(0,nrow=n,ncol=)}
  brownian <- mapply(function(x, x_me) sigma.X.estimate(phy, ta, x, x_me),
                     x = split(t(Xs), colnames(Xs)),
                     x_me = split(t(mv.predictors), colnames(Xs)),
                     SIMPLIFY = FALSE)

  
  sigma_squared <- sapply(brownian, function(x) x$sigma_squared)
  brownian_mean <- sapply(brownian, function(x) x$mean)
  
  combo.data<-as.data.frame(trdata$dat %>%
                       select(any_of(names.traits)))
  names(combo.data)[1]<-"y"

  coef<-lm(y ~ ., data=combo.data)$coef
    
  ols.intercept<-coef[1]
  ols.slope<-coef[-1]

  stan_constraint_data<-list(N = n,
                             Z = dim(Xs)[2],
                             Y = y[,1],
                             mv_response = mv.response[,1],
                             random_cov = matrix(0,nrow = n,ncol=dim(Xs)[2]),
                             mv_random_cov = matrix(0,nrow=n,ncol = dim(Xs)[2]),
                             direct_cov = Xs,
                             mv_direct_cov = mv.predictors,
                             ta = as.matrix(ta),
                             T_term = as.vector(T.term),
                             tia = as.matrix(tia),
                             tja = as.matrix(tja),
                             tij = as.matrix(tij),
                             brownian_mean = t(as.matrix(rev(brownian_mean))),
                             sigma_squared_x = t(as.matrix(rev(sigma_squared))),
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
                           Z = dim(Xs)[2],
                           Y = y[,1],
                           mv_response = mv.response[,1],
                           random_cov = Xs,
                           mv_random_cov = mv.predictors,
                           direct_cov =  matrix(0,nrow = n,ncol=dim(Xs)[2]),
                           mv_direct_cov = matrix(0,nrow = n,ncol=dim(Xs)[2]),
                           ta = as.matrix(ta),
                           T_term = as.vector(T.term),
                           tia = as.matrix(tia),
                           tja = as.matrix(tja),
                           tij = as.matrix(tij),
                           brownian_mean = t(as.matrix(rev(brownian_mean))),
                           sigma_squared_x = t(as.matrix(rev(sigma_squared))),
                           ols_intercept = ols.intercept,
                           ols_slope = ols.slope
                           #mean_log = dist.values$mean_log,
                           #sd_log = dist.values$sd_log,
                           #intercept_sd = dist.values$intercept_sd,
                           #slope_sd = dist.values$slope_sd,
                           #sigma2_y_scale = dist.values$sigma2_y_scale
                           )
  
  stan_OU1_data<-list(N = n,
                          Z = dim(Xs)[2],
                          Y = y[,1],
                          mv_response = mv.response[,1],
                          ta = as.matrix(ta),
                          T_term = as.vector(T.term),
                          tia = as.matrix(tia),
                          tja = as.matrix(tja),
                          tij = as.matrix(tij),
                          Y_mean = mean(y[,1])
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


