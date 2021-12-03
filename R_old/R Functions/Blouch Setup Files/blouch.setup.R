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


blouch.setup<-function(trdata){
  col.names<-names(trdata$dat)
  
  #Define tree parameters  
  phy<-trdata$phy
  #max(branching.times(phy))
  phy$edge.length<-phy$edge.length/max(branching.times(phy)) ## rescale tree to height 1
  
  n<-length(phy$tip.label)
  mrca1 <- ape::mrca(phy)
  times <- ape::node.depth.edgelength(phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(phy$tip.label, phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia
  
  test<-sigma.X.estimate(phy,ta, predictor = (pull(trdata$dat,col.names[2])), mv.predictor = matrix(0,nrow=n,ncol=1))
  
  brownian_mean<-test[1]
  sigma_squared_x<-test[2]
  
  coef<-lm(Y~X,data=data.frame(Y=(pull(trdata$dat,col.names[1])),X=(pull(trdata$dat,col.names[2]))))$coef
  ols.intercept<-coef[1]                         
  ols.slope<-coef[2]
  
  
  stan_constraint_data<-list(N = n,
                             Z = 1,
                             Y = as.vector((pull(trdata$dat,col.names[1]))),
                             mv_response = as.vector(rep(0,n)),
                             random_cov = matrix(0,nrow = n,ncol=1),
                             mv_random_cov = matrix(0,nrow=n,ncol = 1),
                             direct_cov = as.matrix((pull(trdata$dat,col.names[2]))), 
                             mv_direct_cov = matrix(0,nrow = n,ncol=1),
                             ta = as.matrix(ta),
                             T_term = as.vector(T.term),
                             tia = as.matrix(tia),
                             tja = as.matrix(tja),
                             tij = as.matrix(tij),
                             brownian_mean = brownian_mean,
                             sigma_squared_x = sigma_squared_x,
                             ols_intercept = ols.intercept,
                             ols_slope = ols.slope)
  
  #############################################################################################
  #Setting up data - random covariates
  stan_adaptive_data<-list(N = n,
                  Z = 1,
                  Y = as.vector((pull(trdata$dat,col.names[1]))),
                  mv_response = as.vector(rep(0,n)),
                  random_cov = as.matrix((pull(trdata$dat,col.names[2]))),
                  mv_random_cov = matrix(0,nrow=n,ncol = 1),
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
                  ols_slope = ols.slope)
  
  
  return(list(stan_constraint_data,stan_adaptive_data))
}


