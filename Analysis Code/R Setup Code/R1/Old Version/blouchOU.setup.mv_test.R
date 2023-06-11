
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

blouchOU.setup.mv<-function(trdata,name.response.trait,names.direct.traits=NULL,names.random.traits=NULL,names.response.me.trait=NULL,names.direct.me.traits=NULL,names.random.me.traits=NULL){
  
  #Trait names are in the order: trait 1 is adapting to trait 2, 3 ...
  #y ~ b1X2+b2X2 ...
  #Data should already be logged
  #Added for SBR1 - mv code
  
  n<-length(trdata$phy$tip.label)
  if(is.null(names.direct.traits)){
    direct<-matrix(0,nrow=n,ncol=0)}  
  if(is.null(names.random.traits)){
    random<-matrix(0,nrow=n,ncol=0)}  
  
  response<-as.data.frame(trdata$dat %>%
                            select(any_of(name.response.trait)))
  if(!is.null(names.direct.traits)){
    direct<-as.data.frame(trdata$dat %>%
                            select(any_of(names.direct.traits)))
  }
  
  if(!is.null(names.random.traits)){
    random<-as.data.frame(trdata$dat %>%
                            select(any_of(names.random.traits)))
  }
  
  if(!is.null(names.response.me.trait)){
    mv.response<-as.data.frame(trdata$dat %>%select(any_of(names.response.me.trait)))
  }else{
    mv.response<-data.frame(rep(0,n))
  }
  if(!is.null(names.direct.me.traits)){
    mv.direct<-as.data.frame(trdata$dat %>%select(any_of(names.direct.me.traits)))
  }else{
    mv.direct<-matrix(0,nrow=n,ncol=0)
  }
  
  if(!is.null(names.random.me.traits)){
    mv.random<-as.data.frame(trdata$dat %>%select(any_of(names.random.me.traits)))
    N.me<-n
  }else{
    mv.random<-matrix(0,nrow=n,ncol=length(names.random.traits))
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
  if(!is.null(names.random.traits)){
    brownian <- mapply(function(x, x_me) sigma.X.estimate(phy, ta, x, x_me),
                       x = split(t(random), colnames(random)),
                       x_me = split(t(mv.random), colnames(random)),
                       SIMPLIFY = FALSE)
    
    
    sigma_squared <- sapply(brownian, function(x) x$sigma_squared)
    brownian_mean <- sapply(brownian, function(x) x$mean)
  }
  else{
    sigma_squared <- matrix(0,length(names.random.traits),1)
    brownian_mean <- matrix(0,length(names.random.traits),1)
    
  }
  
  
  
  combo.data<-as.data.frame(trdata$dat %>%
                              select(any_of(c(name.response.trait,
                                              names.direct.traits,names.random.traits))))
  names(combo.data)[1]<-"y"
  
  coef<-lm(y ~ ., data=combo.data)$coef
  ols.intercept<-coef[1]
  ols.slope<-as.array(coef[-1]) #Hint: https://groups.google.com/g/stan-users/c/uqxC0Aeg2YY
  
  #if(is.null(names.random.me.traits)){
  #  mv.response<-data.frame(rep(0,n))
  #  }
  
  if(is.null(names.random.me.traits)){
    mv.random<-matrix(0,nrow=n,ncol=0)
  }
  
  #######################################
  stan_data<-list(N = n,
                  Z = length(names.direct.traits)+length(names.random.traits), #Total trait number
                  Z_direct = length(names.direct.traits), #Total direct trait number
                  Z_random = length(names.random.traits), #Total random trait number
                  Z_me_direct = length(names.direct.me.traits),
                  Z_me_random = length(names.random.me.traits),
                  Y = response[,1],
                  mv_response = mv.response[,1],
                  random_cov = random,
                  mv_random_cov = mv.random,
                  direct_cov = direct,
                  mv_direct_cov = mv.direct,
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
                      Z = length(names.direct.traits)+length(names.random.traits),
                      Y = response[,1],
                      mv_response = as.vector(mv.response),
                      ta = as.matrix(ta),
                      T_term = as.vector(T.term),
                      tia = as.matrix(tia),
                      tja = as.matrix(tja),
                      tij = as.matrix(tij),
                      Y_mean = mean(response[,1])
  )
  
  return(list(stan_data,stan_OU1_data))
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


