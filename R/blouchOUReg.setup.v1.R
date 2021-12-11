
# Save this file as `R/blouchOUReg.setup.v1.R`

#' This is the R setup file for blouchOUReg_v1.stan
#' It includes the function sigma.X.estimate, taken from the R Pacakge Slouch (Kopperud et al. 2020).
#'
#' @param trdata Data formatted by make.treedata function from the R Pacakge treeplyr
#' @param names.traits Vector of trait names
#' @return An object of class list for use in blouchOUReg_v1.stan
#' @export

blouchOUReg.setup.v1<-function(trdata,names.traits){
  #Trait names are in the order trait 1 is adapting to trait 2, data should already be logged
  #Names.traits = categorical trait, trait 1, trait 2, mv trait 1, mv trait 2
  
  #Define tree parameters  
  phy<-trdata$phy
  #max(branching.times(phy))
  phy$edge.length<-phy$edge.length/max(branching.times(phy)) ## rescale tree to height 1
  
#############################################################################################  
  #Data formatting drawn from Slouch
  parent <- function(phy, x){
    m <- which(phy$edge[, 2] == x)
    return(phy$edge[m, 1])
  }
  
  lineage.nodes <- function(phy, x){
    k <- x
    N <- length(phy$tip.label)
    while(x != N + 1){
      k <- c(k, parent(phy, x))
      x <- tail(k, n = 1)
    }
    return(k)
  }
  
  #lineages <- lapply(1:n, function(e) lineage.constructor(phy, e, regimes)) #; names(lineages) <- phy$tip.label
  
  
  ###########################################################################################################################
  lineage.constructor <- function(phy, e, regimes){
    nodes <- lineage.nodes(phy, e)
    which.regimes = lapply(levels(regimes), function(x) {res <- match(regimes[nodes], x); res[is.na(res)] <- 0; return(res) })
    names(which.regimes) <- levels(regimes)
    
    nodes_time <-  ape::node.depth.edgelength(phy)[nodes]
    timeflip <- nodes_time[1] - nodes_time ## Time from tip to node(s)
    t_end <- tail(timeflip, n = -1) ## Time from tip to end of segment(s)
    t_beginning <- head(timeflip, n = -1) ## Time from tip to beginning of segment(s)
    regime_time <- c(t_end - t_beginning, 0)
    
    return(list(nodes = nodes, 
                nodes_time = nodes_time,
                t_end = t_end,
                t_beginning = t_beginning,
                regime_time = regime_time,
                which.regimes = which.regimes))
  }
  
  ########################################################################
  #Eq. 3 from Hansen, 1997
  weights_segments <- function(a, lineage){
    res <- c(exp(-a * lineage$t_beginning) - exp(-a * lineage$t_end), 
             exp(-a * lineage$nodes_time[1]))
    return(res)
  }
  #######################################################
  weights_regimes <- function(a, lineage) {#Take in individual lineage
    #nt <- lineage$nodes_time
    res <- weights_segments(a, lineage) #Calculate individual lineage weighting - res is how long lineage spent in each regime *-a
    w <- vapply(lineage$which.regimes, function(e) sum(e*res), FUN.VALUE = 0) ## Sum coefficients for each of the regimes in the which.regimes grouping.
    return(w)
  }
  ###################################################################
  weights_regimes_brown <- function(lineage){
    w <- sapply(lineage$which.regimes, function(e) sum(e*lineage$regime_time))
    return(w)
  }
  
  weight.matrix <- function(phy, a, lineages){
    if(a > 300000000000) a <- 300000000000
    res <- t(vapply(lineages, function(x) weights_regimes(a, x), 
                    FUN.VALUE = numeric(length(lineages[[1]]$which.regimes))) ## Specify type of output
    )
    
    rownames(res) <- phy$tip.label
    return(res)
  }
  #############################################################
  weight.matrix.brown <- function(lineages){
    res <- t(sapply(lineages, weights_regimes_brown))
    return(res)
  }
  
  ## Thanks to user "snaut" at stackoverflow, http://stackoverflow.com/users/1999873/snaut
  ################################################################
  concat.factor <- function(...){
    as.factor(do.call(c, lapply(list(...), as.character)))
  }
  ###################################################################
  
  
  X<-(pull(trdata$dat,names.traits[3]))
  Y<-(pull(trdata$dat,names.traits[2]))
  n<-length(phy$tip.label)
  mrca1 <- ape::mrca(phy)
  times <- ape::node.depth.edgelength(phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(phy$tip.label, phy$tip.label))
  
  #ta[ta==0]<-0.000000000001 #Added because dividing by ta in newest Blouch
  
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia
  
  if(is.na(names.traits[4])!=TRUE){
    mv.response<-as.vector(pull(trdata$dat,names.traits[4]))}
  else{
    mv.response<-as.vector(rep(0,n))}
  
  if(is.na(names.traits[5])!=TRUE){
    mv.pred<-matrix(pull(trdata$dat,names.traits[5]),nrow=n,ncol=1)}
  else{
    mv.pred<-matrix(0,nrow=n,ncol=1)}
  
  test<-sigma.X.estimate(phy,ta, X, matrix(0,nrow=n,ncol = 1))
  
  brownian_mean<-test[1]
  sigma_squared_x<-test[2]
  
  coef<-lm(Y~((pull(trdata$dat,names.traits[1]))+X))$coef
  #print(coef)
  
  ############################################################################
  #Fixed regimes
  regimes_internal <- phy$node.label
  fixed.fact <- (pull(trdata$dat,names.traits[1]))
  regimes_tip <- fixed.fact
  regimes <- concat.factor(regimes_tip, regimes_internal)
  #regimes <- factor(rep("OU1",length(regimes)))# For OU1
  lineages <- lapply(1:n, function(e) lineage.constructor(phy, e, regimes)) #; names(lineages) <- phy$tip.label
  
  store<-NULL
  
  for(i in 1:length(lineages)){
    store<-c(store,length(lineage.nodes(phy,i)))
  }
  max.node.length<-max(store)  
  
  n_regimes<-length(unique(regimes))
  nodes<-matrix(0,length(lineages),max.node.length)
  nodes_time<-matrix(0,length(lineages),max.node.length)
  t_end<-matrix(0,length(lineages),max.node.length)
  t_beginning<-matrix(0,length(lineages),max.node.length)
  regime_time<-matrix(0,length(lineages),max.node.length)
  regimes_matrix<-data.frame(matrix(0,length(lineages),max.node.length))
  
  
  
  for(i in 1:length(lineages)){
    nodes[i,1:length(lineages[[i]]$nodes)]<-lineages[[i]]$nodes
    nodes_time[i,1:length(lineages[[i]]$nodes_time)]<-lineages[[i]]$nodes_time
    t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
    t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
    regime_time[i,1:length(lineages[[i]]$regime_time)]<-lineages[[i]]$regime_time
    regimes_matrix[i,1:length(regimes[lineages[[i]]$nodes])]<-as.numeric(regimes[lineages[[i]]$nodes])
  }
  
  ols.intercept<-coef[1]                         
  #ols.reg<-coef[2:n_regimes]                         
  if(n_regimes<length(coef)){
    ols.slope<-coef[n_regimes+1]
  }
  
  
  ##################################################
  #############################################################################################
  #Setting up data
  stan_constraint_data<-list(N = n,
                  Z = 1,
                  Y = Y,
                  mv_response = mv.response,
                  random_cov = matrix(0,nrow=n,ncol = 1),
                  mv_random_cov = matrix(0,nrow=n,ncol = 1),
                  direct_cov = as.matrix(X),
                  mv_direct_cov = mv.pred,
                  ta = as.matrix(ta),
                  T_term = as.vector(T.term),
                  tia = as.matrix(tia),
                  tja = as.matrix(tja),
                  tij = as.matrix(tij),
                  brownian_mean = brownian_mean,
                  sigma_squared_x = sigma_squared_x,
                  ols_intercept = ols.intercept,
                  ols_slope = ols.slope,
                  #ols_reg = ols.reg,
                  
                  ###########################################################################################
                  #Fixed regimes
                  n_regimes = n_regimes,
                  n_lineages = length(lineages),
                  max_node_length = max.node.length,
                  nodes = as.matrix(nodes),
                  nodes_time = as.matrix(nodes_time),
                  t_end = as.matrix(t_end),
                  t_beginning = as.matrix(t_beginning),
                  regime_time = as.matrix(regime_time),
                  regimes_matrix = data.matrix(regimes_matrix) #Converts to integers - 0's are now the ends of branches
  )
  
  #############################################################################################
  #Setting up data - random covariates
  stan_adaptive_data<-list(N = n,
                             Z = 1,
                             Y = Y,
                             mv_response = mv.response,
                             random_cov = as.matrix(X),
                             mv_random_cov = mv.pred,
                             direct_cov = matrix(0,nrow=n,ncol = 1),
                             mv_direct_cov = matrix(0,nrow=n,ncol = 1),
                             ta = as.matrix(ta),
                             T_term = as.vector(T.term),
                             tia = as.matrix(tia),
                             tja = as.matrix(tja),
                             tij = as.matrix(tij),
                             brownian_mean = brownian_mean,
                             sigma_squared_x = sigma_squared_x,
                             ols_intercept = ols.intercept,
                             ols_slope = ols.slope,
                             #ols_reg = ols.reg,
                           
                             
                             ###########################################################################################
                             #Fixed regimes
                             n_regimes = n_regimes,
                             n_lineages = length(lineages),
                             max_node_length = max.node.length,
                             nodes = as.matrix(nodes),
                             nodes_time = as.matrix(nodes_time),
                             t_end = as.matrix(t_end),
                             t_beginning = as.matrix(t_beginning),
                             regime_time = as.matrix(regime_time),
                             regimes_matrix = data.matrix(regimes_matrix) #Converts to integers - 0's are now the ends of branches
  )
  
  
  return(list(stan_constraint_data,stan_adaptive_data))
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


