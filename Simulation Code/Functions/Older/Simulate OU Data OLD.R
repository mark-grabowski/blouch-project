blouch_simXYOUReg.mvSlouch.setup<-function(trdata,num.direct,num.random,num.ME.direct,num.ME.random,slope,hl,vy,vY0,vX0,theta,reg.values=NULL,reg=NULL,mc=NULL,sim.errors=NULL){
  #Added for SBR1 - mv code
  #Modified to allow for correlated predictors - works ith blouchOUReg_v1_5.stan
  #Simulates Y and Xs following direct effect or adaptive model using mvSlouch
  library(stats)
  n<-length(trdata$phy$tip.label)
  tree<-trdata$phy
  num.X<-num.direct+num.random
  names.direct.traits<-NULL
  names.random.traits<-NULL
  M.error<-NULL
  Sxx <- 1
    
  shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
  edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
  if(sim.errors==TRUE){
    M.error<-rlnorm(1+num.random+num.direct,log(0.1),1) #Want to make a N*M list but cannot right now
    mv.response<-rep(M.error[1],n)
    }
  else{
    M.error<-0
  }

  if(num.direct>0){ #Direct effect model
    A<-log(2)/hl
    B<- -A*slope
    Syy<-sqrt(vy*(2*A))
    times <- ape::node.depth.edgelength(tree)
    T.term <- times[1:n]
    reg.matrix<-data.frame(matrix(reg.values,ncol=length(reg.values),nrow=num.direct))
    names(reg.matrix)<-unique(edge.regimes)
    reg.matrix<-as.matrix(reg.matrix)
    OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                           B=matrix(B,ncol=1,nrow=1),mPsi=reg.matrix,
                           Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                           Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),M.error=M.error,
                           starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),
                           Syy=matrix(Syy,ncol=1,nrow=1),
                           B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
      
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=as.vector(edge.regimes)) #Simulate data
  }
    #expected.slope<-slope*(Syy^2+Sxx)*(1-exp(-A*T.term[1]))/(2*Sxx*T.term[1])
  else{ #Adaptive model simulation
    A<-log(2)/hl
    B<- -A*slope
    Syy<-sqrt(vy*(2*A))
    
    times <- ape::node.depth.edgelength(tree)
    T.term <- times[1:n]
    #Evolutionary sloperession based on optimal slope ression and A  
    rho<- (1-(1 - exp(-A * T.term))/(A * T.term))
    evol.slope<-(rho*slope)[1]
    
    reg.matrix<-data.frame(matrix(reg.values,ncol=length(reg.values),nrow=num.random))
    names(reg.matrix)<-unique(edge.regimes)
    reg.matrix<-as.matrix(reg.matrix)
      
    OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                           B=matrix(B,ncol=1,nrow=1),mPsi=reg.matrix,
                           Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                           Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),#M.error=M.error,
                           starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
      
    
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=as.vector(edge.regimes)) #Simulate data
  }
  response<-OUBMdata[,1]

  if(num.direct>0){ #Only works for all direct or all random traits
    names.direct.traits<-names(OUBMdata)[2:(num.direct+1)]
  }
  if(num.random>0){
    names.random.traits<-names(OUBMdata)[2:(num.random+1)]
  }
  if(num.direct==0){
    direct<-matrix(0,nrow=n,ncol=0)}  
  if(num.random==0){
    random<-matrix(0,nrow=n,ncol=0)}  
  
  if(num.direct!=0){
    direct<-data.frame(OUBMdata[,2:(1+num.direct)])
    if(mc==TRUE){
      direct<-direct-colMeans(direct)}
  }
  
  if(num.random!=0){
    random<-data.frame(OUBMdata[,2:(1+num.random)])
    if(mc==TRUE){
      random<-random-colMeans(random)}
  }
  
  if(sim.errors==TRUE){
    mv.response<-rep(M.error[1],n)
    if(num.ME.direct!=0){
      mv.direct<-matrix(M.error[2:num.ME.direct],nrow=n,ncol=num.ME.direct)
      mv.random<-matrix(0,nrow=n,ncol=0)
      }
    if(num.ME.random!=0){
      mv.random<-matrix(M.error[2:num.ME.random],nrow=n,ncol=num.ME.random)
      mv.direct<-matrix(0,nrow=n,ncol=0)}
  }else{
    mv.response<-rep(0,n)
    mv.random<-matrix(0,nrow=n,ncol=num.random)
    mv.direct<-matrix(0,nrow=n,ncol=0)
  }

  
  #tree$edge.length<-tree$edge.length/max(branching.times(tree)) ## rescale tree to height 1
  ##############################
  n<-length(tree$tip.label)
  mrca1 <- ape::mrca(tree)
  times <- ape::node.depth.edgelength(tree)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(tree$tip.label, tree$tip.label))
  
  #ta[ta==0]<-0.000000000001 #Added because dividing by ta in newest Blouch
  
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia
  
  ########################
  #Added for SBR1 - mv code
  
  #mv.pred<-matrix(0,nrow=n,ncol=)}
  if(num.random!=0){
    if(num.random==1){
      brownian <- list()
      X <- random[,1]
      X_me <- mv.random[,1] 
      brownian[[1]] <-  sigma.X.estimate(tree, ta, X, X_me)
      names(brownian) <- colnames(random)
      sigma_squared <- sapply(brownian, function(x) x$sigma_squared)
      brownian_mean <- sapply(brownian, function(x) x$mean)}
    if(num.random>1){
      rownames(random)<-tree$tip.label
      brownian<-sigma.X.estimate.mv(tree,X)
      sigma_squared<-matrix(brownian[[2]],dim(X)[2],dim(X)[2])
      brownian_mean<-brownian[[1]]
    }
  }else{
    sigma_squared <- matrix(0,0,0)
    brownian_mean <- matrix(0,0,1)
  }
  
  #combo.data<-as.data.frame(trdata$dat %>%
  #                      select(any_of(c(name.response.trait,
  #                                      names.direct.traits,names.random.traits))))
  combo.data<-cbind(response,direct,random)
  names(combo.data)[1]<-"y"
  
  coef<-lm(y ~ ., data=combo.data)$coef
  ols.intercept<-coef[1]
  ols.slope<-as.array(coef[-1])
  

  ############################################################################
  #Fixed regimes
  regimes_internal <-trdata$phy$node.label
  regimes_tip <- trdata$dat$regimes
  regimes <- concat.factor(regimes_tip, regimes_internal)
  lineages <- lapply(1:n, function(e) lineage.constructor(tree, e, regimes)) #; names(lineages) <- tree$tip.label
  
  store<-NULL
  
  for(i in 1:length(lineages)){
    store<-c(store,length(lineage.nodes(tree,i)))
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
  
  if(num.ME.random==0){
    mv.random<-matrix(0,nrow=n,ncol=0)
  }
  ############################################################################################
  #Setting up data
  stan_data<-list(N = n,
                  Z = num.direct+num.random, #Total predictor number
                  Z_direct = num.direct, #Total direct trait number
                  Z_random = num.random, #Total random trait number
                  Z_me_direct = num.ME.direct,
                  Z_me_random = num.ME.random,
                  Y = response,
                  mv_response = mv.response,
                  random_cov = random,
                  mv_random_cov = mv.random,
                  direct_cov = as.matrix(direct,n,num.direct),
                  mv_direct_cov = mv.direct,
                  ta = as.matrix(ta),
                  T_term = as.vector(T.term),
                  tia = as.matrix(tia),
                  tja = as.matrix(tja),
                  tij = as.matrix(tij),
                  brownian_mean = t(as.matrix(brownian_mean)),
                  sigma_squared_x = as.matrix(sigma_squared),
                  ols_intercept = ols.intercept,
                  ols_slope = ols.slope,
                  
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
                  regimes_matrix = data.matrix(regimes_matrix) #Converts to integers - 0's are now the ends of branches)
)
  
  

  return(stan_data)
}


sigma.X.estimate <-
  function (tree, ta, predictor, mv.predictor) {
    predictor <- matrix(predictor, nrow = length(tree$tip.label))
    N <- length(tree$tip.label)
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

sigma.X.estimate.mv<-function(tree,predictors){
  library(geiger)
  return(list(mean = rep(0,length(predictors)),sigma_squared = ratematrix(tree,predictors)))
}

#############################################################################################  
#Data formatting drawn from Slouch
parent <- function(tree, x){
  m <- which(tree$edge[, 2] == x)
  return(tree$edge[m, 1])
}

lineage.nodes <- function(tree, x){
  k <- x
  N <- length(tree$tip.label)
  while(x != N + 1){
    k <- c(k, parent(tree, x))
    x <- tail(k, n = 1)
  }
  return(k)
}

#lineages <- lapply(1:n, function(e) lineage.constructor(tree, e, regimes)) #; names(lineages) <- tree$tip.label

###########################################################################################################################
lineage.constructor <- function(tree, e, regimes){
  nodes <- lineage.nodes(tree, e)
  which.regimes = lapply(levels(regimes), function(x) {res <- match(regimes[nodes], x); res[is.na(res)] <- 0; return(res) })
  names(which.regimes) <- levels(regimes)
  
  nodes_time <-  ape::node.depth.edgelength(tree)[nodes]
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

weight.matrix <- function(tree, a, lineages){
  if(a > 300000000000) a <- 300000000000
  res <- t(vapply(lineages, function(x) weights_regimes(a, x), 
                  FUN.VALUE = numeric(length(lineages[[1]]$which.regimes))) ## Specify type of output
  )
  
  rownames(res) <- tree$tip.label
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



blouch_simXYOU.mvSlouch.setup<-function(trdata,num.direct,num.random,num.ME.direct,num.ME.random,slope,hl,vy,vY0,vX0,theta,mc=NULL,sim.errors=NULL){
  #Added for SBR1 - mv code
  #Modified to allow for correlated predictors - works ith blouchOUReg_v1_5.stan
  #Simulates Y and Xs following direct effect or adaptive model using mvSlouch
  library(stats)
  n<-length(trdata$phy$tip.label)
  tree<-trdata$phy
  num.X<-num.direct+num.random
  names.direct.traits<-NULL
  names.random.traits<-NULL
  M.error<-NULL
  Sxx <- 1
  
  if(sim.errors==TRUE){
    M.error<-rlnorm(1+num.random+num.direct,log(0.1),1) #Want to make a N*M list but cannot right now
    mv.response<-rep(M.error[1],n)
  }
  else{
    M.error<-0
  }
  
  if(num.direct>0){ #Direct effect model
    A<-log(2)/hl
    B<- -A*slope
    Syy<-sqrt(vy*(2*A))
    times <- ape::node.depth.edgelength(tree)
    T.term <- times[1:n]
    OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                         B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                         Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                         Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),M.error=M.error,
                         starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
    
    
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
  }
  else{ #Adaptive model simulation
    A<-log(2)/hl
    B<- -A*slope
    Syy<-sqrt(vy*(2*A))
    
    times <- ape::node.depth.edgelength(tree)
    T.term <- times[1:n]              
    
    #Evolutionary sloperession based on optimal sloperession and A  
    #rho<- (1-(1 - exp(-A * T.term))/(A * T.term))
    #evol.slope<-(rho*opt.slope)[1]
    
    OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                         B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                         Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                         Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
    
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    }
  
  response<-OUBMdata[,1]
  
  if(num.direct>0){ #Only works for all direct or all random traits
    names.direct.traits<-names(OUBMdata)[2:(num.direct+1)]
  }
  if(num.random>0){
    names.random.traits<-names(OUBMdata)[2:(num.random+1)]
  }
  
  if(num.direct==0){
    direct<-matrix(0,nrow=n,ncol=0)}  
  if(num.random==0){
    random<-matrix(0,nrow=n,ncol=0)}  
  
  if(num.direct!=0){
    direct<-data.frame(OUBMdata[,2:(1+num.direct)])
    if(mc==TRUE){
      direct<-direct-colMeans(direct)}
  }
  
  if(num.random!=0){
    random<-data.frame(OUBMdata[,2:(1+num.random)])
    if(mc==TRUE){
      random<-random-colMeans(random)}
  }
  
  if(sim.errors==TRUE){
    mv.response<-rep(M.error[1],n)
    if(num.ME.direct!=0){
      mv.direct<-matrix(M.error[2:num.ME.direct],nrow=n,ncol=num.ME.direct)
      mv.random<-matrix(0,nrow=n,ncol=0)
    }
    if(num.ME.random!=0){
      mv.random<-matrix(M.error[2:num.ME.random],nrow=n,ncol=num.ME.random)
      mv.direct<-matrix(0,nrow=n,ncol=0)}
  }else{
    mv.response<-rep(0,n)
    mv.random<-matrix(0,nrow=n,ncol=num.random)
    mv.direct<-matrix(0,nrow=n,ncol=0)
  }
  n<-length(tree$tip.label)
  mrca1 <- ape::mrca(tree)
  times <- ape::node.depth.edgelength(tree)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(tree$tip.label, tree$tip.label))
  
  #ta[ta==0]<-0.000000000001 #Added because dividing by ta in newest Blouch
  
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia
  
  if(num.random!=0){
    if(num.random==1){
      brownian <- list()
      X <- random[,1]
      X_me <- mv.random[,1] 
      brownian[[1]] <-  sigma.X.estimate(tree, ta, X, X_me)
      names(brownian) <- colnames(random)
      sigma_squared <- sapply(brownian, function(x) x$sigma_squared)
      brownian_mean <- sapply(brownian, function(x) x$mean)}
    if(num.random>1){
      rownames(random)<-tree$tip.label
      brownian<-sigma.X.estimate.mv(tree,X)
      sigma_squared<-matrix(brownian[[2]],dim(X)[2],dim(X)[2])
      brownian_mean<-brownian[[1]]
    }
  }else{
    sigma_squared <- matrix(0,0,0)
    brownian_mean <- matrix(0,0,1)
  }
  
  combo.data<-cbind(response,direct,random)
  names(combo.data)[1]<-"y"
  
  coef<-lm(y ~ ., data=combo.data)$coef
  ols.intercept<-coef[1]
  ols.slope<-as.array(coef[-1])
  
  
  if(num.ME.random==0){
    mv.random<-matrix(0,nrow=n,ncol=0)
  }
  ############################################################################################
  #Setting up data
  stan_data<-list(N = n,
                  Z = num.direct+num.random, #Total predictor number
                  Z_direct = num.direct, #Total direct trait number
                  Z_random = num.random, #Total random trait number
                  Z_me_direct = num.ME.direct,
                  Z_me_random = num.ME.random,
                  Y = response,
                  mv_response = mv.response,
                  random_cov = random,
                  mv_random_cov = mv.random,
                  direct_cov = as.matrix(direct,n,num.direct),
                  mv_direct_cov = mv.direct,
                  ta = as.matrix(ta),
                  T_term = as.vector(T.term),
                  tia = as.matrix(tia),
                  tja = as.matrix(tja),
                  tij = as.matrix(tij),
                  brownian_mean = t(as.matrix(brownian_mean)),
                  sigma_squared_x = as.matrix(sigma_squared),
                  ols_intercept = ols.intercept,
                  ols_slope = ols.slope)
  
  
  
  return(stan_data)
}

