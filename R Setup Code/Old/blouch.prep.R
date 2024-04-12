#' Functions to setup data for Blouch
#'#Note measurement error should be entered as standard deviation, not variance
#' @export
#' @param trdata An object of the class treedata from function treeplyr
#' @param Y  Vector containing name of column in treedata containing response variable
#' @param Y_error  Vector containing name of column in treedata containing error of response variable
#' @param X Vector containing name(s) of column in treedata containing predictor variable(s)
#' @param X_error Vector containing name(s) of column in treedata containing error of predictor variable(s)
#' @param Z_direct Vector containing number of direct effect predictor traits
#' @param Z_adaptive Vector containing number of adaptive predictor traits
#' @param reg_column Vector containing name of regime column in treedata
#' @return A list object ready for Blouch's Stan functions
#'



ts_fxn<-function(phy){ #Calculate t
  n<-length(phy$tip.label)
  mrca1 <- ape::mrca(phy) #Node numbers for MRCA of tips
  times <- ape::node.depth.edgelength(phy) #Time from root to tips of each node, starting with the tips
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(phy$tip.label, phy$tip.label)) #Matrix with time from root to MRCA of pairs of tips - pulls out values of times that correspond with node numbers - integers
  T.term <- times[1:n] #Times from root for tips
  tia <- times[1:n] - ta #Times from root to tips - times from root to MRCA = times from MRCA to tips
  tja <- t(tia) #Transpose of the times from MRCA to tips matrix
  #tij <- tja + tia #Sum of matrix and its transpose - total time separating species
  tij<-cophenetic(phy)
  #return(list(ta,tia,tja,tij,T.term))
  return(list(ta,tij,T.term,tja))
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


#############################################################################################
#Data formatting drawn from Slouch

parent <- function(phy, x){ #Returns parent node of offspring node given node number
  m <- which(phy$edge[, 2] == x)
  return(phy$edge[m, 1])
}

lineage.nodes <- function(phy, x){ #Given a certain node, return the list of all parent nodes back to the root of the tree
  k <- x #Save x in k
  N <- length(phy$tip.label) #total number of tips on tree
  while(x != N + 1){ #while node is not equal to number of tips +1  - starting node - 51 here
    k <- c(k, parent(phy, x)) #Return node at beginning of edge
    x <- tail(k, n = 1) #x is assigned value at end of k, so end of the list of beginning nodes
    #50->99->89->51 0- tracing lineage back by nodes
  }
  return(k)
}

lineage.constructor <- function(phy, e, anc_maps="regimes", regimes, ace){ #Revised 2022 Slouch version
  #e = 50 - tip
  #regimes[50]<-"OU2"
  nodes <- lineage.nodes(phy, e) #Given a certain node, return the list of all parent nodes back to the root of the tree
  #[1] 50 99 87 51
  min_age <- min(ape::node.depth.edgelength(phy)[nodes]) #min root to node time
  #[1] 1.0000000 0.4794736 0.1406307 0.0000000
  if(anc_maps == "regimes"){
    lineage_regimes <- rev(regimes[nodes]) #Reverse order of regimes from that in nodes object
    #[1] OU1 OU1 OU1 OU2
    #Levels: OU1 OU2
    which.regimes <- lapply(levels(regimes), function(x) {res <- match(regimes[nodes], x); res[is.na(res)] <- 0; return(res) })
    #Determine which regimes each node is in
    #x is the levels of the regimes, match takes the regimes at the nodes in the lineage and determines whether which level of regime the node belongs to
    #Any NAs get 0 - this happens when regimes are not assigned on a lineage
    #[[1]]
    #[1] 0 1 1 1 - Tip is OU2, so gets 0 for OU1 here but 1 for OU2 below - reverse order
    #[[2]]
    #[1] 1 0 0 0
    times <-  ape::node.depth.edgelength(phy)[nodes] #Root to node time
    #[1] 1.0000000 0.4794736 0.1406307 0.0000000
    timeflip <- times[1] - times ## Time from tips to node(s)
    #[1] 0.0000000 0.5205264 0.8593693 1.0000000
  }else if(anc_maps == "simmap"){
    ## Simmap splits up each edge into sub-edges, depending on the split. So, we use edges instead of nodes, and introduce sub-edges
    edge_is <- which(phy$edge[,2] %in% nodes) #indices of edges that end with the nodes of lineage
    #[1] 50 99 87 51 - nodes
    #[1] 72 96 98 - edges - 98 has 99/50, 96 has 87/99, 72 has 51/87
    subedges <- unlist(lapply(edge_is, function(i) phy$maps[[i]]))
    #maps = a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.
    #e.g. [[3]]
    #OU1        OU2
    #0.04604509 0.07161615
    #matches indices of edges with the regimes for each edge, and returns named vector with length of time spent in each regime per edge
    # OU1       OU1       OU1
    #0.1406307 0.3388429 0.5205264
    simmap_regimes <- rev(names(subedges))
    #Saves the reversed name of the regimes from the order in subedges
    which.regimes <- lapply(levels(regimes), function(x) {res <- match(simmap_regimes, x); res[is.na(res)] <- 0; return(res)})
    #Returns which regimes are at each edge
    #[[1]] Standard regime scoring
    #[1] 1 1 1

    #[[2]]
    #[1] 0 0 0]
    # Problem. simmap does not provide root estimate. Assuming root estimate is equal to the oldest branch estimate
    root <- lapply(which.regimes, function(e) tail(e, n= 1))
    root_reg <- tail(simmap_regimes, n=1) #Remove first element - tip value

    #returns list with regime score for root - last regime score from which regime to the right
    #[[1]]
    #[1] 1

    #[[2]]
    #[1] 0
    which.regimes <- lapply(seq_along(levels(regimes)), function(x) c(which.regimes[[x]], root[[x]]))
    #Adds root to which regime scoring
    #[[1]]
    #[1] 1 1 2 2
    #[[2]]
    #[1] 0 0 0 0
    timeflip <- cumsum(c(min_age, unname(subedges)))
    #Minimum root to node time + time spent in each regime per edge - added together cimulative sum
    #[1] 0.0000000 0.1406307 0.4794736 1.0000000
    times <- rev(timeflip)
    # save the regimes in this lineage
    lineage_regimes <- as.factor(c(root_reg,names(subedges)))
  }

  #stop()
  names(which.regimes) <- levels(regimes)
  #$OU1
  #[1] 1 1 1 1
  #$OU2
  #[1] 0 0 0 0
  t_end <- tail(timeflip, n = -1) #Remove first element - tip value
  #[1] 0.0000000 0.5205264 0.8593693 1.0000000 - original
  #[1] 0.5205264 0.8593693 1.0000000 -
  t_beginning <- head(timeflip, n = -1) #Remove last element - root value
  #[1] 0.0000000 0.5205264 0.8593693
  regime_time <- c(t_end - t_beginning, min_age) #Calculate time within a regime?
  #Sum(time at end of linege - time at beginning of lineage)
  #[1] 0.5205264 0.3388429 0.1406307 0.0000000
  return(list(nodes = nodes,
              times = times,
              t_end = t_end,
              t_beginning = t_beginning,
              regime_time = regime_time,
              which.regimes = which.regimes,
              lineage_regimes = lineage_regimes))
}

weights_segments <- function(a, lineage){#For individual lineage, determine the weighting of each segment
  #t_beginning and t_end are both vectors, and subtracting them from each other lines up the beginning and end of one segment
  #because if tge tail/head procedure above
  res <- c(exp(-a * lineage$t_beginning) - exp(-a * lineage$t_end),
           exp(-a * lineage$times[1]))
  return(res)
}

weights_regimes <- function(a, lineage) {#For individual lineage, sum up the segments in each regimes
  #nt <- lineage$nodes_time
  res <- weights_segments(a, lineage) ## Rcpp wrapper, equivalent to above commented out code
  w <- vapply(lineage$which.regimes, function(e) sum(e*res), FUN.VALUE = 0) ## Sum coefficients for which the respective regime is equal
  return(w) #Return named vector with regimes weights for individual lineage
}

weight.matrix <- function(phy, a, lineages){ #Wrapper to apply weights_regimes to each lineage
  res <- t(vapply(lineages, function(x) weights_regimes(a, x),
                  FUN.VALUE = numeric(length(lineages[[1]]$which.regimes))) ## Specify type of output
  )

  rownames(res) <- phy$tip.label
  return(res)
}

## Thanks to user "snaut" at stackoverflow, http://stackoverflow.com/users/1999873/snaut
concat.factor <- function(...){
  as.factor(do.call(c, lapply(list(...), as.character)))
}

############################################################################################################
#Direct Effect Model setup
############################################################################################################

blouch.direct.prep<-function(trdata,Y,Y_error,X,X_error){
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  N<-length(trdata$phy$tip.label)
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  datX<-as.matrix(dat[X])
  datXerror<-as.matrix(dat[X_error])
  Z_direct<-dim(X_error)[1]
  #print(paste(N,Z))
  #Direct effect model w/ Statistical Rethinking ME Correction
  dat<-list(N=N,Z_direct=Z,Y_obs=as.vector(t(dat[Y])),X_obs=matrix(datX,nrow=N,ncol=Z),
            Y_error=as.vector(t(dat[Y_error])),X_error=matrix(datXerror,nrow=N,ncol=Z),ta=ta,tij=tij)
  return(dat)
}
############################################################################################################
#Adaptive Model setup
############################################################################################################

blouch.adapt.prep<-function(trdata,Y,Y_error,X,X_error){
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  N<-length(trdata$phy$tip.label)
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  datX<-as.matrix(dat[X])
  datXerror<-as.matrix(dat[X_error])

  if(length(X)==1){Z<-1}
  else{Z<-dim(X)[2]}
  if(is.na(Y_error)!=TRUE){mv.response<-dat[Y_error]}
  else{mv.response<-as.vector(rep(0,n))}
  if(is.na(X_error)!=TRUE){mv.pred<-matrix(datXerror,nrow=N,ncol=Z)}
  else{mv.pred<-matrix(0,nrow=n,ncol=Z)}

  test<-sigma.X.estimate(phy,ta, predictor = datX, mv.predictor = mv.pred)

  brownian_mean<-test[1]
  sigma_squared_x<-test[2]

  #Direct effect model w/ Statistical Rethinking ME Correction
  dat<-list(N=N,Z_adapt=Z,Y_obs=as.vector(t(dat[Y])),X_obs=matrix(datX,nrow=N,ncol=Z),
            Y_error=as.vector(t(dat[Y_error])),X_error=mv.pred,
            ta=ta,tij=tij,tja=tja,T_term=T_term,sigma2_x=sigma_squared_x)
  return(dat)
}

############################################################################################################
#Direct effect Adaptive Model setup
############################################################################################################

blouch.direct.adapt.prep<-function(trdata,Y,Y_error,X,X_error,Z_direct,Z_adaptive){
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  N<-length(trdata$phy$tip.label)
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]]
  T_term<-ts[[3]]
  tja<-ts[[4]]

  datX<-as.matrix(dat[X])
  #return(Z_adapt)
  datX.direct<-datX[,1:Z_direct]
  datX.adapt<-datX[,(Z_direct+1):(Z_adaptive+Z_direct)]

  datXerror<-as.matrix(dat[X_error])

  datXerror.direct<-datXerror[,1:Z_direct]
  datXerror.adapt<-datXerror[,(Z_direct+1):(Z_adaptive+Z_direct)]
  #return(datXerror.direct)
  if(is.na(Y_error)!=TRUE){mv.response<-dat[Y_error]}
  else{mv.response<-as.vector(rep(0,N))}

  if(any(is.na(datXerror.direct)!=TRUE)){mv.pred.direct<-matrix(datXerror.direct,nrow=N,ncol=Z_direct)}
  else{mv.pred.direct<-matrix(0,nrow=N,ncol=Z_direct)}

  if(any(is.na(datXerror.adapt)!=TRUE)){mv.pred.adapt<-matrix(datXerror.adapt,nrow=N,ncol=Z_adaptive)}
  else{mv.pred.adapt<-matrix(0,nrow=N,ncol=Z_adaptive)}

  test<-sigma.X.estimate(phy,ta, predictor = datX.adapt, mv.predictor = mv.pred.adapt)

  brownian_mean<-test[1]
  sigma_squared_x<-test[2]

  Z<-Z_direct+Z_adaptive
  Z_X_error<-Z_direct+Z_adaptive

dat<-list(N=N,Z_direct=Z_direct,Z_adaptive=Z_adaptive,Z_X_error=Z_X_error,
          Y_obs=as.vector(t(dat[Y])),X_obs=matrix(datX,nrow=N,ncol=Z),
          Y_error=as.vector(t(dat[Y_error])),X_error=matrix(datXerror,nrow=N,ncol=Z),
          ta=ta,tij=tij,tja=tja,T_term=T_term,sigma2_x=sigma2_x)
return(dat)
}
############################################################################################################
#Multi-optima model setup
############################################################################################################
blouch.reg.prep<-function(trdata,Y,Y_error,reg.column,anc_maps="regimes"){
  #Data should be send in treeplyr format
  #Only allows for regime shifts at nodes at present, not SIMMAP
  #Get Phylogeny info
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  N<-length(trdata$phy$tip.label)

  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  #Get internal regimes
  regimes_internal <-trdata$phy$node.label
  regimes_tip <- dat[reg.column][,1]

  regimes <- concat.factor(regimes_tip, regimes_internal)
  #anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
  #return(lineages)
############################################################################################################
  nodes<-NULL
  store<-NULL
  reg_num_lineage<-NULL
  for(i in 1:length(lineages)){
    store<-c(store,length(lineage.nodes(trdata$phy,i))) #Calculate max node height
    reg_num_lineage<-c(reg_num_lineage,length(unique(lineages[[i]]$lineage_regimes)))
    nodes<-c(nodes,length(lineages[[i]]$nodes))
  }
  max_node_num<-max(store)
  times<-matrix(0,length(lineages),max_node_num)
  t_end<-matrix(0,length(lineages),max_node_num)
  t_beginning<-matrix(0,length(lineages),max_node_num)
  reg_match<-data.frame(matrix(0,length(lineages),max_node_num))

  for(i in 1:length(lineages)){
    times[i,1:length(lineages[[i]]$times)]<-lineages[[i]]$times
    t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
    t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
    reg_match[i,1:length(lineages[[i]]$lineage_regimes)]<-rev(as.numeric(lineages[[i]]$lineage_regimes))
  }
  ############################################################################################################
  #Y_obs<-trdata$dat$Y
  #Y_error<-trdata$dat$Y_error #Standard error not variance

  reg_tips<-dat[reg.column][,1]
  reg_tips<-as.numeric(as.factor(reg_tips))
  Dmat<-cophenetic(trdata$phy) #Time separating tips, same as tij matrix in Slouch/Blouch code

  ############################################################################################################
  #print(as.vector(t(dat[Y_error])))
  dat<-list(N=N,n_reg=length(unique(regimes)),max_node_num=max_node_num,Y_obs=as.vector(t(dat[Y])),Y_error=as.vector(t(dat[Y_error])),ta=ta,tij=tij,
            t_beginning=t_beginning,t_end=t_end,times=times,reg_match=reg_match,nodes=nodes,Dmat=Dmat)
  return(dat)
}
############################################################################################################
############################################################################################################

blouch.reg.direct.prep<-function(trdata,Y,Y_error,X,X_error,Z_direct,reg.column){
  #Data should be send in treeplyr format
  #Only allows for regime shifts at nodes at present, not SIMMAP
  #Get Phylogeny info

  anc_maps="regimes"
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  #return(dat)
  N<-length(trdata$phy$tip.label)

  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T_term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  #Get internal regimes
  regimes_internal <-trdata$phy$node.label
  regimes_tip <- dat[reg.column][,1]

  regimes <- concat.factor(regimes_tip, regimes_internal)
  #return(regimes)
  #anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
  ############################################################################################################
  nodes<-NULL
  store<-NULL
  reg_num_lineage<-NULL
  for(i in 1:length(lineages)){
    store<-c(store,length(lineage.nodes(trdata$phy,i))) #Calculate max node height
    reg_num_lineage<-c(reg_num_lineage,length(unique(lineages[[i]]$lineage_regimes)))
    nodes<-c(nodes,length(lineages[[i]]$nodes))
  }
  max_node_num<-max(store)
  times<-matrix(0,length(lineages),max_node_num)
  t_end<-matrix(0,length(lineages),max_node_num)
  t_beginning<-matrix(0,length(lineages),max_node_num)
  reg_match<-data.frame(matrix(0,length(lineages),max_node_num))

  for(i in 1:length(lineages)){
    times[i,1:length(lineages[[i]]$times)]<-lineages[[i]]$times
    t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
    t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
    reg_match[i,1:length(lineages[[i]]$lineage_regimes)]<-rev(as.numeric(lineages[[i]]$lineage_regimes))
  }
  ############################################################################################################
  reg_tips<-dat[reg.column][,1]
  reg_tips<-as.numeric(as.factor(reg_tips))
  Dmat<-cophenetic(trdata$phy) #Time separating tips, same as tij matrix in Slouch/Blouch code

  dat<-list(N=N,n_reg=length(unique(regimes)),Z_direct=Z_direct,Z_X_error=Z_direct,max_node_num=max_node_num,
            Y_obs=as.vector(t(dat[Y])),X_obs=data.matrix(dat[X]),#,nrow=N,ncol=Z_direct),
            Y_error=as.vector(t(dat[Y_error])),X_error=data.matrix(dat[X_error]),#matrix(dat[X_error],nrow=N,ncol=Z_direct),
            ta=ta,tij=tij,tja=tja,T_term=T_term,t_beginning=t_beginning,
            t_end=t_end,times=times,reg_match=reg_match,nodes=nodes,reg_tips=reg_tips)

  return(dat)
}

blouch.reg.adapt.prep<-function(trdata,Y,Y_error,X,X_error,Z_adaptive,reg.column){
  #Data should be send in treeplyr format
  #Only allows for regime shifts at nodes at present, not SIMMAP
  #Get Phylogeny info
  anc_maps="regimes"
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  #return(dat)
  N<-length(trdata$phy$tip.label)

  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T_term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  #Get internal regimes
  regimes_internal <-trdata$phy$node.label
  regimes_tip <- dat[reg.column][,1]
  regimes <- concat.factor(regimes_tip, regimes_internal)
  #anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
  ############################################################################################################
  nodes<-NULL
  store<-NULL
  reg_num_lineage<-NULL
  for(i in 1:length(lineages)){
    store<-c(store,length(lineage.nodes(trdata$phy,i))) #Calculate max node height
    reg_num_lineage<-c(reg_num_lineage,length(unique(lineages[[i]]$lineage_regimes)))
    nodes<-c(nodes,length(lineages[[i]]$nodes))
  }
  max_node_num<-max(store)
  times<-matrix(0,length(lineages),max_node_num)
  t_end<-matrix(0,length(lineages),max_node_num)
  t_beginning<-matrix(0,length(lineages),max_node_num)
  reg_match<-data.frame(matrix(0,length(lineages),max_node_num))

  for(i in 1:length(lineages)){
    times[i,1:length(lineages[[i]]$times)]<-lineages[[i]]$times
    t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
    t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
    reg_match[i,1:length(lineages[[i]]$lineage_regimes)]<-rev(as.numeric(lineages[[i]]$lineage_regimes))
  }
  ############################################################################################################

  datX<-as.matrix(dat[X])
  datXerror<-as.matrix(dat[X_error])
  mv.pred.adapt<-matrix(datXerror,nrow=N,ncol=Z_adaptive)
  test<-sigma.X.estimate(phy,ta, predictor = datX, mv.predictor = mv.pred.adapt)

  brownian_mean<-test[1]
  sigma_squared_x<-test[2]

  reg_tips<-dat[reg.column][,1]
  reg_tips<-as.numeric(as.factor(reg_tips))
  Dmat<-cophenetic(trdata$phy) #Time separating tips, same as tij matrix in Slouch/Blouch code

  dat<-list(N=N,n_reg=length(unique(regimes)),Z_adaptive=Z_adaptive,Z_X_error=Z_adaptive,
            max_node_num=max_node_num,
            Y_obs=as.vector(t(dat[Y])),X_obs=data.matrix(dat[X]),
            Y_error=as.vector(t(dat[Y_error])),X_error=data.matrix(dat[X_error]),
            sigma2_x=sigma2_x,ta=ta,tij=tij,tja=tja,
            T_term=T_term,t_beginning=t_beginning,
            t_end=t_end,times=times,reg_match=reg_match,nodes=nodes,reg_tips=reg_tips)

  return(dat)
}

blouch.reg.direct.adapt.prep<-function(trdata,Y,Y_error,X,X_error,Z_direct,Z_adaptive,reg.column){
  #Data should be send in treeplyr format
  #Only allows for regime shifts at nodes at present, not SIMMAP
  #Get Phylogeny info
  anc_maps="regimes"
  phy<-trdata$phy
  dat<-data.frame(trdata$dat)
  #return(dat)
  N<-length(trdata$phy$tip.label)

  n<-length(trdata$phy$tip.label)
  mrca1 <- ape::mrca(trdata$phy)
  times <- ape::node.depth.edgelength(trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
  T_term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia

  #Get internal regimes
  regimes_internal <-trdata$phy$node.label
  regimes_tip <- dat[reg.column][,1]
  regimes <- concat.factor(regimes_tip, regimes_internal)
  #anc_maps<-"regimes"
  lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
  ############################################################################################################
  nodes<-NULL
  store<-NULL
  reg_num_lineage<-NULL
  for(i in 1:length(lineages)){
    store<-c(store,length(lineage.nodes(trdata$phy,i))) #Calculate max node height
    reg_num_lineage<-c(reg_num_lineage,length(unique(lineages[[i]]$lineage_regimes)))
    nodes<-c(nodes,length(lineages[[i]]$nodes))
  }
  max_node_num<-max(store)
  times<-matrix(0,length(lineages),max_node_num)
  t_end<-matrix(0,length(lineages),max_node_num)
  t_beginning<-matrix(0,length(lineages),max_node_num)
  reg_match<-data.frame(matrix(0,length(lineages),max_node_num))

  for(i in 1:length(lineages)){
    times[i,1:length(lineages[[i]]$times)]<-lineages[[i]]$times
    t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
    t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
    reg_match[i,1:length(lineages[[i]]$lineage_regimes)]<-rev(as.numeric(lineages[[i]]$lineage_regimes))
  }
  ############################################################################################################
  Z_X_error<-length(X_error) #Number of X traits with error

  reg_tips<-dat[reg.column][,1]
  reg_tips<-as.numeric(as.factor(reg_tips))
  Dmat<-cophenetic(trdata$phy) #Time separating tips, same as tij matrix in Slouch/Blouch code
  ############################################################################################################

  datX<-as.matrix(dat[X])
  #return(Z_adapt)
  datX.direct<-datX[,1:Z_direct]
  datX.adapt<-datX[,(Z_direct+1):(Z_adaptive+Z_direct)]

  datXerror<-as.matrix(dat[X_error])

  datXerror.direct<-datXerror[,1:Z_direct]
  datXerror.adapt<-datXerror[,(Z_direct+1):(Z_adaptive+Z_direct)]
  #return(datXerror.direct)
  if(is.na(Y_error)!=TRUE){mv.response<-dat[Y_error]}
  else{mv.response<-as.vector(rep(0,N))}

  if(any(is.na(datXerror.direct)!=TRUE)){mv.pred.direct<-matrix(datXerror.direct,nrow=N,ncol=Z_direct)}
  else{mv.pred.direct<-matrix(0,nrow=N,ncol=Z_direct)}

  if(any(is.na(datXerror.adapt)!=TRUE)){mv.pred.adapt<-matrix(datXerror.adapt,nrow=N,ncol=Z_adaptive)}
  else{mv.pred.adapt<-matrix(0,nrow=N,ncol=Z_adaptive)}

  test<-sigma.X.estimate(phy,ta, predictor = datX.adapt, mv.predictor = mv.pred.adapt)

  brownian_mean<-test[1]
  sigma_squared_x<-test[2]

  Z<-Z_direct+Z_adaptive
  Z_X_error<-Z_direct+Z_adaptive

  dat<-list(N=N,n_reg=length(unique(regimes)),Z_direct=Z_direct,Z_adaptive=Z_adaptive,Z_X_error=Z_X_error,max_node_num=max_node_num,
            Y_obs=as.vector(t(dat[Y])),X_obs=data.matrix(dat[X]),Y_error=as.vector(t(dat[Y_error])),
            X_error=data.matrix(dat[X_error]),sigma2_x=sigma2_x,ta=ta,tij=tij,
            tja=tja,T_term=T_term,t_beginning=t_beginning,
            t_end=t_end,times=times,reg_match=reg_match,nodes=nodes,reg_tips=reg_tips)


  return(dat)
}
