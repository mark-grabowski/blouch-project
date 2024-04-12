#Experimental code
#Code to simulate data for regime painting on phylogeny with multiple varying effects for dicrect and adaptive model
#blouchOU_reg_multidirectadaptive_ME_VarEff_nc.stan
#To be used with Blouch SBR1 - Validation Code.R

rm(list=ls())

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

#Script to simulate data to test Blouch OU regimes model
library(devtools)
library(ape)
library(slouch)
library(rstan)
library(treeplyr)
library(ggplot2)
library(ggsci)
library(MASS)

########################################################################################################
#Four regimes with 1 direct and 1 adaptive trait and multiple slopes per optima but single alpha parameter
set.seed(10)

tree.10K<-read.tree("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
#tree.10K<-read.tree("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
N<-50 #Number of species
#set.seed(1) #Set seed to get same random species each time

phy <- keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N]) 
phy<-multi2di(phy)

l.tree<-max(branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree 

#Set regimes - manually - 2 regimes
#Locate nodes
plot(phy,no.margin=TRUE,edge.width=2,cex=0.7)
nodelabels(frame="none",adj=c(1.1,-0.4))
tiplabels()

#Paint Regimes on Tree
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Simulation Code/Functions/set.converge.regimes.R") #Macbook Pro
#source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Simulation Code/Functions/set.converge.regimes.R") #Mac Studio

shifts<-c(84) #Location of nodes with regime shifts
trdata<-data.frame(phy$tip.label)
trdata<-make.treedata(phy,trdata)
trdata<-set.converge.regimes(trdata,shifts)


#Check if manual setting code worked
shifts.total<-c(trdata$dat$regimes,trdata$phy$node.label)
edge.regimes <- factor(shifts.total[trdata$phy$edge[,2]])
print(edge.regimes)
#Get ggplot colors used for plot to make on tree
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

reg.colors<-gg_color_hue(length(unique(trdata$dat$regimes)))

reg_tips<-trdata$dat$regimes
reg_tips<-as.numeric(as.factor(reg_tips))

print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Phylogeny info
n<-length(trdata$phy$tip.label)
mrca1 <- ape::mrca(trdata$phy)
times <- ape::node.depth.edgelength(trdata$phy)
ta <- matrix(times[mrca1], nrow=n, dimnames = list(trdata$phy$tip.label, trdata$phy$tip.label))
T.term <- times[1:n]
tia <- times[1:n] - ta
tja <- t(tia)
tij <- tja + tia

regimes_internal <-trdata$phy$node.label
regimes_tip <- trdata$dat$regimes
regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"regimes"
lineages <- lapply(1:n, function(e) lineage.constructor(trdata$phy, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

###############################################################################################################################################################################
hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.01 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0
Sxx<-10 #Look at effects

ts<-ts_fxn(phy)
ta<-ts[[1]]
tij<-ts[[2]]
T_term<-ts[[3]]
tja<-ts[[4]]

Xa<-fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
Xd<-fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
names(Xa)<-phy$tip.label
names(Xd)<-phy$tip.label
phenogram(phy,Xd,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
phenogram(phy,Xa,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

Xs<-cbind(Xd,Xa)
sigma2_x<-matrix(1,1,1)

Z_adaptive<-1
Z_direct<-1

optima_matrix<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
pred_X<-calc_mixed_dmX(a,T_term,Xs,Z_direct,Z_adaptive)
optima<-c(2,1)
#beta<-c(0.25,0.15,0.35,0.1) #Two Optima/Two Slopes
beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Two traits on columns, two regimes on vertical

mu<-matrix(NA,N,1)
for(i in 1:N){
#  mu[i] = optima_matrix[i,]%*%optima+beta[reg_tips[i]]%*%pred_X[i]
  mu[i] = optima_matrix[i,]%*%optima+pred_X[i,]%*%t(beta[reg_tips[i],])

  }



n_reg<-length(unique(regimes))
V<-calc_adaptive_V(a, sigma2_y, ta,  tij,  tja,  T_term,  beta[,2],  sigma2_x, Z_adaptive, n_reg)
Y<-mvrnorm(n=1,mu,V)

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


dat<-list(N=N,n_reg=length(unique(regimes)),Z_diret=Z_direct,Z_adaptive=Z_adaptive,Y_obs=Y,X_obs=matrix(Xs,nrow=N,ncol=Z_direct+Z_adaptive),
          max_node_num=max_node_num,ta=ta,tij=tij,tja=tja,T_term=T_term,t_beginning=t_beginning,
          t_end=t_end,times=times,reg_match=reg_match,nodes=nodes)

##################################################################################################################
#Simulate errors - original Hansen setup
Z_X_error<-Z_direct+Z_adaptive #Number of X traits with error
X_error<-matrix(0.01,nrow=N,ncol=2)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-Xs+rnorm(N,0,0.01)

dat<-list(N=N,n_reg=length(unique(regimes)),Z_direct=Z_direct,Z_adaptive=Z_adaptive,Z_X_error=Z_X_error,max_node_num=max_node_num,
          Y_obs=Y_with_error,X_obs=matrix(X_with_error,nrow=N,ncol=Z_direct+Z_adaptive),
          Y_error=Y_error,X_error=matrix(X_error,nrow=N,ncol=Z_X_error),
          sigma2_x=sigma2_x,ta=ta,tij=tij,tja=tja,T_term=T_term,t_beginning=t_beginning,
          t_end=t_end,times=times,reg_match=reg_match,nodes=nodes,reg_tips=reg_tips)

########################################################################################################
#Plot of data
df<-data.frame(Y=dat$Y_obs,X=dat$X_obs,Regimes=regimes_tip)


#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.1<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X.1,color=Regimes))+
  
  geom_abline(intercept=optima[1],slope=beta[1,1],alpha=0.5,linetype=2)+
  geom_abline(intercept=optima[2],slope=beta[2,1],alpha=0.5,linetype=2)+
  
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("Y") + xlab("Direct effect trait")+
  scale_color_npg()

slope.plot.1

#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.2<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X.2,color=Regimes))+
  
  geom_abline(intercept=optima[1],slope=beta[1,2],alpha=0.5,linetype=2)+
  geom_abline(intercept=optima[2],slope=beta[2,2],alpha=0.5,linetype=2)+
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("Y") + xlab("Adaptive trait")+
  scale_color_npg()

slope.plot.2



