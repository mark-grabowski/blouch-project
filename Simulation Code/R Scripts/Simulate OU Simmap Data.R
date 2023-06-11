#Code to simulate data for regime painting on phylogeny
#To be used with Blouch SBR1 - Validation Code.R
rm(list=ls())

calc_direct_V<-function(phy, sigma2_y, a){ #Calculate V matrix for direct effect and regime-only models
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  Vt<-sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) * exp(-a * tij)) #ta - time from root to tips, tij  - total time separating spcies
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



################################################################################################
#Make SIMMAP using tip regimes
#50 tips
#Basic Setup
#Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997) Original Stan
#library(parallel)
#num.cores=parallel::detectCores()-2 #For mclapply function
#sim.num<-2 #25
#loops.num<-1 #2
#set.seed(1)

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
#Using node 54

#Paint Regimes on Tree
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Simulation Code/Functions/set.converge.regimes.R") #Macbook Pro
#source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Simulation Code/Functions/set.converge.regimes.R") #Mac Studio

shifts<-c(54) #Location of nodes with regime shifts
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

print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Make simmap tree
#Setup dataset for info on model fits
cl<-parallel::detectCores()-2 #Cl for fitMk function
x<-trdata$dat$regimes
names(x)<-trdata$phy$tip.label
fit.ER<-fitMk(trdata$phy,x,model="ARD")

fit<-fit.ER
fittedQ<-matrix(NA,length(fit$states),length(fit$states))
fittedQ[]<-c(0,fit$rates)[fit$index.matrix+1]
diag(fittedQ)<-0
diag(fittedQ)<--rowSums(fittedQ)
colnames(fittedQ)<-rownames(fittedQ)<-fit$states
tree<-trdata$phy
simmap_trees<-make.simmap(tree,x,Q=fittedQ,nsim=2)
simmap_trees[[1]]$mapped.edge

simmap_tree1<-simmap_trees[[1]]

#Simulate data using simmap
n<-length(simmap_tree1$tip.label)
mrca1 <- ape::mrca(simmap_tree1)
times <- ape::node.depth.edgelength(simmap_tree1)
ta <- matrix(times[mrca1], nrow=n, dimnames = list(simmap_tree1$tip.label, simmap_tree1$tip.label))
T.term <- times[1:n]
tia <- times[1:n] - ta
tja <- t(tia)
tij <- tja + tia

#regimes_internal <-trdata$phy$node.label
#regimes_tip <- trdata$dat$regimes
#regimes <- concat.factor(regimes_tip, regimes_internal)
regimes<-as.factor(simmap_tree1$node.label)
anc_maps<-"simmap"
lineages <- lapply(1:n, function(e) lineage.constructor(simmap_tree1, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

#chr [1:3] "OU1" "OU1" "OU1"
#Simulate Y based on V and incorporating regimes
#Setup parameters

hl<-0.1 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
a<-log(2)/hl
vy<-0.1 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));
optima<-c(0.35,0.25) #Intercepts for two regimes

dmX<-weight.matrix(simmap_tree1, a, lineages) #Slouch approach
mu<-dmX%*%optima #Simulate mu for Y
V<-calc_direct_V(phy, sigma2_y, a)
Y<-mvrnorm(n=1,mu,V)

nodes<-NULL
store<-NULL
times_store<-NULL
reg_num_lineage<-NULL
for(i in 1:length(lineages)){
  store<-c(store,length(lineage.nodes(trdata$phy,i))) #Calculate max node height
  reg_num_lineage<-c(reg_num_lineage,length(unique(lineages[[i]]$lineage_regimes)))
  nodes<-c(nodes,length(lineages[[i]]$nodes))
  times_store<-c(times_store,length(lineages[[i]]$times))
  
}
max_node_num<-max(store)  
max_time_seg<-max(times_store)  

times<-matrix(0,length(lineages),max_time_seg)
t_end<-matrix(0,length(lineages),max_time_seg)
t_beginning<-matrix(0,length(lineages),max_time_seg)
reg_match<-data.frame(matrix(0,length(lineages),max_time_seg))

for(i in 1:length(lineages)){
  #  print(i)
  times[i,1:length(lineages[[i]]$times)]<-lineages[[i]]$times
  t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
  t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
  reg_match[i,1:length(lineages[[i]]$lineage_regimes)]<-rev(as.numeric(lineages[[i]]$lineage_regimes))
}


dat<-list(N=N,n_reg=length(unique(regimes)),max_node_num=max_time_seg,Y_obs=Y,ta=ta,tij=tij,t_beginning=t_beginning,t_end=t_end,times=times,reg_match=reg_match,nodes=times_store)
plot(simmap_tree1)

################################################################################################
#Multi-simmap trees
#Make SIMMAP using tip regimes
#50 tips
#Basic Setup
#Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997) Original Stan
#library(parallel)
#num.cores=parallel::detectCores()-2 #For mclapply function
#sim.num<-2 #25
#loops.num<-1 #2
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
#Using node 54

#Paint Regimes on Tree
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Simulation Code/Functions/set.converge.regimes.R") #Macbook Pro
#source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Simulation Code/Functions/set.converge.regimes.R") #Mac Studio

shifts<-c(54) #Location of nodes with regime shifts
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

print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)

#Make simmap tree
#Setup dataset for info on model fits
cl<-parallel::detectCores()-2 #Cl for fitMk function
x<-trdata$dat$regimes
names(x)<-trdata$phy$tip.label
fit.ER<-fitMk(trdata$phy,x,model="ER")

fit<-fit.ER
fittedQ<-matrix(NA,length(fit$states),length(fit$states))
fittedQ[]<-c(0,fit$rates)[fit$index.matrix+1]
diag(fittedQ)<-0
diag(fittedQ)<--rowSums(fittedQ)
colnames(fittedQ)<-rownames(fittedQ)<-fit$states
tree<-trdata$phy
simmap_trees<-make.simmap(tree,x,Q=fittedQ,nsim=2)
simmap_trees[[1]]$mapped.edge

simmap_tree1<-simmap_trees[[1]]

#Simulate data using simmap
n<-length(simmap_tree1$tip.label)
mrca1 <- ape::mrca(simmap_tree1)
times <- ape::node.depth.edgelength(simmap_tree1)
ta <- matrix(times[mrca1], nrow=n, dimnames = list(simmap_tree1$tip.label, simmap_tree1$tip.label))
T.term <- times[1:n]
tia <- times[1:n] - ta
tja <- t(tia)
tij <- tja + tia

#regimes_internal <-trdata$phy$node.label
#regimes_tip <- trdata$dat$regimes
#regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"simmap"

#chr [1:3] "OU1" "OU1" "OU1"
#Simulate Y based on V and incorporating regimes
#Setup parameters

hl<-0.1
a<-log(2)/hl
sigma2_y<-0.1
#alpha<-4 #Intecept
optima<-c(0.35,0.25) #Intercepts for two regimes
num.trees<-length(simmap_trees)
times_list<-list()
t_end_list<-list()
t_beg_list<-list()
reg_match_list<-list()
times_store_list<-list()

for(k in 1:num.trees){
  simmap_tree<-simmap_trees[[k]]
  regimes<-as.factor(simmap_tree$node.label)
  lineages <- lapply(1:n, function(e) lineage.constructor(simmap_tree, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch
  dmX<-weight.matrix(simmap_tree, a, lineages) #Slouch approach
  mu<-dmX%*%optima #Simulate mu for Y
  V<-calc_direct_V(phy, sigma2_y, a)
  Y<-mvrnorm(n=1,mu,V)
  
  nodes<-NULL
  store<-NULL
  times_store<-NULL
  reg_num_lineage<-NULL
  times_array <- array(0, dim=c(num.trees,2,3))
  
  for(i in 1:length(lineages)){
    store<-c(store,length(lineage.nodes(trdata$phy,i))) #Calculate max node height
    reg_num_lineage<-c(reg_num_lineage,length(unique(lineages[[i]]$lineage_regimes)))
    nodes<-c(nodes,length(lineages[[i]]$nodes))
    times_store<-c(times_store,length(lineages[[i]]$times))
  }
  max_node_num<-max(store)  
  max_time_seg<-max(times_store)  
  
  times<-matrix(0,length(lineages),max_time_seg)
  t_end<-matrix(0,length(lineages),max_time_seg)
  t_beginning<-matrix(0,length(lineages),max_time_seg)
  reg_match<-data.frame(matrix(0,length(lineages),max_time_seg))
  
  for(i in 1:length(lineages)){
    #  print(i)
    times[i,1:length(lineages[[i]]$times)]<-lineages[[i]]$times
    t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
    t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
    reg_match[i,1:length(lineages[[i]]$lineage_regimes)]<-rev(as.numeric(lineages[[i]]$lineage_regimes))
  }
  times_list[[length(times_list) + 1]]<-times
  t_end_list[[length(t_end_list) + 1]]<-t_end
  t_beg_list[[length(t_beg_list) + 1]]<-t_beginning
  reg_match_list[[length(reg_match_list) + 1]]<-reg_match
  times_store_list[[length(times_store_list) + 1]]<-times_store
}



#Turn list into 
library(abind)
times_matrix<-times_list[[1]]
t_end_matrix<-t_end_list[[1]]
t_beg_matrix<-t_beg_list[[1]]
reg_match_matrix<-reg_match_list[[1]]
times_store_matrix<-times_store_list[[1]]

for(i in 2:k){
  times_matrix<-abind(times_matrix,times_list[[i]],along=3)
  t_end_matrix<-abind(t_end_matrix,t_end_list[[i]],along=3)
  t_beg_matrix<-abind(t_beg_matrix,t_beg_list[[i]],along=3)
  reg_match_matrix<-abind(reg_match_matrix,reg_match_list[[i]],along=3)
  times_store_matrix<-cbind(times_store_matrix,times_store_list[[i]])
}
times_matrix<-aperm(times_matrix, c(3,1,2))
t_end_matrix<-aperm(t_end_matrix, c(3,1,2))
t_beg_matrix<-aperm(t_beg_matrix, c(3,1,2))
reg_match_matrix<-aperm(reg_match_matrix, c(3,1,2))

dat<-list(N=N,n_reg=length(unique(regimes)),T=num.trees,max_node_num=max_time_seg,Y_obs=Y,ta=ta,tij=tij,t_beginning=t_beg_matrix,t_end=t_end_matrix,times=times_matrix,reg_match=reg_match_matrix,nodes=times_store_matrix)

