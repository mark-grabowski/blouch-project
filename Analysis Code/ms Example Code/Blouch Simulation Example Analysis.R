#Code for SBR1 Simulation Analysis
#Code to simulate data for regime painting on phylogeny with multiple varying slopes for adaptive model
#Also has correlated varying effects in some cases
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
########################################################################################################

#Script to simulate data to test Blouch OU regimes model
library(devtools)
library(ape)
library(slouch)
library(treeplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(MASS)
library(rstan)
library(rethinking)
#For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
#options(mc.cores = 8)
rstan_options(auto_write = TRUE)

dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
arch <- ifelse(R.version$arch == "aarch64", "arm64", "x86_64")
cat(paste("\nCXX14FLAGS += -O3 -mtune=native -arch", arch, "-ftemplate-depth-256"),
    file = M, sep = "\n", append = FALSE)

########################################################################################################
#Milestone 14/15
#Four regimes with one adaptive trait and multiple slopes per optima but single alpha parameter
set.seed(10)

tree.10K<-read.tree("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
#tree.10K<-read.tree("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
N<-100 #Number of species
#N<-50 #Number of species
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

############################################################################################################
#Paint Regimes on Tree
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Simulation Code/Functions/set.converge.regimes.R") #Macbook Pro
#source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Simulation Code/Functions/set.converge.regimes.R") #Mac Studio

shifts<-c(164,192,104) #Location of nodes with regime shifts #100 species
#shifts<-c(83,72,65) #Location of nodes with regime shifts #50 species
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

#reg.colors<-ggsci::pal_aaas("default",alpha=0.7)(4)
reg.colors<-ggsci::pal_npg(palette=c("nrc"),alpha=1)(4)

print(reg.colors)
plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1,show.tip.label=FALSE)
#tiplabels(pch=19,cex=1,col=reg.colors[factor(trdata$dat$regimes)])


reg_tips<-trdata$dat$regimes
reg_tips<-as.numeric(as.factor(reg_tips))

############################################################################################################
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

#########################
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

X<-fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
sigma2_x<-matrix(1,1,1)
Z_adaptive<-1
names(X)<-phy$tip.label
phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

optima_matrix<-weight.matrix(trdata$phy, a, lineages) #Slouch approach
pred_X<-calc_adaptive_dmX(a,T_term,X)
#optima<-c(2,1.5,1,0.5)
optima<-c(1,2,3,4)
beta<-c(0.75,0.5,0.35,0.25) #Two Optima/Two Slopes
#beta<-c(0.25,0.15,0.35,0.1) #Two Optima/Two Slopes

#beta<-data.frame(matrix(c(0.25,0.15,0.35,0.1),ncol=2,nrow=2)) #Two traits on columns, two regimes on vertical

mu<-matrix(NA,N,1)
for(i in 1:N){
  mu[i] = optima_matrix[i,]%*%optima+beta[reg_tips[i]]%*%pred_X[i]
  #mu[i] = optima_matrix[i,]%*%optima+pred_X[i,]%*%t(beta[reg_tips[i],])
  
}

n_reg<-length(unique(regimes))
V<-calc_adaptive_V(a, sigma2_y, ta,  tij,  tja,  T_term,  beta,  sigma2_x, Z_adaptive, n_reg)
Y<-mvrnorm(n=1,mu,V)

##################################################################################################################
#Simulate errors - original Hansen setup
Z_X_error<-1 #Number of X traits with error
X_error<-matrix(0.01,nrow=N,ncol=Z_X_error)
X_error<-data.frame(X_error)
#names(X_error)<-c("Xd_error","Xa_error")
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-X+rnorm(N,0,0.01)

############################################################################################################
#Code using blouch.prep function
############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_with_error,Y_error,X_with_error,X_error)))
############################################################################################################
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/R Setup Code/blouch.prep.R")
dat<-blouch.reg.adapt.prep(trdata,"Y_with_error","Y_error","X_with_error","X_error","regimes")

############################################################################################################
#Original Code
############################################################################################################
#nodes<-NULL
#store<-NULL
#reg_num_lineage<-NULL
#for(i in 1:length(lineages)){
#  store<-c(store,length(lineage.nodes(trdata$phy,i))) #Calculate max node height
#  reg_num_lineage<-c(reg_num_lineage,length(unique(lineages[[i]]$lineage_regimes)))
#  nodes<-c(nodes,length(lineages[[i]]$nodes))
#}
#max_node_num<-max(store)  
#times<-matrix(0,length(lineages),max_node_num)
#t_end<-matrix(0,length(lineages),max_node_num)
#t_beginning<-matrix(0,length(lineages),max_node_num)
#reg_match<-data.frame(matrix(0,length(lineages),max_node_num))

#for(i in 1:length(lineages)){
#  times[i,1:length(lineages[[i]]$times)]<-lineages[[i]]$times
#  t_end[i,1:length(lineages[[i]]$t_end)]<-lineages[[i]]$t_end
#  t_beginning[i,1:length(lineages[[i]]$t_beginning)]<-lineages[[i]]$t_beginning
#  reg_match[i,1:length(lineages[[i]]$lineage_regimes)]<-rev(as.numeric(lineages[[i]]$lineage_regimes))
#}

##################################################################################################################
#Simulate errors - original Hansen setup
#Z_X_error<-1 #Number of X traits with error
#X_error<-matrix(0.01,nrow=N,ncol=1)
#Y_error<-rep(0.01,N)
#Y_with_error<-Y+rnorm(N,0,0.01)
#X_with_error<-X+rnorm(N,0,0.01)

#2 Regimes with direct effect model with regime info for tips
#dat<-list(N=N,n_reg=length(unique(regimes)),Z_adaptive=Z_adaptive,Z_X_error=Z_X_error,max_node_num=max_node_num,
#          Y_obs=Y_with_error,X_obs=matrix(X_with_error,nrow=N,ncol=Z_adaptive),
#          Y_error=Y_error,X_error=matrix(X_error,nrow=N,ncol=Z_X_error),
#          sigma2_x=sigma2_x,ta=ta,tij=tij,tja=tja,T_term=T_term,t_beginning=t_beginning,
#          t_end=t_end,times=times,reg_match=reg_match,nodes=nodes,reg_tips=reg_tips)
############################################################################################################
#Prior Exploration Plot
lm.allometric<-summary(lm(dat$Y_obs~dat$X_obs))
lm.allometric$coefficients

#Prior vs. Posterior Plot
library(ggsci)
library(rethinking)

alpha.sims<-rnorm(100,lm.allometric$coefficients[1],1.25)
beta.sims<-rnorm(n=100,lm.allometric$coefficients[2],0.25)

df<-data.frame(Y=dat$Y_obs,X=dat$X_obs[,1])
names(df)<-c("Y","X")

slope.plot<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=regimes_tip))+
  geom_abline(intercept=alpha.sims,slope=beta.sims,alpha=0.1)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("Y") + xlab("Adaptive Predictor")+
  scale_color_npg()

slope.plot

########################################################################################################
#Milestone 16 - mlm with varying effects
#Priors
#hl ~ lognormal(log(0.25),0.75);
#vy ~ exponential(20);
#optima_bar ~ normal(2.88,1.5);//Original 4 regimes
#beta_bar ~ normal(0.31,0.25); //Original 4 regimes
#Rho ~ lkj_corr(4);

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_adapt_mlm_ve.stan")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/chatgpt_blouchOU_reg_adapt_mlm_ve.stan")

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve.stan")
fit.reg.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000)
print(fit.reg.adapt.mlm.ve,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
plot(precis(fit.reg.adapt.mlm.ve,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.reg.adapt.mlm.ve)
########################################################################################################
#Milestone 17 - mlm with varying effects - non-centered version
#Regime model with multiadaptive model with measurement error and varying effects - non-centered version
#Priors
#hl ~ lognormal(log(0.25),0.75);
#vy ~ exponential(20);
#L_Rho ~ lkj_corr_cholesky(2);
#sigma ~ normal(0,1);
#optima_bar ~ normal(2.88,1.5);//Original 4 regimes
#beta_bar ~ normal(0.31,0.25); //Original 4 regimes

#setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_adapt_mlm_ve_nc.stan")
#stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
#fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000)
#print(fit.reg.adapt.mlm.ve.nc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
#plot(precis(fit.reg.adapt.mlm.ve.nc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
#post<-extract(fit.reg.adapt.mlm.ve.nc)
########################################################################################################
#Milestone 15 - varying effects model
#Combination of regime model with adaptive model with measurement error and varying slopes
#Priors
#hl ~ lognormal(log(0.25),0.75);
#vy ~ exponential(20);
#optima ~ normal(2.88,1.5);//Original 4 regimes
#for(i in 1:(Z_adaptive)){
#  beta[,i] ~ normal(0.31,0.25);
#}

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_adapt_ve.stan")
stan_model <- stan_model("blouchOU_reg_adapt_ve.stan")
fit.reg.adapt.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000)
print(fit.reg.adapt.ve,pars = c("hl","vy","optima","beta","beta_e"))
plot(precis(fit.reg.adapt.ve,depth=3,pars = c("hl","vy","optima","beta","beta_e")))
post<-extract(fit.reg.adapt.ve)
########################################################################################################
#Milestone 10 - mlm with varying intercepts
#Multilevel model - multilevel optima with adaptive predictor and measurement error
#Priors
#hl ~ lognormal(log(0.25),0.75);
#vy ~ exponential(20);
#optima_bar ~ normal(2.88,0.5);
#beta ~ normal(0.31,0.1);
#sigma ~ normal(0,1);

#setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_adapt_mlm_vi.stan")
#stan_model <- stan_model("blouchOU_reg_adapt_mlm_vi.stan")
#fit.reg.adapt.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =4000,cores=2)
#print(fit.reg.adapt.mlm.vi,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
#plot(precis(fit.reg.adapt.mlm.vi,depth=2,pars = c("hl","vy","optima","optima_bar","beta","sigma")))

#post<-extract(fit.reg.adapt.mlm.vi)

########################################################################################################
#Milestone 6 - basic model
#Regimes with adaptation model with measurement error
#Priors
#hl ~ lognormal(log(0.25),0.75);
#vy ~ exponential(20);
#optima ~ normal(2.88,0.5);
#beta ~ normal(0.31,0.1);

#setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_adapt.stan")
#stan_model <- stan_model("blouchOU_reg_adapt.stan")

#fit.reg.adapt<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =4000,cores=2)

#print(fit.reg.adapt,pars = c("hl","vy","optima","beta","beta_e"))
#post<-extract(fit.reg.adapt)

########################################################################################################
#Hl Plot prior vs. posterior - assume posterior has been extracted using extract(model) and stored in post

hl.sims<-data.frame(rlnorm(n=1000,meanlog=log(0.25),sdlog=0.25))
#hl.sims<-data.frame(hl.sims[hl.sims<3])
names(hl.sims)<-"prior.hl.sims"

hl.post<-data.frame(post$hl)
names(hl.post)<-"post.hl.sims"

hl.plot<-ggplot()+
  geom_density(aes(prior.hl.sims,fill="prior.hl.sims"),alpha=0.2,data=hl.sims)+
  geom_density(aes(post.hl.sims,fill="post.hl.sims"),alpha=0.2,data=hl.post)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #labs(title="Prior vs. Posterior Distribution ",x="Half-life", y = "Density")+
  labs(title="",x="Half-life", y = "Density")+
  
  #scale_fill_manual(labels=c("Posterior","Prior"))+
  geom_vline(xintercept=c(hl),linetype=2)+
  scale_fill_npg(name="",labels=c("Posterior","Prior"))

hl.plot
########################################################################################################
vy.sims<-rexp(n=1000,rate=20)
vy.sims<-data.frame(vy.sims)
names(vy.sims)<-"prior.vy.sims"


vy.post<-data.frame(post$vy)
names(vy.post)<-"post.vy.sims"


vy.plot<-ggplot()+
  geom_density(aes(prior.vy.sims,fill="prior.vy.sims"),alpha=0.2,data=vy.sims)+
  geom_density(aes(post.vy.sims,fill="post.vy.sims"),alpha=0.2,data=vy.post)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #labs(title="Prior vs. Posterior Distribution ",x="vy", y = "Density")+
  labs(title="",x="vy", y = "Density")+
  geom_vline(xintercept=c(vy),linetype=2)+
  
  #scale_fill_manual(labels=c("Posterior","Prior"))+
  scale_fill_npg(name="",labels=c("Posterior","Prior"))

vy.plot
########################################################################################################
#Adaptation model - multiple regimes

a.sims<-log(2)/hl.sims;
sigma2_y.sims<-vy.sims*(2*(log(2)/hl.sims));
#Sigma<-matrix(c(0.25,0,0,0.25),2,2)
#Sigma<-matrix(c(0.25,0,0,0,0,0.25,0,0,0,0,0.25,0,0,0,0,0.25),4,4)
beta.sims<-replicate(length(beta),rnorm(n=1000,0,0.25))

#x<-seq(from=0,to=1,by=0.001)
#V.sim<-calc_direct_V(phy,a.sims,sigma2_y.sims)
library(scales)
mypal <- pal_npg("nrc", alpha = 0.4)(2)

plot( NULL , xlim=c(0,1) , ylim=c(0,0.3) , xlab="Time since MRCA" , ylab="Covariance" ,cex.axis=0.75, mgp=c(1.25,0.25,0),tcl=-0.25)
for (i in 1:30){
  curve(calc_multiadaptive_cov_plot(a.sims[i,],sigma2_y.sims[i,],beta.sims[i,],x,Z_adaptive,n_reg) , add=TRUE , lwd=4 ,col=mypal[2]) #Prior - blue
}

for (i in 1:30){
  curve(calc_multiadaptive_cov_plot(post$a[i],post$sigma2_y[i],as.numeric(data.frame(post$beta)[i,]),x,Z_adaptive,n_reg) , add=TRUE , lwd=4 , col=mypal[1]) #Posterior - red
  #curve(calc_multiadaptive_cov(post$a[i],post$sigma2_y[i],beta,x,Z_adaptive,n_reg) , add=TRUE , lwd=4 , col=col.alpha(2,0.5))
}

#multiple traits


par(mar=c(3,3,0.25,0.25))
covariance.plot <- recordPlot()
dev.off()

#Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)

covariance.plot
########################################################################################################
#Four regimes - Slope plots for adaptive model with measurement error
#For main ms
library(ggsci)
library(rethinking)
X<-X_with_error
Y<-Y_with_error

optima.sims<-rnorm(100,1.89,1.5)
beta.sims<-rnorm(100, 0.5,0.25)

optima.post<-post$optima
beta.post<-data.frame(post$beta)
names(beta.post)<-c("post.beta.1","post.beta.2","post.beta.3","post.beta.4")


mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,2]}

mu.link.21<-function(x.seq){optima.post[,3]+x.seq*beta.post[,3]}
mu.link.22<-function(x.seq){optima.post[,4]+x.seq*beta.post[,4]}

x.seq <- seq(from=min(X), to=max(X) , length.out=100)
mu.11 <- sapply(x.seq , mu.link.11 )
mu.12 <- sapply(x.seq , mu.link.12 )
mu.21 <- sapply(x.seq , mu.link.21 )
mu.22 <- sapply(x.seq , mu.link.22 )


mu.mean.11<-colMeans(mu.11)
mu.mean.12<-colMeans(mu.12)
mu.mean.21<-colMeans(mu.21)
mu.mean.22<-colMeans(mu.22)


mu.mean.11<-data.frame(as.numeric(mu.mean.11))
mu.mean.12<-data.frame(as.numeric(mu.mean.12))
names(mu.mean.11)<-"mu.mean.11"
names(mu.mean.12)<-"mu.mean.12"

mu.mean.21<-data.frame(as.numeric(mu.mean.21))
mu.mean.22<-data.frame(as.numeric(mu.mean.22))
names(mu.mean.21)<-"mu.mean.21"
names(mu.mean.22)<-"mu.mean.22"



mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=PI , prob=0.89 )


mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)

mu.CI.21 <- apply( mu.21 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.22 <- apply( mu.22 , MARGIN=2, FUN=PI , prob=0.89 )

mu.CI.21<-data.frame(t(data.frame(mu.CI.21)),x.seq)
mu.CI.22<-data.frame(t(data.frame(mu.CI.22)),x.seq)


names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.21)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.22)<-c("min.5.5","max.94.5","x.seq")

#df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
df<-data.frame(Y=dat$Y_obs,X=dat$X_obs,Regimes=regimes_tip)
df11<-data.frame(x.seq,mu.mean.11)
df12<-data.frame(x.seq,mu.mean.12)
df21<-data.frame(x.seq,mu.mean.21)
df22<-data.frame(x.seq,mu.mean.22)

mypal <- pal_npg("nrc", alpha = 0.7)(length(beta))

#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.1<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=Regimes))+#,size=0.75)+
  geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.1)+ #Prior
  
  geom_abline(intercept=optima[1],slope=beta[1],alpha=0.5,linetype=2)+
  geom_abline(intercept=optima[2],slope=beta[2],alpha=0.5,linetype=2)+
  geom_abline(intercept=optima[3],slope=beta[3],alpha=0.5,linetype=2)+
  geom_abline(intercept=optima[4],slope=beta[4],alpha=0.5,linetype=2)+
  
  geom_line(data=df11,aes(x=x.seq,y=mu.mean.11),linetype=1)+
  geom_ribbon(data=mu.CI.11,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.2)+
  geom_line(data=df12,aes(x=x.seq,y=mu.mean.12),linetype=1)+
  geom_ribbon(data=mu.CI.12,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.2)+
  
  geom_line(data=df21,aes(x=x.seq,y=mu.mean.21),linetype=1)+
  geom_ribbon(data=mu.CI.21,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.2)+
  geom_line(data=df22,aes(x=x.seq,y=mu.mean.22),linetype=1)+
  geom_ribbon(data=mu.CI.22,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.2)+
  
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  # Right -> inside the plot area
  theme(
    legend.position = c(.8, .3),
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6)
  )+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("Y") + xlab("Adaptive trait")+
  scale_color_npg()

slope.plot.1

########################################################################################################
#For Main ms
fig<-ggarrange(hl.plot, vy.plot, "",slope.plot.1, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")

ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch - not online/For ms/Figures/Fig1.pdf", plot = fig, width=7, height=7 )

pdf(file = "/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch - not online/For ms/Figures/Fig1C.pdf",   # The directory you want to save the file in
    width = 3.6, # The width of the plot in inches
    height = 3.4) # The height of the plot in inches
covariance.plot
dev.off()

########################################################################################################
#Model Comparisons

########################################################################################################
#PSIS-LOO
#library(rethinking)
#compare(fit.reg.direct.ve,fit.reg.direct.vs,func=PSIS)
#compare(fit.reg.adapt.mlm.ve,fit.reg.adapt.mlm.ve.nc,fit.reg.adapt.mlm.vi,fit.reg.adapt.vs,fit.reg.adapt,func=PSIS)
#compare(fit.reg.adapt.mlm.ve,fit.reg.adapt.mlm.ve.nc,fit.reg.adapt.vs,func=PSIS)

library(loo) #Mlm varying effects model
loo_mlm_ve <- loo(fit.reg.adapt.mlm.ve, save_psis = TRUE)
print(loo_mlm_ve)
plot(loo_mlm_ve) #4X6
plot(loo_mlm_ve,label_points=TRUE) #Label outliers

#library(loo) #Varying effects model
loo_ve <- loo(fit.reg.adapt.ve, save_psis = TRUE)
print(loo_ve)
plot(loo_ve) #4X6
plot(loo_ve,label_points=TRUE) #Label outliers

loo_compare(loo_mlm_ve, loo_ve)
#elpd_diff se_diff
#model1  0.0       0.0   
#model2 -1.1       1.3   

-1.1 + c(-1,1)*1.4*1.96


#library(loo) #MLM varying intercepts model
loo_vi <- loo(fit.reg.adapt.mlm.vi, save_psis = TRUE)
print(loo_vi)
plot(loo_vi) #4X6

#library(loo) #Varying effects model - non-centered
#loo_ve_nc <- loo(fit.reg.adapt.mlm.ve.nc, save_psis = TRUE)
#print(loo_ve_nc)
#plot(loo_ve_nc)

#library(loo) #Basic model
loo_basic <- loo(fit.reg.adapt, save_psis = TRUE)
print(loo_basic)
plot(loo_basic)

plot(loo_basic,label_points=TRUE) #Label outliers

loo_compare(loo_mlm_ve, loo_ve, loo_vi, loo_basic)


########################################################################################################
########################################################################################################
#Bayes Factors
library(bridgesampling)
lml.fit.reg.adapt.mlm.ve<-bridge_sampler(fit.reg.adapt.mlm.ve,silent=TRUE)
#lml.fit.reg.adapt.mlm.ve.nc<-bridge_sampler(fit.reg.adapt.mlm.ve.nc,silent=TRUE,maxiter=5000)
lml.fit.reg.adapt.ve<-bridge_sampler(fit.reg.adapt.ve,silent=TRUE,maxiter=5000)
lml.fit.reg.adapt.mlm.vi<-bridge_sampler(fit.reg.adapt.mlm.vi,silent=TRUE)
lml.fit.reg.adapt<-bridge_sampler(fit.reg.adapt,silent=TRUE,maxiter=5000)

#bridgesampling::bf(lml.fit.reg.adapt.mlm.ve, lml.fit.reg.adapt.mlm.ve.nc)
bridgesampling::bf(lml.fit.reg.adapt.ve, lml.fit.reg.adapt.mlm.ve)
bridgesampling::bf(lml.fit.reg.adapt.mlm.ve, lml.fit.reg.adapt.mlm.vi)
bridgesampling::bf(lml.fit.reg.adapt.mlm.ve, lml.fit.reg.adapt)

########################################################################################################
#Traceplots #4X10
traceplot(fit.reg.adapt.mlm.ve,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta"))
traceplot(fit.reg.adapt.vs,pars = c(c("hl","vy","optima","beta","beta_e")))

########################################################################################################
#Prior predictive checks
#Based on Milestone 16 - mlm with varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_adapt_mlm_ve_priorpc.stan")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/chatgpt_blouchOU_reg_adapt_mlm_ve.stan")

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_priorpc.stan")
fit.reg.adapt.mlm.ve.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000, algorithm=c("Fixed_param"))
#print(fit.reg.adapt.mlm.ve.priorpc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
#plot(precis(fit.reg.adapt.mlm.ve.priorpc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.reg.adapt.mlm.ve.priorpc)
mypal <- pal_aaas("default", alpha = 1)(4)

#plot(post$Y_sim_obs[3,],dat$Y_obs)

df<-data.frame(Y=post$Y_sim_obs[3,],X=dat$Y_obs,Regimes=regimes_tip)

priorpc.plot<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=as.factor(Regimes)))+
  geom_abline(intercept=0,slope=1,alpha=0.05)+ #Prior
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  ggtitle("Prior Predictive Check")+
  ylab("Simulated Y") + xlab("True Y")+
  scale_color_manual(name="Regimes",values=mypal,labels=c('OU1', 'OU2', 'OU3', 'OU4'))


priorpc.plot


########################################################################################################
#Posterior predictive checks
#Based on Milestone 16 - mlm with varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_adapt_mlm_ve_postpc.stan")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/chatgpt_blouchOU_reg_adapt_mlm_ve.stan")

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_postpc.stan")
fit.reg.adapt.mlm.ve.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
#print(fit.reg.adapt.mlm.ve.postpc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e","Y_sim_obs"))
#plot(precis(fit.reg.adapt.mlm.ve.postpc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e","Y_sim_obs"))))
post<-extract(fit.reg.adapt.mlm.ve.postpc)

plot(post$Y_sim_obs[1,],dat$Y_obs)

mypal <- pal_aaas("default", alpha = 1)(4)

#plot(post$Y_sim_obs[3,],dat$Y_obs)

df<-data.frame(Y=post$Y_sim_obs[3,],X=dat$Y_obs,Regimes=regimes_tip)

priorpc.plot<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=as.factor(Regimes)))+
  geom_abline(intercept=0,slope=1,alpha=0.05)+ #Prior
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  ggtitle("Posterior Predictive Check")+
  ylab("Simulated Y") + xlab("True Y")+
  scale_color_manual(name="Regimes",values=mypal,labels=c('OU1', 'OU2', 'OU3', 'OU4'))


priorpc.plot


########################################################################################################

