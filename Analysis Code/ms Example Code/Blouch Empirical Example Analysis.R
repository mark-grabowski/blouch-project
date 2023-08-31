#Code for SBR1 Analysis of Cervidae Data
rm(list=ls())

calc_direct_V_plot<-function(phy, sigma2_y, a){
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  Vt<-sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) * exp(-a * tij)) #ta - time from root to tips, tij  - total time separating spcies
  return(Vt)
}

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

########################################################################################################
#Script to Run Empirical Analysis
library(devtools)
library(ape)
library(slouch)
library(rstan)
library(treeplyr)
library(ggplot2)
library(ggsci)
library(MASS)
library(ggpubr)
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
set.seed(10)

########################################################################################################
#For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

########################################################################################################
#Load Data

cervid.tree<-read.nexus("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch - not online/Sharable Data/cervidae_renamed.tre")
cervid.dataset<-read.csv("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch - not online/Sharable Data/updated_mean_dat_no_OL_BGS_MG.csv")

#Remove Muntiacus_atherodes, Elaphodus_cephalophus - rudamentary and female Rangifer. 
cervid.dataset<-filter(cervid.dataset,Genus_Species != "Sinomegaloceros_yabei" & ! Genus_Species =="Alces_alces_gigas"  & ! Genus_Species == "Muntiacus_truongsonensis" & ! Genus_Species ==  "Mazama_temama"& ! Genus_Species ==  "Muntiacus_feae" & ! Genus_Species ==  "Muntiacus_atherodes" & ! Genus_Species ==  "Elaphodus_cephalophus")


cervid.trdata <- make.treedata(cervid.tree, cervid.dataset,name_column="Genus_Species")

cervid.trdata<-filter(cervid.trdata,!(is.na(log_ant_vol)) & !(is.na(log_psl)))
cervid.trdata<-mutate(cervid.trdata,mc.log_psl=cervid.trdata$dat$log_psl-(mean(cervid.trdata$dat$log_psl)))

cervid.trdata<-mutate(cervid.trdata,me.log_ant_vol=cervid.trdata$dat$log_vol_var_est/cervid.trdata$dat$n)
cervid.trdata<-mutate(cervid.trdata,me.log_psl=cervid.trdata$dat$log_psl_var_est/cervid.trdata$dat$n)
cervid.trdata<-filter(cervid.trdata,!(is.na(me.log_ant_vol)) & !(is.na(me.log_psl)))

#Mean standardize predictor
cervid.trdata$dat$ms_log_psl<-cervid.trdata$dat$log_psl-mean(cervid.trdata$dat$log_psl)
cervid.trdata.BGS<-filter(cervid.trdata,!(is.na(BGS)))

########################################################################################################
#Main text analysis adding in BGS effect as adaptive predictor
########################################################################################################

#cervid.trdata.BGS
#Rescale tree height to 1
l.tree<-max(branching.times(cervid.trdata.BGS$phy))
cervid.trdata.BGS$phy$edge.length<-cervid.trdata.BGS$phy$edge.length/l.tree ## rescale tree to height 1
max(branching.times(cervid.trdata.BGS$phy))


tip.label<-cervid.trdata.BGS$phy$tip.label
phy<-cervid.trdata.BGS$phy
Dmat<-cophenetic(phy) #Time separating tips, same as tij matrix in Slouch/Blouch code

########################################################################################################
#Regime reconstruction usign ace from ape
#source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Simulation Code/Functions/set.converge.regimes.R") #Macbook Pro
reconstructed.BGS <- ace(as.factor(cervid.trdata.BGS$dat$BGS), cervid.trdata.BGS$phy, type = "d")

#Frugivore Diet Categories, group size categories, sociality
internal.regimes.BGS<- apply(reconstructed.BGS$lik.anc, 
                                     1, 
                                     function(e) colnames(reconstructed.BGS$lik.anc)[which.max(e)])

#Assign internal regimes to node.label on phylogeny - for use by Blouch prep functions
cervid.trdata.BGS$phy$node.label<-as.factor(internal.regimes.BGS)

#Check if manual setting code worked
shifts.total<-unlist(list(as.factor(cervid.trdata.BGS$dat$BGS),factor(internal.regimes.BGS)))
edge.regimes <- factor(shifts.total[cervid.trdata.BGS$phy$edge[,2]])
print(edge.regimes)

reg.colors<-pal_aaas("default", alpha = 1)(length(unique(edge.regimes)))

#For phylogeny figure - 6X6 portrait
plot(cervid.trdata.BGS$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.5)

########################################################################################################
#Simulate errors - for use with blouchOU_reg_direct
N<-length(cervid.trdata.BGS$phy$tip.label)
Z_direct<-1
Z_X_error<-1
Y_obs<-cervid.trdata.BGS$dat$log_ant_vol
X_obs<-cervid.trdata.BGS$dat$ms_log_psl
Y_error<-sqrt(cervid.trdata.BGS$dat$me.log_ant_vol) #Standard error not variance
X_error<-sqrt(cervid.trdata.BGS$dat$me.log_psl) #Standard error not variance
trdata<-cervid.trdata.BGS
regimes_tip <- cervid.trdata.BGS$dat$BGS

############################################################################################################
#Make trdata file
trdata$dat<-cbind(trdata$dat,data.frame(cbind(Y_obs,Y_error,X_obs,X_error)))
############################################################################################################
#Test Blouch prep code - Regimes + Adaptive Model
source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/R Setup Code/blouch.prep.R")
dat<-blouch.reg.direct.prep(trdata,"Y_obs","Y_error","X_obs","X_error",1,"BGS")

############################################################################################################
#Prior Exploration Plot
lm.allometric<-summary(lm(dat$Y_obs~dat$X_obs))
lm.allometric$coefficients

#Prior vs. Posterior Plot
library(ggsci)
library(rethinking)

alpha.sims<-rnorm(100,lm.allometric$coefficients[1],0.75)
beta.sims<-rnorm(n=100,lm.allometric$coefficients[2],1.75)

df<-data.frame(Y=dat$Y_obs,X=dat$X_obs[,1])
names(df)<-c("Y","X")

slope.plot<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=cervid.trdata.BGS$dat$BGS))+
  geom_abline(intercept=alpha.sims,slope=beta.sims,alpha=0.1)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("log Antler Volume") + xlab("log Posterior Skull Length")+
  scale_color_aaas()

slope.plot

############################################################################################################
#Milestone 14 - mlm with varying effects
#Combination of regime model with direct effect model with measurement error and correlated varying effects

#hl ~ lognormal(log(0.25),0.25);
#vy ~ exponential(5);
#optima_bar ~ normal(-1.179507,0.75);
#beta_bar ~ normal(6.304451,1.75);
#Rho ~ lkj_corr(4);
#sigma ~ normal(0,1);

#setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_direct_mlm_ve.stan")

#stan_model <- stan_model("blouchOU_reg_direct_mlm_ve.stan")
#fit.reg.direct.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000)
#print(fit.reg.direct.mlm.ve,pars = c("hl","vy","optima","beta","sigma","optima_bar","beta_bar"))
#plot(precis(fit.reg.direct.mlm.ve,depth=3,pars = c("hl","vy","optima","beta","Rho","sigma","optima_bar","beta_bar")))
#post<-extract(fit.reg.direct.mlm.ve)
#1: There were 176 divergent transitions after warmup. See

########################################################################################################
#Milestone 14 -mlm with varying effects, non-centered
#Combination of regime model with multitrait direct effect model with measurement error and correlated varying effects - non-centered

#Priors
#hl ~ lognormal(log(0.25),0.25);
#vy ~ exponential(10);
#L_Rho ~ lkj_corr_cholesky(2);
#sigma ~ normal(0,1);
#optima_bar ~ normal(-1.179507,0.75);
#beta_bar ~ normal(6.304451,1.75);
  
#Combination of regime model with multitrait direct effect model with measurement error and correlated varying effects - non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_direct_mlm_ve_nc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_direct_mlm_ve_nc.stan")
#fit.reg.direct.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000)
fit.reg.direct.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000,control = list(adapt_delta = 0.99))
print(fit.reg.direct.mlm.ve.nc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","Z"))
plot(precis(fit.reg.direct.mlm.ve.nc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta")))
post<-extract(fit.reg.direct.mlm.ve.nc)

########################################################################################################
#Milestone 12 - non-mlm varying effects
#Combination of regime model with multiple traits direct effect model with measurement error and varying slopes
#Priors
#hl ~ lognormal(log(0.25),0.25);
#vy ~ exponential(10);
#optima_beta ~ normal(-1.179507,0.75);
#beta ~ normal(6.304451,1.75);

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_direct_ve.stan")

stan_model <- stan_model("blouchOU_reg_direct_ve.stan")
fit.reg.direct.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000)
print(fit.reg.direct.ve,pars = c("hl","vy","optima","beta"))
plot(precis(fit.reg.direct.ve,depth=3,pars = c("hl","vy","optima","beta")))
post<-extract(fit.reg.direct.ve)

########################################################################################################
#Milestone 9 - mlm with varying intercepts
#Multilevel model - multilevel optima with direct effects predictor
#Priors
#hl ~ lognormal(log(0.25),0.25);
#vy ~ exponential(10);
#sigma ~ normal(0,1);
#optima_beta ~ normal(-1.179507,0.75);
#beta ~ normal(6.304451,1.5);

#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_direct_mlm_vi.stan")
##setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")

#stan_model <- stan_model("blouchOU_reg_direct_mlm_vi.stan")
#First run
#fit.reg.direct.mlm.vi<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =4000,cores=2)
#print(fit.reg.direct.mlm.vi,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
#plot(precis(fit.reg.direct.mlm.vi,depth=2,pars = c("hl","vy","optima","optima_bar","sigma")))

#post<-extract(fit.reg.direct.mlm.vi)
########################################################################################################
#Milestone 9 - mlm with varying intercepts - non-centered
#Multilevel model - multilevel optima with direct effects predictor
#Priors
#hl ~ lognormal(log(0.25),0.25);
#vy ~ exponential(10);
#sigma ~ normal(0,1);
#optima_beta ~ normal(-1.179507,0.75);
#beta ~ normal(6.304451,1.5);

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_direct_mlm_vi_nc.stan")

stan_model <- stan_model("blouchOU_reg_direct_mlm_vi_nc.stan")
#First run
fit.reg.direct.mlm.vi.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =4000,cores=2,control = list(adapt_delta = 0.99))
print(fit.reg.direct.mlm.vi.nc,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
plot(precis(fit.reg.direct.mlm.vi.nc,depth=2,pars = c("hl","vy","optima","optima_bar","sigma")))

post<-extract(fit.reg.direct.mlm.vi.nc)
########################################################################################################
#Milestone 5 - basic model with varying intercepts
#Regimes with direct effect model and measurement error
#Priors
#hl ~ lognormal(log(0.25),0.25);
#vy ~ exponential(10);
#optima ~ normal(-1.179507,0.75);
#beta ~ normal(6.304451,1.5);

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_direct.stan")
stan_model <- stan_model("blouchOU_reg_direct.stan")

fit.reg.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =4000,cores=2)

print(fit.reg.direct,pars = c("hl","vy","optima","beta"))
post<-extract(fit.reg.direct)
########################################################################################################
#Compare slopes
df<-data.frame(post$beta)
mean(df[,1]-df[,2])
PI(df[,1]-df[,2],prob=0.95)
sum((df[,1]-df[,2])<0)/4000
sum((df[,1]-df[,2])>=0)/4000
########################################################################################################
#Composite Figures
#Priors vs. posterior plots
hl.sims<-data.frame(rlnorm(n=1000,meanlog=log(0.25),sdlog=0.75))
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
  #geom_vline(xintercept=c(hl),linetype=2)+
  scale_fill_npg(name="",labels=c("Posterior","Prior"))

hl.plot
########################################################################################################

vy.sims<-rexp(n=1000,rate=10)
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
  #geom_vline(xintercept=c(vy),linetype=2)+
  
  #scale_fill_manual(labels=c("Posterior","Prior"))+
  scale_fill_npg(name="",labels=c("Posterior","Prior"))

vy.plot


########################################################################################################
a.sims<-log(2)/hl.sims
sigma2_y.sims<-vy.sims*(2*(log(2)/hl.sims))
#x<-seq(from=0,to=1,by=0.001)
#V.sim<-calc_direct_V(phy,a.sims,sigma2_y.sims)
library(scales)
mypal <- pal_npg("nrc", alpha = 0.4)(2)

#curve(sigma2_y.sims[i,] /(2 * a.sims[i,]) * ((1 - exp(-2 * a.sims[i,] * (1-(x/2)))) * exp(-a.sims[i,] * x)) , add=TRUE , lwd=4 , col=mypal[2]) #Prior - blue

plot( NULL , xlim=c(0,1) , ylim=c(0,0.4) , xlab="Time since MRCA" , ylab="Covariance" ,cex.axis=0.75, mgp=c(1.25,0.25,0),tcl=-0.25)
for (i in 1:30){
  curve(sigma2_y.sims[i,] /(2 * a.sims[i,]) * ((1 - exp(-2 * a.sims[i,] * (1-(x/2)))) * exp(-a.sims[i,] * x)) , add=TRUE , lwd=4 , col=mypal[2]) #Prior - blue
}

for (i in 1:30){
  curve(post$sigma2_y[i] /(2 * post$a[i]) * ((1 - exp(-2 * post$a[i] * (1-(x/2)))) * exp(-post$a[i] * x)) , add=TRUE , lwd=4 , col=mypal[1]) #Posterior - red
}

par(mar=c(3,3,0.25,0.25))
covariance.plot <- recordPlot()
dev.off()

#Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)

covariance.plot

########################################################################################################
#Varying intercepts
library(ggsci)
library(rethinking)
X<-dat$X_obs
Y<-dat$Y_obs

mypal <- pal_aaas("default", alpha = 1)(3)

optima.sims<-rnorm(100,-1.179507,0.75)
beta.sims<-rnorm(n=100,6.304451,1.5)

optima.post<-post$optima
beta.post<-data.frame(post$beta)
names(beta.post)<-c("post.beta.1")


mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,1]}
mu.link.13<-function(x.seq){optima.post[,3]+x.seq*beta.post[,1]}


x.seq <- seq(from=min(X), to=max(X) , length.out=100)
mu.11 <- sapply(x.seq , mu.link.11 )
mu.12 <- sapply(x.seq , mu.link.12 )
mu.13 <- sapply(x.seq , mu.link.13 )


mu.mean.11<-colMeans(mu.11)
mu.mean.12<-colMeans(mu.12)
mu.mean.13<-colMeans(mu.13)


mu.mean.11<-data.frame(as.numeric(mu.mean.11))
mu.mean.12<-data.frame(as.numeric(mu.mean.12))
mu.mean.13<-data.frame(as.numeric(mu.mean.13))

names(mu.mean.11)<-"mu.mean.11"
names(mu.mean.12)<-"mu.mean.12"
names(mu.mean.13)<-"mu.mean.13"


mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.13 <- apply( mu.13 , MARGIN=2, FUN=PI , prob=0.89 )


mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)
mu.CI.13<-data.frame(t(data.frame(mu.CI.13)),x.seq)


names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.13)<-c("min.5.5","max.94.5","x.seq")

#df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
df<-data.frame(Y=dat$Y_obs,X=dat$X_obs,Regimes=regimes_tip)
df11<-data.frame(x.seq,mu.mean.11)
df12<-data.frame(x.seq,mu.mean.12)
df13<-data.frame(x.seq,mu.mean.13)

#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.1<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=as.factor(Regimes)))+
  geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.05)+ #Prior
  
  geom_ribbon(data=mu.CI.11,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df11,aes(x=x.seq,y=mu.mean.11),linetype=1,linewidth=1,alpha=0.3,color=mypal[1])+
  geom_ribbon(data=mu.CI.12,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df12,aes(x=x.seq,y=mu.mean.12),linetype=1,linewidth=1,alpha=0.3,color=mypal[2])+
  geom_ribbon(data=mu.CI.13,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df13,aes(x=x.seq,y=mu.mean.13),linetype=1,linewidth=1,alpha=0.3,color=mypal[3])+
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("log Antler Volume (l)") + xlab("log Posterior Skull Length (cm)")+
  scale_color_manual(name="Breeding Group \nSize",values=mypal,labels=c('1-2', '3-5', '>5'))
#scale_color_discrete(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))+
#scale_color_aaas()

slope.plot.1
#Label points

#5X5.5 Plot
########################################################################################################
#Varying slopes/effects models
library(ggsci)
library(rethinking)
X<-dat$X_obs
Y<-dat$Y_obs

mypal <- pal_aaas("default", alpha = 1)(3)

optima.sims<-rnorm(100,-1.179507,0.75)
beta.sims<-rnorm(n=100,6.304451,1.75)

optima.post<-post$optima
beta.post<-data.frame(post$beta)
names(beta.post)<-c("post.beta.1","post.beta.2","post.beta.3")


mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,2]}
mu.link.13<-function(x.seq){optima.post[,3]+x.seq*beta.post[,3]}


x.seq <- seq(from=min(X), to=max(X) , length.out=100)
mu.11 <- sapply(x.seq , mu.link.11 )
mu.12 <- sapply(x.seq , mu.link.12 )
mu.13 <- sapply(x.seq , mu.link.13 )


mu.mean.11<-colMeans(mu.11)
mu.mean.12<-colMeans(mu.12)
mu.mean.13<-colMeans(mu.13)


mu.mean.11<-data.frame(as.numeric(mu.mean.11))
mu.mean.12<-data.frame(as.numeric(mu.mean.12))
mu.mean.13<-data.frame(as.numeric(mu.mean.13))

names(mu.mean.11)<-"mu.mean.11"
names(mu.mean.12)<-"mu.mean.12"
names(mu.mean.13)<-"mu.mean.13"


mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.13 <- apply( mu.13 , MARGIN=2, FUN=PI , prob=0.89 )


mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)
mu.CI.13<-data.frame(t(data.frame(mu.CI.13)),x.seq)


names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.13)<-c("min.5.5","max.94.5","x.seq")

#df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
df<-data.frame(Y=dat$Y_obs,X=dat$X_obs,Regimes=regimes_tip)
df11<-data.frame(x.seq,mu.mean.11)
df12<-data.frame(x.seq,mu.mean.12)
df13<-data.frame(x.seq,mu.mean.13)

#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.1<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=as.factor(Regimes)))+
  geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.05)+ #Prior
  
  geom_ribbon(data=mu.CI.11,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df11,aes(x=x.seq,y=mu.mean.11),linetype=1,linewidth=1,alpha=0.5,color=mypal[1])+
  geom_ribbon(data=mu.CI.12,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df12,aes(x=x.seq,y=mu.mean.12),linetype=1,linewidth=1,alpha=0.5,color=mypal[2])+
  geom_ribbon(data=mu.CI.13,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df13,aes(x=x.seq,y=mu.mean.13),linetype=1,linewidth=1,alpha=0.5,color=mypal[3])+
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("log Antler Volume (l)") + xlab("log Posterior Skull Length (cm)")+
  scale_color_manual(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))
  #scale_color_discrete(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))+
    #scale_color_aaas()

slope.plot.1
########################################################################################################
#For Fig. 2
#Varying slopes/effects models
library(ggsci)
library(rethinking)
X<-dat$X_obs
Y<-dat$Y_obs

mypal <- pal_aaas("default", alpha = 1)(3)

optima.sims<-rnorm(100,-1.179507,0.75)
beta.sims<-rnorm(n=100,6.304451,1.75)

optima.post<-post$optima
beta.post<-data.frame(post$beta)
names(beta.post)<-c("post.beta.1","post.beta.2","post.beta.3")


mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,2]}
mu.link.13<-function(x.seq){optima.post[,3]+x.seq*beta.post[,3]}


x.seq <- seq(from=min(X), to=max(X) , length.out=100)
mu.11 <- sapply(x.seq , mu.link.11 )
mu.12 <- sapply(x.seq , mu.link.12 )
mu.13 <- sapply(x.seq , mu.link.13 )


mu.mean.11<-colMeans(mu.11)
mu.mean.12<-colMeans(mu.12)
mu.mean.13<-colMeans(mu.13)


mu.mean.11<-data.frame(as.numeric(mu.mean.11))
mu.mean.12<-data.frame(as.numeric(mu.mean.12))
mu.mean.13<-data.frame(as.numeric(mu.mean.13))

names(mu.mean.11)<-"mu.mean.11"
names(mu.mean.12)<-"mu.mean.12"
names(mu.mean.13)<-"mu.mean.13"


mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.13 <- apply( mu.13 , MARGIN=2, FUN=PI , prob=0.89 )


mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)
mu.CI.13<-data.frame(t(data.frame(mu.CI.13)),x.seq)


names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.13)<-c("min.5.5","max.94.5","x.seq")

#df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
df<-data.frame(Y=dat$Y_obs,X=dat$X_obs,Regimes=regimes_tip)
df11<-data.frame(x.seq,mu.mean.11)
df12<-data.frame(x.seq,mu.mean.12)
df13<-data.frame(x.seq,mu.mean.13)

#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.1<-ggplot()+  
  geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.04)+ #Prior
  
  geom_ribbon(data=mu.CI.11,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df11,aes(x=x.seq,y=mu.mean.11),linetype=1,linewidth=1,alpha=0.75,color=mypal[1])+
  geom_ribbon(data=mu.CI.12,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df12,aes(x=x.seq,y=mu.mean.12),linetype=1,linewidth=1,alpha=0.75,color=mypal[2])+
  geom_ribbon(data=mu.CI.13,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.1)+
  geom_line(data=df13,aes(x=x.seq,y=mu.mean.13),linetype=1,linewidth=1,alpha=0.75,color=mypal[3])+
  geom_point(data=df,aes(y=Y,x=X,color=as.factor(Regimes)))+
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("log Antler Volume (l)") + xlab("log Posterior Skull Length (cm)")+
  scale_color_manual(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))
#scale_color_discrete(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))+
#scale_color_aaas()

slope.plot.1

########################################################################################################
#Labels for points
#Varying slopes/effects models
library(ggsci)
library(rethinking)
X<-dat$X_obs
Y<-dat$Y_obs

mypal <- pal_aaas("default", alpha = 1)(3)

optima.sims<-rnorm(100,-1.179507,0.75)
beta.sims<-rnorm(n=100,6.304451,1.75)

optima.post<-post$optima
beta.post<-data.frame(post$beta)
names(beta.post)<-c("post.beta.1","post.beta.2","post.beta.3")


mu.link.11<-function(x.seq){optima.post[,1]+x.seq*beta.post[,1]}
mu.link.12<-function(x.seq){optima.post[,2]+x.seq*beta.post[,2]}
mu.link.13<-function(x.seq){optima.post[,3]+x.seq*beta.post[,3]}


x.seq <- seq(from=min(X), to=max(X) , length.out=100)
mu.11 <- sapply(x.seq , mu.link.11 )
mu.12 <- sapply(x.seq , mu.link.12 )
mu.13 <- sapply(x.seq , mu.link.13 )


mu.mean.11<-colMeans(mu.11)
mu.mean.12<-colMeans(mu.12)
mu.mean.13<-colMeans(mu.13)


mu.mean.11<-data.frame(as.numeric(mu.mean.11))
mu.mean.12<-data.frame(as.numeric(mu.mean.12))
mu.mean.13<-data.frame(as.numeric(mu.mean.13))

names(mu.mean.11)<-"mu.mean.11"
names(mu.mean.12)<-"mu.mean.12"
names(mu.mean.13)<-"mu.mean.13"


mu.CI.11 <- apply( mu.11 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.12 <- apply( mu.12 , MARGIN=2, FUN=PI , prob=0.89 )
mu.CI.13 <- apply( mu.13 , MARGIN=2, FUN=PI , prob=0.89 )


mu.CI.11<-data.frame(t(data.frame(mu.CI.11)),x.seq)
mu.CI.12<-data.frame(t(data.frame(mu.CI.12)),x.seq)
mu.CI.13<-data.frame(t(data.frame(mu.CI.13)),x.seq)


names(mu.CI.11)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.12)<-c("min.5.5","max.94.5","x.seq")
names(mu.CI.13)<-c("min.5.5","max.94.5","x.seq")

#df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
df<-data.frame(Y=dat$Y_obs,X=dat$X_obs,Regimes=regimes_tip,Names=cervid.trdata.BGS$phy$tip.label)
df11<-data.frame(x.seq,mu.mean.11)
df12<-data.frame(x.seq,mu.mean.12)
df13<-data.frame(x.seq,mu.mean.13)



#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot.1<-ggplot()+  
  geom_text(data=df,aes(y=Y,x=X,color=as.factor(Regimes),label=Names))+
  geom_abline(intercept=optima.sims,slope=beta.sims,alpha=0.04)+ #Prior
  
  geom_ribbon(data=mu.CI.11,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.125)+
  geom_line(data=df11,aes(x=x.seq,y=mu.mean.11),linetype=1,linewidth=1,alpha=0.5,color=mypal[1])+
  geom_ribbon(data=mu.CI.12,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.125)+
  geom_line(data=df12,aes(x=x.seq,y=mu.mean.12),linetype=1,linewidth=1,alpha=0.5,color=mypal[2])+
  geom_ribbon(data=mu.CI.13,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.125)+
  geom_line(data=df13,aes(x=x.seq,y=mu.mean.13),linetype=1,linewidth=1,alpha=0.5,color=mypal[3])+
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("log Antler Volume (l)") + xlab("log Posterior Skull Length (cm)")+
  scale_color_manual(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))
#scale_color_discrete(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))+
#scale_color_aaas()

slope.plot.1

########################################################################################################
#Save Plots
fig<-ggarrange(hl.plot, vy.plot, "",slope.plot.1, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")

ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch - not online/For ms/Figures/antler_BGS.pdf", plot = fig, width=7, height=7 )

pdf(file = "/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch - not online/For ms/Figures/antler_BGS_cov.pdf",   # The directory you want to save the file in
    width = 3.57, # The width of the plot in inches
    height = 3.4) # The height of the plot in inches

covariance.plot
dev.off()


############################################################################################################
#Model comparison #PSIS
#loo_compare(loo_mlm_ve, loo_ve, loo_vi, loo_basic)


#library(loo) #Varying effects model
#loo_mlm_ve <- loo(fit.reg.direct.mlm.ve, save_psis = TRUE)
#print(loo_mlm_ve)
#plot(loo_mlm_ve)

library(loo) #MLM Varying effects model - non-centered
loo_mlm_ve_nc <- loo(fit.reg.direct.mlm.ve.nc, save_psis = TRUE)
print(loo_mlm_ve_nc)
plot(loo_mlm_ve_nc)
plot(loo_mlm_ve_nc,label_points=TRUE)

#library(loo) #Non-mlm varying effects model
loo_ve <- loo(fit.reg.direct.ve, save_psis = TRUE)
print(loo_ve)
plot(loo_ve)
plot(loo_ve,label_points=TRUE)

#library(loo) #MLM Varying intercepts model
loo_mlm_vi_nc <- loo(fit.reg.direct.mlm.vi.nc, save_psis = TRUE)
print(loo_mlm_vi_nc)
plot(loo_mlm_vi_nc)
plot(loo_mlm_vi_nc,label_points=TRUE)

#library(loo) #Basic model
loo_basic <- loo(fit.reg.direct, save_psis = TRUE)
print(loo_basic)
plot(loo_basic)
plot(loo_basic,label_points=TRUE)

loo_compare(loo_mlm_ve_nc, loo_ve, loo_mlm_vi_nc, loo_basic)
-0.8 + c(-1,1)*0.9*-1.86

############################################################################################################
#Bayes Factors
library(bridgesampling)
#lml.fit.reg.direct.mlm.ve<-bridge_sampler(fit.reg.direct.mlm.ve,silent=TRUE)
lml.fit.reg.direct.mlm.ve.nc<-bridge_sampler(fit.reg.direct.mlm.ve.nc,silent=TRUE)
#lml.fit.reg.direct.mlm.vi<-bridge_sampler(fit.reg.direct.mlm.vi,silent=TRUE)
lml.fit.reg.direct.mlm.vi.nc<-bridge_sampler(fit.reg.direct.mlm.vi.nc,silent=TRUE)
lml.fit.reg.direct.ve<-bridge_sampler(fit.reg.direct.ve,silent=TRUE)
lml.fit.reg.direct<-bridge_sampler(fit.reg.direct,silent=TRUE)

#bridgesampling::bf(lml.fit.reg.direct.mlm.ve, lml.fit.reg.direct.mlm.ve.nc)
bridgesampling::bf(lml.fit.reg.direct.ve, lml.fit.reg.direct.mlm.ve.nc)
bridgesampling::bf(lml.fit.reg.direct.ve, lml.fit.reg.direct.mlm.vi.nc)
bridgesampling::bf(lml.fit.reg.direct.ve, lml.fit.reg.direct)

############################################################################################################
#Trankplots/Traceplots 4X10
#Best model - multilevel model, non-centered varying effects
traceplot(fit.reg.direct.mlm.ve.nc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta"))
traceplot(fit.reg.direct.ve,pars = c("hl","vy","optima_beta","beta"))
traceplot(fit.reg.direct.mlm.vi.nc,pars = c("hl","vy","optima","optima_bar","beta","sigma"))
traceplot(fit.reg.direct,pars = c("hl","vy","optima","beta"))


trankplot(fit.reg.direct.mlm.ve.nc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta"))
trankplot(fit.reg.direct.vs,pars = c(c("hl","vy","optima_beta","beta")))
############################################################################################################

#Prior predictive checks on best model
########################################################################################################
#Milestone 14 -mlm with varying effects, non-centered - prior predictive check version
#Combination of regime model with multitrait direct effect model with measurement error and correlated varying effects - non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_direct_mlm_ve_nc_priorpc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_direct_mlm_ve_nc_priorpc.stan")
fit.reg.direct.mlm.ve.nc.priorpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000, algorithm=c("Fixed_param"))

#print(fit.reg.direct.mlm.ve.nc.priorpc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta"))
#plot(precis(fit.reg.direct.mlm.ve.nc.priorpc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta")))
post<-extract(fit.reg.direct.mlm.ve.nc.priorpc)
mypal <- pal_aaas("default", alpha = 1)(3)

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
  scale_color_manual(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))


priorpc.plot
############################################################################################################
#Posterior predictive checks on best model
########################################################################################################
#Milestone 14 -mlm with varying effects, non-centered - prior predictive check version
#Combination of regime model with multitrait direct effect model with measurement error and correlated varying effects - non-centered
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Validation Code/Model Checking/blouchOU_reg_direct_mlm_ve_nc_postpc.stan")
#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_reg.stan")

stan_model <- stan_model("blouchOU_reg_direct_mlm_ve_nc_postpc.stan")
fit.reg.direct.mlm.ve.nc.postpc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000)

post<-extract(fit.reg.direct.mlm.ve.nc.postpc)
mypal <- pal_aaas("default", alpha = 1)(3)

#plot(post$Y_sim_obs[3,],dat$Y_obs)

df<-data.frame(Y=post$Y_sim_obs[3,],X=dat$Y_obs,Regimes=regimes_tip)

postpc.plot<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X,color=as.factor(Regimes)))+
  geom_abline(intercept=0,slope=1,alpha=0.05)+ #Prior
  
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  ggtitle("Posterior Predictive Check")+
  ylab("Simulated Y") + xlab("True Y")+
  scale_color_manual(name="Breeding Group Size",values=mypal,labels=c('1-2', '3-5', '>5'))


postpc.plot

