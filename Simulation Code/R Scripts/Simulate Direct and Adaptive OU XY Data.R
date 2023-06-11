#Simulate both direct effect and adaptive model as predictors
#Generative Code

rm(list=ls())
#Script to simulate X an Y OU Data
#04/26/2023
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
calc_adaptive_V<-function(a, sigma2_y, ta,  tij,  tja,  T_term,  beta,  sigma2_x){
  N<-dim(ta)[1];
  Z<-length(beta);
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  T_term<-ts[[3]]
  tja<-ts[[4]]
  ones<-rep(1,Z)
  ti<-matrix(T_term,length(T_term),N);
  if(Z==1){var_opt<-beta^2*sigma2_x[1,1]}
  else{var_opt<-as.numeric(t(beta)%*% sigma2_x %*% ones)}
  term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) * exp(-a * tij)
  term1 = (1 - exp(-a * ti)) / (a * ti)
  term2 = exp(-a * tja) * (1 - exp(-a * ti)) / (a * ti)
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

########################################################################################################
#Basic Setup
library(geiger)
library(MASS)
library(phytools)
library(ggplot2)
#Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997) Original Stan

tree.10K<-read.tree("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
#tree.10K<-read.tree("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
N<-50 #Number of species
set.seed(10) #Set seed to get same random species each time

phy <- keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N]) 
phy<-multi2di(phy)

l.tree<-max(branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree 

tip.label<-phy$tip.label
Dmat<-cophenetic(phy) #Time separating tips, same as tij matrix in Slouch/Blouch code
#Dmat<-Dmat[tip.label,tip.label]/max(Dmat)


########################################################################################################
#Direct effect and Adaptive Model
#Setup parameters
hl<-0.75 #0.1, 0.25, 0.75 - testing options
a<-log(2)/hl
vy<-0.1 #0.25,0.5 - testing options
sigma2_y<-vy*(2*(log(2)/hl));

vX0<-0
vY0 <- 0

ts<-ts_fxn(phy)
ta<-ts[[1]]
tij<-ts[[2]]
T_term<-ts[[3]]
tja<-ts[[4]]
Z_direct<-1
Z_adaptive<-1
Z<-Z_direct+Z_adaptive
sigma2_x<-matrix(1,1,1)

Xd<-rnorm(N,0,1)
names(Xd)<-phy$tip.label
phenogram(phy,Xd,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
Xa<-fastBM(phy,a=vX0,sig2=sigma2_x[1,1],internal=FALSE) #Simulate X BM variable on tree, with scaling 10
phenogram(phy,Xa,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data
#sigma2_x<-ratematrix(phy,Xa) #Calculate evolutionary v/cv matrix
Xs<-cbind(Xd,Xa)

alpha<-2 #Intecept
beta<-c(0.35,0.25) #Slopes
dmX<-calc_mixed_dmX(a,T_term,Xs,Z_direct,Z_adaptive)
mu<-alpha+dmX%*%beta #Simulate mu for Y

V<-calc_adaptive_V(a, sigma2_y, ta,  tij,  tja,  T_term,  beta[(Z_direct+1):(Z_adaptive+Z_direct)],  sigma2_x)
Y<-mvrnorm(n=1,mu,V)

df<-data.frame(Y=Y,X=Xs)

ggplot(data=df,aes(x=X.Xd,y=X.Xa))+
  geom_point()
summary(lm(X.Xa~X.Xd,df))

ggplot(data=df,aes(x=X.Xd,y=Y))+
  geom_point()

summary(lm(Y~X.Xd,df))

ggplot(data=df,aes(x=X.Xa,y=Y))+
  geom_point()

summary(lm(Y~X.Xa,df))

###############################
#Adaptive model w/o ME- blouchOU_adaptive.stan
dat<-list(N=N,Z_direct=Z_direct,Z_adaptive=Z_adaptive,Y_obs=Y,X_obs=matrix(Xs,nrow=N,ncol=Z),ta=ta,tij=tij,tja=tja,T_term=T_term,sigma2_x=sigma2_x)

#Include only direct effect variable
dat<-list(N=N,Z_direct=Z_direct,Z_adaptive=0,Y_obs=Y,X_obs=matrix(Xd,nrow=N,ncol=Z_direct),ta=ta,tij=tij,tja=tja,T_term=T_term,sigma2_x=matrix(0,0,0))

#Include only adaptive variable
dat<-list(N=N,Z_direct=0,Z_adaptive=1,Y_obs=Y,X_obs=matrix(Xa,nrow=N,ncol=Z_adaptive),ta=ta,tij=tij,tja=tja,T_term=T_term,sigma2_x=sigma2_x)


########################################################################################################
#Simulate errors
Z_X_error<-2 #Number of X traits with error
X_error<-matrix(0.01,nrow=N,ncol=Z_X_error)
Y_error<-rep(0.01,N)
Y_with_error<-Y+rnorm(N,0,0.01)
X_with_error<-apply(Xs,2,function(X){X+rnorm(N,0,0.01)})

#Adaptive model w/ Statistical Rethinking ME Correction - blouchOU_adaptive_ME_SR.stan
dat<-list(N=N,Z_direct=Z_direct,Z_adaptive=Z_adaptive,Z_X_error=Z_X_error,Y_obs=Y_with_error,X_obs=matrix(X_with_error,nrow=N,ncol=Z),Y_error=Y_error,X_error=X_error,ta=ta,tij=tij,tja=tja,T_term=T_term,sigma2_x=sigma2_x)
########################################################################################################


