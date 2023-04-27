#Direct effect and random effect simulations

rm(list=ls())
#Script to simulate X an Y OU Data
#04/26/2023
calc_SR_V<-function(dmat, var_y_anc, a){ 
  N<-dim(dmat)[1]
  K<-matrix(NA,N,N)
  for (i in 1:(N-1)) {
    K[i,i]<-var_y_anc
    for (j in (i+1):N) {
      K[i, j]<-var_y_anc * exp(-a * dmat[i,j]) #Based on statistical rethinking Gaussian Process approach
      K[j, i]<-K[i, j]
    }
  }
  return(K)
}

calc_direct_V<-function(phy, sigma2_y, a){
  ts<-ts_fxn(phy)
  ta<-ts[[1]]
  tij<-ts[[2]] #Same as Cophenetic Distance matrix
  Vt<-sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) * exp(-a * tij)) #From Hansen (1997)
  #ta - time from root to tips, tij  - total time separating spcies
  #Variance of y/()
  return(Vt)
  }

ts_fxn<-function(phy){ #Calculate t
  n<-length(phy$tip.label)
  mrca1 <- ape::mrca(phy) #Node numbers for MRCA of tips
  times <- ape::node.depth.edgelength(phy) #Time from root to tips of each node, starting with the tips
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(phy$tip.label, phy$tip.label)) #Matrix with time from root to MRCA of pairs of tips - pulls out values of times that correspond with node numbers - integers
  #T.term <- times[1:n] #Times from root for tips
  #tia <- times[1:n] - ta #Times from root to tips - times from root to MRCA = times from MRCA to tips
  #tja <- t(tia) #Transpose of the times from MRCA to tips matrix
  #tij <- tja + tia #Sum of matrix and its transpose - total time separating species
  tij<-cophenetic(phy)
  #return(list(ta,tia,tja,tij,T.term))
  return(list(ta,tij))
}

  


library(MASS)
library(phytools)
library(ggplot2)
#Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997) Original Stan

tree.10K<-read.tree("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Original Submission/Blouch Testing/Phylogeny/10KPrimateTree.tre")
N<-50 
#set.seed(1) #Set seed to get same random species each time

phy <- keep.tip(tree.10K,sample(tree.10K$tip.label)[1:N]) 

l.tree<-max(branching.times(phy)) ## rescale tree to height 1
phy$edge.length<-phy$edge.length/l.tree 

tip.label<-phy$tip.label
Dmat<-cophenetic(phy) #Time separating tips, same as tij matrix in Slouch/Blouch code
#Dmat<-Dmat[tip.label,tip.label]/max(Dmat)

#Setup parameters
hl<-0.1
a<-log(2)/hl
#a<-6.931472 #Half life of 0.1, fast evoluion
sigma2_y<-1
vX0<-0
vY0 <- 0
vy<-0.1
Sxx<-10 #Look at effects
#var_y_anc<-1
ts<-ts_fxn(phy)
ta<-ts[[1]]
tij<-ts[[2]]

#var_y_anc<-sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)))
#var_y_anc<-var_y_anc[1,1]

#V<-calc_SR_V(Dmat,var_y_anc,a)
#SR_V<-var_y_anc*V

V<-calc_direct_V(phy,sigma2_y,a)

X<-fastBM(phy,a=vX0,sig2=Sxx,internal=FALSE) #Simulate X BM variable on tree, with scaling 10
phenogram(phy,X,spread.labels=TRUE,spread.cost=c(1,0)) #Plot X data

alpha<-4 #Intecept
beta<-0.25 #Slope

mu<-alpha+X*beta #Simulate mu for Y

#Simulate direct effect Y trait
Y<-mvrnorm(n=1,mu,V)

df<-data.frame(Y=Y,X=X)
names(df)<-c("Y","X")

ggplot(data=df,aes(x=X,y=Y))+
  geom_point()


summary(lm(Y~X,df))

#dat<-list(N=N,Z=1,Y=Y,X=X,Dmat=Dmat)
dat<-list(N=N,Z=1,Y=Y,X=X,ta=ta,tij=tij)
