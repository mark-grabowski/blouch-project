########################################################################################################
#Make simmap tree
#Setup dataset for info on model fits
cl<-parallel::detectCores()-2 #Cl for fitMk function
x<-cervid.trdata.BGS$dat$BGS
names(x)<-cervid.trdata.BGS$phy$tip.label
fit.ER<-fitMk(cervid.trdata.BGS$phy,x,model="ER")

fit<-fit.ER
fittedQ<-matrix(NA,length(fit$states),length(fit$states))
fittedQ[]<-c(0,fit$rates)[fit$index.matrix+1]
diag(fittedQ)<-0
diag(fittedQ)<--rowSums(fittedQ)
colnames(fittedQ)<-rownames(fittedQ)<-fit$states
tree<-cervid.trdata.BGS$phy
simmap_trees<-make.simmap(tree,x,Q=fittedQ,nsim=1)
#simmap_trees[[1]]$mapped.edge
simmap_trees$mapped.edge
simmap_tree1<-simmap_trees
#simmap_tree1<-simmap_trees[[1]]

#Simulate data using simmap
n<-length(simmap_tree1$tip.label)
mrca1 <- ape::mrca(simmap_tree1)
times <- ape::node.depth.edgelength(simmap_tree1)
ta <- matrix(times[mrca1], nrow=n, dimnames = list(simmap_tree1$tip.label, simmap_tree1$tip.label))
T.term <- times[1:n]
tia <- times[1:n] - ta
tja <- t(tia)
tij <- tja + tia

#regimes_internal <-cervid.trdata$phy$node.label
regimes_tip <- cervid.trdata.BGS$dat$BGS

#regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"simmap"
regimes<-colnames(simmap_tree1$mapped.edge)

lineages <- lapply(1:n, function(e) lineage.constructor(simmap_tree1, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

############################################################################################################
############################################################################################################
#Create SIMMAPs
#cl<-parallel::detectCores()-2 #Cl for fitMk function
num.cores=parallel::detectCores()-2 #For mclapply function
saved.discrete.models<-data.frame("Model"=NA,"AIC ER"=NA,"AIC SYM"=NA,"AIC ARD"=NA)

#########################################
tree<-cervid.trdata.BGS$phy
x<-as.vector(cervid.trdata.BGS$dat$BGS)
names(x)<-cervid.trdata.BGS$phy$tip.label

## fit model so we don't have to repeatedly recompute Q
fit.ER<-fitMk(cervid.trdata.BGS$phy,x,model="ER")
fit.SYM<-fitMk(cervid.trdata.BGS$phy,x,model="SYM")
fit.ARD<-fitMk(cervid.trdata.BGS$phy,x,model="ARD")

saved.discrete.models<-rbind(saved.discrete.models,c("BGS",AIC(fit.ER),AIC(fit.SYM),AIC(fit.ARD)))
test<-c(AIC(fit.ER),AIC(fit.SYM),AIC(fit.ARD))
num.min.AIC<-match(min(test),test)
fit<-list(fit.ER,fit.SYM,fit.ARD)[[num.min.AIC]]

fittedQ<-matrix(NA,length(fit$states),length(fit$states))
fittedQ[]<-c(0,fit$rates)[fit$index.matrix+1]
diag(fittedQ)<-0
diag(fittedQ)<--rowSums(fittedQ)
colnames(fittedQ)<-rownames(fittedQ)<-fit$states

loops.num<-1
sim.num<-1
trees<-mclapply(1:loops.num,function(n,tree,x,fixedQ) make.simmap(tree,x,Q=fixedQ,nsim=sim.num),
                tree=tree,x=x,fixedQ=fittedQ,mc.cores=num.cores)

trees<-do.call(c,trees) ## combine trees
if(!("multiSimmap"%in%class(trees))) class(trees)<-c("multiSimmap",class(trees))
BGS.trees<-trees
simmap_tree1<-BGS.trees[[1]]

############################################################################################################
#Single SIMMAP Tree
n<-length(simmap_tree1$tip.label)
mrca1 <- ape::mrca(simmap_tree1)
times <- ape::node.depth.edgelength(simmap_tree1)
ta <- matrix(times[mrca1], nrow=n, dimnames = list(simmap_tree1$tip.label, simmap_tree1$tip.label))
T.term <- times[1:n]
tia <- times[1:n] - ta
tja <- t(tia)
tij <- tja + tia

#regimes_internal <-cervid.trdata$phy$node.label
regimes_tip <- cervid.trdata.BGS$dat$BGS
regimes<-as.factor(simmap_tree1$node.label)
lineages <- lapply(1:n, function(e) lineage.constructor(simmap_tree, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

#regimes <- concat.factor(regimes_tip, regimes_internal)
anc_maps<-"simmap"
regimes<-colnames(simmap_tree1$mapped.edge)
lineages <- lapply(1:n, function(e) lineage.constructor(simmap_tree1, e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch



############################################################################################################
#Summarize SIMMAP Discrete Character Mappings based on highest posterior probabilities

summary(BGS.trees)
BGS.tree<-summary(BGS.trees)

regimes<- apply(BGS.tree$ace,1,function(e) colnames(BGS.tree$ace)[which.max(e)])

#shifts.total<-c(cervid.trdata.BGS$dat$BGS,cervid.trdata.BGS$phy$node.label)
edge.regimes <- factor(regimes[cervid.trdata.BGS$phy$edge[,2]])
print(edge.regimes)

reg.colors<-pal_aaas("default", alpha = 0.4)(length(unique(edge.regimes)))

reg_tips<-cervid.trdata.BGS$dat$BGS
reg_tips<-as.numeric(as.factor(reg_tips))

print(reg.colors)
plot(cervid.trdata.BGS$phy,edge.color = reg.colors[edge.regimes], edge.width = 3, cex = 0.2)

############################################################################################################
#Calculate lineages 

lineages <- lapply(1:n, function(e) lineage.constructor(BGS.trees[[1]], e, anc_maps, regimes, ace)) #Trace lineage from tips (n) to root and determine regimes of each node or branch

nodes<-NULL
store<-NULL
reg_num_lineage<-NULL
for(i in 1:length(lineages)){
  store<-c(store,length(lineage.nodes(cervid.trdata.BGS$phy,i))) #Calculate max node height
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

########################################################################################################
#Allometric model first
########################################################################################################
#Rescale tree height to 1
l.tree<-max(branching.times(cervid.trdata$phy))
cervid.trdata$phy$edge.length<-cervid.trdata$phy$edge.length/l.tree ## rescale tree to height 1
max(branching.times(cervid.trdata$phy))


tip.label<-cervid.trdata$phy$tip.label
phy<-cervid.trdata$phy
Dmat<-cophenetic(phy) #Time separating tips, same as tij matrix in Slouch/Blouch code

########################################################################################################
#Phylogeny info
n<-length(cervid.trdata$phy$tip.label)
mrca1 <- ape::mrca(cervid.trdata$phy)
times <- ape::node.depth.edgelength(cervid.trdata$phy)
ta <- matrix(times[mrca1], nrow=n, dimnames = list(cervid.trdata$phy$tip.label, cervid.trdata$phy$tip.label))
T.term <- times[1:n]
tia <- times[1:n] - ta
tja <- t(tia)
tij <- tja + tia

########################################################################################################
N<-length(tip.label)
Z<-1
Y_obs<-cervid.trdata$dat$log_ant_vol
X_obs<-cervid.trdata$dat$ms_log_psl
Y_error<-sqrt(cervid.trdata$dat$me.log_ant_vol) #Standard error not variance
X_error<-sqrt(cervid.trdata$dat$me.log_psl) #Standard error not variance

dat<-list(N=N,Z=Z,Y_obs=Y_obs,X_obs=matrix(X_obs,nrow=N,ncol=Z),Y_error=Y_error,X_error=matrix(X_error,nrow=N,ncol=Z),ta=ta,tij=tij)

########################################################################################################
#Prior Exploration Plot
lm.allometric<-summary(lm(dat$Y_obs~dat$X_obs))
lm.allometric$coefficients

#Prior vs. Posterior Plot
library(ggsci)
library(rethinking)

alpha.sims<-rnorm(100,lm.allometric$coefficients[1],0.75)
beta.sims<-rnorm(n=100,lm.allometric$coefficients[2],1.5)

df<-data.frame(Y=dat$Y_obs,X=dat$X_obs[,1])
names(df)<-c("Y","X")

slope.plot<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X))+
  geom_abline(intercept=alpha.sims,slope=beta.sims,alpha=0.15)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("log Antler Volume") + xlab("log Posterior Skull Length")+
  scale_color_npg()

slope.plot

########################################################################################################
#Allometric analysis using direct effect model with measurement error
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_direct.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_direct.stan")

fit.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter =2000, cores=2)

print(fit.direct,pars = c("hl","vy","alpha","beta"))
post<-extract(fit.direct)

########################################################################################################
#Priors vs. posterior plots
#Half-life Plot
hl.sims<-data.frame(rlnorm(n=1000,meanlog=log(0.25),sdlog=0.75))
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
#Vy plot
vy.sims<-rexp(n=1000,rate=1)
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
#Covariance plot
mypal <- pal_npg("nrc", alpha = 0.4)(2)

a.sims<-log(2)/hl.sims;
sigma2_y.sims<-vy.sims*(2*(log(2)/hl.sims));
#x<-seq(from=0,to=1,by=0.001)
#V.sim<-calc_direct_V(phy,a.sims,sigma2_y.sims)


plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="Time since MRCA" , ylab="Covariance" ,cex.axis=0.75, mgp=c(1.25,0.25,0),tcl=-0.25)
for (i in 1:30){
  curve(sigma2_y.sims[i,] /(2 * a.sims[i,]) * ((1 - exp(-2 * a.sims[i,] * (1-(x/2)))) * exp(-a.sims[i,] * x)) , add=TRUE , lwd=4 , col=mypal[2])
}

for (i in 1:30){
  curve(post$sigma2_y[i] /(2 * post$a[i]) * ((1 - exp(-2 * post$a[i] * (1-(x/2)))) * exp(-post$a[i] * x)) , add=TRUE , lwd=4 , col=mypal[1])
}

par(mar=c(3,3,0.25,0.25))
covariance.plot <- recordPlot()
dev.off()

#Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)

covariance.plot
########################################################################################################
#Prior vs. Posterior Plot

X<-dat$X_obs
Y<-dat$Y_obs

library(ggsci)
library(rethinking)

alpha.sims<-rnorm(100,lm.allometric$coefficients[1],0.75)
beta.sims<-rnorm(n=100,lm.allometric$coefficients[2],1.5)

alpha.post<-data.frame(post$alpha)
names(alpha.post)<-"post.alpha"

beta.post<-data.frame(post$beta)
names(beta.post)<-"post.beta"

mu.link<-function(x.seq){alpha.post+x.seq*beta.post}
x.seq <- seq(from=min(X)-0.1, to=max(X)+0.1 , length.out=100)
mu <- sapply(x.seq , mu.link )
mu.mean <-lapply( mu , mean )
mu.mean<-data.frame(as.numeric(mu.mean))
names(mu.mean)<-"mu.mean"

mu.CI <- lapply( mu , PI , prob=0.89 )
mu.CI<-data.frame(t(data.frame(mu.CI)),x.seq)
names(mu.CI)<-c("min.5.5","max.94.5","x.seq")

#df<-data.frame(Y=stan_sim_data$Y,X=stan_sim_data$direct_cov)
df<-data.frame(Y=Y,X=X)
names(df)<-c("Y","X")

df2<-data.frame(x.seq,mu.mean)

#slope.prior.plot<-ggplot(data=reg.trdata$dat,aes(y=Sim1,x=X))+
slope.plot<-ggplot()+  
  geom_point(data=df,aes(y=Y,x=X))+
  geom_abline(intercept=alpha.sims,slope=beta.sims,alpha=0.05)+
  #geom_abline(intercept=alpha,slope=beta,alpha=0.5,linetype=2)+
  geom_line(data=df2,aes(x=x.seq,y=mu.mean),linetype=1)+
  geom_ribbon(data=mu.CI,aes(x=x.seq,ymin=min.5.5,ymax=max.94.5),linetype=2,alpha=0.25)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #ggtitle("Prior vs. Posterior for Intercept and Slope")+
  ylab("log Antler Volume (l)") + xlab("log Posterior Skull Length (cm)")+
  scale_color_npg()

slope.plot
#5X5.5 Plot
########################################################################################################
#Save Plots
fig<-ggarrange(hl.plot, vy.plot, "",slope.plot, ncol=2,nrow=2, labels = c("A)","B)","C)","D)"),common.legend = TRUE,legend="top")

ggsave("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch - not online/For ms/Figures/antler_allometry.pdf", plot = fig, width=7, height=7 )

pdf(file = "/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch - not online/For ms/Figures/antler_allometry_cov.pdf",   # The directory you want to save the file in
    width = 3.57, # The width of the plot in inches
    height = 3.4) # The height of the plot in inches

covariance.plot
dev.off()

