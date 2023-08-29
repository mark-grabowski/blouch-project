#Blouch testing code
#1. log likelihood testing

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
hl.sims<-data.frame(rlnorm(n=1000,meanlog=log(0.25),sdlog=0.5))
#hl.sims<-data.frame(hl.sims[hl.sims<3])
names(hl.sims)<-"prior.hl.sims"


hl.plot<-ggplot()+
  geom_density(aes(prior.hl.sims,fill="prior.hl.sims"),alpha=0.2,data=hl.sims)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #labs(title="Prior vs. Posterior Distribution ",x="Half-life", y = "Density")+
  labs(title="",x="Half-life", y = "Density")+
  
  #scale_fill_manual(labels=c("Posterior","Prior"))+
  geom_vline(xintercept=c(hl),linetype=2)+
  scale_fill_npg(name="",labels=c("Prior"))

hl.plot
########################################################################################################
vy.sims<-rexp(n=1000,rate=30)
vy.sims<-data.frame(vy.sims)
names(vy.sims)<-"prior.vy.sims"




vy.plot<-ggplot()+
  geom_density(aes(prior.vy.sims,fill="prior.vy.sims"),alpha=0.2,data=vy.sims)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  
  #labs(title="Prior vs. Posterior Distribution ",x="vy", y = "Density")+
  labs(title="",x="vy", y = "Density")+
  geom_vline(xintercept=c(vy),linetype=2)+
  
  #scale_fill_manual(labels=c("Posterior","Prior"))+
  scale_fill_npg(name="",labels=c("Prior"))

vy.plot
#Testing Simulation Example - Adaptive MLM model
########################################################################################################
#Milestone 16 - mlm with varying effects
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_reg_adapt_mlm_ve_logL.stan")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/chatgpt_blouchOU_reg_adapt_mlm_ve.stan")

stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_logL.stan")
fit.reg.adapt.mlm.ve<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =2000)
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

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Finished Versions/blouchOU_reg_adapt_mlm_ve_nc.stan")
stan_model <- stan_model("blouchOU_reg_adapt_mlm_ve_nc.stan")
fit.reg.adapt.mlm.ve.nc<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter =4000)
print(fit.reg.adapt.mlm.ve.nc,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e"))
plot(precis(fit.reg.adapt.mlm.ve.nc,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.reg.adapt.mlm.ve.nc)
########################################################################################################
#Hl Plot prior vs. posterior - assume posterior has been extracted using extract(model) and stored in post

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
#https://mc-stan.org/loo/articles/loo2-with-rstan.html
library(loo) #Varying effects model - non-centered
log_lik_1 <- extract_log_lik(fit.reg.adapt.mlm.ve, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 2) 
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
plot(loo_1)
comp <- loo_compare(loo_1, loo_2)

log_ratios <- -log_lik_1
psis_result <- psis(log_ratios)
plot(psis_result, label_points = TRUE)


########################################################################################################
#Testing Simple direct effect model - no regimes
#Run code from Simulate Direct OU XY Data.R to simulate data
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/")
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/blouch/Stan Models Milestones/Testing Versions/blouchOU_direct_logL.stan")
stan_model <- stan_model("blouchOU_direct_logL.stan")
fit.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,cores=2,iter = 2000)
print(fit.direct,pars = c("hl","vy","alpha","beta"))
#plot(precis(fit.direct,depth=3,pars = c("hl","vy","optima_bar","beta_bar","Rho","sigma","optima","beta","beta_e")))
post<-extract(fit.direct)

########################################################################################################

library(loo) #Varying effects model - non-centered
log_lik_1 <- extract_log_lik(fit.direct, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 2) 
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1)
plot(loo_1)
comp <- loo_compare(loo_1, loo_2)

log_ratios <- -log_lik_1
psis_result <- psis(log_ratios)
plot(psis_result, label_points = TRUE)
