#Blouch Model Validation Code
#Run after Simulate OU XY Data.R
#Using recoded versions of Blouch

library(rstan)
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


#Statistical Rethinking approach to calculating V/CV
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_SR.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_SR.stan")

fit.npi.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter = 4000,control=list(adapt_delta=0.80))

print(fit.npi.direct,pars = c("hl","alpha","beta","var_anc"))


#Using direct effect model from Hansen (1997) for V/CV matrix
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOU_direct.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_direct.stan")

fit.npi.direct<- rstan::sampling(object = stan_model,data = dat,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.npi.direct,pars = c("hl","vy","alpha","beta","sigma2_y"))







#Using Simulated Regime data from: SBR1 - Subsample Primate Regimes

#Basic direct stan model - no regimes
intercept_test<-rnorm(100,stan_sim_data$ols_intercept,0.1)
slope_test<-rnorm(100,stan_sim_data$ols_slope,0.1)

setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/blouchOU_v1_9.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_v1_9.stan")

fit.npi.direct<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.npi.direct,pars = c("hl","beta","beta_e","vy"))

#For downstream analysis and plots




Basic adaptive stan model - no regimes
Priors
intercept_test<-rnorm(100,stan_sim_data$ols_intercept,0.3)
slope_test<-rnorm(100,stan_sim_data$ols_slope,0.1)

```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/blouchOU_v1_9.stan") #Macbook Pro

stan_model <- stan_model("blouchOU_v1_9.stan")

fit.npi.adaptive<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.npi.adaptive,pars = c("hl","beta","beta_e","vy"))

#For downstream analysis and plots
```


Basic no varying effects regimes model
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/1. Standard MV Versions/blouchOUReg_v1_5.stan") #Macbook Pro

stan_model <- stan_model("blouchOUReg_v1_5.stan")

#setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_v1_6.stan")

#stan_model <- stan_model("blouchOUReg_v1_6.stan")

#setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")
#stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_v1_7.stan")

#stan_model <- stan_model("blouchOUReg_v1_7.stan")


fit.npi.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.npi.regimes,pars = c("hl","beta","beta_e","vy"))

#For downstream analysis and plots
```





Plotting
```{r}
library(rethinking)

stan_dens(fit.npi.regimes,pars = c("hl","beta","vy"))
trankplot(fit.npi.regimes,pars = c("hl","beta","vy"))
plot(precis(fit.npi.regimes,depth=1,pars = c("hl","beta","vy")))

#For downstream analysis and plots
ext.fit <- rstan::extract(fit.npi.regimes)


#Plot correlation between posterior distributions
library(ggplot2)
hl<-ext.fit$hl
vy<-ext.fit$vy
beta<-ext.fit$beta

new.data<-data.frame(hl,vy,beta)
ggplot(new.data,aes(x=hl,y=vy))+
  geom_point()

ggplot(new.data,aes(x=hl,y=beta))+
  geom_point()

pairs(fit.npi.regimes,pars=c("hl","beta","vy"))

```



Multilevel model - Varying Intercepts, centered priors
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Macbook Pro

#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")#Mac Studio
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Mac Studio

stan_model <- stan_model("blouchOUReg_MLI.stan")

fit.mli.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.mli.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy"))

```


Plotting
```{r}
stan_dens(fit.mli.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy"))

#For downstream analysis and plots
ext.fit.adaptive.test <- rstan::extract(fit.mli.regimes)

library(ggplot2)
hl<-ext.fit.adaptive.test$hl
vy<-ext.fit.adaptive.test$vy
beta<-ext.fit.adaptive.test$beta

new.data<-data.frame(hl,vy,beta)
ggplot(new.data,aes(x=hl,y=vy))+
  geom_point()

ggplot(new.data,aes(x=hl,y=beta))+
  geom_point()


pairs(fit.mli.regimes,pars=c("hl","beta","vy"))

```

Multilevel model - Varying Intercepts, non-centered priors
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI_nc.stan") #Macbook Pro

#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")#Mac Studio
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Mac Studio

stan_model <- stan_model("blouchOUReg_MLI_nc.stan")

fit.mli.nc.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.9))

print(fit.mli.nc.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy","log_lik"))

#For downstream analysis and plots
```



Multilevel model - Varying Effects, centered priors
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLM.stan") #Macbook Pro

#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")#Mac Studio
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Mac Studio

stan_model <- stan_model("blouchOUReg_MLM.stan")

fit.mli.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.80))

print(fit.mli.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy","log_lik"))

#For downstream analysis and plots
```


Plotting
```{r}
stan_dens(fit.mli.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy"))

#For downstream analysis and plots
ext.fit.adaptive.test <- rstan::extract(fit.mli.regimes)

library(ggplot2)
hl<-ext.fit.adaptive.test$hl
vy<-ext.fit.adaptive.test$vy
beta<-ext.fit.adaptive.test$beta

new.data<-data.frame(hl,vy,beta)
ggplot(new.data,aes(x=hl,y=vy))+
  geom_point()

ggplot(new.data,aes(x=hl,y=beta))+
  geom_point()


pairs(fit.mli.regimes,pars=c("hl","beta","vy"))

```

Multilevel model - Varying Effects, non-centered priors
```{r}
setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/") #Macbook Pro
stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI_nc.stan") #Macbook Pro

#setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/")#Mac Studio
#stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro/Documents/Academic/Research/Current Projects/Blouch project/R1 blouch-testing branch/Stan Models Milestones/Testing Versions/blouchOUReg_MLI.stan") #Mac Studio

stan_model <- stan_model("blouchOUReg_MLI_nc.stan")

fit.mli.nc.regimes<- rstan::sampling(object = stan_model,data = stan_sim_data,chains = 2,iter = 2000,control=list(adapt_delta=0.9))

print(fit.mli.nc.regimes,pars = c("hl","beta_bar","sigma","beta","beta_e","vy","log_lik"))

#For downstream analysis and plots
```

Model Comparisons
```{r}
log_lik_M.mli<-extract_log_lik(fit.mli.regimes, merge_chains = FALSE)
log_lik_M.npi<-extract_log_lik(fit.npi.regimes, merge_chains = FALSE)

log_ratios <- -1 * log_lik_M.mli
r_eff <- relative_eff(exp(-log_ratios))
psis_result <- psis(log_ratios, r_eff = r_eff)
str(psis_result)
plot(psis_result)


r_eff <- relative_eff(exp(log_lik_M.mli), cores = 2) 
loo_1 <- loo(log_lik_M.mli, r_eff = r_eff, cores = 2)

r_eff <- relative_eff(exp(log_lik_M.npi), cores = 2) 
loo_2 <- loo(log_lik_M.npi, r_eff = r_eff, cores = 2)

loo_compare(waic(log_lik_M.mli),waic(log_lik_M.npi))
```



