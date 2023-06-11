#SBR1 OU Simulation Code
#Works with blouchOU_v1_7.stan and blouchOU.setup.mv_v1_3.R


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

stan.OU.data.setup<-function(OUBMparameters,tree,n){
  #source("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOU.setup.mv_v1_3.R")
  #names.traits<-c("trait_1","trait_2",NA,NA)
  name.response.trait<-c("trait_1")
  names.random.traits<-c("trait_2")
  
  stan.data<-list()
  for(i in 1:n){
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X
    #stan.data[[length(stan.data)+1]]<-blouchOU.setup.v1(OUBM.trdata,names.traits)
    #stan.data<-blouchOU.setup.mv(OUBM.trdata,name.response.trait,names.direct.traits=NULL,names.random.traits)
    
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.mv(OUBM.trdata,name.response.trait,names.direct.traits=NULL,names.random.traits)
    
    
    }
  return(stan.data)
}

stan.BM.data.setup<-function(BMparameters,tree,n){
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/Blouch Setup Files/blouchOU.setup.v1.R")
  #n<-50
  names.traits<-c("trait_1","trait_2",NA,NA)
  stan.data<-list()
  for(i in 1:n){
    BMdata<-simulBMProcPhylTree(cervid.trdata$phy,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
    BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
    BM.trdata <- make.treedata(tree, BMdata,name_column="tip.names") #Combine data with 
    #BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.v1(BM.trdata,names.traits)
  }
  return(stan.data)
}

stan.OU.simulation.reps<-function(OUBMparameters,tree,n){
  #setwd("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  #stanc("/Users/markgrabowski/Library/CloudStorage/GoogleDrive-mark.walter.grabowski@gmail.com/Other computers/My MacBook Pro (1)/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_7.stan")
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOU_v1_7.stan")
  stan_model <- stan_model("blouchOU_v1_7.stan")
  
  #######################
  stan_data<-stan.OU.data.setup(OUBMparameters,tree,n)
  #################
  saved.parameters<-data.frame(matrix(NA,n,28))
  parameter.names<-c("optimal.reg.median","optimal.reg.mean","optimal.reg.025","optimal.reg.975",
                     "evol.reg.median","evol.reg.mean","evol.reg.025","evol.reg.975",
                     "t12.median","t12.mean","t12.025","t12.975","vy.median","vy.mean","vy.025","vy.975","sigma2y.median","sigma2y.mean","sigma2y.025","sigma2y.975","Ya.median","Ya.mean","Ya.025","Ya.975","b0.median","b0.mean","b0.025","b0.975")
  names(saved.parameters)<-parameter.names
  #######################
  saved.parameters.mode<-data.frame(matrix(NA,n,7))
  parameter.names<-c("optimal.reg.mode","evol.reg.mode","t12.mode","vy.mode","sigma2y.mode","Ya.mode","b0.mode")
  names(saved.parameters.mode)<-parameter.names
  ###############
  
  for (i in 1:n){
    print(paste("Iteration ",i))
    stan_data_sub<-stan_data[[i]]
    stan_data_sub<-stan_data_sub[[1]]
    #fit.random<- stan(file = stan.model,data = stan_data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    fit<- rstan::sampling(object = stan_model,data = stan_data_sub,chains = 1,iter = 2000,verbose=F,control=list(adapt_delta=0.975))
    par.summary <- summary(fit, pars = c("beta", "beta_e","hl","vy","sigma2_y","Ya","b0"), probs = c(0.025,0.5, 0.975))$summary
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1:4]<-par.summary[2,]
    saved.parameters[i,5:8]<-par.summary[3,]
    saved.parameters[i,9:12]<-par.summary[4,]
    saved.parameters[i,13:16]<-par.summary[5,]
    saved.parameters[i,17:20]<-par.summary[6,]
    saved.parameters[i,21:24]<-par.summary[7,]
    saved.parameters[i,25:28]<-par.summary[8,]
    
        
    fit.ext<-extract(fit)
    saved.parameters.mode[i,1]<-Mode(fit.ext$beta[2])
    saved.parameters.mode[i,2]<-Mode(fit.ext$beta_e)
    saved.parameters.mode[i,3]<-Mode(fit.ext$hl)
    saved.parameters.mode[i,4]<-Mode(fit.ext$vy)
    saved.parameters.mode[i,5]<-Mode(fit.ext$sigma2_y)
    saved.parameters.mode[i,6]<-Mode(fit.ext$Ya)
    saved.parameters.mode[i,7]<-Mode(fit.ext$b0)
    
        }
  return(list(saved.parameters,saved.parameters.mode))
  }

stan.BM.simulation.reps<-function(BMparameters,tree,n){
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/blouchBM_v1.stan")
  stan.model <- "blouchBM_v1.stan"
  #######################
  stan.data<-stan.BM.data.setup(BMparameters,tree,n)
  #######################
  saved.parameters<-data.frame(matrix(NA,n,8))
  parameter.names<-c("evol.reg.median","evol.reg.mean","evol.reg.025","evol.reg.975","sigma2y.median","sigma2y.mean","sigma2y.025","sigma2y.975")
  names(saved.parameters)<-parameter.names
  ###############
  saved.parameters.mode<-data.frame(matrix(NA,n,2))
  parameter.names<-c("evol.reg.mode","sigma2y.mode")
  names(saved.parameters)<-parameter.names
  ###############
  
  for (i in 1:n){
    print(paste("Iteration ",i))
    stan.adaptive.data<-stan.data[[i]]
    stan.adaptive.data<-
    fit.random<- stan(file = stan.model,data = stan.adaptive.data,chains = 1,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    
    par.summary <- summary(fit.random, pars = c("beta","sigma2_y"), probs = c(0.025,0.5, 0.975))$summary
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1:4]<-par.summary[1,]
    saved.parameters[i,5:8]<-par.summary[2,]

    fit.random.ext<-extract(fit.random)
    saved.parameters.mode[i,1]<-Mode(fit.random.ext$beta)
    saved.parameters.mode[i,2]<-Mode(fit.random.ext$sigma2_y)

      }
  return(list(saved.parameters,saved.parameters.mode))
}
