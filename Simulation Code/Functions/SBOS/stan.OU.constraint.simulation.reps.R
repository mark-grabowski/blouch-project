stan.OU.data.setup<-function(OUBMparameters,tree,n){
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/Blouch Setup Files/blouchOU.setup.v1.R")
  #n<-50
  names.traits<-c("trait_1","trait_2",NA,NA)
  stan.data<-list()
  for(i in 1:n){
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.v1(OUBM.trdata,names.traits)
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

stan.OU.constraint.simulation.reps<-function(OUBMparameters,tree,n){
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/blouchOU_test.stan")
  stan.model <- "blouchOU_test.stan"
  #######################
  stan.data<-stan.OU.data.setup(OUBMparameters,tree,n)
  #######################
  saved.parameters<-data.frame(matrix(NA,n,20))
  parameter.names<-c("optimal.reg.median","optimal.reg.mean","optimal.reg.025","optimal.reg.975",
                     "evol.reg.median","evol.reg.mean","evol.reg.025","evol.reg.975",
                     "t12.median","t12.mean","t12.025","t12.975","sigma2y.median","sigma2y.mean","sigma2y.025","sigma2y.975","vy.median","vy.mean","vy.025","vy.975")
  names(saved.parameters)<-parameter.names
  ###############
  
  for (i in 1:n){
    print(paste("Iteration ",i))
    stan.adaptive.data<-stan.data[[i]]
    stan.constraint.data<-stan.adaptive.data[[1]]
    stan.adaptive.data<-stan.adaptive.data[[2]]
    fit.constraint<- stan(file = stan.model,data = stan.constraint.data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    par.summary <- summary(fit.random, pars = c("beta", "beta_evol","hl","sigma2_y","vy"), probs = c(0.025,0.5, 0.975))$summary
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1:4]<-par.summary[1,]
    saved.parameters[i,5:8]<-par.summary[2,]
    saved.parameters[i,9:12]<-par.summary[3,]
    saved.parameters[i,13:16]<-par.summary[4,]
    saved.parameters[i,17:20]<-par.summary[5,]
    
      }
  return(saved.parameters)
  }

stan.BM.constraint.simulation.reps<-function(BMparameters,tree,n){
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
  
  for (i in 1:n){
    print(paste("Iteration ",i))
    stan.adaptive.data<-stan.data[[i]]
    stan.adaptive.data<-stan.adaptive.data[[2]]
    fit.random<- stan(file = stan.model,data = stan.adaptive.data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    par.summary <- summary(fit.random, pars = c("beta","sigma2_y"), probs = c(0.025,0.5, 0.975))$summary
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1:4]<-par.summary[1,]
    saved.parameters[i,5:8]<-par.summary[2,]

      }
  return(saved.parameters)
}
