stan.OU.tree.data.setup<-function(OUBMparameters,tips,n,A,opt.slope,sym){
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/Blouch Setup Files/blouchOU.setup.v1.R")
  #n<-50
  names.traits<-c("trait_1","trait_2",NA,NA)
  stan.data<-list()
  for(i in 1:n){
    tips.evol.slope<-NULL
    if(sym != 1){
      #Following Uyeda and Harmon, 2014
      #tree<-sim.bd.taxa(tips, numbsim=1, lambda=0.1, mu=0, frac = 1, complete = TRUE,
       #               stochsampling = FALSE)[[1]]
      
      #Following bayou tutorial: https://github.com/uyedaj/bayou/blob/master/tutorial.md
      tree <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = tips)
      #tree <- sim.bdtree(b = 1, d = 0.1, stop = "taxa", n = tips)
      
    }
    else{
      #Ultrametric tree
      tree <- sim.taxa(1, n=tips, waitsp=function()2,symmetric=TRUE, complete=TRUE)
      tree<-tree[[1]]
      }
    l.tree<-max(branching.times(tree))
    tree$edge.length<-tree$edge.length/l.tree ## rescale tree to height 1
    
    n<-length(tree$tip.label)
    mrca1 <- ape::mrca(tree)
    times <- ape::node.depth.edgelength(tree)
    T.term <- times[1:n]              
    
    rho<- (1-(1 - exp(-A * T.term))/(A * T.term))
    tips.evol.slope<-c(tips.evol.slope,(rho*opt.slope)[1])
    
    ###
    OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
    OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
    OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
    #OUBM.trdata$dat$trait_2<-OUBM.trdata$dat$trait_2-mean(OUBM.trdata$dat$trait_2) #Mean standardize X
    stan.data[[length(stan.data)+1]]<-c(blouchOU.setup.v1(OUBM.trdata,names.traits),tips.evol.slope)
  }
  return(stan.data)
}

stan.BM.tree.data.setup<-function(BMparameters,tips,n,sym){
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/Blouch Setup Files/blouchOU.setup.v1.R")
  #n<-50
  names.traits<-c("trait_1","trait_2",NA,NA)
  stan.data<-list()
  for(i in 1:n){
    if(sym != 1){
      #Following Uyeda and Harmon, 2014
      #tree<-sim.bd.taxa(tips, numbsim=1, lambda=0.1, mu=0, frac = 1, complete = TRUE,
      #               stochsampling = FALSE)[[1]]
      
      #Following bayou tutorial: https://github.com/uyedaj/bayou/blob/master/tutorial.md
      tree <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = tips)
    }
    else{
      #Ultrametric tree
      tree <- sim.taxa(1, n=tips, waitsp=function()2,symmetric=TRUE, complete=TRUE)
      tree<-tree[[1]]
    }

    l.tree<-max(branching.times(tree))
    tree$edge.length<-tree$edge.length/l.tree ## rescale tree to height 1
    
    #########################
    BMdata<-simulBMProcPhylTree(tree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
    BMdata<-data.frame(BMdata,"tip.names"=row.names(BMdata))
    BM.trdata <- make.treedata(tree, BMdata,name_column="tip.names") #Combine data with 
    #BM.trdata$dat$trait_2<-BM.trdata$dat$trait_2-mean(BM.trdata$dat$trait_2) #Mean standardize X
    stan.data[[length(stan.data)+1]]<-blouchOU.setup.v1(BM.trdata,names.traits)
  }
  return(stan.data)
}

stan.OU.simulation.tips.reps<-function(OUBMparameters,tips,n,opt.slope,sym){
  library(TreeSimGM)
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/blouchOU_test.stan")
  stan.model <- "blouchOU_test.stan"

  #stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/blouchOU_test_v2.stan")
  #stan.model <- "blouchOU_test_v2.stan"
  
  #######################
  stan.data<-stan.OU.tree.data.setup(OUBMparameters,tips,n,as.numeric(OUBMparameters$A),opt.slope,sym)
  #######################
  saved.parameters<-data.frame(matrix(NA,n,20))
  parameter.names<-c("optimal.reg.median","optimal.reg.mean","optimal.reg.025","optimal.reg.975",
                     "evol.reg.median","evol.reg.mean","evol.reg.025","evol.reg.975",
                     "t12.median","t12.mean","t12.025","t12.975","sigma2y.median","sigma2y.mean","sigma2y.025","sigma2y.975","vy.median","vy.mean","vy.025","vy.975")
  names(saved.parameters)<-parameter.names
  evol.slope<-NULL
  ###############
  
  for (i in 1:n){
    print(paste("Iteration ",i))
    #print(stan.data[i])
    stan.ind.data<-stan.data[[i]]
    #print(stan.adaptive.data[2])
    stan.adaptive.data<-stan.ind.data[[2]]
    stan.constraint.data<-stan.ind.data[[1]]
    evol.slope<-c(evol.slope,stan.ind.data[[3]])
    #print(evol.slope)
    
    fit.random<- stan(file = stan.model,data = stan.adaptive.data,chains = 2,iter = 2000,verbose = F, control=list(adapt_delta=0.95))
    par.summary <- summary(fit.random, pars = c("beta", "beta_evol","hl","sigma2_y","vy"), probs = c(0.025,0.5, 0.975))$summary
    par.summary<-par.summary[, c("50%","mean","2.5%","97.5%")]
    saved.parameters[i,1:4]<-par.summary[1,]
    saved.parameters[i,5:8]<-par.summary[2,]
    saved.parameters[i,9:12]<-par.summary[3,]
    saved.parameters[i,13:16]<-par.summary[4,]
    saved.parameters[i,17:20]<-par.summary[5,]
    
      }
  return(list(saved.parameters,evol.slope=evol.slope))
}

stan.BM.simulation.tips.reps<-function(BMparameters,tips,n,sym){
  library(TreeSimGM)
  setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1")
  stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/Stan Functions/blouch v1/blouchBM_v1.stan")
  stan.model <- "blouchBM_v1.stan"
  #######################
  stan.data<-stan.BM.tree.data.setup(BMparameters,tips,n,sym)
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
