repeated_optima_sim<-function(trdata,n,names.traits,optima){
  #names.traits in order (names.direct.traits,names.direct.me.traits,names.fixed.fact)
  source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/R Setup Code/R1/blouchOUReg.setup.mv_v1_3.R")
      setwd("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1")
      names.direct.traits<-names.traits[1]
      names.direct.me.traits<-names.traits[2]
      names.fixed.fact<-names.traits[3]
      mlm.error<-NULL
      npm.error<-NULL
      saved.mlm.optima<-NULL
      saved.npm.optima<-NULL
      mlm.inaccuracy<-NULL
      mlm.optima.inaccuracy<-NULL
      npm.inaccuracy<-NULL
      npm.optima.inaccuracy<-NULL
      
      #stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOUReg_test.stan")
      #stan_model <- stan_model("blouchOUReg_test.stan")

      stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOUReg_opt_test.stan")
      stan_model <- stan_model("blouchOUReg_opt_test.stan")
      
            
      for(i in 1:n){
        #Make n separate stan data files
        name.response.trait<-paste("Sim",i,sep="")
        
        #Sim.data<-trdata$dat[name.response.trait][[1]]#For Mean Centered Y
        #Sim.mean<-mean(Sim.data)#For Mean Centered Y
        #Sim.data.mc<-data.frame(Sim.data-mean(Sim.data))#For Mean Centered Y

        #name.response.trait<-paste("Sim.MC",i,sep="")#For Mean Centered Y
        #names(Sim.data.mc)<-name.response.trait #For Mean Centered Y
        #trdata$dat<-data.frame(trdata$dat,Sim.data.mc)#For Mean Centered Y
        print(name.response.trait)
        
        #One random covariate w/o ME
        #stan_data<-blouchOUReg.setup.mv(ruminant.trdata,names.fixed.fact,name.response.trait,names.direct.traits,names.random.traits=NULL)
      
        #One direct covariate w/o ME
        stan_test_data<-blouchOUReg.setup.mv(trdata,names.fixed.fact,name.response.trait,names.direct.traits,names.random.traits=NULL)
        fit.mlm<- rstan::sampling(object = stan_model,data = stan_test_data,chains = 1,iter = 2000,control=list(adapt_delta=0.9))

        #For downstream analysis and plots
        extract.fit.mlm <- rstan::extract(fit.mlm)
        #return(extract.fit.mlm)
        
        mu.post.beta<-apply(extract.fit.mlm$beta[,1:length(optima)],2,median)
        saved.mlm.optima<-c(saved.mlm.optima,mu.post.beta)
        mlm.error<-c(mlm.error,abs(mu.post.beta-optima))
        mlm.inaccuracy<-c(mlm.inaccuracy,mean((mu.post.beta-optima)^2))
        
        
        #Mean centered Y
        #mu.post.beta<-apply(extract.fit.mlm$beta[,1:length(optima)],2,median)+Sim.mean
        #mlm.error<-c(mlm.error,abs(mu.post.beta-optima))
        #mlm.inaccuracy<-c(mlm.inaccuracy,mean((mu.post.beta-optima)^2))
        
        
        
        }
      stanc("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Blouch project/R1/Code/Stan Code/R1/blouchOUReg_v1_5.stan")
      stan_model <- stan_model("blouchOUReg_v1_5.stan")
      
      for(i in 1:n){
        name.response.trait<-paste("Sim",i,sep="")
        #Sim.data<-trdata$dat[name.response.trait][[1]]#For Mean Centered Y
        #Sim.mean<-mean(Sim.data)#For Mean Centered Y
        
        #name.response.trait<-paste("Sim.MC",i,sep="") #For Mean Centered Y

        print(name.response.trait)
        
        stan_test_data<-blouchOUReg.setup.mv(trdata,names.fixed.fact,name.response.trait,names.direct.traits,names.random.traits=NULL)
        fit.npm<- rstan::sampling(object = stan_model,data = stan_test_data,chains = 1,iter = 2000,control=list(adapt_delta=0.9))
        
        #For downstream analysis and plots
        extract.fit.npm <- rstan::extract(fit.npm)
        
        mu.post.beta<-apply(extract.fit.npm$beta[,1:length(optima)],2,mean)
        saved.npm.optima<-c(saved.npm.optima,mu.post.beta)
        
        npm.error<-c(npm.error,abs(mu.post.beta-optima))
        npm.inaccuracy<-c(npm.inaccuracy,mean((mu.post.beta-optima)^2))
        
        #Mean centered Y
        #mu.post.beta<-apply(extract.fit.npm$beta[,1:length(optima)],2,median)+Sim.mean
        #npm.error<-c(npm.error,abs(mu.post.beta-optima))
        #npm.inaccuracy<-c(npm.inaccuracy,mean((mu.post.beta-optima)^2))
        
        
        
      }
      return(list(rbind(saved.mlm.optima,saved.npm.optima),
                  rbind(mlm.error,npm.error),rbind(mlm.inaccuracy,npm.inaccuracy))) 
  }
  
      

