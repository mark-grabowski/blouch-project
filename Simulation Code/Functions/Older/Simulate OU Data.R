simulate.OU.direct<-function(tree,hl,vy,b1,vX0,vY0,Sxx,theta,num.X){
  #A=(a), Syy=(s), Sxx=(b2), Syx=(b1), mPsi=(r)
  
  A<-log(2)/hl
  Syy<-sqrt(vy*(2*A))
  #Sxx<-10
  Syx<-b1
  mPsi<-theta
  B<- -A*b1
  
  tip.num<-length(tree$tip.label)
  times <- ape::node.depth.edgelength(tree)
  T.term <- times[1:tip.num]              
  expected.slope<-b1*(Syy^2+Sxx)*(1-exp(-A*T.term[1]))/(2*Sxx*T.term[1])

  #OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),B=matrix(0,ncol=1,nrow=1),mPsi=matrix(0,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),Syx=matrix(Syx,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1))
  OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                       B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                       Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                       Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
  
  
  OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
  OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
  OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
  #return(OUBM.trdata)
  stan.data<-stan.data.setup(OUBM.trdata,1,0,0,0)
  return(stan.data)
  }

simulate.OU.adaptive<-function(tree,hl,vy,b1,vX0,vY0,Sxx,theta,num.X){
  #A=(a), Syy=(s), Sxx=(b2), Syx=(b1), mPsi=(r)
  
  A<-log(2)/hl
  Syy<-sqrt(vy*(2*A))
  #Sxx<-10
  Syx<-b1
  mPsi<-theta
  B<- -A*b1
  
  tip.num<-length(tree$tip.label)
  times <- ape::node.depth.edgelength(tree)
  T.term <- times[1:tip.num]              

  OUBMparameters<-list(vY0=matrix(vY0,ncol=1,nrow=1),A=matrix(A,ncol=1,nrow=1),
                       B=matrix(B,ncol=1,nrow=1),mPsi=matrix(theta,ncol=1,nrow=1),
                       Syy=matrix(Syy,ncol=1,nrow=1),vX0=matrix(vX0,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1),
                       Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1),starting_point_for_optim=list(A=matrix(A,ncol=1,nrow=1),Syy=matrix(Syy,ncol=1,nrow=1),B=matrix(B,ncol=1,nrow=1),Sxx=matrix(Sxx,ncol=1,nrow=1)))
  
  
  OUBMdata<-simulMVSLOUCHProcPhylTree(tree,OUBMparameters,regimes=NULL) #Simulate data
  OUBMdata<-data.frame(OUBMdata,"tip.names"=row.names(OUBMdata))
  OUBM.trdata <- make.treedata(tree, OUBMdata,name_column="tip.names") #Combine data with 
  #return(OUBM.trdata)
  stan.data<-stan.data.setup(OUBM.trdata,0,1,0,0)
  return(stan.data)
}


  
stan.data.setup<-function(OUBM.trdata,num.direct,num.random,num.ME.direct,nun.ME.random){
  n<-length(OUBM.trdata$phy$tip.label)
  Z<-num.direct+num.random
  Y<-data.frame(OUBM.trdata$dat[,1])
  #return(Y)
  if(num.random>0){
    random<-data.frame(OUBM.trdata$dat[,-1])
    direct<-matrix(0,nrow=n,ncol=0)
    }
  else{
    direct<-data.frame(OUBM.trdata$dat[,-1])
    random<-matrix(0,nrow=n,ncol=0)
    }
  
  mrca1 <- ape::mrca(OUBM.trdata$phy)
  times <- ape::node.depth.edgelength(OUBM.trdata$phy)
  ta <- matrix(times[mrca1], nrow=n, dimnames = list(OUBM.trdata$phy$tip.label, OUBM.trdata$phy$tip.label))
  
  T.term <- times[1:n]
  tia <- times[1:n] - ta
  tja <- t(tia)
  tij <- tja + tia
  
  mv.response<-rep(0,n)
  mv.random<-matrix(0,nrow=n,ncol=0)
  mv.direct<-matrix(0,nrow=n,ncol=0)
  
  if(num.random!=0){
      rownames(random)<-OUBM.trdata$phy$tip.label
      brownian<-sigma.X.estimate.mv(OUBM.trdata$phy,random)
      sigma_squared<-matrix(brownian[[2]],dim(random)[2],dim(random)[2])
      brownian_mean<-brownian[[1]]
  }else{
    sigma_squared <- matrix(0,0,0)
    brownian_mean <- matrix(0,0,1)
  }
  
  combo.data<-cbind(Y,direct,random)

  names(combo.data)[1]<-"y"
  
  coef<-lm(y ~ ., data=combo.data)$coef
  ols.intercept<-coef[1]
  ols.slope<-as.array(coef[-1])
  #Setting up data
  stan_data<-list(N = n,
                  Z = num.direct+num.random, #Total predictor number
                  Z_direct = num.direct, #Total direct trait number
                  Z_random = num.random, #Total random trait number
                  Z_me_direct = 0,
                  Z_me_random = 0,
                  Y = Y[,1],
                  mv_response = mv.response,
                  random_cov = random,
                  mv_random_cov = mv.random,
                  direct_cov = direct,
                  mv_direct_cov = mv.direct,
                  ta = as.matrix(ta),
                  T_term = as.vector(T.term),
                  tia = as.matrix(tia),
                  tja = as.matrix(tja),
                  tij = as.matrix(tij),
                  brownian_mean = t(as.matrix(brownian_mean)),
                  sigma_squared_x = as.matrix(sigma_squared),
                  ols_intercept = ols.intercept,
                  ols_slope = ols.slope)
  }  
  

sigma.X.estimate.mv<-function(phy,predictors){
  return(list(mean = rep(0,length(predictors)),sigma_squared = geiger::ratematrix(phy,predictors)))
}