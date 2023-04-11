set.regimes<-function(trdata,regimes){
  #returns a new column of the dataset - regimes
  #also returns internal node assignment
  #Make sure you send a merged tree and data file using bayou setup
  
  library(phytools)
  
  getDescendants<-function(tree,node,curr=NULL){
    #plot(tree,font=0.25); nodelabels(bg="white")
    if(is.null(curr)) curr<-vector()
    daughters<-tree$edge[which(tree$edge[,1]==node),2]
    curr<-c(curr,daughters)
    w<-which(daughters>=length(tree$tip))
    if(length(w)>0) for(i in 1:length(w)) 
      curr<-getDescendants(tree,daughters[w[i]],curr)
    return(curr)
  }
  
  
  
  
  #Change color for Plots
  library(ggplot2)
  #source("/Users/markgrabowski/Documents/Academic/Research/Current Projects/Brain Size Diet + Social System/Programs/getDescendants.R")
  
  #Get ggplot colors used for plot to make on tree
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
  }
  
  n.tips<-length(trdata$phy$tip.label)
  num.internal.nodes<-n.tips+(1:trdata$phy$Nnode)
  n.internal.nodes<-length(num.internal.nodes)
  
  #First set all exant tips to OU1
  trdata$dat$regimes<-"OU1"
  #Assign all internal nodes to OU1
  internal.nodes.regimes<-rep("OU1",n.internal.nodes)
  
  #Shifts from Bayou Species + BoM AAPA 03202019 Code - saved on Google Drive
  #OU2
  for(i in 1:(length(regimes))){
    print(i)
    print(regimes[i])
    if(regimes[i]<=n.tips)#shift on single branch
    {
      #Assign external nodes/tips to OUi
      trdata$dat$regimes[regimes[i]]<-paste("OU",sep="",i+1)
      #trdata$phy$tiplabel
      
    }
    else{
      saved.decendants<-getDescendants(trdata$phy,node=regimes[i])
      external.nodes<-saved.decendants[saved.decendants<=n.tips]
      internal.nodes<-saved.decendants[saved.decendants>n.tips]
      
      #Assign external nodes/tips to OUi
      trdata$dat$regimes[external.nodes]<-paste("OU",sep="",i+1)
      
      #Assign internal nodes to OUi
      internal.nodes.regimes[internal.nodes-n.tips]<-paste("OU",sep="",i+1)
      
      #Assign node where shift occurs to OUi
      internal.nodes.regimes[regimes[i]-n.tips]<-paste("OU",sep="",i+1)
    }
  }  
  trdata$phy$node.label<-internal.nodes.regimes
  #Make colors for regimes
  reg.colors<-gg_color_hue(length(regimes)+1)
  
  #Combine external coding and internal coding to plot tree with colored shifts
  regimes.total<-c(trdata$dat$regimes,internal.nodes.regimes)
  edge.regimes <- factor(regimes.total[trdata$phy$edge[,2]])
  print(edge.regimes)
  print(reg.colors)
  plot(trdata$phy,edge.color = reg.colors[edge.regimes], edge.width = 1, cex = 0.2)
  #return(list(trdata,internal.nodes.regimes))
  return(trdata)
}
