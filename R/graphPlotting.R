

plotGraph <- function(graph,myl=NULL)
{
  vdeg <- degree(graph)
  mypalette<-brewer.pal(9,"BuGn")
  vcols <- mypalette[round(rank(vdeg-max(vdeg))/length(vdeg) * 7)+1]
  
  mypalette <- brewer.pal(9,"YlGnBu")
  weights <- E(graph)$weight+1
  ecols <- round(((weights/max(weights)*8)+1))
  range(ecols)
  ecols <- mypalette[ecols]
  
  x <- plot(graph, layout=myl,  vertex.size=3, vertex.frame.color="white",  edge.color= ecols,vertex.color=vcols,vertex.label.font=1,vertex.label.cex=0.7,vertex.label.color=1)
  return(x)
}

plotGraph2 <- function(graph,myl){
  vdeg <- degree(graph)
  mypalette<-brewer.pal(9,"BuGn")
  #vcols <- mypalette[round(rank(vdeg-max(vdeg))/length(vdeg) * 7)+1]
  
  mypalette <- brewer.pal(9,"RdBu")
  weights <- E(graph)$color
  weights[weights>0] <- sqrt(weights[weights>0])
  weights[weights<0] <- -1* sqrt(-1*weights[weights<0])
  weights <- weights+1
  ecols <- round(((weights/2 * 8)+1))
  ecols <- mypalette[ecols]
  
  x <- plot(graph, layout=myl,  vertex.size=3, vertex.frame.color="white",  edge.color= ecols,vertex.label.font=1,vertex.label.cex=0.7,vertex.label.color=1)
  return(x)
}

plotGraphAnn <- function(graph,myl,vfc="white",family=1){
  library(RColorBrewer)
  
  mypalette <- brewer.pal(9,"RdBu")
  weights <- E(graph)$color
  weights[weights>0] <- (weights[weights>0])^0.2
  weights[weights<0] <- -1* (-1*weights[weights<0])^0.2
  weights <- weights+1
  ecols <- round(((weights/2 * 8)+1))
  ecols <- mypalette[ecols]
  
  x <- plot(graph, layout=myl,  vertex.size=6, vertex.frame.color=vfc, edge.width=4, edge.color= ecols,vertex.label=V(graph)$protname ,vertex.label.family=family,vertex.label.font=1 ,vertex.label.dist=0.4, vertex.label.cex=1 , vertex.label.color=1)
  return(x)
}
