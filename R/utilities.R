## 
##
getMolfun <- function(pfilename){
  pgl <- read.csv(pfilename , stringsAsFactors=FALSE)
  gofunction <- pgl$GOMolecularFunction
  dd <- lapply(gofunction , function(x){d <- strsplit(x,";"); return(d[[1]][1])})
  dd <- lapply(dd , function(x){d <- strsplit(x,","); return(d[[1]][1])})
  tmp <- (unlist(dd))
  return(list(ids=pgl$Mapped.IDs, molfun=tmp))
}

##
## Removes rows with nas from matrix
##
rmNarows <- function(plev, thresh=0 )
{
  x <- apply(plev,1,function(x){sum(is.na(x))})
  plev <- plev[-which(x>thresh),]
}



createGraph <- function(adjMatrix,protnam){
  library(igraph)
  graph.speed=graph.adjacency(adjmatrix=adjMatrix, mode="undirected", weighted=TRUE)
  V(graph.speed)$name = protnam
  temp.degree=degree(graph.speed)
  pname <- protnam
  pname[temp.degree==0]=NA
  V(graph.speed)$label=pname
  return(graph.speed)
}



giveCols2names <- function(protnam,tt,ml){
  mcol <- rep("",length(protnam))
  for(i in 1:length(protnam))
  {
    x <- (tt[which(tt$ids==protnam[i]),2])
    if(length(x)==0)
    {
      mcol[i] <- "white" 
    }
    else if(is.na(x)){
      mcol[i] <- "white"
    }
    else{
      #cat(protnam[i]," " ,ml[[x]],"\n")
      if(length(ml[[x]])==0)
      {
        mcol[i] <- "lightgray"
      }
      else
      {
        mcol[i] <- ml[[x]]
      }
    }
  }
  return(mcol)
}

fixSampleNames <- function(cnames)
{
  cnames<-sub("\\.","\\_",cnames)
  return(cnames)
}

splitcolnames <- function(cnames)
{
  cnamessplit <- strsplit(cnames,"\\_")
  colnamestable <- NULL
  for(i in cnamessplit)
  {
    colnamestable <- rbind(colnamestable,i)
  }
  #colnamestable <- as.data.frame(colnamestable)
  #ftable(colnamestable$Time,colnamestable$strain,colnamestable$replicate,colnamestable$condition)
  return(colnamestable)
}


split2table <- function(protnames,split="\\||\\_")
{
  protnames[1:10]
  cnamessplit <- strsplit(protnames,split)
  protnam <- NULL
  for(i in cnamessplit)
  {
    protnam <- rbind(protnam,i)
  }
  return(protnam)
}
