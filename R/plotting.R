#' image plot with labels
#'
#' @param x - matrix
#' @param labels - colnames(x)
#' @param cex.axis - size of axis lables
#' @param cex - size of labels
#' @export
imagelabels = function(x, labels=colnames(x),cex=1,cex.axis=0.5)
{
  image(x, axes = F )
  axis( 1, at=seq(0,1,length=length((labels))) , labels=labels,cex.axis=cex.axis, las=2, cex=cex )
  axis( 2, at=seq(0,1,length=length((labels))) , labels=labels,cex.axis=cex.axis, las=1, cex=cex )
}
#' altman bland
#' @export
altmanbland = function(x,y,main=""){
  idx<-apply(cbind(x,y),1,function(x){sum(x==0)==0})
  mean  = (x+y)/2
  absdiff = abs( x-y )
  #scatter.smooth(mean[idx],absdiff[idx],span=1/3,col=2,pch=".",cex=2,xlab="(y+x)/2",ylab="abs(x-y)",main=main,log="xy")  
  plot(mean,absdiff,log="xy",xlab="(y+x)/2",ylab="abs(x-y)",pch="x",cex=0.5,main=main)
  lines(lowess(mean,absdiff),col=2,lwd=2)
}
#' plot QQ plot
#' @export
pairsQQ = function(dataframesel,main=""){
  pairs(dataframesel, panel = function(x,y){
    r <- qqplot(x , y , plot.it = FALSE)
    points(r,pch=".",cex=2)
    abline(0,1,col=2)
  }
        , lower.panel=NULL
        ,main = main
  )
}
#' plot ratios against retention time
#' @export
pairsRatio = function(dataframesel,RT,main=""){
  pairs(dataframesel, panel = function(x,y,...){
    r <- points(RT , log2(x / y),pch=".",... )
    abline(h = 0,col=2)
  }
        ,xlim=range(RT)
        ,ylim=range(log2(x/y),na.rm=TRUE)
        , lower.panel=NULL
        ,main = main
  )
}
#' normal pairs plot with different pch and plus abline
#' @param pln data matrix or data.frame as normally passed to pairs
#' @param ... params usually passed to pairs
#' @export
mypairs = function(dataframesel,...){
  pairs(dataframesel, panel = function(x,y){
    points(x, y, pch=".")
    abline(0,1,col=2)
  }
        , lower.panel=NULL,...
  )
}
#' heatmap2 to wrapper
#' @export
myheat = function(pln,main="",distf=dist,hclustf=hclust,mypalette,mycol){
  tmp <- heatmap.2( as.matrix(pln) , trace="none" , scale="none" , col=mypalette ,
                    ColSideColors=mycol ,
                    #RowSideColors=protcol,
                    labRow="",
                    cexRow=0.1 + 1/log10(dim(pln)[1]),
                    cexCol=0.1 + 1/log10(dim(pln)[2]),
                    distfun=distf,hclustfun=hclustf,
                    margins=c(5,5),main=main
  )
}
#' volcano plot
#' @param foldchange - fold change values
#' @param pvals - pvalues
myvolcanoplot = function(foldchange, pvals , pthresh = 0.05, ratiothresh = 2, xlab ="log2(T/N)" ,ylab = "-log10(P)"){
  d <- data.frame(ratio = foldchange, pvals = pvals )
  plot(d$ratio,-log10(d$pvals),col="#00000033",pch=19,xlab=xlab, ylab=ylab)
  upsubset<-subset(d,pvals < pthresh & ratio > ratiothresh)
  points(upsubset$ratio,-log10(upsubset$pvals),col=2,pch=19)
  points(upsubset$ratio,-log10(upsubset$pvals),col=1,pch=1)
  abline(h=-log10(pthresh),col="gray")
  downsubset<-subset(d,pvals<pthresh & ratio < -ratiothresh)
  points(downsubset$ratio,-log10(downsubset$pvals),col=3,pch=19)
  points(downsubset$ratio,-log10(downsubset$pvals),col=1,pch=1)
  abline(v=c(-ratiothresh,ratiothresh))
  return(list(upsubset=upsubset,downsubset=downsubset))
}


