#' image plot with labels
#'
#' @param x - matrix
#' @param labels - colnames(x)
#' @param cex.axis - size of axis lables
#' @param cex - size of labels
#' @export
#' @examples
#' x = matrix(rnorm(20*20),ncol=20)
#' imagelabels(x)
imagelabels = function(x, labels=colnames(x),cex=1,cex.axis=0.5,main=NULL,col = heat.colors(12))
{
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(3,1), heights=c(1,1))
  
  image(x, axes = F, main =main, col=col)
  axis( 1, at=seq(0,1,length=length((labels))) , labels=labels,cex.axis=cex.axis, las=2, cex=cex )
  axis( 2, at=seq(0,1,length=length((labels))) , labels=labels,cex.axis=cex.axis, las=1, cex=cex )
  
  colorlevels = seq(min(x,na.rm = TRUE),max(x,na.rm = TRUE),length=length(col))
  image(1, seq(0,1,length=length(colorlevels)),
        matrix(data=colorlevels, nrow=1),
        col=col,xlab="",ylab="",
        axes=FALSE)
  axis( 2, at=seq(0,1,length=length((colorlevels))) , labels=round(colorlevels,digits=2),cex.axis=cex.axis, las=1, cex=cex )
  layout(1)
}
#' if you need an colorscale to you imagelables use this
#' @param data the data matrix
#' @param colors used
#' @export
colorscale = function(data,colors=heat.colors(12)){
  nrc = length(colors)
  z  = seq( min(data) , max(data) , length=nrc)
  image(1, seq(0,1,length=nrc), matrix(z,1,nrc) ,axes=F,ylab="",xlab="")
  axis( 2, at=seq(0,1,length=nrc) , labels=round(z,digits=2), las=2 )
}
#' altman-bland plot 
#' @description
#' plots abs(x-y) against (y+x)/2
#' @param x - input intensities
#' @param y -
#' @param main plotting parameters
#' @param pch -
#' @param log -
#' @export
#' @examples
#' data(SDat)
#' altmanbland(SDat$Intensity[,1],SDat$Intensity[,2],log="xy")
#' 
altmanbland = function(x,y,main="",pch=".",log=""){
  idx<-apply(cbind(x,y),1,function(x){sum(x==0)==0})
  mean  = (x+y)/2
  absdiff = abs( x-y )
  #scatter.smooth(mean[idx],absdiff[idx],span=1/3,col=2,pch=".",cex=2,xlab="(y+x)/2",ylab="abs(x-y)",main=main,log="xy")  
  plot(mean,absdiff,log=log,xlab="(y+x)/2",ylab="abs(x-y)",pch="x",cex=0.5,main=main)
  lines(lowess(mean,absdiff),col=2,lwd=2)
}
#' pairsplot of QQ plots
#' @export
#' @examples
#' data(SDat)
#' pairsQQ( SDat$Intensity )
#' @seealso \code{\link{qqplot}} and  \code{\link{pairs}}
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
#' plot log2(x/y) against retention time
#' @export
#' @examples
#' data(SDat)
#' rtord = orderByRT(SDat)
#' pairsRatio( rtord$Intensity, rtord$RT )
#' @seealso \code{\link{pairs}} and \code{\link{pairsDifference}}
pairsRatio = function(dataframesel,RT,main="",ylim=c(-2,2)){
  pairs(dataframesel, panel = function(x,y,...){
    y = log2(x / y)
    r <- points(RT[!is.na(y)] , y[!is.na(y)],pch=".",... )
    abline(h = 0,col=2)
  }
  ,xlim=range(RT)
  ,ylim=ylim
  , lower.panel=NULL
  ,main = main
  )
}
#' plot x-y against retention time
#' @export
#' @param dataf dataframe
#' @param RT retention time vector
#' @param main
#' @param ylim
#' @param maxPanel maximum number of datasets to display. dim(dataf)[2] > maxPanel, maxPanel columns are sampled
#' @examples
#' data(SDat)
#' rtord = orderByRT(SDat)
#' pairsDifference( asinh(rtord$Intensity), rtord$RT , ylim=NULL)
#' @seealso \code{\link{pairs}} and \code{\link{pairsRatio}}
pairsDifference = function(dataf,RT,main="",ylim=NULL, maxPanel=8 ){
  stopifnot(dim(dataf)[2] > 1)
  if(dim(dataf)[2] > maxPanel){
    idx = sample(1:dim(dataf)[2],maxPanel)
    dataf = dataf[,idx]
  }
  if(is.null(ylim)){
    nr = dim(dataf)[2]
    tmp = NULL
    for(i in 1:nr){
      for(j in 1:nr){
        if(i != j){
          tmp = range(c(tmp,(dataf[,i] - dataf[,j])), na.rm=TRUE)
        }
      }
    }
    ylim = tmp
  }
  print(ylim)
  
  pairs(dataf, panel = function(x,y,...){
    y = x - y
    r <- points(RT[!is.na(y)] , y[!is.na(y)],pch=".",... )
    abline(h = 0,col=2)
  }
  ,xlim=range(RT)
  ,ylim=ylim
  ,lower.panel=NULL
  ,main = main
  )
}
#' normal pairs plot with different pch and plus abline
#' @param pln data matrix or data.frame as normally passed to pairs
#' @param ... params usually passed to pairs
#' @export
#' @examples
#' data(SDat)
#' mypairs(SDat$Intensity,log="xy",main="small data")
#' @seealso also \code{\link{pairs}}
mypairs = function(dataframesel,...){
  pairs(dataframesel, panel = function(x,y){
    points(x, y, pch=".")
    abline(a=0,b=1,v=0,h=0,col=2)
  }
  , lower.panel=NULL,...
  )
}
#' heatmap2 facade
#' @export
#' @param pln or dataframe with numerical values
#' @param distf distance function
#' @param hclustf clustering function
#' @param palette
#' @param ColSideColors
#' @examples
#' data(SDat)
#' simpleheatmap(SDat$Intensity,ColSideColors=c("red","blue","pink"))
simpleheatmap = function(pln,main="",
                         distf=dist,
                         hclustf=hclust,
                         palette=terrain.colors(12),
                         ColSideColors=NULL)
{
  tmp <- heatmap.2( as.matrix(pln) , trace="none" , scale="none" , col=palette ,
                    ColSideColors=ColSideColors ,
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
#' @param pvals pvalues
#' @param pthresh pvalue threshold
#' @param ratiothresh threshold of foldchange
#' @param labals - optional labels
#' @param cex size of labels
#' @export
#' @examples
#' foldchange <- rnorm(1000) 
#' pval <-rexp(1000)
#' volcanoplot(foldchange, pval)
#'
volcanoplot = function(foldchange,
                       pvals ,
                       pthresh = 0.05,
                       ratiothresh = 2,
                       xlab ="log2(T/N)" ,
                       ylab = "-log10(P)",
                       labels = NULL,
                       cex=0.6
){
  d <- data.frame(ratio = foldchange, pvals = pvals )
  rownames(d) = labels 
  
  bla = tryCatch( plot(d$ratio,-log10(d$pvals),col="#00000033",pch=19,xlab=xlab, ylab=ylab),
                  warning=function(bla){ dev.off(); return(1)} 
  )
  if(!is.null(bla)){
    plot(d$ratio,-log10(d$pvals),col=1,pch=19,xlab=xlab, ylab=ylab)
  }
  
  
  upsubset<-subset(d,pvals < pthresh & ratio > ratiothresh)
  points(upsubset$ratio,-log10(upsubset$pvals),col=2,pch=19)
  points(upsubset$ratio,-log10(upsubset$pvals),col=1,pch=1)
  if(length(rownames(upsubset)) > 0){
    text(upsubset$ratio, -log10(upsubset$pvals),rownames(upsubset),cex=cex,pos=4)
  }
  
  abline(h=-log10(pthresh),col="gray")
  downsubset<-subset(d,pvals<pthresh & ratio < -ratiothresh)
  points(downsubset$ratio,-log10(downsubset$pvals),col=3,pch=19)
  points(downsubset$ratio,-log10(downsubset$pvals),col=1,pch=1)
  if(length(rownames(downsubset)) > 0){
    text(downsubset$ratio, -log10(downsubset$pvals),rownames(upsubset),cex=cex,pos=2)
  }
  
  abline(v=c(-ratiothresh,ratiothresh),lty=2)
  abline(v =0,lty=2,lwd=1.5)
  return(list(upsubset=upsubset,downsubset=downsubset))
}

