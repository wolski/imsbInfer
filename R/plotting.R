## image labels
imagelabels = function(x, labels=colnames(x))
{
  image(x, axes = F )
  axis( 1, at=seq(0,1,length=length(colnames(labels))) , labels=labels,cex.axis=0.8 )
  axis( 2, at=seq(0,1,length=length(colnames(labels))) , labels=labels,cex.axis=0.8 )
}


## altman bland
altmanbland = function(x,y,main){
  idx<-apply(cbind(x,y),1,function(x){sum(x==0)==0})
  mean  = (x+y)/2
  absdiff = abs( x-y )
  scatter.smooth(mean[idx],absdiff[idx],span=1/3,col=2,pch=".",cex=2,xlab="(y+x)/2",ylab="abs(x-y)",main=main,log="xy")  
}


##plot QQ plot
pairsQQ = function(dataframesel,main){
  pairs(dataframesel, panel = function(x,y){
    r <- qqplot(x , y , plot.it = FALSE)
    points(r,pch="*")
    abline(0,1,col=2)
  }
        , lower.panel=NULL
        ,main = main
  )
}

mypairs = function(dataframesel){
  pairs(dataframesel, panel = function(x,y){
    points(x, y, pch=".")
    abline(0,1)
  }
        , lower.panel=NULL
  )
}

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


myvolcanoplot = function(foldchange, pvals , pthresh = 0.5, ratiothresh = 2 ){
  d <- data.frame(ratio = foldchange, pvals = pvals )
  plot(d$ratio,-log10(d$pvals),col="#00000033",pch=19,xlab="log2(T/N)", ylab="-log10(P)")
  
  d2<-subset(d,pvals < pthresh & ratio > ratiothresh)
  points(d2$ratio,-log10(d2$pvals),col=2,pch=19)
  points(d2$ratio,-log10(d2$pvals),col=1,pch=1)
  d2<-subset(d,pvals<pthresh & ratio < -ratiothresh)
  points(d2$ratio,-log10(d2$pvals),col=3,pch=19)
  points(d2$ratio,-log10(d2$pvals),col=1,pch=1)
}
