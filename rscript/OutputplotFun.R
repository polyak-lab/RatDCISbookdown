OutputplotFun=function(t2, scaleR="none", main="expression:scaled CD", classN="no", sigMat=NULL){
  ax1=heatmap.2(t2, Rowv = NA, Colv = NA, scale = scaleR, col=RdBu[11:1], trace="none", main=main)
  amat=matrix(0, nrow=ncol(t2), ncol=nrow(t2))
  rownames(amat)=rownames(ax1$carpet)
  colnames(amat)=colnames(ax1$carpet)
  
  if (classN=="wgs"){
    midx2=Cdata$TumorID[match( rownames(tempz2), sapply(strsplit(Cdata$WGS, ".final.bam"), function(x) x[1]))]
    sidx2=sapply(sapply(midx2, function(x) grep(x, rownames(amat))), length)
    sidx=unlist(sapply(midx2, function(x) grep(x, rownames(amat))))
    amat[sidx, ]=tempz2[which(sidx2!=0), match(colnames(amat), colnames(tempz2)) ]
    hs_ax1=which(amat>0, arr.ind = T)
    image(ax1$carpet, col=RdBu[11:1], xaxt="none", yaxt="none", main=main)
    axis(1, at=seq(0, 1, length=ncol(t2)), rownames(amat), las=2, cex=0.7)
    axis(2, at=seq(0, 1, length=nrow(t2)), colnames(amat), las=2, cex=0.7)
    text((hs_ax1[ ,1]-1)/(nrow(amat)-1), na.omit(hs_ax1[ ,2]-1)/(ncol(amat)-1) , "*")
  }
  else{
    image(ax1$carpet, col=RdBu[11:1], xaxt="none", yaxt="none", main=main)
    axis(1, at=seq(0, 1, length=ncol(t2)), rownames(amat), las=2, cex=0.7)
    axis(2, at=seq(0, 1, length=nrow(t2)), colnames(amat), las=2, cex=0.7)
    hs_ax1=which(sigMat[match(colnames(amat), rownames(sigMat)) ,
                        match(rownames(amat), colnames(sigMat))]>0, arr.ind = T)
    text((hs_ax1[ ,2]-1)/(nrow(amat)-1), na.omit(hs_ax1[ ,1]-1)/(ncol(amat)-1) , "*")
  }
  
}