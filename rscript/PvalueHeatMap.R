## Function: Heatmap with P value (asterisks)

PvalHM=function(tab, tabP, title="My plot with P values"){
  library("RColorBrewer")
  par(oma=c(5, 5, 0,0))
  image(tab, col=brewer.pal(11, "RdBu")[11:1], xaxt="n", yaxt="n",
          main=title)
  nr=nrow(tab)
  nc=ncol(tab)
  xval=seq(0, 1, 1/(nr-1))
  yval=seq(0, 1, 1/(nc-1))
  axis(1, at=xval, rownames(tab), las=2)
  axis(2, at=yval, colnames(tab), las=2)
  for(i in 1:nr){
    for (j in 1:nc){
      text(xval[i], yval[j], tabP[i,j], cex=1.5)
    }
  }
}
