ContTable=function(tab, title, chisqtest=F,ylabL="ylab", xlabL="PAM50?"){
  library("RColorBrewer")
  if (chisqtest==T){
    a1=chisq.test(tab)
    tit2=sprintf("Chisq = %g", round(a1$p.value*100)/100)
  } else {
    tit2=" "
  }
  nr=nrow(tab)
  nc=ncol(tab)
  if (chisqtest==T){
    l1=(a1$observed-a1$expected)
    l1[which(l1<0, arr.ind=T)]=0
    image(l1, col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
          xlab=xlabL, ylab=ylabL,
          main=sprintf("%s %s", title, tit2))
  }else{
    image(scale(tab), col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
          xlab=xlabL, ylab=ylabL,
          main=sprintf("%s %s", title, tit2))}
  xval=seq(0, 1, 1/(nr-1))
  yval=seq(0, 1, 1/(nc-1))
  axis(1, at=xval, rownames(tab))
  axis(2, at=yval, colnames(tab))
  for(i in 1:nr){
    for (j in 1:nc){
      text(xval[i], yval[j], tab[i,j], cex=1.5)
    }
  }
}
