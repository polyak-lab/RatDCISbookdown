## EXtra functions

### 1. MergeContigs
### In a list of genomic intervals, merge a contig if its adjacent to one another

Merge_contig=function(x){
  lx1=x$end[1: (nrow(x)-1)]+1
  lx2=x$start[2:nrow(x)]
  
  midx=which(lx1==lx2)
  
  ConseqIdx=diff(midx)
  Conseqidx2=ifelse(ConseqIdx==1, 1, 0)
  Conseqidx3=rle(Conseqidx2)
  Nrep=Conseqidx3$lengths[which(Conseqidx3$values==1)]
  lx3=which(diff(Conseqidx2)==1)+1
  
  sgroups=lapply(1:length(midx), function(x) c(midx[x], midx[x]+1))
  
  sgroups2=list()
  
  count=1
  i=1
 while (i < length(sgroups)){
    if (i%in%lx3){
      x1=which(lx3==i)
      temp=unique(unlist(c(sgroups[i:(i+Nrep[x1])])))
      sgroups2[[count]]=temp
      i=i+Nrep[x1]+1
      count=count+1
    } else {
      sgroups2[[count]]=sgroups[[i]]
      i=i+1
      count=count+1
    }
  }
  
  LocX=data.frame()
  xrep=x[-unique(unlist(sgroups)), ]
  for (i in 1:length(sgroups2)){
    tx=x[sgroups[[i]], ]
    LocB=tx[1, ]
    LocB$end=tx$end[nrow(tx)]
    LocB$gainProportion=mean(tx$gainProportion)
    LocB$lossProportion=mean(tx$lossProportion)
    LocB$gainFrequency=mean(tx$gainFrequency)
    LocB$lossFrequency=mean(tx$lossFrequency)
    LocX=rbind(LocX, LocB)
  }
    
  Output=rbind(xrep, LocX)
  Output=Output[order(Output$chromosome, Output$start), ]
  
  return(Output)
  
}