########
# load in matrix data and evaluate the composition and spatial location
# folder with all data
#######
FilePath="Raw-matrix/Ep-SMA-switch/"
FileList=dir(FilePath, ".txt")
FNames=sapply(strsplit(FileList, "[_|.txt]"	), function(x) x[1])

####
# load in info about Ids that need to be switched.
# Exmplae:
# Sample ID | SMA:EpCAM | SMA
# means to re-label all SMA:EpCAM to SMA in sample SampleID
# Sample ID | EpCAM | SMA : means the EpCAM/SMA fractions have been switched
####
switchTab=read.csv("~/Dropbox (Partners HealthCare)/Scanned slides/Matrix/run_script/rename_id_list.csv", header=F)



library(RColorBrewer)
palette(brewer.pal(9, "Set1"))

#FinalMatrix="~/Dropbox (Partners HealthCare)/Carlos' leukocyte project/Matrix/Locational-infoFinal/"
#FileList2=dir(FinalMatrix, ".csv")
#FFNames2=sapply(strsplit(FileList, "[_|.csv]"	), function(x) x[1])
#idx1=match(FNames, FFNames2)


for (m in 1:length(FileList)){

FileIn1=paste(FilePath, FileList[m], sep="")
#FileIn2=paste(FinalMatrix, FileList2[idx1[m]], sep="")
FName=FNames[m]

print(sprintf('running sample %s', FName))
#D1=read.csv(FileIn2, stringsAsFactors = F)
D1=read.delim(FileIn1, sep="\t", stringsAsFactors=F)

lx1=which(D1$Parent=="Tumor")

if (length(lx1)>1E4){
	
#D1=D1[order(D1$Centroid.X.µm, D1$Centroid.X.µm), ]	
#D2=D2[order(D2$Centroid.X.µm, D2$Centroid.X.µm), ]	

D1=D1[which(D1$Parent=="Tumor"), ]
D1$Class[which(D1$Class=="")]="Unclass"

mx=match(FName, switchTab$V1)

if (!is.na(mx)){
a1=switchTab[ mx, 2]

if (a1=="SMA:EpCAM"){
	x1=which(D1$Class=="EpCAM: SMA")
	D1$Class[x1]=as.character(switchTab$V3[mx])
	FName=paste(FName, "DP-switch", sep="_")
} else {
	x1=which(D1$Class=="EpCAM")
	x2=which(D1$Class=="SMA")
	D1$Class[x2]="EpCAM"
	D1$Class[x1]="SMA"
	FName=paste(FName, "Ep-SMA-switch", sep="_")
}	
}

########################
# Summary boxplots of the cell types and size, EpCam, CD8 and SMA
########################

D1$Class2=D1$Class
D1$Class2[grep("CD8", D1$Class2)]="CD8"
Lx1=table(D1$Class2)
x1=names(Lx1)[which(Lx1==1)]
D1$Class2[which(D1$Class2==x1)]="Unclass"

colA=palette()[factor(D1$Class2)]
termsA=levels(factor(D1$Class2))
colB=colA[match(levels(factor(D1$Class)), D1$Class)]
## Write the first pdf
pdf(sprintf("%s_%s_image_summary_stats_rm_borders.pdf", FName, Sys.Date()), height=8, width=8)
par(mfrow=c(2,3))
boxplot(D1$Nucleus..Area~D1$Class, las=2, ylab="Nucleus Area", col=colB)
boxplot(D1$Cell..Area~D1$Class, las=2, ylab="Cell Area", col=colB)
boxplot(D1$Nucleus..Circularity~D1$Class, las=2, ylab="Nucleus Circularity", col=colB)
boxplot(D1$Cell..Circularity~D1$Class, las=2, ylab="Cell Circularity", col=colB)
boxplot(D1$Cell..Eccentricity~D1$Class, las=2, ylab="Cell Eccentricity", col=colB)
boxplot(D1$Nucleus.Cell.area.ratio~D1$Class, las=2, ylab="Nucleus Cell Area Ratio", col=colB)
## Staining Patterns
par(mfrow=c(2,2))
boxplot(D1$Cell..CD8.mean~D1$Class, las=2, ylab="CD8mean", col=colB)
boxplot(D1$Cell..DAPI.mean~D1$Class, las=2, ylab="DAPI mean", col=colB)
boxplot(D1$Cell..EpCAM.mean~D1$Class, las=2, ylab="EpCAM mean", col=colB)
boxplot(D1$Cell..SMA.mean~D1$Class, las=2, ylab="SMA mean", col=colB)

par(mfrow=c(2,2))
barplot(Lx1, las=2, ylab="number of cells", col=colB)
barplot(Lx1/sum(Lx1), las=2, ylab="fraction cells", col=colB)
barplot(Lx1[-length(Lx1)]/sum(Lx1[-length(Lx1)]), las=2, ylab="fraction cells rm Unclass", col=colB)
dev.off()

#########################
## Overview of the image
#########################

# Label any CD8 as one class:

pdf(sprintf("%s_%s_image_overview_rm_borders.pdf", FName, Sys.Date()), height=10, width=10)
## also plot a zoom in: 
xyvals=round(c(mean(D1$Centroid.X.µm), mean(D1$Centroid.Y.µm)))
par(mar = c(6, 3, 2, 2),xpd = TRUE)
plot(D1$Centroid.X.µm, D1$Centroid.Y.µm, col=factor(D1$Class2), pch=".", xlab="X co-ord", ylab="Y co-ord", main=FName)
legend("bottomright", inset=c(0.1, -0.05), levels(factor(D1$Class2)), pch=19, col=c(1:length(levels(factor(D1$Class2)))), horiz=T)
par(mar = c(2, 3, 2, 2),xpd = F)
plot(D1$Centroid.X.µm, D1$Centroid.Y.µm, col=factor(D1$Class2), pch=20, xlab="X co-ord", ylab="Y co-ord", 
     xlim=c(xyvals[1]-500, xyvals[1]+500), ylim=c(xyvals[2]-500, xyvals[2]+500))
legend("bottomright", termsA, col=c(1:length(termsA)), pch=19, lwd=2)
dev.off()

#########################
## Spatial statistics
#########################

library(spatstat)
library(abind)
library(reshape2)

windowOut=ripras(D1$Centroid.X.µm, D1$Centroid.Y.µm)
ppOutR=ppp(D1$Centroid.X.µm, D1$Centroid.Y.µm, marks=factor(D1$Class2), poly=windowOut$bdry)

#####################################
## Part A: Nearest Neighbour analysis
#####################################
knnOverview=nndist(ppOutR, by=marks(ppOutR), k=c(1:5))

termsB=paste("^",termsA, ".dist", sep="")

knn1=knnOverview[ ,grep("dist.1", colnames(knnOverview))]
knn3=sapply(termsB, function(x) rowMeans(knnOverview[ ,grep(x, colnames(knnOverview))[1:3]]))
knn5=sapply(termsB, function(x) rowMeans(knnOverview[ ,grep(x, colnames(knnOverview))]))
knnDistMatrix=abind(knn1, knn3, knn5, along=3)

cutoffdist=seq(10, 30, by=5)

ProportionMat=list()
KStestpvalues=list()  
MeanDist=list()
MeanDist2=list()

 pdf(sprintf("%s_%s_knn_dist_rm_borders.pdf", FName, Sys.Date()), width=9, height=5)
for (i in 1:length(termsA)){
  xind=which(marks(ppOutR)==termsA[i])
  par(mfrow=c(1, 3))
  ## create a confidence interval
  
  ptest=matrix(NA, ncol=length(termsA), nrow=3)
  colnames(ptest)=termsA
  rownames(ptest)=c("knn1", "knn3", "knn5")
  mdist=ptest
  mdist2=ptest
  prop=matrix(NA, ncol=length(termsA), nrow=3*length(cutoffdist))
  colnames(prop)=termsA
  rownames(prop)=paste(rep(rownames(ptest), each=length(cutoffdist)), cutoffdist, sep="-")
  for (k in 1:3){
   plot(NA, xlim=c(0, 200), ylim=c(0, 1), main=paste(termsA[i], rownames(ptest)[k]), xlab="distance")
    legend("bottomright", termsA, col=c(1:length(termsA)), pch=19, lwd=1)
  for (j in 1:length(termsA)){
    ptest[k, j]=ks.test(knnDistMatrix[xind, 1, k], knnDistMatrix[xind, j, k])$p.value
    prop[(5*(k-1)+1):(5*k), j]=sapply(cutoffdist, function(x) length(which(knnDistMatrix[xind, j, k]<x)))
    mdist[k,j]=mean(knnDistMatrix[xind, j, k])
    mdist2[k,j]=median(knnDistMatrix[xind, j, k])
    lines(ecdf(knnDistMatrix[xind, j, k]), lwd=1, col=j)
      }
  }
  KStestpvalues[[i]]=ptest
  ProportionMat[[i]]=prop/length(xind)
  MeanDist[[i]]=mdist
  MeanDist2[[i]]=mdist2
  
  }
dev.off()

names(KStestpvalues)=termsA
names(ProportionMat)=termsA
names(MeanDist)=termsA
names(MeanDist2)=termsA

ax1=melt(KStestpvalues)
ax2=melt(ProportionMat)
ax3=melt(MeanDist)
ax4=melt(MeanDist2)

Ax3=cbind(ax1, ax3[ ,3], ax4[ ,3])
colnames(Ax3)=c("knn", "NearestCellType", "KS.pval", "CellType", "MeanDistance", "MedianDistance")


write.csv(Ax3, file=sprintf("%s_%s_knn_Metrics_rm_borders.csv", FName, Sys.Date()))
write.csv(ax2, file=sprintf("%s_%s_ProportionInteracting_rm_borders.csv", FName, Sys.Date()))

##############################
## Part B: M-H testing
##############################
library(divo)

sizetest=seq(50, 500, by=50)

Xout=data.frame()
SizeSum=data.frame()

for (gridsize in sizetest){
#gridsize=100 ## can change this gradually

xlimV=c(min(D1$Centroid.X.µm), max(D1$Centroid.X.µm))
ylimV=c(min(D1$Centroid.Y.µm), max(D1$Centroid.Y.µm))

xsize=diff(xlimV)/1000
ysize=diff(ylimV)/1000

xcut=seq(xlimV[1], xlimV[2], length=round(diff(xlimV)/gridsize)+1)
ycut=seq(ylimV[1], ylimV[2], length=round(diff(ylimV)/gridsize)+1)

ppTessR=tess(ppOutR, xgrid=xcut, ygrid=ycut)
nx=lapply(termsA, function(x) quadratcount(ppOutR[ppOutR$marks==x], tess=ppTessR))
NXlist=sapply(nx, function(x) as.vector(x))
colnames(NXlist)=termsA

## metrics to report out
Nexpected=colMeans(NXlist)
Ntiles=ppTessR$n
SizeSum=rbind(SizeSum, c(gridsize, Ntiles, Nexpected))
##
xres<-mh(NXlist, resample=1000)
xres2=melt(xres)
xres3=cbind(xres2[ xres2$L1=="Mean", 1:3], xres2[ xres2$L1=="Lower.Quantile", 3], xres2[ xres2$L1=="Upper.Quantile", 3], size=gridsize)
xres3=xres3[-which(xres3[ ,3]==0 | xres3[ ,3]==1), ]
Xout=rbind(Xout, xres3)
}

colnames(SizeSum)=c("gridsize", "Ntiles", paste("Nexpected", termsA, sep="."))
write.csv(SizeSum, file=sprintf("%s_%s_MH-setup-summary_rm_borders.csv", FName, Sys.Date()), row.names = F)
colnames(Xout)=c("Var1", "Var2", "MH.mean", "MH.lower", "MH.upper", "gridsize")
write.csv(Xout, file=sprintf("%s_%s_MH-values_rm_borders.csv", FName, Sys.Date()), row.names = F)
Xout$comp=paste(Xout$Var1, Xout$Var2)

#Xout$SE=(Xout$MH.upper-Xout$MH.lower)/(2*1.96)
#Xout$z=Xout$MH.mean/Xout$SE

library(ggplot2)

pdf(sprintf("%s_%s_MH_summary_rm_borders.pdf", FName, Sys.Date()), width=10, height=8)
p<-ggplot(Xout, aes(x=gridsize, y=MH.mean, col=comp))+geom_point()+geom_errorbar(aes(ymin=MH.lower, ymax=MH.upper))
print(p)
p2<-ggplot(Xout, aes(x=gridsize, y=MH.mean, col=Var1))+geom_point()+geom_errorbar(aes(ymin=MH.lower, ymax=MH.upper))+facet_wrap(~Var2, scale="free")
print(p2)
dev.off()


write.csv(D1[ ,c("Class", "Centroid.X.µm", "Centroid.Y.µm", "Class2")], file=sprintf("%s_%s_edited_raw_matrix_rm_borders.csv", FName, Sys.Date()))

}}