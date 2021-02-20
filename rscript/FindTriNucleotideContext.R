FindTriNucleotideContext=function(fName){
Input=read.delim(fName, sep=" ", header=F)
idx=c(grep("N", Input$V6), grep("-", Input$V5))
if (length(idx)>0){
  Input=Input[-idx, ]
}
Input$V7=substr(Input$V6, 1,1)
Input$V8=substr(Input$V6, 3, 3)
## change the mutation type:
# want C --> X, T--> X
Input$change=paste(Input$V4, ">",Input$V5, sep="")

Input$Var1B=ifelse(Input$V4=="G", "C", ifelse(Input$V4=="A", "T", as.character(Input$V4)))
Input$Var1C=ifelse(Input$V5=="G", "C", ifelse(Input$V5=="A", "T", ifelse(Input$V5=="T", "A", "G")))

Input$change[grep("G>", Input$change)]=paste(paste(Input$Var1B, ">", Input$Var1C, sep=""))[grep("G>", Input$change)]
Input$change[grep("A>", Input$change)]=paste(paste(Input$Var1B, ">", Input$Var1C, sep=""))[grep("A>", Input$change)]

## change the codon if the 
Input$codon=paste(ifelse(Input$V8=="G", "C", ifelse(Input$V8=="A", "T", ifelse(Input$V8=="T", "A", "G"))),
                  Input$Var1B,
                  ifelse(Input$V7=="G", "C", ifelse(Input$V7=="A", "T", ifelse(Input$V7=="T", "A", "G"))), sep="")
Input$codon[which(Input$Var1B==Input$V4)]=as.character(Input$V6[which(Input$Var1B==Input$V4)])

Ax2=table(Input$change, Input$codon)
Ax3=melt(Ax2)
Ax3=Ax3[-which(Ax3$value==0), ]
Ax3
}
