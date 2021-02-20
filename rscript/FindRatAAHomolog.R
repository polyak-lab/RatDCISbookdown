# library(biomaRt)
# #library(refGenome)
# 
# TS1=unique(MeltCosmic$Gene_Symbol)
# TS2=unique(c(as.character(ChangList$Gene), as.character(BaileyList$Gene)))
# TS2=setdiff(TS2, unique(Rat2HumProt$HGNC.symbol))
# Rat2HumProt = getLDS(attributes = c("rgd_symbol", "ensembl_peptide_id"), filters = "rgd_symbol",  values = TS1 , mart = rat, attributesL = c("hgnc_symbol","ensembl_peptide_id"), martL = human, uniqueRows=T)
# 
# library(ensembldb)
# library(EnsDb.Hsapiens.v86)
# Hsedb=EnsDb.Hsapiens.v86
# 
# GeneHum=BaileyList$Gene
# AA1=substr(BaileyList$Mutation, 3, 3)
# AA2=substr(BaileyList$Mutation, nchar(BaileyList$Mutation), nchar(BaileyList$Mutation))
# Loc=substr(BaileyList$Mutation, 4, nchar(BaileyList$Mutation)-1)
# 
# AA1=C2List$AA1  
# AA2=C2List$AA2
# Loc=C2List$Loc
# GeneHum=C2List$Gene

FindRatAAHomolog=function(GeneHum, AA1, Loc, AA2){

rat = useEnsembl("ensembl", mirror="asia", data="rnorvegicus_gene_ensembl")
human = useEnsembl("ensembl", mirror="asia", dataset = "hsapiens_gene_ensembl")
Hum2RatProt = getLDS(attributes = c("hgnc_symbol", "ensembl_peptide_id", "ensembl_transcript_id"), filters = "hgnc_symbol",  
                     values = GeneHum , mart = human, attributesL = c("rgd_symbol","ensembl_peptide_id", "ensembl_transcript_id"), martL = rat, uniqueRows=T)
colnames(Hum2RatProt)[c(2,4)]=c("Protein.stable.ID.1", "Protein.stable.ID")
RatGene=SymHum2Rat$RGD.symbol[match(GeneHum, SymHum2Rat$HGNC.symbol)]

naid=which(is.na(RatGene)| RatGene=="")
RatGene[naid]=tolower(GeneHum[naid])
RatGene[naid]=firstup(RatGene[naid])

# AA1=substr(BaileyList$Mutation, 3, 3)
# AA2=substr(BaileyList$Mutation, nchar(as.character(BaileyList$Mutation)), nchar(as.character(BaileyList$Mutation)))

ListN=data.frame(Gene=GeneHum, AAno=as.numeric(Loc), AA1=AA1, Variant=AA2, RatGene=RatGene)
#rbind(ChangList[ ,c("Gene", "Residue", "Variants", "RatGene", "AA1", "AAno")], B2[ ,c("Gene", "Residue", "Variants", "RatGene", "AA1", "AAno")])

## also figure out whether the human mutation sites are at the right place?
txHuman=transcripts(Hsedb, filter=GeneNameFilter(na.omit(ListN$Gene)), columns = c("protein_id","gene_id", "tx_id"))
humProt<-proteins(Hsedb,filter=GeneNameFilter(na.omit(ListN$Gene)), return.type="AAStringSet")

ListN$Sequence=NA
ListN$HumProt=NA #names(idx1)[]

## Find the Human sequence and Protein

for (i in 1:nrow(ListN)){
  idVal=na.omit(txHuman$protein_id[which(txHuman$gene_name==ListN$Gene[i])])
  tempA=humProt[which(names(humProt)%in%idVal)]
  tempA=tempA[which(width(tempA)>as.numeric(ListN$AAno[i]))]
  AA=as.character(subseq(tempA, as.numeric(ListN$AAno[i]), as.numeric(ListN$AAno[i])))
  idx1=which(AA==ListN$AA1[i])
  #SubStrSeq=subseq(tempA[idx1[1]], as.numeric(ChangList$AAno[i])-5, as.numeric(ChangList$AAno[i])+5)
  if (length(idx1)>0){
    ax=as.numeric(ListN$AAno[i])-5
    aend=as.numeric(ListN$AAno[i])+5
    if (ax<1){
      ax=1
    }
    if (aend>width(tempA)[idx1[1]]){
      aend=width(tempA)[idx1[1]]
    }
    ListN$Sequence[i]=subseq(tempA[idx1[1]], ax, aend)
    ListN$HumProt[i]=names(AA[idx1[1]])
  }
}


## Find the corresponding rat sequence and peotein
#rmidx=which(is.na(ListN$RatGene))
ratTx=transcripts(RnorV87,filter=GeneNameFilter(ListN$RatGene), columns=c("tx_id","protein_id", "gene_name"))
ratProt<-proteins(RnorV87,filter=GeneNameFilter(ListN$RatGene), return.type="AAStringSet")

nAA=5
# rematch the cases where the length is not as long as the width - may need to double check?
# go through and search for the amino acid and check the values match

ListN$RatProt=NA
ListN$RatAAno=NA
ListN$RatSequence=NA

for (i in 1:nrow(ListN)){
  if (!is.na(ListN$RatGene[i]) & !is.na(ListN$HumProt[i])){
  idVal=na.omit(ratTx$protein_id[which(ratTx$gene_name%in%ListN$RatGene[i])])
  tempA=ratProt[which(names(ratProt)%in%idVal)]
  testAlignment=lapply(tempA, function(x) matchPattern(ListN$Sequence[i],x, min.mismatch =0, max.mismatch = 2))
  ## up to here
   
  lx1=sapply(testAlignment, length)
  
  m1=which(lx1>1)
  if (length(m1)>0){
    testAlignment=lapply(tempA, function(x) matchPattern(ListN$Sequence[i],x, min.mismatch =0, max.mismatch = 1))
  }
  
  lx1=sapply(testAlignment, length)
  idx=which(lx1>0)
  
  if (length(idx)>0){
  testAlignment=testAlignment[idx]
  tempA=tempA[idx]
  
  checkidx=sapply(testAlignment, width)
  startidx=sapply(testAlignment, start)
  endidx=sapply(testAlignment, end)
  
  ## find example where more than 1 is possible
  # 207: skip if no sequence present
  # tempA=tempA[which(width(tempA)>=ListN$AAno[i])]
  # AA=as.character(subseq(tempA, MeltCosmic$AAno[i], MeltCosmic$AAno[i]))
  # useIdx=which(AA==MeltCosmic$AAwt[i])
  # if (length(useIdx)>0){
  #   MeltCosmic$ProtID[i]=idVal[useIdx[1]]
  #   ax1=MeltCosmic$AAno[i]-5
  #   if (ax1<1){
  #     ax1=1
  #   }
  #   MeltCosmic$Seq[i]=as.character(subseq(tempA[useIdx[1]],ax1 , MeltCosmic$AAno[i]+5))
  # }
  ##
  bx=startidx-ListN$AAno[i]
  a1=dim(bx)
  ax=which.min(abs(bx))

  if (length(ax)>0){
      ListN$RatProt[i]=names(tempA)[ax]
      ListN$RatAAno[i]=mean(c(unlist(endidx[ax]), unlist(startidx[ax])))
      ListN$RatSequence[i]=subseq(tempA[ax], unlist(startidx[ax]), unlist(endidx[ax]))
        if (checkidx[ax]!=11){
            shiftAA=(11-checkidx[ax])
            if (ListN$RatAAno[i]<6){
            ListN$RatAAno[i]=round(mean(c(unlist(endidx[ax]), unlist(startidx[ax]-shiftAA))))
            } else {
            ListN$RatAAno[i]=round(mean(c(unlist(endidx[ax]+nAA), unlist(startidx[ax]))))
            }
    }
  }
}}
#naidx=(which(is.na(MeltCosmic$Seq)))
}
ListN
}

## do we actually need this:


FindHumanAAHomolog=function(GeneRat, AA1, Loc, AA2){
  
  human = useEnsembl("ensembl", mirror="asia",  dataset = "hsapiens_gene_ensembl")
  rat = useEnsembl("ensembl",mirror="asia",   data="rnorvegicus_gene_ensembl")
  Rat2HumProt = getLDS(attributesL = c("hgnc_symbol", "ensembl_peptide_id", "ensembl_transcript_id"), filters = "rgd_symbol",  
                       values = unique(GeneRat), martL = human, attributes = c("rgd_symbol","ensembl_peptide_id", "ensembl_transcript_id"), mart = rat, uniqueRows=T)
  #colnames(Hum2RatProt)[c(2,4)]=c("Protein.stable.ID.1", "Protein.stable.ID")
  
  RatGene=Rat2HumProt$HGNC.symbol[match(GeneRat, Rat2HumProt$RGD.symbol)]
  # 
  naid=which(is.na(RatGene)| RatGene=="")
  RatGene[naid]=toupper(GeneRat[naid])
 # RatGene[naid]=firstup(RatGene[naid])
  
  # AA1=substr(BaileyList$Mutation, 3, 3)
  # AA2=substr(BaileyList$Mutation, nchar(as.character(BaileyList$Mutation)), nchar(as.character(BaileyList$Mutation)))
  
  ListN=data.frame(Gene=GeneRat, AAno=as.numeric(Loc), AA1=AA1, Variant=AA2, HumGene=RatGene)
  #rbind(ChangList[ ,c("Gene", "Residue", "Variants", "RatGene", "AA1", "AAno")], B2[ ,c("Gene", "Residue", "Variants", "RatGene", "AA1", "AAno")])
  
  ## also figure out whether the human mutation sites are at the right place?
  txHuman=transcripts(Hsedb, filter=GeneNameFilter(na.omit(ListN$HumGene)), columns = c("protein_id","gene_id", "tx_id", "gene_name"))
  humProt<-proteins(Hsedb,filter=GeneNameFilter(na.omit(ListN$HumGene)), return.type="AAStringSet")
  ratTx=transcripts(RnorV87,filter=GeneNameFilter(ListN$Gene), columns=c("tx_id","protein_id", "gene_name"))
  ratProt<-proteins(RnorV87,filter=GeneNameFilter(ListN$Gene), return.type="AAStringSet")
  
  ListN$RatSequence=NA
  ListN$RatProt=NA
  ListN$HumSequence=NA
  ListN$HumProt=NA #names(idx1)[]
  ListN$HumAAno=NA
  nAA=5
  ## Find the Human sequence and Protein
  
  for (i in  1:nrow(ListN)){
    idVal=na.omit(ratTx$protein_id[which( ratTx$gene_name==ListN$Gene[i])])
    tempA=ratProt[which(names(ratProt)%in%idVal)]
    tempA=tempA[which(width(tempA)>as.numeric(ListN$AAno[i]))]
    AA=as.character(subseq(tempA, as.numeric(ListN$AAno[i]), as.numeric(ListN$AAno[i])))
    idx1=which(AA==ListN$AA1[i])
    
    #browser()
    #SubStrSeq=subseq(tempA[idx1[1]], as.numeric(ChangList$AAno[i])-5, as.numeric(ChangList$AAno[i])+5)
    if (length(idx1)>0){
      ax=as.numeric(ListN$AAno[i])-nAA
      aend=as.numeric(ListN$AAno[i])+nAA
      if (ax<1){
        ax=1
      }
      if (aend>width(tempA)[idx1[1]]){
        aend=width(tempA)[idx1[1]]
      }
      ListN$RatSequence[i]=subseq(tempA[idx1[1]], ax, aend)
      ListN$RatProt[i]=names(AA[idx1[1]])
      
      ## look for the human homolog
        idVal=na.omit(txHuman$protein_id[which(txHuman$gene_name%in%ListN$HumGene[i])])
        tempA=humProt[which(names(humProt)%in%idVal)]
        testAlignment=lapply(tempA, function(x) matchPattern(ListN$RatSequence[i],x, min.mismatch =0, max.mismatch = 2))
      ## up to here
        lx1=sapply(testAlignment, length)
        idx=which(lx1==1)
        
        if (length(idx)>0){
          testAlignment=testAlignment[idx]
          tempA=tempA[idx]
          
          checkidx=sapply(testAlignment, width)
          startidx=sapply(testAlignment, start)
          endidx=sapply(testAlignment, end)
          bx=sapply(startidx, function(x) x-ListN$AAno[i])
          ax=which.min(abs(unlist(bx)))
          
          if (sum(sapply(checkidx, length))!=length(checkidx)){
               ax1=substr(names(ax), 1, nchar(names(ax))-1)
          } else {
            ax1=names(ax)
          }
          
          if (length(ax)>0){
            ListN$HumProt[i]=ax1
            ListN$HumAAno[i]=mean(c(unlist(endidx)[ax], unlist(startidx)[ax]))
            ListN$HumSequence[i]=subseq(tempA[grep(ax1, names(tempA))], unlist(startidx)[ax], unlist(endidx)[ax])
            if (unlist(checkidx)[ax]!=11){
              shiftAA=(11-unlist(checkidx)[ax])
              if (ListN$HumAAno[i]<6){
                ListN$HumAAno[i]=round(mean(c(unlist(endidx[ax]), unlist(startidx[ax]-shiftAA))))
              } else {
                ListN$HumAAno[i]=round(mean(c(unlist(endidx[ax]+nAA), unlist(startidx[ax]))))
              }
            }
          }
        }}
      
  }
  ListN
}
  
  ## Find the corresponding rat sequence and peotein
  #rmidx=which(is.na(ListN$RatGene))
   
  

#####
## Annotate with hotspot data: have any of these mutations been published?
#####
#Rat2HumProt=rbind(Rat2HumProt, Hum2RatProt)

# published lists: Bailey and Chang

# BaileyList$AA1=substr(BaileyList$Mutation, 3, 3)
# BaileyList$AA2=substr(BaileyList$Mutation, nchar(as.character(BaileyList$Mutation)), nchar(as.character(BaileyList$Mutation)))
# BaileyList$AAno=substr(BaileyList$Mutation, 4, nchar(as.character(BaileyList$Mutation))-1)
# BaileyList$RatGene=SymHum2Rat$RGD.symbol[match(BaileyList$Gene, SymHum2Rat$HGNC.symbol)]
# 
# ChangList$RatGene=SymHum2Rat$RGD.symbol[match(ChangList$Gene, SymHum2Rat$HGNC.symbol)]
# ChangList$AA1=substr(ChangList$Residue, 1, 1)
# ChangList$AAno=substr(ChangList$Residue, 2, 7)
# 
# ## figure out the variants
# ChangComb=paste(ChangList$Gene, ChangList$Residue, sep="")
# BaileyComb=paste(BaileyList$Gene, BaileyList$AA1, BaileyList$AAno, sep="")
# Ax=match(BaileyComb, ChangComb)
# B2=BaileyList[which(is.na(Ax)), ]
# B2$Residue=paste(B2$AA1, B2$AAno, sep="")
# B2$Variants=B2$AA2
# 
# ChangListN=rbind(ChangList[ ,c("Gene", "Residue", "Variants", "RatGene", "AA1", "AAno")], B2[ ,c("Gene", "Residue", "Variants", "RatGene", "AA1", "AAno")])
# 
# 
# ## also figure out whether the human mutation sites are at the right place?
# txHuman=transcripts(HSedb, filter=GeneNameFilter(na.omit(ChangListN$Gene)))
# humProt<-proteins(HSedb,filter=GeneNameFilter(na.omit(ChangListN$Gene)), return.type="AAStringSet")
# 
# ChangListN$Sequence=NA
# 
# for (i in 1:nrow(ChangListN)){
#   idVal=na.omit(Rat2HumProt$Protein.stable.ID.1[which(Rat2HumProt$HGNC.symbol==ChangListN$Gene[i])])
#   tempA=humProt[which(names(humProt)%in%idVal)]
#   tempA=tempA[which(width(tempA)>as.numeric(ChangListN$AAno[i]))]
#   AA=as.character(subseq(tempA, as.numeric(ChangListN$AAno[i]), as.numeric(ChangListN$AAno[i])))
#   idx1=which(AA==ChangListN$AA1[i])
#   #SubStrSeq=subseq(tempA[idx1[1]], as.numeric(ChangList$AAno[i])-5, as.numeric(ChangList$AAno[i])+5)
#   if (length(idx1)>0){
#     ax=as.numeric(ChangListN$AAno[i])-5
#     aend=as.numeric(ChangListN$AAno[i])+5
#     if (ax<1){
#       ax=1
#     }
#     if (aend>width(tempA)[idx1[1]]){
#       aend=width(tempA)[idx1[1]]
#     }
#     ChangListN$Sequence[i]=subseq(tempA[idx1[1]], ax, aend)
#   }
# }
# 
# ##########
# ## collect the AA number and value in each
# ## Obtain the amino acid sequence
# #########
# AAout=strsplit(as.character(MeltCosmic$HGVSp_Short), "[0-9]+")
# AAwt=sapply(AAout, function(x) x[1])
# AAmt=sapply(AAout, function(x) x[2])
# AAno=regmatches(MeltCosmic$HGVSp_Short, gregexpr("[[:digit:]]+", MeltCosmic$HGVSp_Short))
# 
# MeltCosmic$AAwt=AAwt
# MeltCosmic$AAmt=AAmt
# MeltCosmic$AAno=NA
# MeltCosmic$AAno[which(sapply(AAno, length)==1)]=as.numeric(as.character(unlist(AAno)))
# 
# ratTx<-transcripts(RnorV87,filter=GeneNameFilter(MeltCosmic$Gene_Symbol), columns=c("tx_id","protein_id", "gene_name"))
# ratProt<-proteins(RnorV87,filter=TxNameFilter(ratTx$tx_id), return.type="AAStringSet")
# 
# nAA=5
# # rematch the cases where the length is not as long as the width - may need to double check?
# # go through and search for the amino acid and check the values match
# 
# for (i in 1:length(MeltCosmic$Gene_Symbol)){
#   idVal=na.omit(ratTx$protein_id[which(ratTx$gene_name==MeltCosmic$Gene_Symbol[i])])
#   tempA=ratProt[which(names(ratProt)%in%idVal)]
#   tempA=tempA[which(width(tempA)>=MeltCosmic$AAno[i])]
#   AA=as.character(subseq(tempA, MeltCosmic$AAno[i], MeltCosmic$AAno[i]))
#   useIdx=which(AA==MeltCosmic$AAwt[i])
#   if (length(useIdx)>0){
#     MeltCosmic$ProtID[i]=idVal[useIdx[1]]
#     ax1=MeltCosmic$AAno[i]-5
#     if (ax1<1){
#       ax1=1
#     }
#     MeltCosmic$Seq[i]=as.character(subseq(tempA[useIdx[1]],ax1 , MeltCosmic$AAno[i]+5))
#   }}
# 
# naidx=(which(is.na(MeltCosmic$Seq)))
# 
# ########
# # Look at the amino acid sequences: 
# # 1. convert to one letter symbol
# # 2. Figure out overlap with published lists
# #######
# 
# ## Get Human protein sequences
# HumProt=proteins(HSedb,filter=ProteinIdFilter(unique(Rat2HumProt$Protein.stable.ID.1)), return.type="AAStringSet")
# ## match each protein sequence rat to human and search for matching peptide:
# 
# MeltCosmic$HumSeq=NA
# MeltCosmic$HumProt=NA
# MeltCosmic$HumStart=NA
# MeltCosmic$HumEnd=NA
# MeltCosmic$HumAAno=NA
# MeltCosmic$HumAA=NA
# 
# MeltCosmic$PossHuman=NA
# MeltCosmic$CompSeq=NA
# 
# for (i in 1:nrow(MeltCosmic)){
#   if (!is.na(MeltCosmic$Seq[i])){
#     HumID=Rat2HumProt$Protein.stable.ID.1[which(Rat2HumProt$Protein.stable.ID==MeltCosmic$ProtID[i])]
#     tempA=HumProt[names(HumProt)%in%HumID]
#     
#     if (length(tempA)>0){
#       t2=NULL
#       count=0
#       
#       while (length(t2)<1 & count<4 ){
#         testAlignment=vmatchPattern(MeltCosmic$Seq[i], tempA, min.mismatch =0, max.mismatch = count) # first do it with 0 mismatch
#         checkidx=abs(start(testAlignment)-MeltCosmic$AAno[i]+5)
#         t2=which.min(unlist(checkidx))
#         checkidxLength=sapply(checkidx, length)
#         count=count+1
#       }
#       if (length(t2)>0){
#         if (length(which(checkidxLength>1))>0){
#           testAlignment=testAlignment[[names(t2)]][t2]
#         } else {
#           testAlignment=testAlignment[[names(t2)]]
#         }
#         tempB=subseq(tempA[names(t2)], start=start(testAlignment), end=end(testAlignment))
#         MeltCosmic$HumSeq[i]=as.character(tempB)
#         MeltCosmic$HumProt[i]=names(t2)
#         MeltCosmic$HumStart[i]=start(testAlignment)
#         MeltCosmic$HumEnd[i]=end(testAlignment)
#         MeltCosmic$HumAA[i]=substr(MeltCosmic$HumSeq[i], nchar(MeltCosmic$HumSeq[i])-5, nchar(MeltCosmic$HumSeq[i])-5)
#         MeltCosmic$HumAAno[i]=MeltCosmic$HumEnd[i]-5
#       }}
#     # check the ChangList
#     Ctemp=ChangListN[which(ChangListN$RatGene==MeltCosmic$Gene_Symbol[i] & !is.na(ChangListN$Sequence)), ]
#     if (nrow(Ctemp)>0){
#       #  Ctemp=Ctemp[-which(is.na(Ctemp$Sequence)), ]
#       testAlignment=vmatchPattern(MeltCosmic$Seq[i], Ctemp$Sequence, min.mismatch =0, max.mismatch = 3)
#       checkidx=sapply(testAlignment, length)
#       t2=which(checkidx>0)
#       if (length(t2)>0){
#         MeltCosmic$PossHuman[i]=paste(Ctemp$Residue[t2], collapse=" ")
#         MeltCosmic$CompSeq[i]=paste(Ctemp$Sequence[t2], collapse=" ")
#       }
#     }
#   }
# }
# 
# MeltCosmic$HumanProteinName=Rat2HumProt$HGNC.symbol[match(MeltCosmic$HumProt, Rat2HumProt$Protein.stable.ID.1)]
# 
# MCcomb=paste(MeltCosmic$HumanProteinName, MeltCosmic$HumAA, MeltCosmic$HumAAno, sep="")
# 
# MeltCosmic$HumanCancerMut=NA
# MeltCosmic$HumanCancerMut[MCcomb%in%ChangComb]=1
# MeltCosmic$HumanCancerMut[MCcomb%in%BaileyComb]=1
# 
