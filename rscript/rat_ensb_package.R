# make rat ensb package
library(ensembldb)
library(EnsDb.Hsapiens.v86)
HSedb<-EnsDb.Hsapiens.v86
library(Biostrings)
library(AnnotationHub)
ah <- AnnotationHub()

## Query all available files for Ensembl release 77 for
## Mus musculus.
query(ah, c("Rattus norvegicus")) # max 92, which annotations 81
RnorV87 <- ah[["AH53239"]]

#save(RnorV87, file="Rnor6_annotation_ensb_87.RData")

# For each of the proteins listed before, query the amino acid sequence in the human homolog as well as the 
## need to load Blist

ratTx<-transcripts(RnorV87,filter=GeneIdFilter(BList$GeneID), columns=c("tx_id","protein_id", "gene_name"))
ratProt<-proteins(RnorV87,filter=GeneIdFilter(BList$GeneID), return.type="AAStringSet")

nAA=5
BList$ProtID=ratTx$protein_id[match(BList$GeneID, ratTx$gene_id)]
Midx=match(BList$ProtID, names(ratProt))

tempProt=ratProt[Midx]
lx=which(!is.na(BList$ProtID))



BList$Seq<-XVector::subseq(ratProt[Midx], BList$AAno-nAA, BList$AAno+nAA)

