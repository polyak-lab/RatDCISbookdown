# make rat ensb package
library(ensembldb)
library(Biostrings)
library(AnnotationHub)
ah <- AnnotationHub()

## Query all available files for Ensembl release 77 for
## Mus musculus.
query(ah, c("Rattus norvegicus")) # max 92, which annotations 81
RnorV87 <- ah[["AH56709"]] #89 #87: AH53239"]]
#save(RnorV87, file="annotation/Rnor6_annotation_ensb_87.RData")

