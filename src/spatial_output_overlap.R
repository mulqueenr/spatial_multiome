library(Seurat)
library(stringi)


rc <- function(nucSeq){
  stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq))
}

reverse_dna <- function(nucSeq){
  stri_reverse(nucSeq)
}


#count of overlap with cells passing WGS (from complexity metrics for now)
spatial<-"/home/rmulqueen/projects/spatial_wgs/data/250129_First_Experiment/curio_pipeline/output/CuioWGS_DCIS41T_RNA_ConfPositioned_seurat_spatial.rds"
wgs<-"/home/rmulqueen/projects/spatial_wgs/data/250129_First_Experiment/dna_cellranger/dna_cell_metrics.csv"
clones<-"/home/rmulqueen/projects/spatial_wgs/data/250129_First_Experiment/dna_cellranger/copykit/dcis41t.scCNA.tsv"

spatial<-readRDS(spatial)
#spatial in spatial@assays$SPATIAL

wgs<-read.table(wgs,sep=",")
colnames(wgs)<-c("cellid","dna_readcount","dna_duplicate_rate","dna_estimated_size")
row.names(wgs)<-wgs$cellid

clones<-read.table(clones,sep="\t")
row.names(clones)<-unlist(lapply(strsplit(row.names(clones),"[.]"),"[[",1))

clones <- clones[row.names(clones) %in% Cells(spatial),]
wgs<-wgs[row.names(wgs) %in% Cells(spatial)]


clones <- clones[row.names(clones) %in% Cells(spatial),]

#trying rev comp of barcode?
clones$revcomp<-paste0(unlist(lapply(unlist(lapply(strsplit(row.names(clones),"-"),"[[",1)),rc)),"-1")

clones$rev_dna<-paste0(unlist(lapply(unlist(lapply(strsplit(row.names(clones),"-"),"[[",1)),reverse_dna)),"-1")

spatial<-AddMetaData(spatial,clones)