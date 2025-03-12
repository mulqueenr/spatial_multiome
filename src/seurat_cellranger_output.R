
library(Seurat)
library(ggplot2)
library(patchwork)
library(stringr)
library(reshape2)
library(plyr)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(optparse)
set.seed(123)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="Sample directory from cellranger output.", metavar="character"),
  make_option(c("-n", "--out_name"), type="character", default="sample", 
              help="Prefix outname for sample.", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

indir = opt$input_directory
prefix = opt$out_name
system(paste0("mkdir -p ",outdir))

#read in data
counts <- Read10X_h5(paste0("./outs/filtered_feature_bc_matrix.h5")) #count data

# create a Seurat object containing the RNA data
dat <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA")

dat$sample<-outname

#measure MT
dat <- PercentageFeatureSet(dat, 
    pattern = "^MT-", 
    col.name = "percent.mt")

# run sctransform
dat <- SCTransform(dat, 
    vars.to.regress = "percent.mt", 
    verbose = FALSE)

#run dim reduction
dat<- RunPCA(dat, 
    verbose = FALSE)

dat <- RunUMAP(dat, 
    dims = 1:30, 
    verbose = FALSE)

dat <- FindNeighbors(dat, 
    dims = 1:30,
    verbose = FALSE)

dat <- FindClusters(dat, 
    verbose = FALSE)

plt<-DimPlot(dat, label = TRUE)
ggsave(plt,file=paste0(prefix,".seurat_cluster.pdf"))

saveRDS(dat,file=paste0(prefix,".seuratObj.rds"))