library(copykit)
library(BiocParallel)
library(optparse)
#library(colorRamp2)
BiocParallel::bpparam()

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="List of single-cell bam files", metavar="character"),
  make_option(c("-p", "--output_prefix"), type="character", default="dcis41t", 
              help="Prefix of output"),
  make_option(c("-c", "--task_cpus"), type="integer", default=300, 
              help="Integer number of cpus")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
cpu_count=opt$task_cpus
prefix=opt$output_prefix

register(MulticoreParam(progressbar = T, workers = cpu_count), default = T)

genome <- "hg38"
resolution <- "220kb" #running at real low res for now
min_bincount = 10

dat<-runVarbin(".", remove_Y=TRUE,is_paired_end=TRUE,resolution="220kb")
dat  <- runMetrics(dat)

pdf(paste0(prefix,".qc_metrics.pdf"))
plotMetrics(dat, metric = c("overdispersion", 
                              "breakpoint_count",
                              "reads_total",
                              "reads_duplicates",
                              "reads_assigned_bins",
                              "percentage_duplicates"),label = "reads_total")
dev.off()

#dat <- dat[,colData(dat)$reads_assigned_bins >= 100000]
colData(dat)$sample_name<-prefix
dat <- findAneuploidCells(dat)
dat <- findOutliers(dat)

pdf(paste0(prefix,".subclone.heatmap.outlier.pdf"))
plotHeatmap(dat, row_split='outlier',n_threads=50)
dev.off()

# kNN smooth profiles
dat <- knnSmooth(dat)

# Create a umap embedding 
dat <- runUmap(dat)
k_clones<-findSuggestedK(dat) #10

dat  <- findClusters(dat, k_superclones=k_clones@metadata$suggestedK-3, k_subclones=k_clones@metadata$suggestedK+3)#output from k_clones
pdf(paste0(prefix,".subclone.umap.pdf"))
plotUmap(dat, label = 'subclones')
dev.off()

pdf(paste0(prefix,".superclone.umap.pdf"))
plotUmap(dat, label = 'superclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
dat <- calcConsensus(dat)
dat <- runConsensusPhylo(dat)
dat <- runPhylo(dat, metric = 'manhattan')
dat <- calcInteger(dat,assay="segment_ratios",method="scquantum")

#col_fun = colorRamp2(c(-0.5, 0, 5), c("blue", "white", "red"))
# Plot a copy number heatmap with clustering annotation
pdf(paste0(prefix,".subclone.heatmap.segment_ratios.pdf"))
plotHeatmap(dat, label = c('reads_assigned_bins','sample_name'),assay="segment_ratios",order='hclust',n_threads=50)
dev.off()

pdf(paste0(prefix,".subclone.heatmap.integer.pdf"))
plotHeatmap(dat, label = c('reads_assigned_bins','sample_name'),assay="integer",order='hclust',n_threads=50)
dev.off()

pdf(paste0(prefix,".subclone.phylo.pdf"))
plotPhylo(dat, label = 'subclones')
dev.off()

dat <- calcInteger(dat, method = 'scquantum', assay = 'smoothed_bincounts')
pdf(paste0(prefix,".subclone.heatmap.pdf"))
plotHeatmap(dat, label = c('reads_total'),order='hclust',n_threads=50) #col=col_fun
dev.off()

saveRDS(dat,file=paste0(prefix,".scCNA.rds"))
write.table(as.data.frame(dat@colData),file=paste0(prefix,".scCNA.tsv"),sep="\t",col.names=T,row.names=T)

