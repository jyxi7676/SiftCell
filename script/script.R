######Scripts for SiftCell paper#######
#(A): Tutorials for running SiftCell framework can be found at colab notebook: https://colab.research.google.com/drive/1ebV42lohhkTWGjkHHjUs1Ma92IXAxx2t
#(B): CellBender is run on colab notebook using GPU, the scripts can be found at : https://colab.research.google.com/drive/1wh6WTUbVNsiwW-V4F7xAtW1BL_qtsm9U https://colab.research.google.com/drive/1fJW6Cygt9QsT8ULSyep9k5JUIq3pXA3v?usp=sharing https://colab.research.google.com/drive/1x33VlEaZHmcr4JeenRJD4b6tKr-3sS5i?usp=sharing
#(C): The following code is for DIEM, Emptydrops, DecontX and SiftCell framwork.


#install packages:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DropletUtils")
library(DropletUtils)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("celda")
library(celda)

devtools::install_github("tpq/peakRAM")
library(peakRAM)

library(devtools)
devtools::install_github("marcalva/diem")
library(diem)

install.packages('Seurat')
library(Seurat)
library(Matrix)



emptyD=function(workingdir)
{
  sce <- read10xCounts(paste0(workingdir,'/DGE'))
  set.seed(100)
  e.out <- emptyDrops(sce)
  return(e.out)
}

diem=function(workingdir)
{

  counts <- read_10x(paste0(workingdir,'/DGE')) # Read 10X data into sparse matrix
  sce <- create_SCE(counts) # Create SCE object from counts

  #Add MT% and MALAT1% (optional, xx seconds)
  mt_genes <- grep(pattern="^mt-", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
  sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
  malat <- grep(pattern="^malat1$", x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)
  sce <- get_gene_pct(x = sce, genes=malat, name="MALAT1")

  # DIEM steps
  sce <- set_debris_test_set(sce)
  sce <- filter_genes(sce)
  sce <- get_pcs(sce)
  sce <- init(sce)
  sce <- run_em(sce)
  sce <- assign_clusters(sce)
  sce <- estimate_dbr_score(sce)

  # Evaluate debris scores
  # sm <- summarize_clusters(sce)
  # plot_clust(sce, feat_x = "n_genes", feat_y = "score.debris",
  #            log_x = TRUE, log_y = FALSE)
  # plot_clust(sce, feat_x = "pct.mt", feat_y = "score.debris",
  #            log_x = TRUE, log_y = FALSE)

  # Call targets using debris score for single-nucleus data
  sce <- call_targets(sce, thresh_score = 0.5)

  # Call targets by removing droplets in debris cluster(s) for single-cell data
  sce <- call_targets(sce, clusters = "debris", thresh = NULL)

  seur <- convert_to_seurat(sce)
  return(seur)
}

decontX=function(workingdir,UMIthreshold)
{
  m=readMM(paste0(workingdir,'/DGE/matrix.mtx'))
  c=read.table(paste0(workingdir,'/DGE/barcodes.tsv'))$V1
  r=read.table(paste0(workingdir,'/DGE/genes'))$V2
  colnames(m)=c
  rownames(m)=r
  m=m[,colSums(m)>UMIthreshold]
  m=m[rowSums(m)>0,]
  orgDGE=as.matrix(m)
  decont=decontX(x=orgDGE)
  return(decont)
}


#brain data:
workingdir='~/brain_nuclei_e18_1k/'
celltype=read.csv('~/Brain_SiftCellMix_celltype_XG_res0.05.csv',row.names=1)

peakRAM(o1=emptyD(workingdir))
peakRAM(o2=diem(workingdir))
peakRAM(SiftCellShuffle(workingdir))
peakRAM(SiftCellBoost(workingdir,threshold=100,expectedN=1000,dataName='Brain',mitoGenes = FALSE))
peakRAM(SiftCellMix(workingdir,celltype))
peakRAM(o3=decontX(workingdir,100))



#Colon data:
workingdir='~/colon/'
celltype=read.csv('~/Colon_SiftCellMix_celltype_XG_res0.1.csv',row.names=1)
emptyD=function(workingdir)
peakRAM(o1=emptyD(workingdir))
peakRAM(o2=diem(workingdir))
peakRAM(SiftCellShuffle(workingdir))
peakRAM(SiftCellBoost(workingdir,threshold=200,expectedN=800,dataName='Colon',mitoGenes = FALSE))
peakRAM(SiftCellMix(workingdir,celltype))
peakRAM(o3=decontX(workingdir,200))


workingdir='~/pbmc10k_v3/'
celltype=read.csv('~/PBMC_SiftCellMix_celltype_XG_res0.05.csv',row.names = 1)
peakRAM(o1=emptyD(workingdir))
peakRAM(o2=diem(workingdir))
peakRAM(SiftCellShuffle(workingdir))
peakRAM(SiftCellBoost(workingdir,threshold=100,expectedN=10000,dataName='PBMC',mitoGenes = FALSE))
peakRAM(SiftCellMix(workingdir,celltype))
peakRAM(o3=decontX(workingdir,100))


