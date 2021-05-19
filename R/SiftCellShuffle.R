#library(Matrix)
#library(reshape2)
#' The function generate shuffleDGE
#' @param workingdir the working directory
#' @import reshape2
#' @import Matrix
#' @import tictoc
#' @importFrom utils read.table write.table
#' @export
SiftCellShuffle = function(workingdir)
{
  #workingdir = "C:/Users/jyxi/Documents/rescue/112385/"
  #sanity checking
  matrixDir = paste0(workingdir,"DGE/matrix.mtx")
  barcodeDir = paste0(workingdir,"DGE/barcodes.tsv")
  geneDir = paste0(workingdir,"DGE/genes.tsv")
  barcodes = read.table(barcodeDir,sep = '\t', header = F)
  genes = read.table(geneDir,sep = '\t', header = F)
  
  #replace dot with underscore in gene names
  genes$V2 = gsub(".", "-", genes$V2, fixed = TRUE)



  m=readLines(matrixDir)
  nline = length(m)
  mline = m[3:nline]
  #spm=colsplit(string=mline[1:100], pattern=" ", names=c("Part1", "Part2","Part3"))
  spm=colsplit(string=mline, pattern=" ", names=c("Part1", "Part2","Part3"))
  n_genes = spm$Part1[1]
  n_bcds = spm$Part2[1]
  num_nonzeros = spm$Part3[1]

  if(nrow(genes) != n_genes)
  {
    stop("gene dimensiton does not match")
  }

  if(nrow(barcodes) != n_bcds)
  {
    stop("barcode dimensiton does not match")
  }


  #sum_umis = num_nonzeros
  numis = spm$Part3[-1]
  geneInd = spm$Part1[-1]
  barcodeInd = spm$Part2[-1]
  totalumi = sum(numis)
  uniq_geneInd = unique(geneInd)
  uniq_barcInd = unique(barcodeInd)

  tic()
  shuffle=dgeShuffle(length(uniq_geneInd),length(uniq_barcInd),sum(numis!=0),numis,geneInd,barcodeInd,totalumi)
  toc()  #121second good!
  shuffleDGE = sparseMatrix(
    i = shuffle$igenes,
    j = shuffle$ibcds,
    x = shuffle$umi
  )

  shuffleDGE = shuffleDGE[rowSums(shuffleDGE)>0,]
  shuffleDGE = shuffleDGE[,colSums(shuffleDGE)>0]
  subDir = "shuffleDGE"
  dir.create(file.path(workingdir, subDir), showWarnings = FALSE)
  setwd(file.path(workingdir, subDir))
  writeMM(shuffleDGE,file = "matrix.mtx")
  unig = unique(shuffle$igenes)
  unig = unig[order(unig,decreasing = F)]
  unib = unique(shuffle$ibcds)



  write.table(genes[unig,], file = "genes.tsv", row.names=FALSE,col.names = FALSE, sep="\t",quote = F)
  write.table(barcodes[unib,], file = "barcodes.tsv", row.names=FALSE,col.names = FALSE, sep="\t",quote=F)
  write.table(genes, file = geneDir, row.names=FALSE,col.names = FALSE, sep="\t",quote = F)

}



#
#
# m=readLines(matrixDir)
# nline = length(m)
#
#
#
# mline = m[3:nline]
# spm=colsplit(string=mline[1:100], pattern=" ", names=c("Part1", "Part2","Part3"))
# spm=colsplit(string=mline, pattern=" ", names=c("Part1", "Part2","Part3"))
#
#
#
# n_genes = spm$Part1[1]
# n_bcds = spm$Part2[1]
# num_nonzeros = spm$Part3[1]
# #sum_umis = num_nonzeros
# numis = spm$Part3[-1]
# geneInd = spm$Part1[-1]
# barcodeInd = spm$Part2[-1]
# totalumi = sum(numis)
# tic()
# shuffle=dgeShuffle(length(unique(geneInd)),length(unique(barcodeInd)),sum(numis!=0),numis,geneInd,barcodeInd,totalumi)
# toc()  #121second good!
# shuffleDGE = sparseMatrix(
#   i = shuffle$igenes,
#   j = shuffle$ibcds,
#   x = shuffle$umi
# )
# writeMM(shuffleDGE,file = "C:/Users/jyxi/Downloads/demo.mtx")
