 install.packages("devtools")
 library(devtools)
 install_github("jyxi7676/SiftCell",force = T)
 library(SiftCell)
#colon data as example
system("gdown 'https://drive.google.com/uc?id=1QMWvkVFUrVvT9CzJ7_Tb-WaHvXifn74G'")
system("tar xzvf DGE.tar.gz")



#SiftCellShuffle
workingdir = "/content/"
setwd(workingdir)
SiftCellShuffle(workingdir)


#Visualization
install.packages("Seurat")
library(Seurat)
library(ggplot2)
library(Matrix)
dir1 = "/content/DGE/"
dir2 = "/content/shuffleDGE/"
dataName = "Sample"
threshold = 200
df.raw = data.frame(BATCH=character(),N.CELL=integer(),BATCH=character(),N.UMI=integer())
genes = read.table(file = paste0(dir1,"/genes.tsv"), sep = "\t")
genes$V2 = gsub('[_ ]', '-', genes$V2)
write.table(genes, file = paste0(dir1,"/genes.tsv"), quote = FALSE, sep = '\t', row.names = F, col.names = F)
obj.org = Read10X(dir1)
obj.org = obj.org[rowSums(obj.org)>0,]
obj.org = obj.org[,colSums(obj.org)>0]
n.umis = colSums(obj.org)
batch1 = paste0(dataName,"Org")
df.raw = rbind(df.raw,data.frame(BATCH=rep(batch1,length(n.umis)),N.CELL=1:length(n.umis),N.UMI=sort(n.umis,decreasing=TRUE)))
obj.org = obj.org[,colSums(obj.org)>threshold]
obj.Org = CreateSeuratObject(obj.org, min.cells =0, min.features = 0, project = batch1)

genes = read.table(file = paste0(dir2,"/genes.tsv"), sep = "\t")
genes$V2 = gsub('[_ ]', '-', genes$V2)
write.table(genes, file = paste0(dir2,"/genes.tsv"), quote = FALSE, sep = '\t', row.names = F, col.names = F)
obj.shf = Read10X(dir2)
obj.shf = obj.shf[rowSums(obj.shf)>0,]
obj.shf = obj.shf[,colSums(obj.shf)>0]
n.umis = colSums(obj.shf)
batch2 = paste0(dataName,"Shf")
df.raw = rbind(df.raw,data.frame(BATCH=rep(batch2,length(n.umis)),N.CELL=1:length(n.umis),N.UMI=sort(n.umis,decreasing=TRUE)))
obj.shf = obj.shf[,colSums(obj.shf)>threshold]
obj.Shf = CreateSeuratObject(obj.shf, min.cells = 0, min.features = 0, project = batch2)

n.org = dim(obj.Org)[2]
n.shf = dim(obj.Shf)[2]
obj.srt = merge(obj.Org, y=c(obj.Shf),add.cell.ids = c(batch1,batch2), project = 'merged')
obj.srt = subset(x=obj.srt,subset= nCount_RNA > threshold)
obj.srt = NormalizeData(object = obj.srt)
obj.srt = FindVariableFeatures(object = obj.srt, selection.method = 'vst')
obj.srt[["percent.mt"]] = PercentageFeatureSet(obj.srt, pattern = "^MT-")
obj.srt = ScaleData(object = obj.srt)
obj.srt = RunPCA(object = obj.srt, npcs = 30)  #depends on data: if too few above UMI cutoff. 124792
obj.srt = FindNeighbors(object = obj.srt, dims = 1:10)
obj.srt = FindClusters(object = obj.srt, resolution = 0.1)
obj.srt = RunTSNE(object = obj.srt, dims = 1:10,check_duplicates = F)
obj.srt = RunUMAP(object = obj.srt, dims = 1:10)
df.cell = cbind(obj.srt@meta.data, obj.srt@reductions$tsne@cell.embeddings, obj.srt@reductions$umap@cell.embeddings, obj.srt@reductions$pca@cell.embeddings)

install.packages("cowplot")
library(cowplot)
library(ggplot2)
p1 <- ggplot(df.cell, aes(tSNE_1,tSNE_2,colour=orig.ident)) + geom_point(size=0.5,alpha=0.3) + theme(legend.position='bottom')
p2 <- ggplot(df.cell, aes(tSNE_1,tSNE_2,colour=as.factor(RNA_snn_res.0.1))) + geom_point(size=0.5,alpha=0.3) + theme(legend.position='bottom')
p3 <- ggplot(df.cell, aes(tSNE_1,tSNE_2,colour=nCount_RNA)) + geom_point(size=0.5,alpha=0.3) +  theme(legend.position='bottom',legend.key.width=unit(0.9,'in')) + scale_color_gradientn(colours=rainbow(7),trans='log10')
p4 <- ggplot(df.cell, aes(UMAP_1,UMAP_2,colour=orig.ident)) + geom_point(size=0.5,alpha=0.3) + theme(legend.position='bottom')
p5 <- ggplot(df.cell, aes(UMAP_1,UMAP_2,colour=as.factor(RNA_snn_res.0.1))) + geom_point(size=0.5,alpha=0.3) + theme(legend.position='bottom')
p6 <- ggplot(df.cell, aes(UMAP_1,UMAP_2,colour=nCount_RNA)) + geom_point(size=0.5,alpha=0.3) + theme(legend.position='bottom',legend.key.width=unit(0.9,'in')) + scale_color_gradientn(colours=rainbow(7),trans='log10')
p=plot_grid(p1,p2,p3,p4,p5,p6,nrow=2,ncol=3)
plot(p)


#SiftCellBoost and visualization
SiftCellBoost("/content/",threshold = 200)  #
cell = read.table("/content/CellContaining.txt",header = F)[,1]
df.cell.org = df.cell[df.cell$orig.ident=='SampleOrg',c("nCount_RNA","tSNE_1","tSNE_2","UMAP_1","UMAP_2")]
df.cell.org$SiftCellBoost = "Soup"
df.cell.org[match(paste0("SampleOrg_",cell),rownames(df.cell.org)),"SiftCellBoost"] = "Cell"
ggplot(df.cell.org,aes(tSNE_1,tSNE_2,color = SiftCellBoost))+geom_point(size=0.5,alpha=0.3)+theme_bw()+theme(legend.position = "bottom")


#SiftCellMix:
system("gdown 'https://drive.google.com/uc?id=14veqr3kA6GKoGKrG3EWcheE37-b1MNJH'")
celltype = read.csv("/content/external_celltype.csv",row.names=1,header = T)
workingdir="/content/"
setwd(workingdir)
SiftCellMix(workingdir,celltype,alpha=0.01,threshold=200)
frac=read.csv('/content/frac.csv',row.names=1)
frac$BARCODE = paste0("SampleOrg_",rownames(frac))
df.cell.org = df.cell[df.cell$orig.ident=='SampleOrg',c("nCount_RNA","tSNE_1","tSNE_2","UMAP_1","UMAP_2")]
df.cell.org$BARCODE = rownames(df.cell.org)
frac_df = merge(df.cell.org,frac,by="BARCODE")
head(frac_df)
ggplot(frac_df,aes(tSNE_1,tSNE_2,color = 1-p.AMB))+geom_point(size=0.5,alpha=0.3)+theme_bw()+scale_color_gradientn(colors = rainbow(7))+theme(legend.position ="bottom")