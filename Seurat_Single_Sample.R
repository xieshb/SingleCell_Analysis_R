suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(DoubletFinder))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)

if(len < 3){
	cat("\n Rscript Seurat_Single_Sample.R 10X_matrix_dir output_dir  Sample_name \n\n")
	quit()
}
path<-getwd()
setwd(path)
input_dir<-args[1]
out_dir<-args[2]
sample<-args[3]


dir.create(out_dir,recursive=T)

data_10X<-Read10X(data.dir = input_dir,cell.column=1,gene.column=1)
colnames(data_10X)<-paste0(colnames(data_10X),"_",sample)

data_seurat<- CreateSeuratObject(counts = data_10X, project = "Seurat_P", min.cells = 3, min.features = 200)

data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")

Q1_pdf<-paste0(out_dir,"/nFeature_nCount_PercentMt.pdf")
pdf(Q1_pdf)
VlnPlot(data_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

Q2_pdf<-paste0(out_dir,"/nCount_nFeature_And_nCount_PercentMt.pdf")
pdf(Q2_pdf)
plot1 <- FeatureScatter(data_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

max_genes<-nrow(data_seurat@assays$RNA)*0.95
data<-subset(data_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < max_genes & percent.mt < 5)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
#data <- NormalizeData(data)


### DoubletFinder
#sweep.res.list <- paramSweep_v3(data, PCs = 1:pcSelect, sct = F)
#使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
#sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
#bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
#pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
#DoubletRate = ncol(data)*8*1e-6
#homotypic.prop <- modelHomotypic(data$seurat_clusters)
#nExp_poi <- round(DoubletRate*ncol(data)) 
# 使用同源双细胞比例对计算的双细胞比例进行校正 
#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 使用确定好的参数鉴定doublets
#data <- doubletFinder_v3(data, PCs = 1:pcSelect, pN = 0.25, pK = pK_bcmvn,nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

## 结果展示，分类结果在pbmc@meta.data中
#DoubletF<-paste(out_dir,"/DoubletFinder_DimPlot.pdf")
#pdf(DoubletF)
#DimPlot(data, reduction = "umap", group.by = "DF.classifications_0.25_0.3_171")
#dev.off()




data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)

Top10_hvgs<-paste0(out_dir,"/Top10_HVGS.pdf")
pdf(Top10_hvgs)
plot1 <- VariableFeaturePlot(data)+theme(axis.text.x=element_text(size=5))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()






all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
#data<-SCTransform(data)
data <- RunPCA(data, features = VariableFeatures(object = data))

PCA12<-paste0(out_dir,"/PCA12.pdf")
pdf(PCA12)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
dev.off()

PCA_15<-paste0(out_dir,"/PCA_15_DimHeatmap.pdf")
pdf(PCA_15)
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
p1<-JackStrawPlot(data, dims = 1:15)
p2<-ElbowPlot(data)
PCA_JS_EI<-paste0(out_dir,"/PCA_JackStrawPlot_ElbowPlot.pdf")
pdf(PCA_JS_EI)
p1+p2
dev.off()

data <- FindNeighbors(data, dims = 1:10)
data <- RunUMAP(data, dims = 1:10)
data <- RunTSNE(data, dims = 1:10,check_duplicates = FALSE)
data <- FindClusters(data, resolution = 0.5)

p<-DimPlot(data, reduction = "umap")
umap_pdf<-paste0(out_dir,"/DimPlot_umap.pdf")
pdf(umap_pdf)
p
dev.off()

data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers.top2 <- data.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
Fvlnplot<-paste0(out_dir,"/Markers_VlnPlot.pdf")
pdf(Fvlnplot)
VlnPlot(data, features = data.markers.top2$gene,pt.size=0)+xlab(" ")+ theme(text=element_text(size=8))
dev.off()

FeatureP<-paste0(out_dir,"/FeaturePlot.pdf")
pdf(FeatureP)
FeaturePlot(data, features = data.markers.top2$gene) + theme(text=element_text(size=8))
dev.off()

top10<- data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmapP<-paste0(out_dir,"/DoHeatmap_top10.pdf")
pdf(DoHeatmapP)
DoHeatmap(data, features = top10$gene) + NoLegend()+theme(axis.text.y=element_text(size=5))
dev.off()

data@meta.data$Sample<-sample

#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(data)
#data <- RenameIdents(data, new.cluster.ids)
#DimP<-paste(out_dir,"/DimPlot_umap_celltype.pdf")
#DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#dev.off()

RDS<-paste0(out_dir,"/",sample,"_Seurat_data.rds")
saveRDS(data, file = RDS)

















