suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)

if(len < 2 ){
	cat("\nRscript Seurat_multi_Sample.R A.rds,B.rds,C.rds,D.rds  output_dir\n\n")
	quit()
}

outdir<-args[2]
path<-getwd()
setwd(path)

dir.create(outdir,recursive=T)

RDSLIST<-strsplit(args[1],',')
print(RDSLIST[[1]])
len<-length(RDSLIST[[1]])
AA2<-c()
tem<-strsplit(RDSLIST[[1]],"/")
for(i in 1:len){AA2<-c(AA2,tem[[i]][length(tem[[i]])])}
RDSNAME<-gsub(".rds","",AA2)
print(AA2)
#for(i in 1:len){assign(RDSNAME[i],value=readRDS(AA2[i]))}
for(i in 1:len){assign(RDSNAME[i],value=readRDS(RDSLIST[[1]][i]))}
Slist<-c()
for(i in 1:len){Slist<-c(Slist,get(RDSNAME[i]))}
print(Slist)

features <- SelectIntegrationFeatures(object.list = Slist)
anchors <- FindIntegrationAnchors(object.list = Slist, anchor.features = features)
#combined <- IntegrateData(anchorset = anchors,normalization.method = c("SCT"))
combined <- IntegrateData(anchorset = anchors,normalization.method = c("LogNormalize"))

DefaultAssay(combined)<-"integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 20, verbose = FALSE)
combined <- RunUMAP(combined,reduction="pca",dims=1:20)
combined <- RunTSNE(combined,reduction="pca",dims=1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined)

DimPlot_p<-paste(outdir,"/DimPlot_combined.pdf")
pdf(DimPlot_p)
p1 <- DimPlot(combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()
combine.markers<-FindAllMarkers(combined)
markers<-combine.markers[c("gene","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster")]
marker<-paste0(outdir,"/combined_markers.xls")
write.table(markers,file=marker,quote=F,sep="\t",row.names=F)

outRDS<-paste0(outdir,"/combined.rds")
saveRDS(combined,file=outRDS)


