.rs.restartR() 

library(Seurat)
library(SeuratDisk)
library(CellChat)
library(Matrix)
library(dplyr)
library(future)
library(patchwork)


data_dir <- "~/ava/data/GSE118127"  
files <- list.files(data_dir, pattern = "_bc_matrices_h5\\.h5$", full.names = TRUE)
stopifnot(length(files) > 0)


objs <- lapply(files, function(f) {
  mat <- Read10X_h5(f)  # 10x HDF5
  sample <- sub("_bc_matrices_h5\\.h5$", "", basename(f))
  so <- CreateSeuratObject(counts = mat, project = sample, min.cells = 3, min.features = 200)
  so$sample <- sample
  so <- RenameCells(so, add.cell.id = sample)  # cell barcode 前缀化：sample_barcode
  so
})


seu <- Reduce(function(x, y) merge(x, y), objs)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

seu <- RunUMAP(seu, dims = 1:30)
UMAPPlot(seu)
Idents(seu) <- seu$orig.ident
UMAPPlot(seu)
# R: assign_batch
assign_batch <- function(sample_name) {
  s <- as.character(sample_name)
  res <- rep("unKnown", length(s))
  
  res[grepl("sample_[13]-\\d+($|_)", s)] <- "Batch1"
  
  res[res == "unKnown" & grepl("(B1|B2|C1|145)", s)] <- "Batch2"
  res
}

seu$Batch <- assign_batch(seu$sample)
Idents(seu) <- seu$Batch
UMAPPlot(seu)
library(harmony)
seu <- RunHarmony(seu,group.by.vars = 'orig.ident')
seu <- RunUMAP(seu,dims = 1:20,reduction = 'harmony')
DimPlot(seu,reduction = 'umap')+xlab('UMAP 1')+ylab('UMAP 2')+ theme_bw()
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
seu <- FindClusters(seu, resolution = 1)  # 


seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)
seu <- RunTSNE(seu, reduction = "harmony", dims = 1:30)
DimPlot(seu,reduction = 'tsne',label = T,repel = T)

VlnPlot(seu,features = c('CD53','CXCR4')) # cluster 10,15 are Immune
VlnPlot(seu,features = c('HSD17B1','SERPINE2','GSTA1'))  # Granulosa 17
VlnPlot(seu,features = c('IGFBP2','HTRA1','INHBB'))  # cumulus 6
VlnPlot(seu,features = c('KRT18','CITED2','LIPH'))  # mural 13
FeaturePlot(seu,c('KRT18','CITED2','LIPH','AKIRIN1'),reduction = 'tsne')

Idents(seu) <- seu$seurat_clusters
VlnPlot(seu,features = c('GSTA1'))  # Granulosa 17
DimPlot(seu,reduction = 'tsne',label = T,repel = T)
colnames(seu)
VlnPlot(seu,features = c('TAGLN','RGS5')) # Smooth muscle 7,14
VlnPlot(seu,features = c('VWF','CLDN5')) # Endothelial 2,5,18
VlnPlot(seu,features = c('DCN','LUM')) # Theca & Stroma 0,1,3,4,8,9,11,12,16,
seu <- JoinLayers(seu, assay = "RNA", layers = "data")
DE = FindMarkers(seu,ident.1 = '15')
# cluster 15 is macrophage

Idents(seu) <- seu$seurat_clusters
u <- as.character(seu$seurat_clusters)

ct_map <- c(
  "0"  = "Theca_Stroma",
  "1"  = "Theca_Stroma",
  "2"  = "Endothelial",
  "3"  = "Theca_Stroma",
  "4"  = "Theca_Stroma",
  "5"  = "Endothelial",
  "6"  = "cumulus_GC",     
  "7"  = "Smooth_muscle",
  "8"  = "Theca_Stroma",
  "9"  = "Theca_Stroma",
  "10" = "Immune",           
  "11" = "Theca_Stroma",
  "12" = "Theca_Stroma",
  "13" = "Mural_GC",         
  "14" = "Smooth_muscle",
  "15" = "Macrophage",       
  "16" = "Theca_Stroma",
  "17" = "Granulosa",
  "18" = "Endothelial"
)


anno <- unname(ct_map[as.character(seu$seurat_clusters)])
anno[is.na(anno)] <- "Unknown"
seu$celltype <- anno
Idents(seu) <- "celltype"

DimPlot(seu, reduction = "tsne", label = TRUE, repel = TRUE, group.by = "celltype")

anno <- unname(ct_map[u])         #
anno[is.na(anno)] <- "Unknown"    # 

seu$celltype <- anno
Idents(seu) <- "celltype"
DimPlot(seu, reduction = "tsne", label = TRUE, repel = TRUE)
Idents(seu) <- seu$seurat_clusters
DimPlot(seu, reduction = "tsne", label = TRUE, repel = TRUE)
FeaturePlot(seu,'DCN')
Idents(seu) <- "celltype"


print(table(seu$celltype))
# DimPlot(seu, reduction = "umap", label = TRUE, group.by = "celltype")
# cellchat analysis
seu <- JoinLayers(seu, assay = "RNA", layers = "data")
data.input <- GetAssayData(seu,assay = 'RNA',layer = 'data')
labels = seu$celltype
cells <- colnames(data.input)
meta <- data.frame(
  group = as.character(Idents(seu)[cells]),  
  row.names = cells,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- setIdent(cellchat,ident.use = 'group')
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat,raw.use = T)
cellchat <- filterCommunication(cellchat,min.cells = 2)




df.net1 <- subsetCommunication(cellchat)



# saveRDS(seu,file = 'seu.rds')
seu <- readRDS('seu.rds')
samples <- unique(seu$sample)
cellchat_list <- list()
for (s in samples) {
  seu_sub <- subset(seu, subset = sample == s)
  
  data.input <- GetAssayData(seu_sub, assay = "RNA", slot = "data")
  meta <- data.frame(group = seu_sub$celltype, row.names = colnames(seu_sub))
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat_list[[s]] <- cellchat
}

keep_groups <- c("Granulosa","cumulus_GC","Mural_GC")
mode <- "touching"  
extract_lr_from_existing <- function(cc, keep_groups, mode = c("within","touching")) {
  mode <- match.arg(mode)
  
  comm <- subsetCommunication(cc)
  
  if (mode == "within") {
    comm <- comm[comm$source %in% keep_groups & comm$target %in% keep_groups, ]
  } else { # touching
    comm <- comm[comm$source %in% keep_groups | comm$target %in% keep_groups, ]
  }
  
  lr_cols <- intersect(c("ligand","receptor","interaction_name","pathway_name","source","target"), colnames(comm))
  unique_lr <- unique(comm[, lr_cols, drop = FALSE])
  rownames(unique_lr) <- NULL
  return(unique_lr)
}

lr_list <- lapply(cellchat_list, extract_lr_from_existing, keep_groups = keep_groups, mode = mode)


unique_lr_all <- unique(do.call(rbind, lr_list))

target_list_human <- c(
  "BMP2","BMP4","BMP6","BMP7","GDF11","AMH","INHBB",
  "WNT10A","WNT10B","WNT2","WNT4","WNT6","WNT9A","WNT5A",
  "TGFA","AREG","BTC","HBEGF","EREG",
  "FGF1","FGF2","FGF18","FGF9","FGF16","FGF21",
  "IGF1","DHH","IHH","CXCL12","CTF1","LIF","OSM",
  "TNF","TNFSF12","SPP1","RETN","NAMPT",
  "ANGPTL2","ANGPTL4","ANGPT2","MDK",
  "EDN1","EDN2","NPPC","GAS6","PROS1"
)

df.target <- unique_lr_all[unique_lr_all$ligand %in% target_list_human |unique_lr_all$receptor %in% target_list_human, ]

print(unique(df.target$ligand))
print(unique(df.target$receptor))

write.csv(unique_lr_all, "~/ava/data/GSE118127/unique_LR_touching.csv", row.names = FALSE)
saveRDS(cellchat_list,file = "~/ava/data/GSE118127/cellchat_list.rds")
