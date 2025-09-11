require(Seurat)

samp.names <- paste0('YMK',1:6)
dirs <- paste0('data/',samp.names,'/outs/filtered_feature_bc_matrix/')
names(dirs) <- samp.names

tall <- Read10X(dirs[1:3])
tall <- CreateSeuratObject(counts = tall)
ball <- Read10X(dirs[4:6])
ball <- CreateSeuratObject(counts = ball)

rm(dirs, samp.names)

# %mito
tall <- PercentageFeatureSet(tall, pattern = "^MT-", col.name = "percent.mt")
ball <- PercentageFeatureSet(ball, pattern = "^MT-", col.name = "percent.mt")

# tall$lognCount <- log(tall$nCount_RNA)
# ball$lognCount <- log(ball$nCount_RNA)
# 
# VlnPlot(ball, features = c("nFeature_RNA", "lognCount", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(ball, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(ball, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# 
# boxplot(ball$lognCount ~ ball$orig.ident)
# boxplot(ball$nFeature_RNA ~ ball$orig.ident)

# filtering
tall <- subset(tall, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 5)
ball <- subset(ball, subset = nFeature_RNA > 750 & nFeature_RNA < 6000 & nCount_RNA > 1800 & nCount_RNA < 22000 & percent.mt < 5)


# REMOVE CONTAMINANT CELL TYPE #
################################

# run sctransform
tall <- SCTransform(tall, vars.to.regress = "orig.ident", verbose = TRUE)
ball <- SCTransform(ball, vars.to.regress = "orig.ident", verbose = TRUE)
#tall <- SCTransform(tall, verbose = TRUE)

tall <- RunPCA(tall, verbose = FALSE)
ball <- RunPCA(ball, verbose = FALSE)
tall <- RunUMAP(tall, dims = 1:30, verbose = FALSE)
ball <- RunUMAP(ball, dims = 1:30, verbose = FALSE)

tall <- FindNeighbors(tall, dims = 1:30, verbose = FALSE)
ball <- FindNeighbors(ball, dims = 1:30, verbose = FALSE)
tall <- FindClusters(tall, verbose = FALSE)
ball <- FindClusters(ball, verbose = FALSE)

# by cluster
DimPlot(tall, label = TRUE)
# by sample
DimPlot(tall, group.by = 'orig.ident')

# by cluster
DimPlot(ball, label = TRUE)
# by sample
DimPlot(ball, group.by = 'orig.ident')


# tall - remove if umap_1 > 10
tall <- tall[, which(tall@reductions$umap@cell.embeddings[,1] < 10)]

# ball - remove if umap_2 > 7
ball <- ball[, which(ball@reductions$umap@cell.embeddings[,2] < 7)]


# RE-RUN NORMALIZATION #
########################

# run sctransform
tall <- SCTransform(tall, vars.to.regress = "orig.ident", verbose = TRUE)
ball <- SCTransform(ball, vars.to.regress = "orig.ident", verbose = TRUE)
#tall <- SCTransform(tall, verbose = TRUE)

tall <- RunPCA(tall, verbose = FALSE)
ball <- RunPCA(ball, verbose = FALSE)
tall <- RunUMAP(tall, dims = 1:30, verbose = FALSE)
ball <- RunUMAP(ball, dims = 1:30, verbose = FALSE)

tall <- FindNeighbors(tall, dims = 1:30, verbose = FALSE)
ball <- FindNeighbors(ball, dims = 1:30, verbose = FALSE)
tall <- FindClusters(tall, verbose = FALSE)
ball <- FindClusters(ball, verbose = FALSE)

# by cluster
DimPlot(tall, label = TRUE)
# by sample
DimPlot(tall, group.by = 'orig.ident')

# by cluster
DimPlot(ball, label = TRUE)
# by sample
DimPlot(ball, group.by = 'orig.ident')


saveRDS(tall, file='data/tall.rds')
saveRDS(ball, file='data/ball.rds')



