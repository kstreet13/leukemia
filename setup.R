require(Seurat)

samp.names <- paste0('YMK',1:6)
dirs <- paste0('data/',samp.names,'/outs/filtered_feature_bc_matrix/')
names(dirs) <- samp.names

seu <- Read10X(dirs)
seu <- CreateSeuratObject(counts = seu)

rm(dirs, samp.names)

# %mito
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")

# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

# filtering
seu <- subset(seu, subset = nFeature_RNA > 1000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 5)

# doublets?

# run sctransform
seu <- SCTransform(seu, vars.to.regress = "orig.ident", verbose = TRUE)
#seu <- SCTransform(seu, verbose = TRUE)

seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, verbose = FALSE)

# by cluster
DimPlot(seu, label = TRUE)
# by sample
DimPlot(seu, group.by = 'orig.ident')

# clearly separate by leukemia type (different cell lines)
# split by this, remove outlier cell type, analyze separately






tallset <- c(2,5,6,7,8,9,12,13)
ballset <- c(0,1,3,4,10,11)

tallsub <- seu[, seu$seurat_clusters %in% tallset]
ballsub <- seu[, seu$seurat_clusters %in% ballset]

tall <- CreateSeuratObject(counts = tallsub@assays$RNA@layers$counts)
ball <- CreateSeuratObject(counts = ballsub@assays$RNA@layers$counts)

rm(seu, tallset, ballset, tallsub, ballsub)

# %mito
tall <- PercentageFeatureSet(tall, pattern = "^MT-", col.name = "percent.mt")
ball <- PercentageFeatureSet(ball, pattern = "^MT-", col.name = "percent.mt")

# VlnPlot(tall, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

# filtering
seu <- subset(seu, subset = nFeature_RNA > 1000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 5)


