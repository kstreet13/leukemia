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

# VlnPlot(tall, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

# filtering
seu <- subset(seu, subset = nFeature_RNA > 1000 & nCount_RNA > 2000 & nCount_RNA < 30000 & percent.mt < 5)







