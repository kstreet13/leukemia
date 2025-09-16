require(Seurat)

ball <- readRDS('data/ball.rds')

d1 <- DimPlot(ball, label = TRUE)
d2 <- DimPlot(ball, group.by = 'orig.ident')
d1 + d2


# dots

markers <- unique(c('CTNNB1','CD10','CD34','CD38','CD19','CD133','CD19','CD34','CD20','CD24','TCF7','NOTCH1','LEF1','RUNX1','PROM1','CKIT','LGR5','HOX','LMO2','TCF','TCR','BCL11B','GATA3','E2A','EBF1','RAG','PAX5','EZH2','DOT1L','MLL','SUV420H1','SUV420H2','DNMT3A','TET2','BCL6','IL7R','FOXP3','PTCRA',
                    'TBL1X','TBLR1','MYC','BIRC5','CCND1','AXIN2','TCF','CBP','P300','DKK1','BCL2','P53','SIAH1','CDK','JUN','CDH1'))



DotPlot(ball, markers) + RotatedAxis()




# clusters 1,4,8,7 are interesting

ball$myclus <- rep(0, ncol(ball))
ball$myclus[ball$seurat_clusters == 1] <- 1
ball$myclus[ball$seurat_clusters == 4] <- 2
ball$myclus[ball$seurat_clusters == 8] <- 3
ball$myclus[ball$seurat_clusters == 7] <- 4

DimPlot(ball, group.by = 'myclus')

mosaicplot(table(ball$myclus, ball$orig.ident), main='Clusters by Sample', col = 2:4)
mosaicplot(table(ball$seurat_clusters, ball$orig.ident), main='Clusters by Sample', col = 2:4)
# col=c('salmon','goldenrod','lightgreen','dodgerblue','purple')

# Differential Expression
markers <- FindAllMarkers(ball, group.by = 'myclus')
# markers1 <- FindMarkers(ball, ident.1 = 1, ident.2 = 0, group.by = 'myclus')
# markers2 <- FindMarkers(ball, ident.1 = 2, ident.2 = 0, group.by = 'myclus')
# markers3 <- FindMarkers(ball, ident.1 = 3, ident.2 = 0, group.by = 'myclus')
# markers4 <- FindMarkers(ball, ident.1 = 4, ident.2 = 0, group.by = 'myclus')


# violin
VlnPlot(ball, rownames(markers)[1:12], group.by = 'myclus', slot = 'data', log=FALSE, alpha=0.1)


# DE BETWEEN SAMPLES #
######################

de1 <- FindMarkers(ball, ident.1 = 'YMK5', ident.2 = 'YMK4', group.by = 'orig.ident')
# volcano
plot(de1$avg_log2FC, -log10(de1$p_val_adj), col = rgb(0,0,0,.5))

de2 <- FindMarkers(ball, ident.1 = 'YMK6', ident.2 = 'YMK4', group.by = 'orig.ident')
# volcano
plot(de2$avg_log2FC, -log10(de2$p_val_adj), col = rgb(0,0,0,.5))


write.csv(de1, file='~/Desktop/ball_6h-ctrl.csv', quote = FALSE)
write.csv(de2, file='~/Desktop/ball_24h-ctrl.csv', quote = FALSE)

# add more markers to Dot Plot

# DE between samples within groups (4&8, 1&7)

