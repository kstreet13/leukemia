require(Seurat)

tall <- readRDS('data/tall.rds')

d1 <- DimPlot(tall, label = TRUE)
d2 <- DimPlot(tall, group.by = 'orig.ident')
d1 + d2


# dots (use actual known genes)

markers <- unique(c('CTNNB1','CD7','CD1A','CD34','CD38','CD4','CD44','CD8','CD3','TCF7','NOTCH1','LEF1','RUNX1','PROM1','CKIT','LGR5','HOX','LMO2','TCF','TCR','BCL11B','GATA3','E2A','EBF1','RAG','PAX5','EZH2','DOT1L','MLL','SUV420H1','SUV420H2','DNMT3A','TET2','BCL6','IL7R','FOXP3','PTCRA'))

DotPlot(tall, markers) + RotatedAxis()



# clusters 11 and maybe 14 are interesting

tall$myclus <- rep(0, ncol(tall))
tall$myclus[tall$seurat_clusters == 14] <- 1
tall$myclus[tall$seurat_clusters == 11] <- 2

DimPlot(tall, group.by = 'myclus')

mosaicplot(table(tall$myclus, tall$orig.ident), main='Clusters by Sample', col = 2:4)
mosaicplot(table(tall$seurat_clusters, tall$orig.ident), main='Clusters by Sample', col = 2:4)
# col=c('salmon',3,'dodgerblue')

markers2 <- FindMarkers(tall, ident.1 = 2, ident.2 = 0, group.by = 'myclus')
# volcano
plot(markers2$avg_log2FC, -log10(markers2$p_val_adj), col = rgb(0,0,0,.5))

markers1 <- FindMarkers(tall, ident.1 = 1, ident.2 = 0, group.by = 'myclus')
# volcano
plot(markers1$avg_log2FC, -log10(markers1$p_val_adj), col = rgb(0,0,0,.5))


# violin
VlnPlot(tall, markers[1:10], slot = 'data', log=FALSE, alpha=0.1)


# DE BETWEEN SAMPLES #
######################

de1 <- FindMarkers(tall, ident.1 = 'YMK2', ident.2 = 'YMK1', group.by = 'orig.ident')
# volcano
plot(de1$avg_log2FC, -log10(de1$p_val_adj), col = rgb(0,0,0,.5))

de2 <- FindMarkers(tall, ident.1 = 'YMK3', ident.2 = 'YMK1', group.by = 'orig.ident')
# volcano
plot(de2$avg_log2FC, -log10(de2$p_val_adj), col = rgb(0,0,0,.5))


write.csv(de1, file='~/Desktop/tall_6h-ctrl.csv', quote = FALSE)
write.csv(de2, file='~/Desktop/tall_24h-ctrl.csv', quote = FALSE)



# DE between samples within cluster 11 (t-cell leukemia stem cells)


