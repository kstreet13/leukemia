require(Seurat)

tall <- readRDS('data/tall.rds')

d1 <- DimPlot(tall, label = TRUE)
d2 <- DimPlot(tall, group.by = 'orig.ident')
d1 + d2


# dots (use actual known genes)

ballmarkers <- unique(c('CTNNB1','CD10','CD34','CD38','CD19','CD133','CD19','CD34','CD20','CD24','TCF7','NOTCH1','LEF1','RUNX1','PROM1','CKIT','LGR5','HOX','LMO2','TCF','TCR','BCL11B','GATA3','E2A','EBF1','RAG','PAX5','EZH2','DOT1L','MLL','SUV420H1','SUV420H2','DNMT3A','TET2','BCL6','IL7R','FOXP3','PTCRA',
                        'TBL1X','TBLR1','MYC','BIRC5','CCND1','AXIN2','TCF','CBP','P300','DKK1','BCL2','P53','SIAH1','CDK','JUN','CDH1'))
tallmarkers <- unique(c('CTNNB1','CD7','CD1A','CD34','CD38','CD4','CD44','CD8','CD3','TCF7','NOTCH1','LEF1','RUNX1','PROM1','CKIT','LGR5','HOX','LMO2','TCF','TCR','BCL11B','GATA3','E2A','EBF1','RAG','PAX5','EZH2','DOT1L','MLL','SUV420H1','SUV420H2','DNMT3A','TET2','BCL6','IL7R','FOXP3','PTCRA','FZD','LRP5','LRP6','DVL','AXIN','APC','GSK3B','CK1','TBL1X','TBLR1','MYC','BIRC5','CCND1','AXIN2','TCF','LEF1','CBP','P300','EP300','DKK1','BCL2','P53','SIAH1','CDK','JUN','CDH1'))

othermarkers <- unique(c('NFKB','P65','RELA','SMRT','NCOR','HDAC3','GPS2','SCF','CUL1','SKP1','RING','RBX1','ROC1','PI3K','AKT','PTEN','MTOR','MTORC1','MAPK','ERK','JNK','P38','ERK5','TNF','IL','IFN','NOTCH','JAKSTAT','VEGF','PARP1','PPP1R15A','RAD51','BRCA1','BRCA2','BAX','BIM','PUMA','MCL1','ABC','ABCB1','MDR1','ABCC1','ABCG2'))

allmarkers <- unique(c(ballmarkers,tallmarkers,othermarkers)) 


DotPlot(tall, allmarkers, cluster.idents = TRUE) + RotatedAxis()


interestingGenes <- c('CTNNB1','CD7','CD38','CD44','TCF7','LEF1','RUNX1','BCL11B','BIRC5','LMO2','EZH2','TET2','BCL6','TBL1X','JUN','GSK3B','SKP1','RBX1','PTEN','PARP1','PPP1R15A','BRCA1','ABCC1')

DotPlot(tall, interestingGenes, cluster.idents = TRUE) + RotatedAxis()


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
de2 <- FindMarkers(tall, ident.1 = 'YMK3', ident.2 = 'YMK1', group.by = 'orig.ident')

# save tables
write.csv(de1, file='~/Desktop/tall_allCells_6h-ctrl.csv', quote = FALSE)
write.csv(de2, file='~/Desktop/tall_allCells_24h-ctrl.csv', quote = FALSE)

# volcano plots with genes labelled
png(filename = '~/Desktop/T-ALL/tall_allCells_6h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de1$p_val_adj < .99)
plot(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'T-ALL All Cells: 6h - CTRL')
text(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), labels = rownames(de1)[filt])
dev.off()

png(filename = '~/Desktop/T-ALL/tall_allCells_24h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de2$p_val_adj < .99)
plot(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'T-ALL All Cells: 24h - CTRL')
text(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), labels = rownames(de2)[filt])
dev.off()


# DE between samples within cluster 11 (t-cell leukemia stem cells)
# DE: JUST CLUSTER 11 #
#######################

de1 <- FindMarkers(tall, ident.1 = 'YMK2', ident.2 = 'YMK1', group.by = 'orig.ident', subset.ident = '11')
de2 <- FindMarkers(tall, ident.1 = 'YMK3', ident.2 = 'YMK1', group.by = 'orig.ident', subset.ident = '11')

# save tables
write.csv(de1, file='~/Desktop/tall_clus11_6h-ctrl.csv', quote = FALSE)
write.csv(de2, file='~/Desktop/tall_cluss11_24h-ctrl.csv', quote = FALSE)

# volcano plots with genes labelled
png(filename = '~/Desktop/T-ALL/tall_clus11_6h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de1$p_val_adj < .99)
plot(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'T-ALL Cluster 11: 6h - CTRL')
text(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), labels = rownames(de1)[filt])
dev.off()

png(filename = '~/Desktop/T-ALL/tall_clus11_24h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de2$p_val_adj < .99)
plot(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'T-ALL Cluster 11: 24h - CTRL')
text(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), labels = rownames(de2)[filt])
dev.off()


# DE cluster 11 markers
# DE: CLUSTER 11 VS. OTHERS #
#############################

de1 <- FindMarkers(tall, ident.1 = '11')

# save tables
write.csv(de1, file='~/Desktop/T-ALL/tall_clus11-other.csv', quote = FALSE)

# volcano plots with genes labelled
png(filename = '~/Desktop/T-ALL/tall_clus11-other_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de1$p_val_adj < .99)
plot(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'T-ALL All Cells: Cluster 11 - Other')
text(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), labels = rownames(de1)[filt])
dev.off()




