require(Seurat)

ball <- readRDS('data/ball.rds')

d1 <- DimPlot(ball, label = TRUE)
d2 <- DimPlot(ball, group.by = 'orig.ident')
d1 + d2


# dots

ballmarkers <- unique(c('CTNNB1','CD10','CD34','CD38','CD19','CD133','CD19','CD34','CD20','CD24','TCF7','NOTCH1','LEF1','RUNX1','PROM1','CKIT','LGR5','HOX','LMO2','TCF','TCR','BCL11B','GATA3','E2A','EBF1','RAG','PAX5','EZH2','DOT1L','MLL','SUV420H1','SUV420H2','DNMT3A','TET2','BCL6','IL7R','FOXP3','PTCRA',
                    'TBL1X','TBLR1','MYC','BIRC5','CCND1','AXIN2','TCF','CBP','P300','DKK1','BCL2','P53','SIAH1','CDK','JUN','CDH1'))
tallmarkers <- unique(c('CTNNB1','CD7','CD1A','CD34','CD38','CD4','CD44','CD8','CD3','TCF7','NOTCH1','LEF1','RUNX1','PROM1','CKIT','LGR5','HOX','LMO2','TCF','TCR','BCL11B','GATA3','E2A','EBF1','RAG','PAX5','EZH2','DOT1L','MLL','SUV420H1','SUV420H2','DNMT3A','TET2','BCL6','IL7R','FOXP3','PTCRA','FZD','LRP5','LRP6','DVL','AXIN','APC','GSK3B','CK1','TBL1X','TBLR1','MYC','BIRC5','CCND1','AXIN2','TCF','LEF1','CBP','P300','EP300','DKK1','BCL2','P53','SIAH1','CDK','JUN','CDH1'))

othermarkers <- unique(c('NFKB','P65','RELA','SMRT','NCOR','HDAC3','GPS2','SCF','CUL1','SKP1','RING','RBX1','ROC1','PI3K','AKT','PTEN','MTOR','MTORC1','MAPK','ERK','JNK','P38','ERK5','TNF','IL','IFN','NOTCH','JAKSTAT','VEGF','PARP1','PPP1R15A','RAD51','BRCA1','BRCA2','BAX','BIM','PUMA','MCL1','ABC','ABCB1','MDR1','ABCC1','ABCG2'))

allmarkers <- unique(c(ballmarkers,tallmarkers,othermarkers)) 


DotPlot(ball, allmarkers, cluster.idents = TRUE) + RotatedAxis()


interestingGenes <- c('CTNNB1','CD38','CD24','EZH2','TBL1X','BIRC5','BCL2','JUN','CD44','GSK3B','GPS2','PTEN','PARP1','PPP1R15A','BRCA1','BRCA2','BAX')

DotPlot(ball, interestingGenes, cluster.idents = TRUE) + RotatedAxis()




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
de2 <- FindMarkers(ball, ident.1 = 'YMK6', ident.2 = 'YMK4', group.by = 'orig.ident')

# save tables
write.csv(de1, file='~/Desktop/ball_allCells_6h-ctrl.csv', quote = FALSE)
write.csv(de2, file='~/Desktop/ball_allCells_24h-ctrl.csv', quote = FALSE)

# volcano plots with genes labelled
png(filename = '~/Desktop/B-ALL/ball_allCells_6h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de1$p_val_adj < .99)
plot(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'B-ALL All Cells: 6h - CTRL')
text(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), labels = rownames(de1)[filt])
dev.off()

png(filename = '~/Desktop/B-ALL/ball_allCells_24h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de2$p_val_adj < .99)
plot(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'B-ALL All Cells: 24h - CTRL')
text(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), labels = rownames(de2)[filt])
dev.off()




# DE between samples within groups (4&8, 1&7)


ball$group <- rep('other', ncol(ball))
ball$group[ball$seurat_clusters %in% c(4,8)] <- '4_8'
ball$group[ball$seurat_clusters %in% c(1,7)] <- '1_7'
Idents(ball) <- 'group'

# DE: JUST CLUSTERS 4&8 #
#########################

de1 <- FindMarkers(ball, ident.1 = 'YMK5', ident.2 = 'YMK4', group.by = 'orig.ident', subset.ident = '4_8')
de2 <- FindMarkers(ball, ident.1 = 'YMK6', ident.2 = 'YMK4', group.by = 'orig.ident', subset.ident = '4_8')

# save tables
write.csv(de1, file='~/Desktop/ball_clus4and8_6h-ctrl.csv', quote = FALSE)
write.csv(de2, file='~/Desktop/ball_clus4and8_24h-ctrl.csv', quote = FALSE)

# volcano plots with genes labelled
png(filename = '~/Desktop/B-ALL/ball_clus4and8_6h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de1$p_val_adj < .99)
plot(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'B-ALL Clusters 4 & 8: 6h - CTRL')
text(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), labels = rownames(de1)[filt])
dev.off()

png(filename = '~/Desktop/B-ALL/ball_clus4and8_24h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de2$p_val_adj < .99)
plot(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'B-ALL Clusters 4 & 8: 24h - CTRL')
text(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), labels = rownames(de2)[filt])
dev.off()


# DE: JUST CLUSTERS 1&7 #
#########################

de1 <- FindMarkers(ball, ident.1 = 'YMK5', ident.2 = 'YMK4', group.by = 'orig.ident', subset.ident = '1_7')
de2 <- FindMarkers(ball, ident.1 = 'YMK6', ident.2 = 'YMK4', group.by = 'orig.ident', subset.ident = '1_7')

# save tables
write.csv(de1, file='~/Desktop/ball_clus1and7_6h-ctrl.csv', quote = FALSE)
write.csv(de2, file='~/Desktop/ball_clus1and7_24h-ctrl.csv', quote = FALSE)

# volcano plots with genes labelled
png(filename = '~/Desktop/B-ALL/ball_clus1and7_6h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de1$p_val_adj < .99)
plot(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'B-ALL Clusters 1 & 7: 6h - CTRL')
text(de1$avg_log2FC[filt], -log10(de1$p_val_adj[filt]), labels = rownames(de1)[filt])
dev.off()

png(filename = '~/Desktop/B-ALL/ball_clus1and7_24h-ctrl_volcano.png', width = 12, height = 8, units = "in", res = 200)
filt <- which(de2$p_val_adj < .99)
plot(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), col='white',
     xlab = 'log2 Fold Change', ylab = '-log10 p-value',
     main = 'B-ALL Clusters 1 & 7: 24h - CTRL')
text(de2$avg_log2FC[filt], -log10(de2$p_val_adj[filt]), labels = rownames(de2)[filt])
dev.off()




