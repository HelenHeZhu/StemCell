###
### Step02 clustering using Seurat
###
library(Seurat)
library(dplyr)
library(Matrix)
#
data   = read.table('./Intermediate.Data.Files/Data1.Step01.Seurat.Epi.Expression.Matrix.raw.txt', header=TRUE, sep='\t')
symbol = read.table('./Intermediate.Data.Files/Data1.Step01.Seurat.Epi.Symbols.txt', header=TRUE, sep='\t')
#
PG_all       <- CreateSeuratObject(raw.data = data, min.cells = 0, min.genes = 0, project = "PG")
mito.symbol  <- grep(pattern='^mt-', x=symbol$Symbol, value=TRUE)
mito.genes   <- as.character(subset(symbol, Symbol %in% mito.symbol)$Gene_ID)
mito.genes
percent.mito <- Matrix::colSums(PG_all@raw.data[mito.genes,])/Matrix::colSums(PG_all@raw.data)
PG_all       <- AddMetaData(PG_all, metadata=percent.mito, col.name='percent.mito')
VlnPlot(PG_all, features.plot=c('nGene','nUMI','percent.mito'), nCol=3)
#
PG_all       <- NormalizeData(PG_all, normalization.method='LogNormalize', scale.factor=1e4)
PG_all       <- FindVariableGenes(PG_all, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5)
length(PG_all@var.genes)
PG_all <- ScaleData(PG_all, vars.to.regress=c('nUMI','percent.mito'))
PG_all <- RunPCA(PG_all, pc.genes=PG_all@var.genes, do.print=TRUE, pcs.print=1:5, genes.print=5)
PCAPlot(PG_all, dim.1=1, dim.2=2)
PG_all <- ProjectPCA(PG_all, do.print=FALSE)
PCHeatmap(PG_all, pc.use=1:9,   cells.use=500, do.balanced=TRUE, label.columns=FALSE, use.full=FALSE)
PCHeatmap(PG_all, pc.use=10:18, cells.use=500, do.balanced=TRUE, label.columns=FALSE, use.full=FALSE)
#
PG_all <- FindClusters(PG_all, reduction.type='pca', dims.use=1:12, resolution=c(0.6,0.7,0.8,0.9,1,0.5), print.output=0, save.SNN=TRUE)
PG_all <- RunTSNE(PG_all, dims.use=1:12, do.fast=TRUE)
TSNEPlot(PG_all, do.label=T)
PCAPlot(PG_all, dim.1=1, dim.2=2)
#

save(PG_all, symbol, file=file.path('./Intermediate.Data.Files/Data1.Step02.Seurat.Image.Rda'))
write.table(as.matrix(PG_all@raw.data),                        quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.Expression.Matrix.raw.txt' )
write.table(expm1(as.matrix(PG_all@data)),                     quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.Expression.Matrix.nor.txt' )
write.table(PG_all@meta.data,                                  quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.MetaData.Matrix.txt')
write.table(symbol, row.names=T,                               quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.Symbols.txt')
write.table(PG_all@dr$pca@cell.embeddings,                     quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step02.Seurat.PCAs.Matrix.txt')
write.table(PG_all@dr$tsne@cell.embeddings,                    quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step02.Seurat.tSNE.Matrix.txt')
varGenes <- subset(symbol, Gene_ID %in% PG_all@var.genes)
write.table(varGenes,                                          quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step02.Seurat.varGenes.Matrix.txt')
#
PG.markers <- FindAllMarkers(PG_all, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
PG.markers$symbol  <- symbol[PG.markers$gene,'Symbol']
PG.markers$cluster <- as.numeric(PG.markers$cluster)
write.table(PG.markers, quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.Clusters.DGEs.Markers.txt')


##ploting the tsne mappig
tsneplot <- function(PG_all, cluster, shift=0, main='') {
    tsne_matrix       = data.frame(PG_all@dr$tsne@cell.embeddings)
    tsne_matrix$ident = as.numeric(PG_all@meta.data[, cluster])
	tsne_matrix$ident = as.numeric(tsne_matrix$ident) + shift
    col_pal           = rainbow(max(tsne_matrix$ident))
    tsne_matrix %>% dplyr::group_by(ident) %>% summarize(median(tSNE_1), median(tSNE_2)) -> centers
    centers = data.frame(centers)
    plot(tsne_matrix[,1:2], pch=16, cex=0.5, col=col_pal[as.numeric(tsne_matrix$ident)], main=main)
    text(centers[,c(2,3)],as.character(centers[,c(1)]))
}
##
gene_tsne_plot <- function(PG_all, x){
    print (subset(symbol,Symbol %in% x ))
    FeaturePlot(PG_all, features.plot=c( as.character(subset(symbol,Symbol %in% x )$Gene_ID) ), cols.use=c('grey','blue'), reduction.use='tsne',pt.size = 0.7)
}
gene_tsne_plot(PG_all, c('Pecam1','Cd34','Icam1','Vwf'))
gene_tsne_plot(PG_all, c('Lyve1', 'Tek','Vcam1','Cdh5'))
gene_tsne_plot(PG_all, c('Eng','S1pr1','Emcn'))
#####


### 
### heatmaps on the tsne maps
### 
library('RColorBrewer') 
library('dplyr')
#
norm_exp <- function(PG_all){
	fPlot <- t(PG_all@scale.data)
	fPlot <- cbind(fPlot, PG_all@dr$tsne@cell.embeddings)
	return( data.frame(fPlot) )
}
fPlot <- norm_exp(PG_all)
pal   <- colorRampPalette(rev(brewer.pal(n=7, name='RdYlBu')))(200)
#
GeneExpressionTrendsPlot <- function(fPlot, gene){
	GeneID <- rownames(subset(symbol,Symbol %in% gene ))
	ggplot(data.frame(fPlot), aes_string('tSNE_1','tSNE_2', color=GeneID))  + 
		geom_point(size=1,pch=20) + scale_color_gradientn(colours=pal) + ggtitle( gene ) + theme_void() + 
		theme(legend.position='bottom', legend.direction='horizontal', legend.title=element_blank(), plot.title=element_text(face='italic'))
}
#
listsA <- c('Cdh1','Epcam','Cldn1','Ocln','Vim','Fn1','S100a4','Zeb1','Zeb2','Twist1','Twist2','Snai1','Snai2','Prrx1','Prrx2','Foxc2','Cdh2','Acta2')
genes  <- listsA[listsA %in% symbol$Symbol]
for (i in genes){
	p <- GeneExpressionTrendsPlot(fPlot, i)
	ggsave(paste('SupFigures.004.Epithelial.tsne.','heatmap_',i,'.pdf',sep=''), p, width=10,height=10 )
}






















