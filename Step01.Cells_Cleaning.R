###
### Step01 removal of poor quality cells and contaminated non-epithelial cells 
###

library(Seurat)
library(dplyr)
library(Matrix)
#
data             <- read.table('./Raw.Data.Files/gene_cell_exprs_table.xls', header=TRUE, sep='\t')
symbol           <- data[,1:2]
rownames(symbol) <- symbol$Gene_ID
data             <- data[,-c(1,2)]
rownames(data)   <- rownames(symbol)
#
PG_all       <- CreateSeuratObject(raw.data = data, min.cells = 0, min.genes = 0, project = "PG")
mito.symbol  <- grep(pattern='^mt-', x=symbol$Symbol, value=TRUE)
mito.genes   <- as.character(subset(symbol, Symbol %in% mito.symbol)$Gene_ID)
mito.genes
#
percent.mito <- Matrix::colSums(PG_all@raw.data[mito.genes,])/Matrix::colSums(PG_all@raw.data)
PG_all       <- AddMetaData(PG_all, metadata=percent.mito, col.name='percent.mito')
VlnPlot(PG_all, features.plot=c('nGene','nUMI','percent.mito'), nCol=3)

# Poor-quality cells with less than 1000 genes detected, less than 5000 UMIs or more than 5% UMI mapped to mitochondria genes were removed. 
PG_all <- FilterCells(object = PG_all, subset.names = c("nGene", "nUMI", "percent.mito"), low.thresholds = c(1000, 5000, -Inf), high.thresholds = c(Inf, Inf,0.05))
VlnPlot(PG_all, features.plot=c('nGene','nUMI','percent.mito'), nCol=3)
#
PG_all <- NormalizeData(PG_all, normalization.method='LogNormalize', scale.factor=1e4)
PG_all <- FindVariableGenes(PG_all, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5)
length(PG_all@var.genes)
#
PG_all <- ScaleData(PG_all, vars.to.regress=c('nUMI','percent.mito'))
PG_all <- RunPCA(PG_all, pc.genes=PG_all@var.genes, do.print=TRUE, pcs.print=1:5, genes.print=5)
PCAPlot(PG_all, dim.1=1, dim.2=2)
# 
PG_all <- ProjectPCA(PG_all, do.print=FALSE)
PCHeatmap(PG_all, pc.use=1:9,   cells.use=500, do.balanced=TRUE, label.columns=FALSE, use.full=FALSE)
PCHeatmap(PG_all, pc.use=10:18, cells.use=500, do.balanced=TRUE, label.columns=FALSE, use.full=FALSE)
#
PG_all <- FindClusters(PG_all, reduction.type='pca', dims.use=1:9, resolution=c(0.6,0.7,0.8,0.9,1,0.5), print.output=0, save.SNN=TRUE)
PG_all <- RunTSNE(PG_all, dims.use=1:9, do.fast=TRUE)
TSNEPlot(PG_all, do.label=T)
PCAPlot(PG_all, dim.1=1, dim.2=2)
####

####
#### retrieve non-epithelial cells
FeaturePlot(PG_all, features.plot=c( as.character(subset(symbol,Symbol %in% c('Cd74','Cd72','Cd54') )$Gene_ID) ), cols.use=c('yellow','red'), reduction.use='tsne',pt.size = 0.7,nCol=3)
FeaturePlot(PG_all, features.plot=c( as.character(subset(symbol,Symbol %in% c('Eng','S1pr1','Emcn') )$Gene_ID) ), cols.use=c('yellow','red'), reduction.use='tsne',pt.size = 0.7,nCol=3)
#
Immu.Cluster <- FindMarkers(PG_all, ident.1='9', only.pos=TRUE, min.pct=0.25)
Endo.Cluster <- FindMarkers(PG_all, ident.1='4', only.pos=TRUE, min.pct=0.25)
Immu.Cluster$Symbol <- symbol[rownames(Immu.Cluster),'Symbol']
Endo.Cluster$Symbol <- symbol[rownames(Endo.Cluster),'Symbol']
write.table(Immu.Cluster, quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step01.ImmuCluster.DEGs.txt')
write.table(Endo.Cluster, quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step01.EndoCluster.DEGs.txt')
##
save(PG_all, symbol, file=file.path('./Intermediate.Data.Files/Data1.Step01.Seurat.Image.Rda'))
##
epi_cells <- rownames(subset(PG_all@meta.data, res.0.5!=4 & res.0.5!=9))
#
write.table(as.matrix(PG_all@raw.data[,epi_cells]),            quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step01.Seurat.Epi.Expression.Matrix.raw.txt' )
write.table(expm1(as.matrix(PG_all@data[,epi_cells])),         quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step01.Seurat.Epi.Expression.Matrix.nor.txt' )
write.table(subset(PG_all@meta.data, res.0.5!=4 & res.0.5!=9), quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step01.Seurat.Epi.MetaData.Matrix.txt')
write.table(symbol, row.names=T,                               quote=F, sep='\t', file='./Intermediate.Data.Files/Data1.Step01.Seurat.Epi.Symbols.txt')

###
### heatmap for non-epithelial cells markers
### 
library(ComplexHeatmap)
library(circlize)
#
geneSets                 <- subset(symbol, Symbol %in% c('Cd74','Cd72','Cd54','Eng','S1pr1','Emcn'))
geneSets                 <- geneSets[with(geneSets, order(Symbol)),]
cellRanks                <- PG_all@meta.data[with(PG_all@meta.data, order(res.0.5)),]
PartialMatrix            <- PG_all@scale.data[rownames(geneSets), rownames(cellRanks)]
rownames(PartialMatrix)  <- geneSets$Symbol
cellRanks$col            <- rainbow(max(as.numeric(cellRanks$res.0.5))+1)[as.numeric(cellRanks$res.0.5)+1]
#
ha_column <- HeatmapAnnotation(
				df  = data.frame(
					ClusterID = as.numeric(cellRanks$res.0.5)
				),
				col = list(
					ClusterID = colorRamp2(unique(as.numeric(cellRanks$res.0.5)), rainbow(max(as.numeric(cellRanks$res.0.5))+1))
				)
)
ht1 = Heatmap(PartialMatrix, name = "Scaled\nUMI", column_title = "non-epithelial cells markers", top_annotation = ha_column, show_column_names=FALSE, cluster_columns = FALSE,
row_names_gp = gpar(fontsize = 10))
draw(ht1)


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
listsA <- c('Cd74', 'Cd72','Cd54','Eng','S1pr1','Emcn')
genes  <- listsA[listsA %in% symbol$Symbol]
for (i in genes){
	p <- GeneExpressionTrendsPlot(fPlot, i)
	ggsave(paste('SupFigures.004.non-Epithelial.tsne.','heatmap_',i,'.pdf',sep=''), p, width=10,height=10 )
}
