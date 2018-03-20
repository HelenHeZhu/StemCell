###
### clustering and pseudotime calculating
###
library('monocle')
library(dplyr)
#
data     = read.table('./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.Expression.Matrix.raw.txt', header=TRUE, sep='\t')
labels   = read.table('./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.MetaData.Matrix.txt',header=TRUE, sep='\t')
labels[,5:10]           <- labels[,5:10] + 1
labels[,5:ncol(labels)] <- lapply(labels[,5:ncol(labels)], as.character)
labels                  <- labels[,1:10]
symbol           = read.table('./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.Symbols.txt', header=TRUE, sep='\t')
colnames(symbol) = c('Gene_ID', 'gene_short_name')
#
pd <- new('AnnotatedDataFrame', data=labels)
fd <- new('AnnotatedDataFrame', data=symbol[rownames(data),])
HSMM <- newCellDataSet(as.matrix(data), phenoData = pd, featureData = fd, lowerDetectionLimit = 1, expressionFamily=negbinomial.size())
## main estimate_size_and_dispersion
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
## 
HSMM <- detectGenes(HSMM, min_expr = 0.1)
#
disp_table             <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM                   <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F)   
#
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 12, reduction_method = 'tSNE', verbose = T) 
#
HSMM <- clusterCells(HSMM, num_clusters=11, k=12, method='densityPeak')
#
plot_tsne <- function(HSMM, cluster, shift=0, main='') {
	tsne_matrix           = data.frame(t(HSMM@reducedDimA))
	colnames(tsne_matrix) = c('tSNE_1','tSNE_2')
	tsne_matrix$ident = as.numeric(pData(HSMM)[,cluster])
	tsne_matrix$ident = tsne_matrix$ident + shift
	colpal            = rainbow(max(as.numeric(tsne_matrix$ident)))	
	tsne_matrix %>% dplyr::group_by(ident) %>% summarize(median(tSNE_1), median(tSNE_2)) -> centers
	centers = data.frame(centers)
	plot(tsne_matrix[,1:2], pch=16, cex=0.5, col=colpal[as.numeric(tsne_matrix[,'ident'])], main=main)
	text(centers[,c(2,3)],as.character(centers[,c(1)]))
}
#
plot_tsne(HSMM, 'Cluster')
plot_tsne(HSMM, 'res.0.5')
##
##
HSMM2           <- HSMM
expressed_genes <- row.names(subset(fData(HSMM2), num_cells_expressed >= 10))
diff_test_res   <- differentialGeneTest(HSMM2[expressed_genes,], fullModelFormulaStr = "~Cluster", cores=detectCores())
ordering_genes  <- row.names (subset(diff_test_res, qval < 0.01))
HSMM2           <- setOrderingFilter(HSMM2, ordering_genes)
plot_ordering_genes(HSMM2)
HSMM2 <- reduceDimension(HSMM2, max_components = 2, method = 'DDRTree')
HSMM2 <- orderCells(HSMM2,root_state = 5)
plot_cell_trajectory(HSMM2, color_by = "res.0.5")
q()




###
### DiffusionMap, tips & pseudotime
###
library('scran')
library('destiny')
library('plyr')
library('dplyr')
library('ggplot2')
#
data           <- read.table('./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.Expression.Matrix.nor.txt', header=TRUE, sep='\t')
labels         <- read.table('./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.MetaData.Matrix.txt', header=TRUE, sep='\t')
labels         <- labels[colnames(data),]
labels[,5:10]  <- labels[,5:10] + 1
labels         <- labels[,1:10]
#
varGenes <- read.table('./Intermediate.Data.Files/Data1.Step02.Seurat.varGenes.Matrix.txt',header=TRUE, sep='\t')
#
m.vp3 <- data[rownames(varGenes),]
m.vp3 <- t(log(m.vp3+1))
#
dm3  <- DiffusionMap(m.vp3, n_eigs=20, rotate=TRUE)
dpt3 <- DPT(dm3)
#
colpal = rainbow(max(labels$res.0.5))
plot(dpt3, pch=16, cex=0.5, col=colpal[labels$res.0.5] )
q()







#######
####### Go bar plot (go analysis was done on the DAVID website )
####### 
David.Goa <- read.table('Cluster7.DGEs.all.go',  header=TRUE, sep='\t')
David.Gob <- read.table('Cluster7.DGEs.Up.go',   header=TRUE, sep='\t', quote='')
David.Goc <- read.table('Cluster7.DGEs.Down.go', header=TRUE, sep='\t', quote='')
### plot the log10(Pvalue)
gobarplot <- function(davidgo, Categorys, tops=25, width=0.5, widthy=1, ... ) {
	if (Categorys == 'GOTERM_BP_DIRECT' ) { lists <- tidyr::separate(subset(davidgo, Category==Categorys), Term, c('ID', 'Term'), sep='~')	}
	if (Categorys == 'KEGG_PATHWAY' ) {     lists <- tidyr::separate(subset(davidgo, Category==Categorys), Term, c('ID', 'Term'), sep=':')	}
	lists      <- lists[,c('Category','ID','Term','Count','PValue','FDR')]
	lists$log2PValue <- -log10(lists$PValue)
	lists      <- lists[order(lists$PValue),]
	if (tops>nrow(lists)) {tops=nrow(lists)}
	lists      <- lists[1:tops,]
	rownames(lists) <- as.character(lists$Term)
	#
	listA         <- lists[,c('Term','log2PValue')]
	listA$orderID <- 1:nrow(listA)
	if (Categorys == 'GOTERM_BP_DIRECT' ){ fillcol='red' }else{fillcol='cornflowerblue'}
	if (Categorys == 'GOTERM_BP_DIRECT' ){ ggtitles='Cluster7 vs OtherClusters (GO analysis)'}else{ ggtitles='Cluster7 vs OtherClusters (KEGG analysis)'}
	p1 <- ggplot(listA, aes(orderID, log2PValue)) + ggtitle( ggtitles )
	p1 <- p1 + geom_bar(stat="identity" , width = width, fill = fillcol)
	p1 + scale_x_continuous(breaks=listA$orderID, labels=listA$Term) + 
	theme(axis.title.y=element_blank(),axis.text.x = element_text(size = rel(1.2), angle=0),axis.text.y = element_text(size = rel(widthy), angle=0)) + 
	theme(panel.background = element_blank(),panel.border = element_blank()) + coord_flip() 
}
gobarplot(David.Goa, 'KEGG_PATHWAY',     tops=25   )
gobarplot(David.Goa, 'GOTERM_BP_DIRECT', tops=125, widthy=0.6)
gobarplot(David.Goa, 'GOTERM_BP_DIRECT', tops=60, widthy=0.8)
#
gobarplot(David.Gob, 'KEGG_PATHWAY',     tops=25,  )
gobarplot(David.Gob, 'GOTERM_BP_DIRECT', tops=125, widthy=0.6)
gobarplot(David.Goc, 'GOTERM_BP_DIRECT', tops=60, widthy=0.8)
#
gobarplot(David.Goc, 'KEGG_PATHWAY',     tops=25,  )
gobarplot(David.Goc, 'GOTERM_BP_DIRECT', tops=125, widthy=0.6)
gobarplot(David.Goc, 'GOTERM_BP_DIRECT', tops=60, widthy=0.8)

