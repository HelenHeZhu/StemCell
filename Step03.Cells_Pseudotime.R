###
### Step03 slingshot, calculating the pseudotime 
### 
m.cluster  <- read.table('./Intermediate.Data.Files/Data1.Step02.Seurat.Epi.MetaData.Matrix.txt', header=TRUE, sep='\t')
PCA.matrix <- read.table('./Intermediate.Data.Files/Data1.Step02.Seurat.PCAs.Matrix.txt', header=TRUE, sep='\t')
m.cluster  <- m.cluster[rownames(PCA.matrix),]
###
###
library(slingshot)
pca_matrix        <- PCA.matrix
clus.pca          <- m.cluster$res.0.5
sds.pca           <- slingshot(pca_matrix[,1:12], clus.pca, start.clus=6)
save(sds.pca, file=file.path('./Intermediate.Data.Files/Data1.Step03.Slingshot.RImage.Rdata'))
#
Dims.pca          <- sds.pca@reducedDim                                    # pca dimension
clusterLabels.pca <- sds.pca@clusterLabels                                 # cluster labels
connectivity.pca  <- sds.pca@connectivity                                  # 
clusters.pca      <- rownames(connectivity.pca)                            #
nclus.pca         <- nrow(connectivity.pca)                                #
centers.pca       <- t(sapply(clusters.pca,function(clID){
                           x.sub <- Dims.pca[clusterLabels.pca == clID,]
                           return(colMeans(x.sub))
                 }))                                                               # cluster‘s centers in 3 dimension
rownames(centers.pca) <- clusters.pca                                              # add row names
Dims.pca              <- Dims.pca[ clusterLabels.pca %in% clusters.pca, ]          # 
clusterLabels.pca     <- clusterLabels.pca[clusterLabels.pca %in% clusters.pca]
#
sds.pca@lineages
xs <- c(NULL, centers.pca[, 1 ])
ys <- c(NULL, centers.pca[, 2 ])
zs <- c(NULL, centers.pca[, 3 ])
xs <- c( xs, as.numeric(sapply( sds.pca@curves, function(c){     c$s[, 1 ] }) ))
ys <- c( ys, as.numeric(sapply( sds.pca@curves, function(c){     c$s[, 2 ] }) ))
zs <- c( zs, as.numeric(sapply( sds.pca@curves, function(c){     c$s[, 3 ] }) ))
#
rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = 'iso', xlim = range(xs), ylim = range(ys), zlim = range(zs), 
	box=FALSE, axes=FALSE, xlab = '', ylab = '', zlab = '' )
colpal = rainbow(9)
#rgl::plot3d(Dims.pca, col=colpal[as.numeric(clusterLabels.pca)+1], add=TRUE, type='p', size=4, pch=16,alpha=I(1/20))  #, box=FALSE, axes=FALSE)
rgl::abclines3d(max(Dims.pca[,1]),max(Dims.pca[,2]),max(Dims.pca[,3]), a = diag(3), col = "black", lwd=2)
rgl::plot3d(centers.pca, size = 10, add = TRUE, pch=16, col = colpal[as.numeric(rownames(centers.pca))+1], alpha=1)
for (i in 1:(nclus.pca-1)){
	for (j in (i+1):nclus.pca){
		if (connectivity.pca[i,j]==1){
			rgl::lines3d(x=centers.pca[c(i,j),1], y=centers.pca[c(i,j),2], z=centers.pca[c(i,j),3], col='black', lwd=2)
		}
	}
}
rgl::rgl.snapshot('Slingshot.Cluster9.branchs.straightLine.png', fmt='png',top=TRUE)
# 7-2-1-5-6-8   7-2-4-9   7-2-1-3  

rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = 'iso', xlim = range(xs), ylim = range(ys), zlim = range(zs), 
	box=FALSE, axes=FALSE, xlab = '', ylab = '', zlab = '' )
colpal = rainbow(9)
rgl::plot3d(Dims.pca, col=colpal[as.numeric(clusterLabels.pca)+1], add=TRUE, type='p', size=4, pch=16,alpha=I(1/10))  #, box=FALSE, axes=FALSE)
rgl::abclines3d(max(Dims.pca[,1]),min(Dims.pca[,2]),min(Dims.pca[,3]), a = diag(3), col = "black", lwd=2)
rgl::plot3d(centers.pca, size = 10, add = TRUE, pch=16, col = colpal[as.numeric(rownames(centers.pca))+1], alpha=1)
coloursDict = c('red','blue','green','magenta','cyan','yellow','gray','black')
ids = 0
for(c in sds.pca@curves){ 
   ids = ids+1
   rgl::lines3d(c$s[c$tag,1:3], col= coloursDict[ids],lwd=2) 
}
# rgl::rgl.postscript('Corp.Slingshot.Cluster05.branchs.pdf', 'pdf')
# rgl::rgl.postscript('2018.01.05.Corp.Slingshot.Cluster5.branchs.svg', 'svg')
rgl::rgl.snapshot('Slingshot.Cluster9.branchs.fittedCurve.png', fmt='png',top=TRUE)



## order the cells on each branch
for (i in 1:ncol(pseudotime(sds.pca))){
	linedf <- data.frame(pseudotime=pseudotime(sds.pca)[,i], clus.labels = clusterLabels.pca, samples=rownames(Dims.pca))
	linedf <- linedf[with(linedf, order(pseudotime)),]
	medoids <- sapply(levels(linedf$clus.labels), function(clID){
		x.sub <- linedf$pseudotime[linedf$clus.labels == clID]
		col <- colpal[as.numeric(as.character(linedf$clus.labels))+1][which.max(linedf$clus.labels == clID)]
		return( list(means=mean(x.sub, na.rm=TRUE), sdev=sd(x.sub, na.rm=TRUE),col=col) )
	})
	means = unlist(medoids['means',])
	sdev  = unlist(medoids['sdev',])
	col   = unlist(medoids['col',])  
	pdf(paste('2018.02.18.Slingshot.Cluster09.Pseudo.branch_',i,'.pdf',sep=''), width=6,height=3)
	plot(linedf$pseudotime,rep(0, length(linedf$pseudotime)),cex=3,axes=F, pch=16, xlab='', ylab='', col=colpal[as.numeric(as.character(linedf$clus.labels))+1], ylim=c(-0.1, 0.1), xlim = range(linedf$pseudotime, na.rm=TRUE)); 
	abline(h=0, col="black")
	points(x=means,y=rep(0.07, length(means)), col=col, pch=19)
	arrows(means-sdev, rep(0.07, length(means)), means+sdev, rep(0.07, length(means)), length=0.05, angle=90, code=3, col=col)
	dev.off()
}
#
q()

##
## Plotting the heatmaps( ranked by cell clusters )
## 
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(slingshot)
load('./Intermediate.Data.Files/Data1.Step03.Slingshot.RImage.Rdata')
load('./Intermediate.Data.Files/Data1.Step02.Seurat.Image.Rda')  
#

KrtGenes    <- c('Krt5','Krt14','Krt8','Krt18','Krt19','Pbsn')
StemGenes   <- c('Trp63','Tpm2','Bcl2','Nkx3-1','Procr','Cd44',
                'Kit','Itga6','Ly6a','Prom1','Tert','Bmi1','Sox2','Alcam')
ArGenes     <- c('Mcm3','Iqgap2','Nkx3-1','Spdef','Mcm5','Dtl','Birc5','Endod1', 
              'Rfc3','Ncapd3','Atp1b1','Pbsn','Mccc2','Ccng2', 'Cdc25a','Taok3',
              'Acat2','Mcm4','Kcnn2', 'Ncapd3','Fen1')
CycleGenes  <- c('Aurkb','Bub1','Brca1','Brca2','Birc5','Bub1b','Cdc6','Cdkn2d','Ccnb1','Ccng2','Cdkn3','Cdkn2c','Cdca3','Ccnb1','Cdca2','Cdkn1a','Ccne1','Ccne2','Ccnf','Ccnb2','Cenpa','Cenpf','Cdc20','Cdc25b','Cdc25c','Cks1b','Cks2','Ccna2','Cdc25a','Cdca8','Cdk1','Dhfr','E2f1','Ece1','Gmnn','Msh2','Mki67','Mcm2','Nasp','Npat','Ndc80','Pcna','Pttg1','Rrm2','Rad51','Rad21','Rad51ap1','Racgap1','Rrm2','Rpa2','Rad51','Rrm1','Slbp','Tyms','Top2a','Tacc3','Usp1')
emtGenes    <- c('Cdh1','Epcam','Cldn1','Ocln','Vim','Fn1',
               'S100a4','Zeb1','Zeb2','Twist1','Twist2',
               'Snai1','Snai2','Prrx1','Prrx2','Foxc2',
               'Cdh2','Acta2')
#
## Manual running those codes for each gene sets
GeneMatrix <- subset(symbol, Symbol %in% KrtGenes)
#GeneMatrix <- subset(symbol, Symbol %in% StemGenes)
#GeneMatrix <- subset(symbol, Symbol %in% ArGenes)
#GeneMatrix <- subset(symbol, Symbol %in% CycleGenes)
#GeneMatrix <- subset(symbol, Symbol %in% emtGenes)

cellMetadata <- PG_all@meta.data
cellMetadata <- cellMetadata[order(cellMetadata$res.0.5),]
clusterID    <- paste('C', as.numeric(cellMetadata$res.0.5)+1,sep='')
#
scaledMatrix <- PG_all@scale.data[rownames(GeneMatrix),rownames(cellMetadata)]
rownames(scaledMatrix) <- GeneMatrix$Symbol
#
CountsMatrix <- as.matrix(PG_all@data)[rownames(GeneMatrix),rownames(cellMetadata)]
rownames(CountsMatrix) <- GeneMatrix$Symbol
#
htcol  <- rainbow(9)
ha_column = HeatmapAnnotation(df = data.frame(
					cluster = clusterID),
					col     = list(cluster=c(
								'C1'=htcol[1], 'C2'=htcol[2], 'C3'=htcol[3], 'C4'=htcol[4], 
								'C5'=htcol[5], 'C6'=htcol[6], 'C7'=htcol[7], 'C8'=htcol[8], 'C9'=htcol[9] )
					))
ht1 = Heatmap(scaledMatrix, name = "ht1", column_title = "ScaledCounts", top_annotation = ha_column, show_column_names=FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 10))
ht2 = Heatmap(CountsMatrix, name = "ht1", column_title = "RawCounts", top_annotation = ha_column, show_column_names=FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 10))
draw(ht1)
draw(ht2)
#
q()

## -------------------------------------
## plotting heatmaps of different lineages
## Complexheatmaps, ranked by pseudotime
##
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(slingshot)
load('./Intermediate.Data.Files/Data1.Step03.Slingshot.RImage.Rdata')
load('./Intermediate.Data.Files/Data1.Step02.Seurat.Image.Rda')  
##
plotPseodoHeatmap    <- function(x, symbollists, slingshotObj){
	i            <- x                 # branch id
	symbolList   <- symbollists       # input symbol list
	slingshotObj <- slingshotObj      # slingshot object
	##
	PseudoMax            <- max(pseudotime(slingshotObj), na.rm=TRUE)
    cellPseudoRank       <- data.frame(pseudotime  = pseudotime(slingshotObj)[,i], 
									   clus.labels = sds.pca@clusterLabels, 
									   samples=rownames(sds.pca@reducedDim) )
	cellPseudoRank       <- subset(cellPseudoRank, pseudotime!='NA')               
	cellPseudoRank       <- cellPseudoRank[with(cellPseudoRank, order(pseudotime)),]                                      # rank the cells by pseudotime value
	cellPseudoRank$cols  <- rainbow(max(as.numeric(sds.pca@clusterLabels))+1)[as.numeric(cellPseudoRank$clus.labels)]     # cluster color labels
	# 热图注释信息
	annoA     <- data.frame(cellPseudoRank[!duplicated(cellPseudoRank$clus.labels),c(2,4)])
	# 热图注释
	HA.list <- c('C1'=rainbow(9)[1], 'C2'=rainbow(9)[2], 'C3'=rainbow(9)[3], 'C4'=rainbow(9)[4],
				'C5'=rainbow(9)[5], 'C6'=rainbow(9)[6], 'C7'=rainbow(9)[7], 'C8'=rainbow(9)[8],
				'C9'=rainbow(9)[9])
	ha_column <- HeatmapAnnotation( df = data.frame(
					Pseudotime = cellPseudoRank$pseudotime,
					Cluster    = paste('C',as.numeric(cellPseudoRank$clus.labels),sep='')
					),                     # cellPseudoRank$clus.labels -> factors; auto add 1 when convert with as.numeric
					col=list(
					Pseudotime = colorRamp2( c(0, PseudoMax/2, PseudoMax), c('blue','yellow','red')),
					Cluster    = HA.list[as.numeric(annoA$clus.labels)]
					),
					show_annotation_name = TRUE
					)
	#
	genesList               <- subset(symbol, Symbol %in% symbolList)$Gene_ID      
	PartialMatrix           <- PG_all@scale.data[genesList,rownames(cellPseudoRank)]             
	rownames(PartialMatrix) <- symbol[rownames(PartialMatrix),'Symbol']
	ht1 = Heatmap(PartialMatrix, name = "Scaled\nUMI", column_title = paste('Branch',i, sep='_'), top_annotation = ha_column, show_column_names=FALSE, cluster_columns = FALSE,
	row_names_gp = gpar(fontsize = 9))
	draw(ht1)
}
#
KrtGenes    <- c('Krt5','Krt14','Krt8','Krt18','Krt19','Pbsn')
ArGenes     <- c('Mcm3','Iqgap2','Nkx3-1','Spdef','Mcm5','Dtl','Birc5','Endod1', 
                 'Rfc3','Ncapd3','Atp1b1','Pbsn','Mccc2','Ccng2', 'Cdc25a','Taok3',
                 'Acat2','Mcm4','Kcnn2', 'Ncapd3','Fen1')
CycleGenes  <- c('Aurkb','Bub1','Brca1','Brca2','Birc5','Bub1b','Cdc6','Cdkn2d','Ccnb1','Ccng2',
				 'Cdkn3','Cdkn2c','Cdca3','Ccnb1','Cdca2','Cdkn1a','Ccne1','Ccne2','Ccnf','Ccnb2',
				 'Cenpa','Cenpf','Cdc20','Cdc25b','Cdc25c','Cks1b','Cks2','Ccna2','Cdc25a','Cdca8',
				 'Cdk1','Dhfr','E2f1','Ece1','Gmnn','Msh2','Mki67','Mcm2','Nasp','Npat','Ndc80','Pcna',
				 'Pttg1','Rrm2','Rad51','Rad21','Rad51ap1','Racgap1','Rrm2','Rpa2','Rad51','Rrm1','Slbp','Tyms','Top2a','Tacc3','Usp1')
emtGenes    <- c('Cdh1','Epcam','Cldn1','Ocln','Vim','Fn1','S100a4','Zeb1','Zeb2','Twist1','Twist2',
               'Snai1','Snai2','Prrx1','Prrx2','Foxc2','Cdh2','Acta2')
#
Branches.num  <-  ncol(pseudotime(sds.pca))  
plotPseodoHeatmap(1, KrtGenes, sds.pca)
plotPseodoHeatmap(2, KrtGenes, sds.pca)
plotPseodoHeatmap(3, KrtGenes, sds.pca)
plotPseodoHeatmap(1, ArGenes, sds.pca)
plotPseodoHeatmap(2, ArGenes, sds.pca)
plotPseodoHeatmap(3, ArGenes, sds.pca)
plotPseodoHeatmap(1, CycleGenes, sds.pca)
plotPseodoHeatmap(2, CycleGenes, sds.pca)
plotPseodoHeatmap(3, CycleGenes, sds.pca)
plotPseodoHeatmap(1, emtGenes, sds.pca)
plotPseodoHeatmap(2, emtGenes, sds.pca)
plotPseodoHeatmap(3, emtGenes, sds.pca)
#
q()

## -------------------------------------
## plotting heatmaps of different lineages
## Complexheatmaps, ranked by brach clusters
##
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(slingshot)
load('./Intermediate.Data.Files/Data1.Step03.Slingshot.RImage.Rdata')
load('./Intermediate.Data.Files/Data1.Step02.Seurat.Image.Rda')  
##
sds.pca@lineages$Lineage
#
plotBranchClusterHeatmap    <- function(x, symbollists, slingshotObj){
	i            <- x                 # branch id
	symbolList   <- symbollists       # input symbol list
	slingshotObj <- slingshotObj      # slingshot object
	##
	PseudoMax            <- max(pseudotime(slingshotObj), na.rm=TRUE)
    cellPseudoRank       <- data.frame(pseudotime  = pseudotime(slingshotObj)[,i], 
									   clus.labels = sds.pca@clusterLabels, 
									   samples=rownames(sds.pca@reducedDim) )
	cellPseudoRank       <- subset(cellPseudoRank, pseudotime!='NA')               
	cellPseudoRank       <- cellPseudoRank[with(cellPseudoRank, order(pseudotime)),]                                      # rank the cells by pseudotime value
	cellPseudoRank$cols  <- rainbow(max(as.numeric(sds.pca@clusterLabels))+1)[as.numeric(cellPseudoRank$clus.labels)]     # cluster color labels
	clusterRanks <- as.numeric(slingshotObj@lineages[[i]])
	blankFrame   <- cellPseudoRank[c(1,2),]
	rownames(blankFrame) <- c('a','b') 
	for (m in clusterRanks) {                           # 
		blankFrame <- rbind(blankFrame, subset(cellPseudoRank, clus.labels==m))
	}
	cellPseudoRank <- blankFrame[-c(1,2),]
	#
	annoA     <- data.frame(cellPseudoRank[!duplicated(cellPseudoRank$clus.labels),c(2,4)])
	#
	HA.list <- c('C1'=rainbow(9)[1], 'C2'=rainbow(9)[2], 'C3'=rainbow(9)[3], 'C4'=rainbow(9)[4],
				'C5'=rainbow(9)[5], 'C6'=rainbow(9)[6], 'C7'=rainbow(9)[7], 'C8'=rainbow(9)[8],
				'C9'=rainbow(9)[9])
	ha_column <- HeatmapAnnotation( df = data.frame(
#					Pseudotime = cellPseudoRank$pseudotime,
					Cluster    = paste('C',as.numeric(cellPseudoRank$clus.labels),sep='')
					),                     # cellPseudoRank$clus.labels -> factors; auto add 1 when convert with as.numeric
					col=list(
#					Pseudotime = colorRamp2( c(0, PseudoMax/2, PseudoMax), c('blue','yellow','red')),
					Cluster    = HA.list[as.numeric(annoA$clus.labels)]
					),
					show_annotation_name = TRUE
					)
	#
	genesList               <- subset(symbol, Symbol %in% symbolList)$Gene_ID      
	PartialMatrix           <- PG_all@scale.data[genesList,rownames(cellPseudoRank)]             
	rownames(PartialMatrix) <- symbol[rownames(PartialMatrix),'Symbol']
	ht1 = Heatmap(PartialMatrix, name = "Scaled\nUMI", column_title = paste('Branch',i, sep='_'), top_annotation = ha_column, show_column_names=FALSE, cluster_columns = FALSE,
	row_names_gp = gpar(fontsize = 9))
	draw(ht1)
}
#
KrtGenes    <- c('Krt5','Krt14','Krt8','Krt18','Krt19','Pbsn')
ArGenes     <- c('Mcm3','Iqgap2','Nkx3-1','Spdef','Mcm5','Dtl','Birc5','Endod1', 
                 'Rfc3','Ncapd3','Atp1b1','Pbsn','Mccc2','Ccng2', 'Cdc25a','Taok3',
                 'Acat2','Mcm4','Kcnn2', 'Ncapd3','Fen1')
CycleGenes  <- c('Aurkb','Bub1','Brca1','Brca2','Birc5','Bub1b','Cdc6','Cdkn2d','Ccnb1','Ccng2',
				 'Cdkn3','Cdkn2c','Cdca3','Ccnb1','Cdca2','Cdkn1a','Ccne1','Ccne2','Ccnf','Ccnb2',
				 'Cenpa','Cenpf','Cdc20','Cdc25b','Cdc25c','Cks1b','Cks2','Ccna2','Cdc25a','Cdca8',
				 'Cdk1','Dhfr','E2f1','Ece1','Gmnn','Msh2','Mki67','Mcm2','Nasp','Npat','Ndc80','Pcna',
				 'Pttg1','Rrm2','Rad51','Rad21','Rad51ap1','Racgap1','Rrm2','Rpa2','Rad51','Rrm1','Slbp','Tyms','Top2a','Tacc3','Usp1')
emtGenes    <- c('Cdh1','Epcam','Cldn1','Ocln','Vim','Fn1','S100a4','Zeb1','Zeb2','Twist1','Twist2',
               'Snai1','Snai2','Prrx1','Prrx2','Foxc2','Cdh2','Acta2')
#
Branches.num  <-  ncol(pseudotime(sds.pca))  
plotBranchClusterHeatmap(1, KrtGenes, sds.pca)
plotBranchClusterHeatmap(2, KrtGenes, sds.pca)
plotBranchClusterHeatmap(3, KrtGenes, sds.pca)
plotBranchClusterHeatmap(1, ArGenes, sds.pca)
plotBranchClusterHeatmap(2, ArGenes, sds.pca)
plotBranchClusterHeatmap(3, ArGenes, sds.pca)
plotBranchClusterHeatmap(1, CycleGenes, sds.pca)
plotBranchClusterHeatmap(2, CycleGenes, sds.pca)
plotBranchClusterHeatmap(3, CycleGenes, sds.pca)
plotBranchClusterHeatmap(1, emtGenes, sds.pca)
plotBranchClusterHeatmap(2, emtGenes, sds.pca)
plotBranchClusterHeatmap(3, emtGenes, sds.pca)
#
q()


############ Abandoned
#=============================================================================================
# plotting heatmaps of different lineages 
#dim(scaledMatrix)   # 54 9278
#clusterID[1:3]      # "C1" "C1" "C1"
#length(clusterID)   #  9278
#old_Cluster <- as.numeric(cellMetadata$res.0.5)+1
#scaledMatrix_T <- t(scaledMatrix)
#dim(scaledMatrix_T)
#scaledMatrix_T <- data.frame(scaledMatrix_T)
#scaledMatrix_T$old_Cluster <- as.numeric(old_Cluster)
#head(scaledMatrix_T)
## write.csv(scaledMatrix_T,'scaledMatrix_T.csv',row.names = F)
## 7 2 1 5 6 8 
## 7 2 4 9 
## 7 2 1 3
#Old_Index <- scaledMatrix_T$old_Cluster
#Old_Index[which(Old_Index == 7)] <- 11
#Old_Index[which(Old_Index == 2)] <- 12
#Old_Index[which(Old_Index == 1)] <- 13
#Old_Index[which(Old_Index == 5)] <- 14
#Old_Index[which(Old_Index == 8)] <- 15
#Old_Index[which(Old_Index == 4)] <- 16
#Old_Index[which(Old_Index == 9)] <- 17
#Old_Index[which(Old_Index == 3)] <- 18
#Old_Index[which(Old_Index == 6)] <- 19
#scaledMatrix_T$Index <- Old_Index
#Target_Data <- scaledMatrix_T[order(scaledMatrix_T$Index),]
#dim(Target_Data)
#head(Target_Data)
##=================================
## (lineage 1)7 2 1 5 6 8
#New_Target_Data <- Target_Data[c(which(Target_Data$Index == 11),
#                                 which(Target_Data$Index == 12),
#                                 which(Target_Data$Index == 13),
#                                 which(Target_Data$Index == 14),
#                                 which(Target_Data$Index == 19),
#                                 which(Target_Data$Index == 15)),]
#dim(New_Target_Data)
#head(New_Target_Data)
#size_New_Target_Data <- dim(New_Target_Data)
#New_Target_Data_2 <- New_Target_Data[,c(-(size_New_Target_Data[2]-1),-size_New_Target_Data[2])]
#Cor_Data <- t(New_Target_Data_2)
#Cor_Data[1:5,1:3]
#dim(Cor_Data)
#Indexs <- New_Target_Data$old_Cluster
#clusterID_zw <- rep('C1',1)
#for (K in  1:length(Indexs))
#{
#  clusterID_zw[K] <- paste('C',Indexs[K],sep = '')
#}
#ha_column = HeatmapAnnotation(df = data.frame(
#  cluster = clusterID_zw),
#  col     = list(cluster=c(
#    'C1'=htcol[1], 'C2'=htcol[2], 'C6'=htcol[6],
#    'C5'=htcol[5], 'C7'=htcol[7], 'C8'=htcol[8] )
#  ))
#ht2 = Heatmap(Cor_Data, 
#              name = "ht1", 
#              column_title = "ScaledCounts", 
#              top_annotation = ha_column, 
#              show_column_names=FALSE, 
#              cluster_columns = FALSE, 
#              row_names_gp = gpar(fontsize = 10))
#draw(ht2)
#
##=================================
## (lineage 2)7 2 4 9  
#New_Target_Data <- Target_Data[c(which(Target_Data$Index == 11),
#                                 which(Target_Data$Index == 12),
#                                 which(Target_Data$Index == 16),
#                                 which(Target_Data$Index == 17)),]
#dim(New_Target_Data)
#head(New_Target_Data)
#size_New_Target_Data <- dim(New_Target_Data)
#New_Target_Data_2 <- New_Target_Data[,c(-(size_New_Target_Data[2]-1),-size_New_Target_Data[2])]
#Cor_Data <- t(New_Target_Data_2)
#Cor_Data[1:5,1:3]
#dim(Cor_Data)
#Indexs <- New_Target_Data$old_Cluster
#clusterID_zw <- rep('C1',1)
#for (K in  1:length(Indexs))
#{
#  clusterID_zw[K] <- paste('C',Indexs[K],sep = '')
#}
#ha_column = HeatmapAnnotation(df = data.frame(
#  cluster = clusterID_zw),
#  col     = list(cluster=c(
#    'C4'=htcol[4], 'C2'=htcol[2], 
#    'C9'=htcol[9], 'C7'=htcol[7])
#  ))
#ht2 = Heatmap(Cor_Data, 
#              name = "ht1", 
#              column_title = "ScaledCounts", 
#              top_annotation = ha_column, 
#              show_column_names=FALSE, 
#              cluster_columns = FALSE, 
#              row_names_gp = gpar(fontsize = 10))
#draw(ht2)
#
##=================================
## (lineage 3)7 2 1 3
#New_Target_Data <- Target_Data[c(which(Target_Data$Index == 11),
#                                 which(Target_Data$Index == 12),
#                                 which(Target_Data$Index == 13),
#                                 which(Target_Data$Index == 18)),]
#dim(New_Target_Data)
#head(New_Target_Data)
#size_New_Target_Data <- dim(New_Target_Data)
#New_Target_Data_2 <- New_Target_Data[,c(-(size_New_Target_Data[2]-1),-size_New_Target_Data[2])]
#Cor_Data <- t(New_Target_Data_2)
#Cor_Data[1:5,1:3]
#dim(Cor_Data)
#Indexs <- New_Target_Data$old_Cluster
#clusterID_zw <- rep('C1',1)
#for (K in  1:length(Indexs))
#{
#  clusterID_zw[K] <- paste('C',Indexs[K],sep = '')
#}
#ha_column = HeatmapAnnotation(df = data.frame(
#  cluster = clusterID_zw),
#  col     = list(cluster=c(
#    'C1'=htcol[1], 'C2'=htcol[2], 
#    'C7'=htcol[7], 'C3'=htcol[3] )
#  ))
#ht2 = Heatmap(Cor_Data, 
#              name = "ht1", 
#              column_title = "ScaledCounts", 
#              top_annotation = ha_column, 
#              show_column_names=FALSE, 
#              cluster_columns = FALSE, 
#              row_names_gp = gpar(fontsize = 10))
#draw(ht2)



