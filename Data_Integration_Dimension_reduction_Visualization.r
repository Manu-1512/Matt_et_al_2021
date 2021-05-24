###############

# ====================================
# Author: Manu Singh
# Date: 2021-05-18
# Title: colon-seq data analysis for NIBD and CD patients
# ====================================

# =====================================
# Enviornment variables
# =====================================
setwd("~/Collab/Praveen/seurat_obj/Doublet_removed_Annotated")
getwd()
Sys.Date()
main_dir <- getwd()
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()
options(future.globals.maxSize = 6096*1524^2 )
set.seed(786)

# =====================================
# Load libraries
# =====================================
library(Seurat)
library(ggplot2)
library(ggpubr)
library(RCurl)
library(gplots)
library(BSDA)
library(reshape)
library(sctransform)
library(scales)
library(monocle)
##############################


seurat_files <- list.files()

colon_list <- NULL
for(j in 1:length(seurat_files)) {
colon_file <- seurat_files[j]
colon_file_name <- basename(colon_file) %>% gsub(".rds", "", .)
colon_list <- append(colon_list, readRDS(colon_file))
}


integrate <- function(colon_list){
colon_features <- SelectIntegrationFeatures(object.list = colon_list, nfeatures = 3500)
colon_list <- PrepSCTIntegration(object.list = colon_list, anchor.features = colon_features)
colon_anchors <- FindIntegrationAnchors(object.list = colon_list,
normalization.method = "SCT",
anchor.features = colon_features,
 dims = 1:30,
verbose = TRUE)

colon_integrated <- IntegrateData(anchorset = colon_anchors,
normalization.method = "SCT",
dims = 1:30,
verbose = TRUE)
colon_integrated <- FindVariableFeatures(colon_integrated, selection.method = "vst", nfeatures = 3500)
colon_integrated <- ScaleData(colon_integrated, vars.to.regress = c("S.Score", "G2M.Score", "nUMI", "ribo.genes"), features = row.names(colon))

colon_integrated <- RunPCA(colon_integrated, features = VariableFeatures(object = colon_integrated), npcs = 30, verbose = FALSE)
colon_integrated <- FindNeighbors(colon_integrated, dims = 1:30)
colon_integrated <- RunUMAP(colon_integrated, dims = 1:30)
return(colon_integrated)
}

colon_integrated <- integrate(colon_list)

DefaultAssay(colon_integrated) <- "RNA"

colon_integrated@meta.data$celltypes <- Idents(colon_integrated)
colon_integrated@meta.data$patients <- gsub("._$", "", colon_integrated@meta.data$orig.ident)


Idents(colon_integrated) <- factor(Identscolon_integrated), levels = c("EEC", "G2_M_G1_TA", "S_phase_TA",  "Stem", "SPIB+_cluster", "Tuft_cells",
    "CP", "IC", "CA1+_early_colonocyte", "CA1+_late_colonocyte", "CEACAM7+_colonocyte", "secretory_progenitor",
    "immature_goblet", "goblet"))
    

## Plot the UMAP
png("Praveen_UMAP_reannotated.png",  width = 28.7, height = 20.8, units = "cm", res = 600, pointsize = 12)
DimPlot(colon_integrated, reduction = "umap", cols=c("black", "orange","darkorange", "forestgreen", "gold", "yellowgreen", "royalblue", "blue", "cyan", "slateblue","midnightblue","magenta","red","maroon"), pt.size=1)
dev.off()


my_genes <- c("LGR5", "ASCL2", "SMOC2", "ZG16", "MUC2", "CA1", "FCGBP",
  "MEIS1", "BEST4", "SPIB", "CA7", "CEACAM1", "CEACAM7", "GUCA2A",
  "ME1", "CA2", "KLF2", "MDM2", "DNAJB1", "HSPA1B", "MKI67", "UBE2C", "STMN1", "ESRRG", "RP11-420N3.2", "BCL2", 
  "PTGS1", "SOX8", "DRGX", "NKX2-2", "SCGN", "SCG5")


my_genes <- my_genes[c(4:7,12:18,10,27:29,8,9,11,1:3,19:23,30:32)]

png("Praveen_DotPlot_reannotated_first_clustering_genes.png",  width = 38.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
DotPlot(colon_integrated, features=my_genes, scale.min=60, scale.max=100, dot.min=0) + scale_colour_gradient2(low = "white", mid = "steelblue", high = "red2") + RotatedAxis()
dev.off()


############################ Most variable TEs belong to certain families


mvgs  <- VariableFeatures(colon_integrated)
mvte <- mvgs[grep("TELONG", mvgs)]

man  <- (sapply(strsplit(ma,"-CHR"),"[[",1))

dat <- as.data.frame(table(man))


dat <- dat[order(dat$Freq),]

da <- tail(dat, n=36)

row.names(da) <- da[,1]

png("MVTEs_barplots_Praveen_Integrated.png",  width = 26.7, height = 16.8, units = "cm", res = 300, pointsize = 12)
 par(mar=c(15.1,4.1,4.1,2.1))
 my_bar <- barplot(as.matrix(da$Freq), border=T , beside=T, names.arg=row.names(da), ylim=c(0,250), las=2 , space = 0.2,
                  col=c(rep("forestgreen", 9) , rep("gold", 1) , rep("forestgreen", 3), 
                  rep("violet", 1), rep("forestgreen", 1),
                  rep("steelblue1", 3), rep("violet", 1),
                  rep("forestgreen", 1), rep("gold", 1), rep("violet", 3),rep("steelblue1", 1),
                   rep("violet", 3), rep("gold", 1),rep("forestgreen", 1), rep("violet", 6)))
legend("topleft", legend = c("LTR/ERVs","Old L1","Young L1", "THE1/MST") , 
     col = c("forestgreen", "violet", "gold", "steelblue1") , 
     bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0.05))
dev.off()


#################################  

sce <- as.SingleCellExperiment(colon_integrated)
exprMat <- as.matrix(counts(sce))
cluster.averages <- AverageExpression(colon_integrated) 
cellInfo <- data.frame(seuratCluster=Idents(colon_integrated))

my_genes = man

counts <- exprMat[which(row.names(exprMat) %in% my_genes),]
dim(counts)

################ counts was used for Heatmap (using ComplexHeatmap package) of most variable TE loci as presented
################ Available in Singh et al., Cell Reports (2020) manuscript


#### Hunting for genes correlated with TE loci.

var.cor <- correlatePairs(sce, subset.row=my_genes)
ZG16.cor <- var.cor[grep("ZG16", var.cor$gene1),]
ZG16.cor <- ZG16.cor[order(ZG16.cor$FDR),]
head(ZG16.cor)
write.table(ZG16.cor, "HTR2B_significant_top_100_pairs.tsv", row.names=T, col.names=T, sep="\t", dec=".", quote=F) 




################################ Finally the pseudoTime analysis for all the celltypes using Monocle



Epithelial <- colon_integrated

data <- as(as.matrix(Epithelial@assays$RNA@data), 'sparseMatrix')

DefaultAssay(Epithelial) = "RNA"

Epithelial@meta.data$celltypes <- Idents(Epithelial)

pd <- new('AnnotatedDataFrame', data = Epithelial@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

pseudo_Epithelial <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 1,
                              expressionFamily = uninormal())
                              
                              
pseudo_Epithelial <- reduceDimension(pseudo_Epithelial,norm_method="none", 
                        reduction_method="DDRTree",
                        max_components=4,
                        scaling=TRUE,
                        verbose=TRUE,
                        pseudo_expr=1)
                        

expressed_genes <- VariableFeatures(object = Epithelial)

clustering_DEG_genes <- differentialGeneTest(pseudo_Epithelial[expressed_genes,], fullModelFormulaStr='~celltypes', cores=4)


sig_genes <- subset(clustering_DEG_genes, qval < 0.01)


sig_genes[,c("gene_short_name", "pval", "qval")]

my_pseudotime_de <- differentialGeneTest(pseudo_Epithelial,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1)

png("Trajectory_PseudoTime_all_colonocytes_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(pseudo_Epithelial, cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50,  color_by = "Pseudotime") +  scale_colour_gradient2(low = "lightblue", mid = "gold", high = "blue")
dev.off()


png("Trajectory_Celltypes_all_colonocytes_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(pseudo_Epithelial, color_by = "celltypes", cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50)+ 
 scale_color_manual(values=c("forestgreen", "orange", "yellowgreen", "cyan", "purple","red", "maroon"))
dev.off()


png("Trajectory_Patients_all_colonocytes_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(pseudo_Epithelial, color_by = "patients", cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50) + 
 scale_color_manual(values=c("tomato", "midnightblue"))
dev.off()


##################### Only for colonocytes
pbmc.sub1 <- subset(pbmc, idents = c("Stem", "SPIB+_cluster","CP", "IC", "CA1+_early_colonocyte", "CA1+_late_colonocyte", "CEACAM7+_colonocyte"))
data <- as(as.matrix(pbmc.sub1@assays$RNA@data), 'sparseMatrix')
DefaultAssay(pbmc.sub1) = "RNA"

pbmc.sub1@meta.data$celltypes <- Idents(pbmc.sub1)

pd <- new('AnnotatedDataFrame', data = pbmc.sub1@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

fd <- new('AnnotatedDataFrame', data = fData)

fd <- new('AnnotatedDataFrame', data = fData)

pd <- new('AnnotatedDataFrame', data = pbmc.sub1@meta.data)

colonocytes <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 1,
                              expressionFamily = uninormal())

colonocytes <- reduceDimension(colonocytes,norm_method="none", reduction_method="DDRTree",  max_components=4, scaling=TRUE,    verbose=TRUE,  pseudo_expr=1)

colonocytes <- orderCells(colonocytes)

my_pseudotime_de <- differentialGeneTest(colonocytes,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 1)
                                         
                                         

png("Trajectory_PseudoTime_all_colonocytes_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(colonocytes, cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50,  color_by = "Pseudotime") +  scale_colour_gradient2(low = "lightblue", mid = "gold", high = "blue")
dev.off()


png("Trajectory_Celltypes_all_colonocytes_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(colonocytes, color_by = "celltypes", cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50)+ 
 scale_color_manual(values=c("forestgreen", "orange", "yellowgreen", "cyan", "purple","red", "maroon"))
dev.off()


png("Trajectory_Patients_all_colonocytes_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(colonocytes, color_by = "patients", cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50) + 
 scale_color_manual(values=c("tomato", "midnightblue"))


#################### Only for Goblet cells

pbmc.sub1 <- subset(pbmc, idents = c("Stem", "secretory_progenitor", "immature_goblet", "goblet"))
data <- as(as.matrix(pbmc.sub1@assays$RNA@data), 'sparseMatrix')
DefaultAssay(pbmc.sub1) = "RNA"
pbmc.sub1@meta.data$celltypes <- Idents(pbmc.sub1)
pd <- new('AnnotatedDataFrame', data = pbmc.sub1@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))


fd <- new('AnnotatedDataFrame', data = fData)
goblet <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 1,
                              expressionFamily = uninormal())
                              
goblet <- reduceDimension(goblet,norm_method="none", reduction_method="DDRTree",  max_components=4, scaling=TRUE,    verbose=TRUE,  pseudo_expr=1)
goblet <- orderCells(goblet)

png("Trajectory_PseudoTime_all_Goblet_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(goblet, cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50,  color_by = "Pseudotime") +  scale_colour_gradient2(low = "lightblue", mid = "gold", high = "blue")
dev.off()



png("Trajectory_Celltypes_all_goblet_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(goblet, color_by = "celltypes", cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50)+ 
 scale_color_manual(values=c("forestgreen","yellowgreen", "purple", "maroon"))
dev.off()


png("Trajectory_Patients_all_goblet_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(goblet, color_by = "patients", cell_size=0.8, show_backbone=F,  show_branch_points = FALSE, show_tree = TRUE, theta=0.50) + 
 scale_color_manual(values=c("tomato", "midnightblue"))
dev.off()

saveRDS(colon_integrated, file = "colon_integrated_individual_data.rds")
saveRDS(colonocytes, file = "colonocytes_Monocle.rds")
saveRDS(goblet, file = "goblet_Monocle.rds")

