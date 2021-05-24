
###############

# ====================================
# Author: Manu Singh
# Date: 2021-05-18
# Title: colon-seq data analysis for NIBD and CD patients
# ====================================

# =====================================
# Enviornment variables
# =====================================
## Continuing the DATA ANALYSIS PART 
# =====================================
# Load libraries
# =====================================

library(Seurat)

library(reshape)
library(ggplot2)
library(cowplot)
library(reticulate)
library(dplyr)
library(data.table)


counts <- as.data.frame(colon_integrated[["RNA"]]@counts)

TE_fams <- counts[(grep("TELONG",rownames(counts))),]

TosumV1 <- (sapply(strsplit(rownames(TE_fams),"-DUP"),"[[",1))

TE_fams$Tosum <- TosumV1

TE_famsV2 <- data.table(TE_fams)

#TE_famsV2_sum <- TE_famsV2[,lapply(.SD,sum),by=Tosum, .SDcols=(colnames(TE_famsV2)[-(which(colnames(TE_famsV2)=="Tosum"))])]

# The above code is to sum the all reads mapping to TE loci, and assigning to a family

TE_famsV2_sum <- TE_famsV2[,lapply(.SD,mean),by=Tosum, .SDcols=(colnames(TE_famsV2)[-(which(colnames(TE_famsV2)=="Tosum"))])]


rownames(TE_famsV2_sum) <- (unlist(TE_famsV2_sum[,1]))
TE_famsV2_sumMatrix <- as.matrix(TE_famsV2_sum[,-1])
rownames(TE_famsV2_sumMatrix) <- (unlist(TE_famsV2_sum[,1]))


###### Now getting the genes dataframe
All_joined_countMatrixRemaining <- counts[(-grep("TELONG",rownames(counts))),]

All_joined_countMatrixFamilyWise <- (rbind(All_joined_countMatrixRemaining,as.matrix(TE_famsV2_sumMatrix)))

celltypeCondition <- colon_integrated$celltypeCondition
celltypes <- colon_integrated$celltypes
history(max.show=50))
history(max.show=50)


##########
colon_integrated2 <- CreateSeuratObject(counts = All_joined_countMatrixFamilyWise , min.cells = 3, min.features = 2000, project = "TE_Alevin", meta.data = metaDataOursv2)

colon_integrated2@metadata$celltype <- colon_integrated@metadata$celltype
colon_integrated2@metadata$orig.ident <- colon_integrated@metadata$orig.ident
colon_integrated2@metadata$patients <- colon_integrated@metadata$patients

### 
colon_integrated <- subset(colon_integrated, subset = nFeature_RNA > 3000 & nCount_RNA < 80000)

VlnPlot(colon_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
colon_integrated <- NormalizeData(colon_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
colon_integrated <- FindVariableFeatures(colon_integrated, selection.method = "vst", nfeatures = 2500)


###### From here on, we followed all steps same as in the DATA ANALYSIS codes, till the end.

colon.list <- SplitObject(colon_integrated, split.by = "patients")
colon.list <- lapply(X = colon.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2500)
})
features <- SelectIntegrationFeatures(object.list = colon.list)
immune.anchors <- FindIntegrationAnchors(object.list = colon.list, anchor.features = features)
colon_integrated <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(colon_integrated) <- "integrated"
colon_integrated <- ScaleData(colon_integrated, verbose = FALSE, features = all.genes)
colon_integrated <- RunPCA(colon_integrated, npcs = 30, verbose = FALSE)
colon_integrated <- RunUMAP(colon_integrated, reduction = "pca", dims = 1:30)
colon_integrated <- FindNeighbors(colon_integrated, reduction = "pca", dims = 1:30)

colon_integrated@misc$celltype_markers <- FindAllMarkers(colon_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.15)


my_genes <- (colon_integrated@misc$celltype_markers[grep("ERV|MLT|^L1|LTR|^MER", cmark$gene),]$gene)


png("Praveen_DotPlot_TE_families_markers.png",  width = 38.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
DotPlot(colon_integrated, features=my_genes, scale.min=0, scale.max=100, dot.min=0, idents=c("goblet","immature_goblet", "secretory_progenitor", "CA1+_early_colonocyte", "CA1+_late_colonocyte", "CEACAM7+_colonocyte","G2_M_G1_TA", "S_phase_TA",  "Stem", "SPIB+_cluster", "Tuft_cells")) + scale_colour_gradient2(low = "white", mid = "steelblue", high = "red2") + RotatedAxis()
dev.off()



pbmc.sub1 <- subset(pbmc, idents = c("Stem", "SPIB+_cluster","CP", "IC", "CA1+_early_colonocyte", "CA1+_late_colonocyte", "CEACAM7+_colonocyte"))
data <- as(as.matrix(pbmc.sub1@assays$RNA@data), 'sparseMatrix')
DefaultAssay(pbmc.sub1) = "RNA"

pbmc.sub1@meta.data$celltypes <- Idents(pbmc.sub1)

pd <- new('AnnotatedDataFrame', data = pbmc.sub1@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
library(monocle)
fd <- new('AnnotatedDataFrame', data = fData)


colonocytes <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 1,
                              expressionFamily = uninormal())
colonocytes <- reduceDimension(colonocytes,norm_method="none", reduction_method="DDRTree",  max_components=4, scaling=TRUE,    verbose=TRUE,  pseudo_expr=1)
colonocytes <- orderCells(colonocytes)


expressed_genes <- VariableFeatures(pbmc.sub1)

clustering_DEG_TEs <- differentialGeneTest(colonocytes[expressed_genes,], fullModelFormulaStr='~celltypes', cores=4)
sig_genes <- subset(clustering_DEG_TEs, qval < 0.01)

my_genes <- sig_genes[grep("ERV|MLT|^L1|LTR|^MER", sig_genes$gene),]$gene


my_pseudotime_de_G %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_gene

png("Trajectory_Transposable_elements_all_goblet_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(colonocytes, markers = my_genes, use_color_gradient = TRUE, markers_linear=F, cell_size=0.1, show_backbone=F, show_branch_points=F) +  scale_colour_gradient2(low = "lightblue", mid = "red", high = "black")
dev.off()

png("Trajectory_Young_Old_Transposable_elements_all_colonocytes_V2.png",  width = 18.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
plot_cell_trajectory(colonocytes, markers = c("L1PA2","L1PA3","L1PA14","L1PA10"), use_color_gradient = TRUE, markers_linear=F, cell_size=0.6, show_backbone=F, show_branch_points=F) +  scale_colour_gradient2(low = "gold", mid = "purple", high = "blue")
dev.off()

sig_gene_names <- row.names(head(sig_genes, n=39))

png("Trajectory_Heatmap_Transposable_elements_all_goblet_V2.png",  width = 15.8, height = 19.8, units = "cm", res = 600, pointsize = 12)
plot_pseudotime_heatmap(colonocytes[sig_gene_names,],  num_clusters = 3, hmcols = colorRampPalette(c("red","white","blue"))(256),    cores = 1,  show_rownames = F)
dev.off()


png("Trajectory_Heatmap_labelled_Transposable_elements_all_goblet_V2.png",  width = 15.8, height = 19.8, units = "cm", res = 600, pointsize = 12)
plot_pseudotime_heatmap(colonocytes[sig_gene_names,],  num_clusters = 3, hmcols = colorRampPalette(c("red","white","blue"))(256),    cores = 1,  show_rownames = T)
dev.off()




#################### The exact same code was then repeated for goblet cells

###### While analyzing the genes only, every step was repeated as it is (obviously except the importing the data)

colon_integrated <- CreateSeuratObject(counts = All_joined_countMatrixRemaining, min.cells = 5, min.features = 2000, project = "TE_Alevin", meta.data = metaDataOursv2)



