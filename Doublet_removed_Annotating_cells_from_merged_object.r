
###############

# ====================================
# Author: Manu Singh
# Date: 2021-05-18
# Title: colon-seq data analysis for NIBD and CD patients
# ====================================

# =====================================
# Enviornment variables
# =====================================
setwd("~/Collab/Praveen/seurat_obj")
getwd()
Sys.Date()
main_dir <- getwd()
date <- gsub("-", "", Sys.Date())

dir.create(file.path(main_dir, date), showWarnings = FALSE)
setwd(file.path(main_dir, date))

getwd()
options(future.globals.maxSize = 6096*1524^2 )
set.seed(2811)

# =====================================
# Load libraries
# =====================================

library(Seurat)
library(ggplot2)
library(cowplot)
library(reticulate)
library(dplyr)
##############################


######## This is an initial clustering without normalizing for batch effect, which would be corrected later
#### The purpose of this clustering to get the initial annotations which we could later replaced by the original annotations of cell clusters
##### we are now dealing with the merged datasets .see line 96 of previous code

colon <- readRDS("Merged_seurat_objects_Praveen_NIBD_CD_samples")

# from the results of 95th quantiles 
colon <- subset(colon, subset = nFeature_RNA > 3000 & nCount_RNA > 2.5e+04 & nCount_RNA < 1e+05)

# Removing the potential doublets


install.packages("remotes")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

library(DoubletFinder)
library(remotes)

pre_doublets <- colon

Idents(pre_doublets) <- pre_doublets$celltypes


colon_normalized <- SCTransform(pre_doublets)


# parameters as per developer's recommendations
# results were consistent when we played from first 10 PCs till 30 PCs

sweep.res.list_colon <- paramSweep_v3(colon_normalized, PCs = 1:10, sct = FALSE)
sweep.stats_colon <- summarizeSweep(sweep.res.list_colon, GT = FALSE)
bcmvn_colon <- find.pK(sweep.stats_colon)
dev.off()



annotations <- colon_normalized@meta.data$orig.ident
homotypic.prop <- modelHomotypic(annotations)   
nExp_poi <- round(0.075*nrow(colon_normalized@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
colon_normalized <- doubletFinder_v3(colon_normalized, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)

Idents(colon_normalized) <- colon_normalized@meta.data$DF.classifications_0.25_0.09_899

length(WhichCells(colon_normalized, ident="Doublet"))

# less than 1% cells were doublets

jj <- WhichCells(colon_normalized, ident="Doublet")

doublets_removed <- pre_doublets[,!colnames(pre_doublets) %in% jj]

colon_normalized_doublet_removed <- SCTransform(doublets_removed)


colon <- colon_normalized_doublet_removed



### Normalizing in usual way
colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 20000)

# Here we took 3500 (combination of TEs and genes).
#In the case of genes or TEs alone, we only take 2000 rows as most variable

colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 3500)

# Here we are regressing out the impact of ribosomal, nUMI, and cell cycle genes

colon <- ScaleData(colon, vars.to.regress = c("S.Score", "G2M.Score", "nUMI", "ribo.genes"), features = row.names(colon))

colon <- RunPCA(colon, features = VariableFeatures(object = colon))

colon <- JackStraw(colon, num.replicate = 100)
colon <- ScoreJackStraw(colon, dims = 1:20)

colon <- FindNeighbors(colon, dims = 1:20)

colon <- FindClusters(colon, resolution = 0.8)

colon@misc$seurat_markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


metadata <- as.data.frame(cbind(colnames(colon), as.character(colon@meta.data$cell_type),  as.character(colon@meta.data$patients),colon@meta.data$nCount_RNA,colon@meta.data$nFeature_RNA))

colnames(metadata) <- c("coloumns", "cell_types", "Patients", "Total_RNA", "Total_features")

write.table(metadata, "MetaData_NIBD_CD_scRNAseq_for_Matt.tsv", row.names=F, col.names=T, sep="\t", dec=".", quote=F)

matt_is <- <- read.delim("~/Downloads/NIBD-CD_epithelial-cells_metadata.csv", stringsAsFactors=F, sep=",", row.names=1)

row.names(matt_is) <- gsub("_.$", "", row.names(matt_is))

head(row.names(matt_is))

#"NIBD_AAACCCAAGATCCCGC" 
#"NIBD_AAACCCAAGCCTCCAG"
#"NIBD_AAACCCACAAGGCCTC" 
#"NIBD_AAACCCACATGTTCAG"
#"NIBD_AAAGGATCACTTCATT" 
#"NIBD_AAAGGATCAGCAGTCC"


unique(matt_is$cluster_named)

#"immature colonocyte"   
#"CA1+ late colonocyte"  
#"colonocyte progenitor"
#"CEACAM7+ colonocyte"   
#"S phase TA"            
#"immature goblet"      
#"G2-M-G1 TA"            
#"S phase TA-2"          
#"goblet"               
#"SPIB+ cells"
#"stem"
#"CA1+ early colonocyte"
#"Secretory progenitor"  
#"EEC"   


# Reannotating the cells based on the Epithelial cells classification used during the first submission of this paper

Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "immature colonocyte")))))       <- "IC"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "CA1+ late colonocyte")))))      <- "CA1+_late_colonocyte"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "colonocyte progenitor")))))     <- "CP"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "CEACAM7+ colonocyte")))))       <- "CEACAM7+_colonocyte"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "S phase TA")))))                <- "S_phase_TA"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "immature goblet")))))           <- "immature_goblet"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "G2-M-G1 TA")))))                <- "G2_M_G1_TA"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "S phase TA-2")))))              <- "S_phase_TA"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "goblet")))))                    <- "goblet"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "SPIB+ cells")))))               <- "SPIB+_cells"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "stem")))))                      <- "Stem"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "CA1+ early colonocyte")))))     <- "CA1+_early_colonocyte"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "Secretory progenitor")))))      <- "secretory_progenitor"
Idents(object = colon, cells = WhichCells(subset(colon, cells=row.names(subset(matt_is, cluster_named == "EEC")))))                       <- "EEC"


#Here, checking the status of seurat clusters after modifying the annotations.
DimPlot(colon)
dev.off()

#This visualization suggested that there are slightly higher number of cells in this analysis, 
#so we are annotating the rest of cells carefully and only if they share the same clusters as classified above 

Idents(object = colon, cells = WhichCells(colon, idents=18)) <- "Tuft_cells"
Idents(object = colon, cells = WhichCells(colon, idents=c(0,7))) <- "CA1+_late_colonocyte"
Idents(object = colon, cells = WhichCells(colon, idents=1)) <- "G2_M_G1_TA"
Idents(object = colon, cells = WhichCells(colon, idents=c(2,9))) <- "CEACAM7+_colonocyte"
Idents(object = colon, cells = WhichCells(colon, idents=3)) <- "secretory_progenitor"
Idents(object = colon, cells = WhichCells(colon, idents=4)) <- "S_phase_TA"
Idents(object = colon, cells = WhichCells(colon, idents=c(5,10))) <- "CA1+_early_colonocyte"
Idents(object = colon, cells = WhichCells(colon, idents=6)) <- "goblet"
Idents(object = colon, cells = WhichCells(colon, idents=8)) <- "immature_goblet"
Idents(object = colon, cells = WhichCells(colon, idents=11)) <- "SPIB+_cells"
Idents(object = colon, cells = WhichCells(colon, idents=12)) <- "CP"
Idents(object = colon, cells = WhichCells(colon, idents=13)) <- "IC"
Idents(object = colon, cells = WhichCells(colon, idents=14)) <- "Stem"
Idents(object = colon, cells = WhichCells(colon, idents=17)) <- "EEC"


# Note: cluster number 15 and 16 were immune cells or macrophages so were not annotated, and all the unannotated cells will be removed in next step

unique(Idents(colon))
#0 CA1+_late_colonocyte
#1 G2_M_G1_TA
#2 CEACAM7+_colonocyte
#3 secretory progenitor
#4 S_phase_TA, S_phase_TA-2
#5 CA1+_early_colonocyte
#6 goblet
#7 CA1+_late_colonocyte
#8 immature goblet
#9 CEACAM7+_colonocyte
#10 CA1+_early_colonocyte


#13 SPIB+_cells
#14 stem
#17 EEC
#### 18 Tuft cells



DefaultAssay(colon) <- "RNA"

colon <- subset(colon, idents = c("secretory_progenitor","goblet","CP","immature_goblet","CA1+_late_colonocyte","IC","G2_M_G1_TA","CA1+_late_colonocyte","SPIB+_cluster","Stem","S_phase_TA","Tuft_cells","CEACAM7+_colonocyte","EEC","CA1+_early_colonocyte"))
DimPlot(colon)
dev.off()

obj.list <- SplitObject(pbmc, split.by = "orig.ident")

### This is saving the individual samples, Note that these merged datasets were not corrected for batch effect
for(i in 1:length(obj.list)) {
colon_file <- obj.list[i]
colon_file_name <- basename(colon_file) %>% gsub(".rds","",.)
message(paste0("Processing Sample ",j," of ",length(Praveen_files)," - ", colon_file_name))

saveRDS(colon_file, file = paste0("Annotated_", colon_file_name, ".rds"))
}

saveRDS(colon, "Merged_Annotated_seurat_objects_Praveen_NIBD_CD_samples.rds")




