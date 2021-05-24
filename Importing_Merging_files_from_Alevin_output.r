###############

# ====================================
# Author: Manu Singh
# Date: 2021-05-18
# Title: colon-seq data analysis for NIBD and CD patients
# ====================================

# =====================================
# Enviornment variables
# =====================================
setwd("~/Collab/Praveen")
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
library(ggpubr)
library(RCurl)
library(gplots)
library(BSDA)
library(reshape)
library(sctransform)
library(tximport, lib.loc= "/home/ms3559/R/x86_64-pc-linux-gnu-library/3.5")
sessionInfo()
# =====================================
# Load all the objects 
# These objects were processed with Seurat v3 on 2020/05
# =====================================



files <- file.path("/workdir/Manu/Data/colon_seq/Praveen/CD_AC/CD_1/alevin/quants_mat.gz")

pbmc <- CreateSeuratObject(counts = txi , min.cells = 3, min.features = 200, project = "CD_3")



Praveen_files <- list.files()
for(j in 1:length(Praveen_files)) {
colon_file <- Praveen_files[j]
colon_file_name <- basename(colon_file) %>% gsub("_quants_mat.gz","",.)
message(paste0("Processing Sample ",j," of ",length(Praveen_files)," - ", colon_file_name))

#Read Alevin for 10X outputs

colon_dge <- tximport(colon_file, type="alevin")

colon_objs <- CreateSeuratObject(counts = colon_dge$counts, project = colon_file_name, min.cells = 5, min.features = 2500)

#mito reads
#normalize
#scale


colon_objs <- SCTransform(colon_objs, vars.to.regress = c("percent.mt"), verbose = T, variable.features.n = 3000)

saveRDS(colon_objs, file = paste0("/seurat_obj/Seu_", colon_file_name, ".rds"))
}


############################


colon <- merge(x=colon_objs[[1]], y=c(colon_objs[[2]],colon_objs[[3]],colon_objs[[4]],colon_objs[[5]],colon_objs[[6]],colon_objs[[7]]), add.cell.ids = c("CD_1", "NIBD_1", "NIBD_2","NIBD_3","NIBD_4", "CD_2", "CD_3"), project="Praveen");
rm(colon_objs); # save some space
str(colon@meta.data) # examine the structure of the Seurat object meta data
saveRDS(colon, file = sprintf("%s/MergedSeuratObject.rds", outdir));

mito.genes <- grep(pattern = "^MT-", x = rownames(x = colon), value = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = colon, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = colon, slot = 'counts'));
colon[['percent.mito']] <- percent.mito;

ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = colon), value = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = colon, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = colon, slot = 'counts'));
colon[['percent.ribo']] <- percent.ribo;


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
colon <- CellCycleScoring(colon, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

saveRDS(colon, file = "Merged_seurat_objects_Praveen_NIBD_CD_samples", colon_file_name, ".rds"))


# These following commands are to check if data needs to be subset , we will store it for now and use it later
min <- min(colon@meta.data$nFeature_RNA);
m <- median(colon@meta.data$nFeature_RNA)
max <- max(colon@meta.data$nFeature_RNA)    
s <- sd(colon@meta.data$nFeature_RNA)
min1 <- min(colon@meta.data$nCount_RNA)
max1 <- max(colon@meta.data$nCount_RNA)
m1 <- mean(colon@meta.data$nCount_RNA)
s1 <- sd(colon@meta.data$nCount_RNA)
quant_95_counts <- quantile(colon@meta.data$nCount_RNA, 0.95) # calculate value in the 95th percentile
print(paste("Feature stats:",min,m,max,s));
print(paste("UMI stats:",min1,m1,max1,s1,Count93));

