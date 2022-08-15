########################Immune Repertoire Preprocess############################
# rm(list = ls())
# 
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/ImmuneRepertoire/")
# 
# PreprocessVDJ <- function(contig = NULL, type = NULL, export = NULL) {
#   contig_tab <- read.csv(contig)
# 
#   # Select valid consensus
#   contig_tab <- contig_tab[contig_tab$productive %in% c("True", "true") &
#                              contig_tab$full_length %in% c("True", "true") &
#                              contig_tab$high_confidence %in% c("True", "true") &
#                              contig_tab$is_cell %in% c("True", "true"),]
# 
#   # Eliminate Multi/None chain type
#   if (type == "BCR") {
#     contig_tab <- contig_tab[contig_tab$chain %in% c("IGH", "IGK", "IGL"),]
# 
#     # For convenience, IGK and IGL are defined as IGL for data processing
#     # This does not alter the exact names of the corresponding v/j segments
#     contig_tab$chain[contig_tab$chain %in% c("IGK", "IGL")] <- "IGL"
#   }
# 
#   if (type == "TCR") {
#     contig_tab <- contig_tab[contig_tab$chain %in% c("TRA", "TRB"),]
#   }
#   
#   # Eliminate orphan chain (a sole TRA/TRB/IGH/IGL/IGK chain)
#   unpaired_barcode <- names(table(contig_tab$barcode))[table(contig_tab$barcode) == 1]
#   contig_tab <- contig_tab[!contig_tab$barcode %in% unpaired_barcode,]
#   
#   print("unpaired_barcode:")
#   print(unpaired_barcode)
#   
#   # if more than one TCR/BCR or IGH/IGL&IGK paired chains were identified in one cell,
#   # we only kept the dominant paired chain (supported by the largest number of UMIs) for it
#   multi_chain_barcode <- names(table(contig_tab$barcode))[table(contig_tab$barcode) > 2]
#   eliminate <- c()
#   for (i in multi_chain_barcode) {
#     multi_chain_contig <- contig_tab[contig_tab$barcode == i, ]
#     retain_index <- tapply(multi_chain_contig$umis, INDEX = multi_chain_contig$chain, FUN = max)
# 
#     # define which row should be retained through iteration
#     retain <- c()
# 
#     for (j in 1:length(names(retain_index))) {
#       # retain_index is a series of number (umi) with names (specific chain)
#       # names(retain_index)[j] corresponds to specific chain type that is duplicated in a cell, e.g. "TRA" or "IGH"
#       # retain_index[j] corresponds to the chain supported by the maxim umi
#       retainlogic <- multi_chain_contig$chain == names(retain_index)[j] & multi_chain_contig$umis == retain_index[j]
#       retain <- c(retain, rownames(multi_chain_contig)[retainlogic])
#     }
# 
#     # define which row should be eliminated through iteration
#     eliminate <- c(eliminate, rownames(multi_chain_contig)[!rownames(multi_chain_contig) %in% retain])
#   }
# 
#   contig_tab <- contig_tab[!rownames(contig_tab) %in% eliminate,]
# 
#   # On rare occasions duplicated chains were supported by same umi, we retain those
#   # supported by more reads
#   multi_chain_barcode <- names(table(contig_tab$barcode))[table(contig_tab$barcode) > 2]
# 
#   eliminate <- c()
#   for (i in multi_chain_barcode) {
#     multi_chain_contig <- contig_tab[contig_tab$barcode == i, ]
#     retain_index <- tapply(multi_chain_contig$reads, INDEX = multi_chain_contig$chain, FUN = max)
# 
#     retain <- c()
# 
#     for (j in 1:length(names(retain_index))) {
#       # names(retain_index)[j] correspond to specific chain type that is still duplicated in a cell, e.g. "TRA" or "IGH"
#       # retain_index[j] correspond to the chain supported by the maximum reads
#       retainlogic <- multi_chain_contig$chain == names(retain_index)[j] & multi_chain_contig$umis == retain_index[j]
#       retain <- c(retain, rownames(multi_chain_contig)[retainlogic])
#     }
# 
#     # define which row should be eliminated through iteration
#     eliminate <- c(eliminate, rownames(multi_chain_contig)[!rownames(multi_chain_contig) %in% retain])
#   }
#   contig_tab <-  contig_tab[!rownames(contig_tab) %in% eliminate,]
# 
#   filename <- paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/ImmuneRepertoire/ProcessedData/", type, "/", export, sep = "")
#   write.csv(contig_tab, file = filename, quote = F, row.names = F)
# }
# 
# PreprocessVDJ(contig = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/VDJ/Chromium_X/10k_PBMC_5pv2_nextgem_Chromium_X_10k_PBMC_5pv2_nextgem_Chromium_X_vdj_b_filtered_contig_annotations.csv", type = "BCR", export = "Control3.csv")
# PreprocessVDJ(contig = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/VDJ/SC5v2_humanPBMCs_5Kcells/SC5v2_humanPBMCs_5Kcells_Connect_single_channel_SC5v2_humanPBMCs_5Kcells_Connect_single_channel_vdj_b_filtered_contig_annotations.csv", type = "BCR", export = "Control4.csv")
# PreprocessVDJ(contig = "/mnt/TRAINING/Zhongjiacheng/XRY/PBMC-XRY-BCR/filtered_contig_annotations.csv", type = "BCR", export = "Patient.csv")
# 
# PreprocessVDJ(contig = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/VDJ/Chromium_X/10k_PBMC_5pv2_nextgem_Chromium_X_10k_PBMC_5pv2_nextgem_Chromium_X_vdj_t_filtered_contig_annotations.csv", type = "TCR", export = "Control3.csv")
# PreprocessVDJ(contig = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/VDJ/SC5v2_humanPBMCs_5Kcells/SC5v2_humanPBMCs_5Kcells_Connect_single_channel_SC5v2_humanPBMCs_5Kcells_Connect_single_channel_vdj_t_filtered_contig_annotations.csv", type = "TCR", export = "Control4.csv")
# PreprocessVDJ(contig = "/mnt/TRAINING/Zhongjiacheng/XRY/PBMC-XRY-TCR/filtered_contig_annotations.csv", type = "TCR", export = "Patient.csv")

#######################scRNAseq processing for integration######################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/QualityControl/")
# 
# library(Seurat)
# library(DoubletFinder)
# library(dplyr)
# library(ggplot2)
# 
# SeuratPatientPBMC <- CreateSeuratObject(counts =  Read10X(data.dir = "/mnt/TRAINING/Zhongjiacheng/XRY/PBMC-XRY-scRNA/filtered_feature_bc_matrix/"), min.cells = 3)
# SeuratControl1PBMC <- CreateSeuratObject(counts = Read10X(data.dir = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/ATAC+RNA/10k_PBMC_Multiome_nextgem_Chromium_X/filtered_feature_bc_matrix/")$`Gene Expression`, min.cells = 3)
# SeuratControl2PBMC <- CreateSeuratObject(counts = Read10X(data.dir = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/ATAC+RNA/pbmc_granulocyte_sorted_10k/filtered_feature_bc_matrix/")$`Gene Expression`, min.cells = 3)
# SeuratControl3PBMC <- CreateSeuratObject(counts = Read10X(data.dir = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/VDJ/Chromium_X/sample_feature_bc_matrix/"), min.cells = 3)
# SeuratControl4PBMC <- CreateSeuratObject(counts = Read10X(data.dir = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/VDJ/SC5v2_humanPBMCs_5Kcells/sample_feature_bc_matrix/"), min.cells = 3)
# SeuratControl5PBMC <- CreateSeuratObject(counts = Read10X(data.dir = "/mnt/TRAINING/Zhongjiacheng/10xgenomics/RNA/10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix/filtered_feature_bc_matrix/"), min.cells = 3)
# 
# SeuratPatientBALF <- CreateSeuratObject(counts = Read10X(data.dir = "/mnt/TRAINING/Zhongjiacheng/XRY/BALF-XRY-scRNA/filtered_feature_bc_matrix/"), min.cells = 3)
# SeuratControl6BALF <- CreateSeuratObject(counts = read.csv("/mnt/TRAINING/Zhongjiacheng/Controls/BALF/GSM4593890_sample_3_UMI_counts.csv.gz", row.names = 1), min.cells = 3)
# SeuratControl7BALF <- CreateSeuratObject(counts = read.csv("/mnt/TRAINING/Zhongjiacheng/Controls/BALF/GSM4593891_sample_4_UMI_counts.csv.gz", row.names = 1), min.cells = 3)
# SeuratControl8BALF <- CreateSeuratObject(counts = read.csv("/mnt/TRAINING/Zhongjiacheng/Controls/BALF/GSM4593893_sample_6_UMI_counts.csv.gz", row.names = 1), min.cells = 3)
# SeuratControl9BALF <- CreateSeuratObject(counts = read.csv("/mnt/TRAINING/Zhongjiacheng/Controls/BALF/GSM4593896_sample_9_UMI_counts.csv.gz", row.names = 1), min.cells = 3)
# SeuratControl10BALF <- CreateSeuratObject(counts = read.csv("/mnt/TRAINING/Zhongjiacheng/Controls/BALF/GSM4593897_sample_10_UMI_counts.csv.gz", row.names = 1), min.cells = 3)
# 
# PreprocessSeurat <- function(dataset = NULL) {
#   seu <- get(dataset)
#   seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
# 
#   nCount_RNA_mean <- mean(unlist(seu[["nCount_RNA"]]))
#   nCount_RNA_sd <- sd(unlist(seu[["nCount_RNA"]]))
#   nCount_RNA_threshold <- nCount_RNA_mean + 3*nCount_RNA_sd
# 
#   for(i in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
#     Vln <- as.data.frame(cbind(rep(i, nrow(seu@meta.data)), seu@meta.data[i]), stringsAsFactors = F)
#     colnames(Vln) <- c("Identity", "Value")
#     Vln$Value <- as.numeric(Vln$Value )
#     tmp <- ggplot(Vln, aes(x = Identity, y = Value, fill = "red")) +
#       geom_violin(show.legend = F) +
#       theme_light() +
#       theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#     assign(paste(i, "_Vlnplot", sep = ""), tmp)
#   }
# 
#   nFeature_RNA_Vlnplot <- nFeature_RNA_Vlnplot + geom_hline(yintercept = 500, color = "blue", size = 2)
# 
#   nCount_RNA_Vlnplot <- nCount_RNA_Vlnplot + geom_hline(yintercept = nCount_RNA_threshold, linetype = "dashed",
#                                                         color = "red", size = 2)
# 
#   percent.mt_Vlnplot <- percent.mt_Vlnplot + geom_hline(yintercept = 20, linetype = "dashed",
#                                                         color = "red", size = 2)
# 
#   svg(filename = paste(dataset, "/",dataset, "Violin.svg", sep = ""), width = 5, height = 3)
#   print(nFeature_RNA_Vlnplot + nCount_RNA_Vlnplot + percent.mt_Vlnplot)
#   dev.off()
# 
#   pc.use <- 30
# 
#   seu <- subset(seu, subset = nFeature_RNA > 500 & nCount_RNA < nCount_RNA_threshold & percent.mt < 20) %>%
#     NormalizeData() %>%
#     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#     ScaleData() %>%
#     RunPCA() %>%
#     FindNeighbors(dims = 1:pc.use) %>%
#     RunUMAP(dims = 1:pc.use) %>%
#     FindClusters()
# 
#   sweep.res.list <- paramSweep_v3(seu, PCs = 1:pc.use, sct = F)
#   sweep.stats <- summarizeSweep(sweep.res.list, GT = F)
#   bcmvn <- find.pK(sweep.stats)
#   pk_optimal <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
# 
#   if (500 < ncol(seu) & ncol(seu) <= 1000) {
#     DoubletFormationRate <- 0.004
#   }
# 
#   if (1000 < ncol(seu) & ncol(seu) <= 2000) {
#     DoubletFormationRate <- 0.008
#   }
# 
#   if (2000 < ncol(seu) & ncol(seu) <= 3000) {
#     DoubletFormationRate <- 0.016
#   }
# 
#   if (3000 < ncol(seu) & ncol(seu) <= 4000) {
#     DoubletFormationRate <- 0.023
#   }
# 
#   if (4000 < ncol(seu) & ncol(seu) <= 5000) {
#     DoubletFormationRate <- 0.031
#   }
# 
#   if (5000 < ncol(seu) & ncol(seu) <= 6000) {
#     DoubletFormationRate <- 0.039
#   }
# 
#   if (6000 < ncol(seu) & ncol(seu) <=  7000) {
#     DoubletFormationRate <- 0.046
#   }
# 
#   if (7000 < ncol(seu) & ncol(seu) <=  8000) {
#     DoubletFormationRate <- 0.054
#   }
# 
#   if (8000 < ncol(seu) & ncol(seu) <=  9000) {
#     DoubletFormationRate <- 0.061
#   }
# 
#   if (9000 < ncol(seu) & ncol(seu) <  10000) {
#     DoubletFormationRate <- 0.069
#   }
# 
#   if (ncol(seu) > 10000) {
#     DoubletFormationRate <- 0.076
#   }
# 
#   # Exclude potential homotypic doublets (which are difficult to identify) from the inference
#   # nExp_poi corresponds to the total number of doublet predictions produced
#   # nExp_poi.adj is nExp_poi adjusted according to the estimated proportion of homotypic doublets
# 
#   homotypic.prop <- modelHomotypic(seu$seurat_clusters)
#   nExp_poi <- round(DoubletFormationRate*ncol(seu))
#   nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
# 
#   seu <-  doubletFinder_v3(seu, PCs = 1:pc.use, pN = 0.25, pK = pk_optimal, nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
# 
#   seu@meta.data <- seu@meta.data[, -grep(colnames(seu@meta.data), pattern = "pANN_")]
#   colnames(seu@meta.data)[ncol(seu@meta.data)] <- "DFInference"
# 
#   svg(filename = paste(dataset, "/",dataset, "doublet_UMAP.svg", sep = ""), width = 4, height = 3)
#   print(DimPlot(seu,reduction = "umap", group.by = "DFInference") + theme(plot.title = element_blank()))
#   dev.off()
# 
#   svg(filename = paste(dataset, "/",dataset, "separate_UMAP.svg", sep = ""), width = 6.5, height = 3)
#   print(DimPlot(seu,reduction = "umap", split.by = "DFInference") + theme(plot.title = element_blank()))
#   dev.off()
# 
#   # subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 200 & percent.mt < 20)
# 
#   seu <- seu[, which(seu@meta.data$DFInference == "Singlet")]
#   seu@meta.data$orig.ident <- strsplit(dataset, split = "Seurat")[[1]][2]
#   return(seu)
# }
# 
# SeuratPatientPBMC <- PreprocessSeurat("SeuratPatientPBMC")
# SeuratControl1PBMC <- PreprocessSeurat("SeuratControl1PBMC")
# SeuratControl2PBMC <- PreprocessSeurat("SeuratControl2PBMC")
# SeuratControl3PBMC <- PreprocessSeurat("SeuratControl3PBMC")
# SeuratControl4PBMC <- PreprocessSeurat("SeuratControl4PBMC")
# SeuratControl5PBMC <- PreprocessSeurat("SeuratControl5PBMC")
# 
# save(SeuratPatientPBMC, SeuratControl1PBMC, SeuratControl2PBMC, SeuratControl3PBMC,
#      SeuratControl4PBMC, SeuratControl5PBMC, file = "PreprocessedPBMCSamples.Rdata")
# 
# SeuratPatientBALF <-  PreprocessSeurat("SeuratPatientBALF")
# SeuratControl6BALF <- PreprocessSeurat("SeuratControl6BALF")
# SeuratControl7BALF <- PreprocessSeurat("SeuratControl7BALF")
# SeuratControl8BALF <- PreprocessSeurat("SeuratControl8BALF")
# SeuratControl9BALF <- PreprocessSeurat("SeuratControl9BALF")
# SeuratControl10BALF <- PreprocessSeurat("SeuratControl10BALF")
# 
# save(SeuratPatientBALF, SeuratControl6BALF, SeuratControl7BALF, SeuratControl8BALF,
#      SeuratControl9BALF, SeuratControl10BALF, file = "PreprocessedBALFSamples.Rdata")

##########################MergeDatasets-PBMC####################################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/MergeDatasets/")
# 
# library(ggplot2)
# library(Seurat)
# library(dplyr)
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/QualityControl/PreprocessedPBMCSamples.Rdata")
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/cycle.rda")
# 
# pc.use <- 30
# 
# RawSeuratPBMC <- merge(SeuratPatientPBMC,
#                        y = c(SeuratControl1PBMC, SeuratControl2PBMC,
#                              SeuratControl3PBMC, SeuratControl4PBMC,
#                              SeuratControl5PBMC)
# )
# 
# # Cell cycle scoring
# 
# RawSeuratPBMC <- NormalizeData(RawSeuratPBMC) %>%
#   CellCycleScoring(g2m.features = g2m_genes, s.features = s_genes) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# RawSeuratPBMC <- ScaleData(RawSeuratPBMC)
# RawSeuratPBMC <- RunPCA(RawSeuratPBMC)
# 
# # We do not see large differences due to cell cycle phase. Based on this plot, we would not regress out the variation due to cell cycle.
# 
# svg(filename = paste("Cell_cycle_scoring.svg", sep = ""), width = 10, height = 4)
# print(DimPlot(RawSeuratPBMC,
#               reduction = "pca",
#               group.by = "Phase",
#               split.by = "Phase") +
#         theme(plot.title = element_blank()))
# dev.off()
# 
# RawSeuratPBMC <-RunUMAP(RawSeuratPBMC, dims = 1:pc.use)
# 
# svg(filename = paste("SeuratPBMC_Before_batch_correction.svg", sep = ""), width = 6, height = 5)
# print(DimPlot(RawSeuratPBMC, reduction = "umap", group.by = "orig.ident",
#               cols = c("#6EBDF8",
#                        "#6BC567",
#                        "#E680C4",
#                        "#E1B385",
#                        "#A687DF",
#                        "#F98077")) +
#         theme(plot.title = element_blank()))
# dev.off()
# 
# DefaultAssay(RawSeuratPBMC) <- "RNA"
# RawSeuratPBMC <- DietSeurat(RawSeuratPBMC, assays = "RNA")
# SeuratPBMC.list <- SplitObject(RawSeuratPBMC, split.by = "orig.ident")
# SeuratPBMC.list <- lapply(SeuratPBMC.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = T) # This step was intended to display the original normalized counts in subsequent visualization (SCTransform does not normalize the original RNA assay)
#   x <- SCTransform(x, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = T)
# })
# 
# integ_features <- SelectIntegrationFeatures(object.list = SeuratPBMC.list,
#                                             nfeatures = 3000)
# 
# SeuratPBMC.list <- PrepSCTIntegration(object.list = SeuratPBMC.list,
#                                    anchor.features = integ_features)
# 
# integ_anchors <- FindIntegrationAnchors(object.list = SeuratPBMC.list,
#                                         normalization.method = "SCT",
#                                         anchor.features = integ_features)
# 
# SeuratPBMCIntegrated <- IntegrateData(anchorset = integ_anchors,
#                                    normalization.method = "SCT")
# 
# SeuratPBMCIntegrated <- RunPCA(object = SeuratPBMCIntegrated)
# 
# save(SeuratPBMCIntegrated, file = "SeuratPBMCIntegrated.Rdata")

###################PBMC Clustering and Main Cell Type Annotation################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/ReannotationAndElimination/")
# 
# library(ggplot2)
# library(Seurat)
# library(dplyr)
# suppressMessages(library(ArchR))
# 
# source("/mnt/TRAINING/Zhongjiacheng/Customized Codes.R")
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/MergeDatasets/SeuratPBMCIntegrated.Rdata")
# 
# ExactMatch <- function(x) {
#   paste("\\b", x, "\\b", sep = "")
# }
# 
# svg(filename = paste("Elbow_Plot.svg", sep = ""), width = 7, height = 5)
# ElbowPlot(SeuratPBMCIntegrated, ndims = 50) # to determine number of dimensions for clustering
# dev.off()
# 
# pc.use <- 37
# 
# SeuratPBMCIntegrated <- FindNeighbors(SeuratPBMCIntegrated, dims = 1:pc.use, verbose = T)
# SeuratPBMCIntegrated <- FindClusters(SeuratPBMCIntegrated, resolution = 1.2)
# SeuratPBMCIntegrated <- RunUMAP(SeuratPBMCIntegrated, dims = 1:pc.use)
# 
# svg(filename = paste("SeuratPBMC_After_batch_correction.svg", sep = ""), width = 6, height = 5)
# print(DimPlot(SeuratPBMCIntegrated, reduction = "umap",
#               group.by = "orig.ident", cols = c("#6EBDF8",
#                                                 "#6BC567",
#                                                 "#E680C4",
#                                                 "#E1B385",
#                                                 "#A687DF",
#                                                 "#F98077")) +
#         theme(plot.title = element_blank()))
# dev.off()
# 
# # Switch back to RNA assay to perform down-stream differential analyses
# DefaultAssay(SeuratPBMCIntegrated) <- "RNA"
# 
# Reanno.cluster.ids <- levels(SeuratPBMCIntegrated@active.ident)
# 
# Reanno.cluster.ids[c(0, 1, 4, 6, 8, 9, 10, 11, 12, 15, 16, 17, 18, 20, 21, 22, 25)+1] <- "TNK"
# Reanno.cluster.ids[c(2, 3, 5, 14, 19, 23, 26)+1] <- "Myeloid"
# Reanno.cluster.ids[c(7, 13, 29)+1] <- "B"
# Reanno.cluster.ids[c(24)+1] <- "pDC"
# Reanno.cluster.ids[c(27)+1] <- "Platelets"
# Reanno.cluster.ids[c(28)+1] <- "Stem"
# 
# names(Reanno.cluster.ids) <- levels(SeuratPBMCIntegrated@active.ident)
# SeuratPBMCIntegrated <- RenameIdents(SeuratPBMCIntegrated, Reanno.cluster.ids)
# 
# DimPlot(SeuratPBMCIntegrated, label = T)
# 
# SeuratPBMCIntegrated$Group <- ""
# SeuratPBMCIntegrated$Group[grep(SeuratPBMCIntegrated$orig.ident, pattern = "Control")] <- "Control"
# SeuratPBMCIntegrated$Group[grep(SeuratPBMCIntegrated$orig.ident, pattern = "Patient")] <- "Patient"
# 
# SeuratPBMCIntegrated$MainCellType <-as.character(SeuratPBMCIntegrated@active.ident)
# 
# SeuratPBMCIntegrated$GroupMainCellType <- ""
# SeuratPBMCIntegrated$GroupMainCellType <- paste(SeuratPBMCIntegrated$Group, SeuratPBMCIntegrated$MainCellType)
# 
# DEG <- list()
# 
# for (CellType in c("Myeloid", "TNK", "B")) {
#   FirstGroup <- paste("Patient", CellType, sep = " ")
#   SecondGroup <- paste("Control", CellType, sep = " ")
# 
#   DEGPv1To2 <- FindMarkers(object = SeuratPBMCIntegrated[,SeuratPBMCIntegrated$orig.ident %in% c("PatientPBMC",
#                                                                        "Control1PBMC",
#                                                                        "Control2PBMC")],
#                            ident.1 = FirstGroup,
#                            ident.2 = SecondGroup,
#                            group.by = "GroupMainCellType",
#                            logfc.threshold = 0.2,
#                            test.use = "MAST")
# 
#   DEGPv3To5 <- FindMarkers(object = SeuratPBMCIntegrated[,SeuratPBMCIntegrated$orig.ident %in% c("PatientPBMC",
#                                                                        "Control3PBMC",
#                                                                        "Control4PBMC",
#                                                                        "Control5PBMC")],
#                            ident.1 = FirstGroup,
#                            ident.2 = SecondGroup,
#                            group.by = "GroupMainCellType",
#                            logfc.threshold = 0.2,
#                            test.use = "MAST")
# 
#   CommonUp <- intersect(rownames(DEGPv1To2[DEGPv1To2$avg_log2FC > 0,]),
#                         rownames(DEGPv3To5[DEGPv3To5$avg_log2FC > 0,]))
#   CommonDown <- intersect(rownames(DEGPv1To2[DEGPv1To2$avg_log2FC < 0,]),
#                         rownames(DEGPv3To5[DEGPv3To5$avg_log2FC < 0,]))
# 
#   DEG[[CellType]] <- DEGPv1To2[rownames(DEGPv1To2) %in% c(CommonUp, CommonDown),]
# 
#   Retain <- DEG[[CellType]]$p_val_adj < 0.05
# 
#   DEG[[CellType]] <- DEG[[CellType]][Retain, ]
# }
# 
# save(DEG, file = "MainCellTypeDEG.Rdata")
# 
# PotentialBioMarkers <- c("IL13",
# "IL6",
# "IFNG",
# 
# "CXCL3", # GROγ
# "CCL5", # RANTES
# "LEP", # Leptin
# 
# "TNF",
# "IL4",
# "IL10",
# "IL12",
# "IL7",
# "IL18",
# "IL17A",
# "IL1B",
# "IL8",
# "CXCL10", # IP-10
# "CXCL12" # DF-1alpha
# )
# 
# ImmuneInhibitors <- c("LAG3", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "SLAMF6", "BTLA", "TIM3")
# 
# ImmuneMarkers <- grep(paste(c(ImmuneInhibitors, "^CCL", "^CXCL", "^CXCR", "^CCR",
#                               "^TNF", "^TGF", "^IFN", "^GZM", "^STAT", "^PRF", "^DEFB", "^DEFA",
#                               "^LYZ", "S100A8", "S100A9", "^CAMP", "^LCN2", "^LTF",
#                               paste("^IL", c(1:50), sep = "")), collapse = "|"),
#                       rownames(SeuratPBMCIntegrated), value = T)
# 
# DEG$TNK[rownames(DEG$TNK) %in% c(PotentialBioMarkers, ImmuneMarkers),]
# DEG$B[rownames(DEG$B) %in% c(PotentialBioMarkers, ImmuneMarkers),]
# DEG$Myeloid[rownames(DEG$Myeloid) %in% c(PotentialBioMarkers, ImmuneMarkers),]
# 
# DiffImmuneMarkers <- Reduce(union, list(rownames(DEG$TNK[rownames(DEG$TNK) %in% c(PotentialBioMarkers, ImmuneMarkers),]),
#                        rownames(DEG$B[rownames(DEG$B) %in% c(PotentialBioMarkers, ImmuneMarkers),]),
#                        rownames(DEG$Myeloid[rownames(DEG$Myeloid) %in% c(PotentialBioMarkers, ImmuneMarkers),])))
# 
# levels(DiffImmuneMarkers) <- rev(c("LAG3", "TIGIT",
#                                    "CCL4", "CCL4L2", "CCL5",
#                                    "CXCR3", "CXCR4", "CXCR5", "CCR7",
#                                    "TNFRSF1A", "TNFRSF1B", "TNFSF10", "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF25",
#                                    "TNFAIP3", "TNFAIP2", "TNFAIP8",
#                                    "IFNAR2", "IFNGR1",
#                                    "IL1B", "IL15", "IL16", "IL32", "TGFB1",
#                                    "IL2RG", "IL4R", "IL7R", "IL6R", "IL13RA1", "IL17RA", "TGFBR2",
#                                    "IL6ST",
#                                    "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
#                                    "GZMA", "GZMB", "GZMH", "GZMM", "PRF1", "TNF",
#                                    "S100A8", "S100A9", "LYZ"
#                                    ))
# 
# AnalyseSubTypes <- c("Myeloid", "TNK", "B")
# levels(AnalyseSubTypes) <- c("Myeloid", "TNK", "B")
# 
# svg(filename = "MainCellTypesBubble.svg", width = 6, height = 8)
# BubblePlot(SeuratPBMCIntegrated[,SeuratPBMCIntegrated$orig.ident %in% c("Control1PBMC", "Control2PBMC", "PatientPBMC")],
#            Features = DiffImmuneMarkers, CellTypes = AnalyseSubTypes, ColData = "MainCellType", DiffRes = DEG, colorlow = "#A36DBB", colorhigh = "#EBCDE6")
# dev.off()
# 
# svg(filename = "CCL5.svg", width = 8, height = 4)
# FeaturePlot(SeuratPBMCIntegrated, features = "CCL5", split.by = "Group")
# dev.off()
# 
# svg(filename = "IL1B.svg", width = 8, height = 4)
# FeaturePlot(SeuratPBMCIntegrated, features = "IL1B", split.by = "Group")
# dev.off()
# 
# MainCellPalette <- c(B = "#FF7979",
#                      pDC = "#BBF785",
#                      Stem = "#E2DC6C",
#                      Myeloid = "#ACD624",
#                      TNK = "#77C6DD",
#                      Platelets = "#D691EB"
#                      )
# 
# svg(filename = "Cluster_distribution.svg", width = 6, height = 6)
# tmp <- FetchData(SeuratPBMCIntegrated, vars = c("UMAP_1", "UMAP_2", "orig.ident"))
# tmp <- cbind(tmp, SeuratPBMCIntegrated@active.ident)
# colnames(tmp)[ncol(tmp)] <- "MainCellType"
# 
# tmp$Labels <- match(tmp$MainCellType, names(MainCellPalette))
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = MainCellType)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = MainCellPalette[names(MainCellPalette) %in% tmp$MainCellType]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = names(MainCellPalette),
#            color = LabelPositions$color, size = 4,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "none",
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank())
# dev.off()
# 
# #---------------------> Seurat-Style stacked violin plot <----------------------
# svg(filename = paste("All_Cell_Violin.svg", sep = ""), width = 6, height = 4)
# Features <- c("MS4A1", "LILRA4", "MEIS1", "LYZ", "CD3E", "PPBP")
# 
# long <- FetchData(SeuratPBMCIntegrated, vars = Features, slot = "data")
# long$MainType <- SeuratPBMCIntegrated@active.ident
# long$Cell <- rownames(long)
# long <- reshape2::melt(long, id.vars = c("Cell","MainType"), measure.vars = Features,
#                        variable.name = "Feat", value.name = "Expr")
# 
# long$MainType <- factor(x = long$MainType, levels = rev(names(MainCellPalette)))
# long$Color <- MainCellPalette[match(long$MainType, table = names(MainCellPalette))]
# 
# ggplot(long, aes(Expr, factor(MainType))) +
#   geom_violin(scale = "width", adjust = 1, trim = TRUE, aes(fill = MainType, color = MainType)) +
#   scale_fill_manual("legend", values = MainCellPalette) +
#   scale_color_manual("legend", values = MainCellPalette) +
#   scale_x_continuous(expand = c(0, 0), labels = function(x)
#     c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
#   facet_grid(cols = vars(Feat), scales = "free")  +
#   cowplot::theme_cowplot(font_size = 15) +
#   theme(legend.position = "none", panel.spacing = unit(0, "lines"),
#         panel.background = element_rect(fill = NA, color = "black"),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold"),
#         strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())
# dev.off()
# 
# save(SeuratPBMCIntegrated, file = "SeuratPBMC_Annotated_MainCellType.Rdata")
# 
# TNK.RNA <- SeuratPBMCIntegrated[,SeuratPBMCIntegrated@active.ident == "TNK"]
# B.RNA <- SeuratPBMCIntegrated[,SeuratPBMCIntegrated@active.ident == "B"]
# Myeloid.RNA <- SeuratPBMCIntegrated[,SeuratPBMCIntegrated@active.ident %in% c("Myeloid", "pDC")]
# 
# save(TNK.RNA, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/ReannotationAndElimination/SeuratPBMC_Annotated_TNK.Rdata")
# save(B.RNA, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/ReannotationAndElimination/SeuratPBMC_Annotated_B.Rdata")
# save(Myeloid.RNA, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/ReannotationAndElimination/SeuratPBMC_Annotated_Myeloid.Rdata")

################################## Cell chat ###################################
rm(list = ls())
setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/CellChat/")

library(ggplot2)
library(Seurat)
library(dplyr)
library(CellChat)
library(emojifont)

source("/mnt/TRAINING/Zhongjiacheng/Customized Codes.R")

load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/MergeDatasets/SeuratPBMCIntegrated.Rdata")
svg(filename = paste("Elbow_Plot.svg", sep = ""), width = 7, height = 5)
ElbowPlot(SeuratPBMCIntegrated, ndims = 50) # to determine number of dimensions for clustering
dev.off()

pc.use <- 37

SeuratPBMCIntegrated <- FindNeighbors(SeuratPBMCIntegrated, dims = 1:pc.use, verbose = T)
SeuratPBMCIntegrated <- FindClusters(SeuratPBMCIntegrated, resolution = 1.2)
SeuratPBMCIntegrated <- RunUMAP(SeuratPBMCIntegrated, dims = 1:pc.use)

svg(filename = paste("SeuratPBMC_After_batch_correction.svg", sep = ""), width = 6, height = 5)
print(DimPlot(SeuratPBMCIntegrated, reduction = "umap",
              group.by = "orig.ident", cols = c("#6EBDF8",
                                                "#6BC567",
                                                "#E680C4",
                                                "#E1B385",
                                                "#A687DF",
                                                "#F98077")) +
        theme(plot.title = element_blank()))
dev.off()

SeuratPBMCIntegrated$Group <- ""
SeuratPBMCIntegrated$Group[grep(SeuratPBMCIntegrated$orig.ident, pattern = "Control")] <- "Control"
SeuratPBMCIntegrated$Group[grep(SeuratPBMCIntegrated$orig.ident, pattern = "Patient")] <- "Patient"

DefaultAssay(SeuratPBMCIntegrated) <- "RNA"

Reanno.cluster.ids <- levels(SeuratPBMCIntegrated@active.ident)

Reanno.cluster.ids[c(9, 25)+1] <- "NK"
Reanno.cluster.ids[c(17)+1] <- "NKT"
Reanno.cluster.ids[c(4, 8, 10, 20)+1] <- "CD4"
Reanno.cluster.ids[c(6, 12, 15)+1] <- "CD8"
Reanno.cluster.ids[c(1, 11, 16, 22)+1] <- "Naive CD4"
Reanno.cluster.ids[c(0)+1] <- "Naive CD8"
Reanno.cluster.ids[c(18)+1] <- "Treg"
Reanno.cluster.ids[c(7, 13, 29)+1] <- "B"
Reanno.cluster.ids[c(14)+1] <- "CD16+ mono"
Reanno.cluster.ids[c(2, 3, 5, 19, 26)+1] <- "CD14+ mono"
Reanno.cluster.ids[c(23)+1] <- "Dendritic"
Reanno.cluster.ids[c(27)+1] <- "Platelets"
Reanno.cluster.ids[c(24)+1] <- "pDC"
Reanno.cluster.ids[c(28)+1] <- "Stem"
Reanno.cluster.ids[c(21)+1] <- "MAIT/γδ"

names(Reanno.cluster.ids) <- levels(SeuratPBMCIntegrated@active.ident)
SeuratPBMCIntegrated <- RenameIdents(SeuratPBMCIntegrated, Reanno.cluster.ids)

# Discard Platelets and Stem

CellchatSeu <- SeuratPBMCIntegrated[,SeuratPBMCIntegrated@active.ident %in% c("NK", "NKT", "CD4", "CD8",
                                                                              "Naive CD4", "Naive CD8", "Treg",
                                                                              "B", "pDC", "Dendritic",
                                                                              "CD14+ mono", "CD16+ mono")]

CellchatSeu$CellType <- CellchatSeu@active.ident

CellChatPalette <- c('Naive CD8' = "#A8ABEA",
                   'Naive CD4' = "#B2E6EC",
                   'CD4' = "#77C6DD",
                   'CD8' = "#0D97D5",
                   'Treg' = "#4AB4B1",
                   'NK' = "#8672D0",
                   'NKT' = "#6F9DB9",
                   'CD14+ mono' = "#C9EF75",
                   'CD16+ mono' = "#73CF9F",
                   'Dendritic' = "#93E9B6",
                   'B' = "#FF7979",
                   'pDC' = "#BBF785"
)

Cluster2Unicode <- readxl::read_xlsx(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Dec2UnicodeSymbol.xlsx")

svg(filename = "CellChatDim.svg", width = 9.5, height = 7.5)
tmp <- FetchData(CellchatSeu, vars = c("UMAP_1", "UMAP_2", "orig.ident", "CellType"))

tmp$Labels <- match(tmp$CellType, names(CellChatPalette))

LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
                             UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
LabelPositions$color <- "black"

ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = CellType)) +
  geom_point(size = 0.2, alpha = 1) +
  scale_color_manual(values = CellChatPalette[names(CellChatPalette) %in% tmp$CellType]) +
  annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rank(match(rownames(LabelPositions), names(CellChatPalette))),
           color = LabelPositions$color, size = 7,  fontface = "bold") +
  theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
                        legend.key.size = unit(1.3, 'cm'),
                        legend.text = element_text(size = 17),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 9,
                                                  shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$CellType))]
  )))
dev.off()

# ######################### Cell Chat preprocessing ##############################
# 
# object.list <- list(Control = CellChatProcessing(CellchatSeu[,grep(CellchatSeu$orig.ident, pattern = "Control")], Grouping = "CellType"),
#                     Patient = CellChatProcessing(CellchatSeu[,grep(CellchatSeu$orig.ident, pattern = "Patient")], Grouping = "CellType"))
# 
# #### Note that the second dataset is patient !!!!!!!
# 
# cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# 
# UseCellNames <- levels(unlist(cellchat@idents))
# names(UseCellNames) <- CellChatPalette[match(UseCellNames, names(CellChatPalette))]
# 
# gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
# gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
# gg1 + gg2
# 
# # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
# pos.dataset <-"Patient"
# # define a char name used for storing the results of differential expression analysis
# features.name <- pos.dataset
# 
# # perform differential expression analysis
# cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset,
#                                        features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# 
# net <- netMappingDEG(cellchat, features.name = features.name)
# net.up <- subsetCommunication(cellchat, net = net, datasets = "Patient", ligand.logFC = 0.25, receptor.logFC = 0.25)
# net.down <- subsetCommunication(cellchat, net = net, datasets = "Patient", ligand.logFC = -0.25, receptor.logFC = -0.25)
# 
# save(cellchat, net.up, net.down, UseCellNames, CellChatPalette, object.list, file = "cellchat.Rdata")

load(file = "cellchat.Rdata")

# Since the signaling genes in the net.up and net.down might be complex with multi-subunits,
# we can do further deconvolution to obtain the individual signaling genes.

net.up[net.up$source == "CD8" | net.up$target == "CD8",]
net.down[net.down$source == "CD8" | net.down$target == "CD8",]

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

svg(filename = "heatmap interaction number.svg", width = 6, height = 6)
netVisual_diffInteraction(cellchat, weight.scale = T, color.use = names(UseCellNames))
dev.off()

svg(filename = "heatmap interaction strength.svg", width = 6, height = 6)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", color.use = names(UseCellNames))
dev.off()

gg1 <- netVisual_heatmap(cellchat, color.use = names(UseCellNames))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.use = names(UseCellNames))
#> Do heatmap based on a merged object

svg(filename = "heatmap interaction heatmap.svg", width = 7.5, height = 3.5)
print(gg1 + gg2)
dev.off()

# compare all the interactions sending from antigen presenting cells to T cells
svg(filename = "Antigen presenting cells to Patient T Up.svg", width = 6, height = 6)
netVisual_chord_gene(object.list[[2]], sources.use = grep(paste(c("B", "Dendritic"), collapse = "|"), UseCellNames), net = net.up,
                     targets.use = grep(paste(c("CD8", "NK", "CD4", "Treg"), collapse = "|"), UseCellNames), color.use = names(UseCellNames), lab.cex = 0.45,
                     title.name = "Increased signaling from patient's antigen presenting cells to T cells", show.legend = FALSE)
dev.off()

svg(filename = "Antigen presenting cells to Patient T Down.svg", width = 6, height = 6)
netVisual_chord_gene(object.list[[2]], sources.use = grep(paste(c("B", "Dendritic"), collapse = "|"), UseCellNames), net = net.down,
                     targets.use = grep(paste(c("CD8", "NK", "CD4", "Treg"), collapse = "|"), UseCellNames), color.use = names(UseCellNames), lab.cex = 0.45,
                     title.name = "Reduced signaling from patient's antigen presenting cells to T cells", show.legend = FALSE)
dev.off()

svg(filename = "Legend.svg", width = 9, height = 9)
netVisual_chord_gene(object.list[[2]], sources.use = grep(paste(c("B", "Dendritic", "CD8", "NK", "CD4", "Treg"), collapse = "|"), UseCellNames), net = net.up,
                     targets.use = grep(paste(c("Treg"), collapse = "|"), UseCellNames), color.use = names(UseCellNames), lab.cex = 0.5,
                     title.name = "Increased signaling from patient's antigen presenting cells to T cells", show.legend = TRUE)
dev.off()

svg(filename = "upgulated and down-regulated signaling ligand-receptor pairs.svg", width = 10, height = 10)
netVisual_bubble(cellchat, sources.use = grep(paste(c("B", "Dendritic"), collapse = "|"), UseCellNames),
                 targets.use = grep(paste(c("CD8", "NK", "CD4", "Treg"), collapse = "|"), UseCellNames),
                 comparison = c(2, 1), angle.x = 45, remove.isolate = FALSE)
dev.off()

###########################scATAC Preprocessing#################################
# rm(list = ls())
# setwd(dir = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/")
# 
# suppressMessages(library(ArchR))
# addArchRThreads(threads = 16)
# addArchRGenome("hg38")
# 
# minTSS <- 9
# minFrags <- 1500
# 
# ArrowFiles <- createArrowFiles(
#   inputFiles = c("/mnt/TRAINING/Zhongjiacheng/XRY/PBMC-XRY-scATAC/fragments.tsv.gz",
#                  "/mnt/TRAINING/Zhongjiacheng/10xgenomics/ATAC+RNA/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_atac_fragments.tsv.gz",
#                  "/mnt/TRAINING/Zhongjiacheng/10xgenomics/ATAC+RNA/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"),
#   sampleNames = c("Patient-scATAC", "Control1-scATAC", "Control2-scATAC"),
#   minTSS = minTSS,
#   minFrags = minFrags,
#   addTileMat = T,
#   addGeneScoreMat = T, force = T
# )
# 
# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search with doublet projection.
#   LSIMethod = 1, outDir = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/Doublet/"
# )
# 
# projectPBMC <- ArchRProject(
#   ArrowFiles = ArrowFiles,
#   outputDirectory = "projectPBMC/",
#   copyArrows = T #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
# )
# 
# getAvailableMatrices(projectPBMC)
# 
# head(projectPBMC$cellNames)
# quantile(projectPBMC$TSSEnrichment)
# 
# getCellColData(projectPBMC, select = "nFrags")
# getCellColData(projectPBMC, select = c("log10(nFrags)", "TSSEnrichment"))
# 
# svg(filename = "QualityControl/Ridges of TSSEnrichment.svg", width = 8, height = 8)
# plotGroups(
#   ArchRProj = projectPBMC,
#   groupBy = "Sample",
#   colorBy = "cellColData",
#   name = "TSSEnrichment",
#   plotAs = "ridges",
#   baseSize = 18,
#   pal = c("Control1-scATAC" = "#6EBDF8",
#                              "Control2-scATAC" = "#6BC567",
#                              "Patient-scATAC" = "#F98077")
# )
# dev.off()
# 
# projectPBMC <- filterDoublets(projectPBMC)
# 
# projectPBMC <- addIterativeLSI(ArchRProj = projectPBMC)
# 
# saveArchRProject(ArchRProj = projectPBMC, outputDirectory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/PreprocessedATAC/", load = F, overwrite = T)
# 
# ##
# 
# rm(list = ls())
# setwd(dir = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/")
# 
# projectPBMC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/PreprocessedATAC/"))
# 
# projectPBMC <- addUMAP(
#   ArchRProj = projectPBMC,
#   reducedDims = "IterativeLSI",
#   name = "UMAP",
#   nNeighbors = 30,
#   minDist = 0.5,
#   metric = "cosine",
#   seed = 2
# )
# 
# svg(filename = "BatchCorrection/UmapATAC before batch correction.svg", width = 8, height = 8.5)
# tmp <- getEmbedding(ArchRProj = projectPBMC, embedding = "UMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$Sample <- projectPBMC$Sample
# 
# SamplePalette <- c("Control1-scATAC" = "#6EBDF8",
#                    "Control2-scATAC" = "#6BC567",
#                    "Patient-scATAC" = "#F98077")
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = Sample)) +
#   geom_point(size = 0.2, alpha = 1) + scale_color_manual(values = SamplePalette) +
#   theme_ArchR() + theme(legend.position = "bottom", legend.title = element_blank(),
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(override.aes = list(size = 10)))
# dev.off()
# 
# projectPBMC <- addHarmony(
#   ArchRProj = projectPBMC,
#   reducedDims = "IterativeLSI",
#   name = "Harmony",
#   groupBy = "Sample"
# )
# 
# projectPBMC <- addUMAP(
#   ArchRProj = projectPBMC,
#   reducedDims = "Harmony",
#   name = "HarmonyUMAP",
#   nNeighbors = 30,
#   minDist = 0.5,
#   metric = "cosine",
#   seed = 2
# )
# 
# svg(filename = "BatchCorrection/UmapATAC after batch correction.svg", width = 8, height = 8.5)
# tmp <- getEmbedding(ArchRProj = projectPBMC, embedding = "HarmonyUMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$Sample <- projectPBMC$Sample
# 
# SamplePalette <- c("Control1-scATAC" = "#6EBDF8",
#                    "Control2-scATAC" = "#6BC567",
#                    "Patient-scATAC" = "#F98077")
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = Sample)) +
#   geom_point(size = 0.2, alpha = 1) + scale_color_manual(values = SamplePalette) +
#   theme_ArchR() + theme(legend.position = "bottom", legend.title = element_blank(),
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(override.aes = list(size = 10)))
# dev.off()
# 
# saveArchRProject(ArchRProj = projectPBMC, outputDirectory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/BatchCorrectedATAC/", load = F, overwrite = T)

####################scATAC-scRNA Integration-Main Cell Types####################
# rm(list = ls())
# 
# library(ggplot2)
# library(Seurat)
# library(dplyr)
# suppressMessages(library(ArchR))
# addArchRThreads(threads = 16)
# 
# setwd(dir = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Main cell types/")
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/ReannotationAndElimination/SeuratPBMC_Annotated_MainCellType.Rdata")
# projectPBMC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/BatchCorrectedATAC/"))
# 
# SeuratPBMCIntegrated <- SeuratPBMCIntegrated[,SeuratPBMCIntegrated$orig.ident %in% c("Control1PBMC", "Control2PBMC", "PatientPBMC")]
# SeuratPBMCIntegrated <- DietSeurat(SeuratPBMCIntegrated, assays = "RNA")
# 
# SeuratPBMCIntegrated$MainCellTypes <- SeuratPBMCIntegrated@active.ident
# 
# getAvailableMatrices(projectPBMC)
# 
# groupList <- SimpleList()
# for (i in c("Control1", "Control2", "Patient")) {
#   RNA.cells <- colnames(SeuratPBMCIntegrated[,grep(SeuratPBMCIntegrated$orig.ident, pattern = i)])
#   ATAC.cells <- grep(projectPBMC$cellNames, pattern = i, value = T)
# 
#   groupList[[i]] <- SimpleList(
#     ATAC = ATAC.cells,
#     RNA = RNA.cells
#   )
# }
# 
# projectPBMC <- addGeneIntegrationMatrix(
#   ArchRProj = projectPBMC,
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = SeuratPBMCIntegrated,
#   addToArrow = F,
#   groupList = groupList,
#   groupRNA = "MainCellTypes",
#   nameCell = "predictedCell_Co",
#   nameGroup = "predictedGroup_Co",
#   nameScore = "predictedScore_Co",
#   force = T
# )
# 
# print(plotEmbedding(
#   projectPBMC,
#   colorBy = "cellColData",
#   name = "predictedGroup_Co",
#   embedding = "HarmonyUMAP",
# ))
# 
# svg(filename =  "Constrained_Main_Cell_Type_Transfer_Score.svg", width = 6, height = 6)
# print(plotEmbedding(
#   projectPBMC,
#   colorBy = "cellColData",
#   name = "predictedScore_Co",
#   embedding = "HarmonyUMAP"
# ))
# dev.off()
# 
# svg(filename = "Constrained_Main_Cell_Type_Transfer_Group.svg", width = 6, height = 6)
# MainCellPalette <- c(B = "#FF7979",
#                      Stem = "#E2DC6C",
#                      pDC = "#BBF785",
#                      Myeloid = "#ACD624",
#                      TNK = "#77C6DD",
#                      Platelets = "#D691EB"
#                      )
# 
# tmp <- getEmbedding(ArchRProj = projectPBMC, embedding = "HarmonyUMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$predictedGroup_Co <- projectPBMC$predictedGroup_Co
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[grep(rownames(tmp), pattern = "Control"),"UMAP_1"]), INDEX = tmp[grep(rownames(tmp), pattern = "Control"),"predictedGroup_Co"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[grep(rownames(tmp), pattern = "Control"),"UMAP_2"]), INDEX = tmp[grep(rownames(tmp), pattern = "Control"),"predictedGroup_Co"], FUN = mean)))
# 
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = predictedGroup_Co)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = MainCellPalette[names(MainCellPalette) %in% tmp$predictedGroup_Co]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = names(MainCellPalette)[match(rownames(LabelPositions), names(MainCellPalette))],
#            color = LabelPositions$color, size = 4,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "none",
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank())
# dev.off()
# 
# TNK <- rownames(projectPBMC)[projectPBMC$predictedGroup_Co == "TNK"]
# B <- rownames(projectPBMC)[projectPBMC$predictedGroup_Co == "B"]
# Myeloid <- rownames(projectPBMC)[projectPBMC$predictedGroup_Co %in% c("Myeloid", "pDC")]
# 
# saveArchRProject(ArchRProj = subsetCells(projectPBMC, TNK), outputDirectory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/AnnotatedATAC/TNK/", load = F, overwrite = T)
# saveArchRProject(ArchRProj = subsetCells(projectPBMC, B), outputDirectory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/AnnotatedATAC/B/", load = F, overwrite = T)
# saveArchRProject(ArchRProj = subsetCells(projectPBMC, Myeloid), outputDirectory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/AnnotatedATAC/Myeloid/", load = F, overwrite = T)