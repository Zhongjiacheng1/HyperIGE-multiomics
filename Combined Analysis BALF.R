########################## MergeDatasets-BALF ##################################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratBALF/MergeDatasets/")
# 
# library(ggplot2)
# library(Seurat)
# library(dplyr)
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/QualityControl/PreprocessedBALFSamples.Rdata")
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/cycle.rda")
# 
# pc.use <- 30
# 
# RawSeuratBALF <- merge(SeuratPatientBALF,
#                        y = c(SeuratControl6BALF, SeuratControl7BALF,
#                              SeuratControl8BALF, SeuratControl9BALF,
#                              SeuratControl10BALF)
# )
# 
# # Cell cycle scoring
# 
# RawSeuratBALF <- NormalizeData(RawSeuratBALF) %>%
#   CellCycleScoring(g2m.features = g2m_genes, s.features = s_genes) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# RawSeuratBALF <- ScaleData(RawSeuratBALF)
# RawSeuratBALF <- RunPCA(RawSeuratBALF)
# 
# # Considerable variations due to cell cycle phase were observed. Based on this plot, we would regress out the variation due to cell cycle.
# 
# svg(filename = paste("Cell_cycle_scoring.svg", sep = ""), width = 10, height = 4)
# print(DimPlot(RawSeuratBALF,
#               reduction = "pca",
#               group.by = "Phase",
#               split.by = "Phase") +
#         theme(plot.title = element_blank()))
# dev.off()
# 
# RawSeuratBALF <-RunUMAP(RawSeuratBALF, dims = 1:pc.use)
# 
# svg(filename = paste("SeuratBALF_Before_batch_correction.svg", sep = ""), width = 6, height = 5)
# print(DimPlot(RawSeuratBALF, reduction = "umap", group.by = "orig.ident",
#               cols = c("#6EBDF8",
#                        "#6BC567",
#                        "#E680C4",
#                        "#E1B385",
#                        "#A687DF",
#                        "#F98077")) +
#         theme(plot.title = element_blank()))
# dev.off()
# 
# DefaultAssay(RawSeuratBALF) <- "RNA"
# RawSeuratBALF <- DietSeurat(RawSeuratBALF, assays = "RNA")
# SeuratBALF.list <- SplitObject(RawSeuratBALF, split.by = "orig.ident")
# SeuratBALF.list <- lapply(SeuratBALF.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = T) # This step was intended to display the original normalized counts in subsequent visualization (SCTransform does not normalize the original RNA assay)
#   x <- SCTransform(x, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = T)
# })
# 
# integ_features <- SelectIntegrationFeatures(object.list = SeuratBALF.list,
#                                             nfeatures = 3000)
# 
# SeuratBALF.list <- PrepSCTIntegration(object.list = SeuratBALF.list,
#                                    anchor.features = integ_features)
# 
# integ_anchors <- FindIntegrationAnchors(object.list = SeuratBALF.list,
#                                         normalization.method = "SCT",
#                                         anchor.features = integ_features)
# 
# SeuratBALFIntegrated <- IntegrateData(anchorset = integ_anchors,
#                                    normalization.method = "SCT")
# 
# SeuratBALFIntegrated <- RunPCA(object = SeuratBALFIntegrated)
# 
# save(SeuratBALFIntegrated, file = "SeuratBALFIntegrated.Rdata")
# 
######################### BALF Initial Integration #############################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratBALF/CellTypeIdentification/")
# 
# library(ggplot2)
# library(Seurat)
# library(dplyr)
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratBALF/MergeDatasets/SeuratBALFIntegrated.Rdata")
# 
# svg(filename = paste("Elbow_Plot.svg", sep = ""), width = 7, height = 5)
# ElbowPlot(SeuratBALFIntegrated, ndims = 50) # to determine number of dimensions for clustering
# dev.off()
# 
# pc.use <- 42
# SeuratBALFIntegrated <- FindNeighbors(SeuratBALFIntegrated, dims = 1:pc.use, verbose = T)
# SeuratBALFIntegrated <- FindClusters(SeuratBALFIntegrated, resolution = 3)
# SeuratBALFIntegrated <- RunUMAP(SeuratBALFIntegrated, dims = 1:pc.use)
# 
# svg(filename = paste("SeuratBALF_After_batch_correction.svg", sep = ""), width = 6, height = 5)
# DimPlot(SeuratBALFIntegrated, reduction = "umap", group.by = "orig.ident", cols = c("#6EBDF8",
#                                                                                     "#6BC567",
#                                                                                     "#E680C4",
#                                                                                     "#E1B385",
#                                                                                     "#A687DF",
#                                                                                     "#F98077")) +
#   theme(plot.title = element_blank())
# dev.off()
# 
# DimPlot(SeuratBALFIntegrated, reduction = "umap", label = T)
# 
# # Switch back to RNA assay to perform down-stream differential analyses
# DefaultAssay(SeuratBALFIntegrated) <- "RNA"
# 
# BALFPBMCIntegrated.Markers <- FindAllMarkers(SeuratBALFIntegrated, only.pos = F)
# 
# save(BALFPBMCIntegrated.Markers, SeuratBALFIntegrated, file = "SeuratBALFIntegratedandMarkers.Rdata")

######################## Combined Analysis MainCellTypes #######################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/MainCellTypes")
# 
# library(ggplot2)
# library(Seurat)
# library(dplyr)
# 
# source("/mnt/TRAINING/Zhongjiacheng/Customized Codes.R")
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratBALF/CellTypeIdentification/SeuratBALFIntegratedandMarkers.Rdata")
# 
# BALFPBMCIntegrated.Markers$Delta <- BALFPBMCIntegrated.Markers$pct.1 - BALFPBMCIntegrated.Markers$pct.2
# top10 <- BALFPBMCIntegrated.Markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_log2FC)
# 
# Reanno.cluster.ids <- levels(SeuratBALFIntegrated@active.ident)
# 
# Reanno.cluster.ids <- rep("Myeloid", length(Reanno.cluster.ids))
# Reanno.cluster.ids[c(46)+1] <- "Goblet" # TFF3
# Reanno.cluster.ids[c(44)+1] <- "Hillock" # KRT13
# Reanno.cluster.ids[c(38)+1] <- "Ciliated" # DNAH5
# Reanno.cluster.ids[c(4, 13, 20, 24, 42)+1] <- "TNK"
# Reanno.cluster.ids[c(43)+1] <- "B"
# 
# names(Reanno.cluster.ids) <- levels(SeuratBALFIntegrated@active.ident)
# SeuratBALFIntegrated <- RenameIdents(SeuratBALFIntegrated, Reanno.cluster.ids)
# 
# SeuratBALFIntegrated$CellType <- SeuratBALFIntegrated@active.ident
# 
# MainCellPalette <- c("Myeloid" = "#ACD624",
#                      "TNK" = "#77C6DD",
#                      "B" = "#FF7979",
#                      "Hillock" = "#D5D5D5",
#                      "Ciliated" = "#B6B6B6",
#                      "Goblet" = "#868686"
# )
# 
# SeuratBALFIntegrated$Group <- ""
# SeuratBALFIntegrated$Group[grep(SeuratBALFIntegrated$orig.ident, pattern = "Control")] <- "Control"
# SeuratBALFIntegrated$Group[grep(SeuratBALFIntegrated$orig.ident, pattern = "Patient")] <- "Patient"
# 
# SeuratBALFIntegrated$MainCellType <- as.character(SeuratBALFIntegrated@active.ident)
# 
# SeuratBALFIntegrated$MainCellType[SeuratBALFIntegrated$MainCellType %in% c("φ", "Cycling φ", "Dendritic")] <-  "Myeloid"
# 
# SeuratBALFIntegrated$GroupMainCellType <- paste(SeuratBALFIntegrated$Group, SeuratBALFIntegrated$MainCellType)
# 
# DEG <- list()
# 
# for (CellType in c("Myeloid", "TNK", "B")) {
#   FirstGroup <- paste("Patient", CellType, sep = " ")
#   SecondGroup <- paste("Control", CellType, sep = " ")
# 
#   tmp <- FindMarkers(object = SeuratBALFIntegrated,
#                            ident.1 = FirstGroup,
#                            ident.2 = SecondGroup,
#                            group.by = "GroupMainCellType",
#                            logfc.threshold = 0.2,
#                            test.use = "MAST")
# 
#   Retain <- tmp$p_val_adj < 0.05
# 
#   DEG[[CellType]] <- tmp[Retain, ]
# }
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
# ImmuneMarkers <- grep(paste(c("^CCL", "^CXCL", "^CXCR", "^CCR",
#                               "^TNF", "^TGF", "IFN",
#                               paste("^IL", c(1:50),
#                                     sep = "")), collapse = "|"), rownames(SeuratBALFIntegrated), value = T)
# 
# DiffImmuneMarkers <- Reduce(union, list(rownames(DEG$TNK[rownames(DEG$TNK) %in% c(PotentialBioMarkers, ImmuneMarkers),]),
#                        rownames(DEG$B[rownames(DEG$B) %in% c(PotentialBioMarkers, ImmuneMarkers),]),
#                        rownames(DEG$Myeloid[rownames(DEG$Myeloid) %in% c(PotentialBioMarkers, ImmuneMarkers),])))
# 
# PBMCImmuneMarkers <- c("LAG3", "TIGIT",
#                        "CCL4", "CCL4L2", "CCL5",
#                        "CXCR3", "CXCR4", "CXCR5", "CCR7",
#                        "TNFRSF1A", "TNFRSF1B", "TNFSF10", "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF25",
#                        "TNFAIP3", "TNFAIP2", "TNFAIP8",
#                        "IFNAR2", "IFNGR1",
#                        "IL1B", "IL15", "IL16", "IL32", "TGFB1",
#                        "IL2RG", "IL4R", "IL7R", "IL6R", "IL13RA1", "IL17RA", "TGFBR2",
#                        "IL6ST",
#                        "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
#                        "GZMA", "GZMH", "GZMH", "GZMM", "PRF1", "TNF",
#                        "S100A8", "S100A9", "LYZ"
# )
# 
# levels(DiffImmuneMarkers) <- rev(c("CCL3", "CCL4", "CCL4L2", "CCL3L1", "CCL18", "CXCL5", "CXCL8","CXCL10", "CXCL16",
#                                    "CXCR2", "CXCR4", "CXCR6", "CCR1", "CCRL2", "CCR4", "CCR5",
#                                    "TNFRSF1B", "TNFRSF9", "TNFSF12", "TNFRSF18", "TNFRSF4", "TNFSF13", "TNFAIP3", "TNFAIP2", "TNFAIP6", "TNFRSF14",
#                                    "IFNAR1", "IFNGR1", "IFNGR2",
#                                    "TGFBI",
#                                    "IL1B", "IL4I1", "IL16",
#                                    "IL2RG", "IL2RB", "IL3RA", "IL6R", "IL7R", "IL10RB", "TGFBR1", "TGFBR2",
#                                    "IL1RN", "IL6ST"))
# 
# AnalyseSubTypes <- c("Myeloid", "TNK", "B")
# levels(AnalyseSubTypes) <- c("Myeloid", "TNK", "B")
# 
# svg(filename = "MainCellTypesBubble.svg", width = 6, height = 8)
# BubblePlot(SeuratBALFIntegrated, Features = DiffImmuneMarkers, CellTypes = AnalyseSubTypes, DiffRes = DEG, ColData = "MainCellType", colorlow = "#BC9800", colorhigh = "#F4ECC8")
# dev.off()
# 
# intersect(PBMCImmuneMarkers, DiffImmuneMarkers)
# 
# svg(filename = "CCL4.svg", width = 8, height = 4)
# FeaturePlot(SeuratBALFIntegrated, features = "CCL4", split.by = "Group")
# dev.off()
# 
# svg(filename = "CCL4L2.svg", width = 8, height = 4)
# FeaturePlot(SeuratBALFIntegrated, features = "CCL4L2", split.by = "Group")
# dev.off()
# 
# svg(filename = "Cluster_distribution.svg", width = 6, height = 6)
# tmp <- FetchData(SeuratBALFIntegrated, vars = c("UMAP_1", "UMAP_2", "orig.ident"))
# tmp <- cbind(tmp, SeuratBALFIntegrated@active.ident)
# colnames(tmp)[ncol(tmp)] <- "MainCellType"
# 
# tmp$Labels <- match(tmp$MainCellType, names(MainCellPalette))
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
# LabelPositions$color <- "black"
# 
# CyCLabeladjust <- names(MainCellPalette)
# CyCLabeladjust[which(CyCLabeladjust == "Cycling φ")] <- "Cycling \n φ"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = MainCellType)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = MainCellPalette) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = CyCLabeladjust,
#            color = LabelPositions$color, size = 4,  fontface = "bold") +
#   ArchR::theme_ArchR() + theme(legend.position = "none",
#                                legend.text = element_text(size = 17),
#                                axis.title = element_blank(),
#                                axis.text = element_blank(),
#                                axis.ticks = element_blank())
# dev.off()
# 
# #---------------------> Seurat-Style stacked violin plot <----------------------
# 
# svg(filename = paste("All_Cell_Violin.svg", sep = ""), width = 6, height = 4.5)
# Features <- c("MARCO", "CD3E", "MS4A1", "KRT13", "DNAH5", "TFF3")
# 
# long <- FetchData(SeuratBALFIntegrated, vars = Features, slot = "data")
# long$MainType <- SeuratBALFIntegrated@active.ident
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
# Myeloid.BALF <- SeuratBALFIntegrated[,SeuratBALFIntegrated@active.ident %in% c("Myeloid")]
# TNK.BALF <- SeuratBALFIntegrated[,SeuratBALFIntegrated@active.ident %in% c("TNK")]
# B.BALF <- SeuratBALFIntegrated[,SeuratBALFIntegrated@active.ident %in% c("B")]
# 
# save(Myeloid.BALF, file = "Myeloid.BALF.Rdata")
# save(TNK.BALF, file = "TNK.BALF.Rdata")
# save(B.BALF, file = "B.BALF.Rdata")

################## Combined BALF Myeloid cells Reclustering ####################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/Myeloid/")
# 
# library(ggplot2)
# library(Seurat)
# library(dplyr)
# 
# load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/MainCellTypes/Myeloid.BALF.Rdata")
# 
# pc.use <- 30
# 
# Myeloid.BALF <- DietSeurat(Myeloid.BALF, assays = "RNA")
# 
# # Cell cycle scoring
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/cycle.rda")
# Myeloid.BALF <- NormalizeData(Myeloid.BALF) %>%
#   CellCycleScoring(g2m.features = g2m_genes, s.features = s_genes) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# Myeloid.BALF <- ScaleData(Myeloid.BALF)
# Myeloid.BALF <- RunPCA(Myeloid.BALF)
# 
# #Considerable differences due to cell cycle phase were observed. Based on this plot, we would regress out the variation due to cell cycle.
# 
# png(filename = paste("Cell_cycle_scoring.png", sep = ""), width = 6000, height = 2000, res = 300)
# DimPlot(Myeloid.BALF,
#         reduction = "pca",
#         group.by = "Phase",
#         split.by = "Phase")
# dev.off()
# 
# Myeloid.BALF <-RunUMAP(Myeloid.BALF, dims = 1:pc.use)
# 
# png(filename = paste("RawBALFPBMC_Myeloid_Before_batch_correction.png", sep = ""), width = 2500, height = 2000, res = 300)
# DimPlot(Myeloid.BALF, reduction = "umap", group.by = "orig.ident")
# dev.off()
# 
# DefaultAssay(Myeloid.BALF) <- "RNA"
# Myeloid.BALF <- DietSeurat(Myeloid.BALF, assays = "RNA")
# Myeloid.BALF.list <- SplitObject(Myeloid.BALF, split.by = "orig.ident")
# Myeloid.BALF.list <- lapply(Myeloid.BALF.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = T) # This step was intended to display the original normalized counts in subsequent visualization (SCTransform does not normalize the original RNA assay)
#   x <- SCTransform(x, vars.to.regress = c("percent.mt", "nCount_RNA", "S.Score", "G2M.Score"), verbose = T)
# })
# 
# integ_features <- SelectIntegrationFeatures(object.list = Myeloid.BALF.list,
#                                             nfeatures = 3000)
# 
# Myeloid.BALF.list <- PrepSCTIntegration(object.list = Myeloid.BALF.list,
#                                    anchor.features = integ_features)
# 
# integ_anchors <- FindIntegrationAnchors(object.list = Myeloid.BALF.list,
#                                         normalization.method = "SCT",
#                                         anchor.features = integ_features)
# 
# Myeloid.BALF.integrated <- IntegrateData(anchorset = integ_anchors,
#                                    normalization.method = "SCT")
# 
# Myeloid.BALF.integrated <- RunPCA(object = Myeloid.BALF.integrated)
# save(Myeloid.BALF.integrated, file = "Myeloid.BALF.integrated.Rdata")
# 
# load("Myeloid.BALF.integrated.Rdata")
# 
# png(filename = paste("Elbow_Plot.png", sep = ""), width = 2500, height = 2000, res = 300)
# ElbowPlot(Myeloid.BALF.integrated, ndims = 50) # to determine number of dimensions for clustering
# dev.off()
# 
# pc.use <- 45
# Myeloid.BALF.integrated <- FindNeighbors(Myeloid.BALF.integrated, dims = 1:pc.use, verbose = T)
# Myeloid.BALF.integrated <- FindClusters(Myeloid.BALF.integrated, resolution = 3)
# Myeloid.BALF.integrated <- RunUMAP(Myeloid.BALF.integrated, dims = 1:pc.use)
# 
# DefaultAssay(Myeloid.BALF.integrated) <- "RNA"
# 
# Markers.Myeloid.BALF.integrated <- FindAllMarkers(Myeloid.BALF.integrated, only.pos = F)
# 
# save(Myeloid.BALF.integrated, Markers.Myeloid.BALF.integrated, file = "Myeloid_BALF_reclustering.Rdata")

############### Combined Analysis BALF Myeloid cells Annotation ################
rm(list = ls())
setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/Myeloid/Analysis/")

library(Seurat)
library(dplyr)
library(ggplot2)
library(emojifont)
library(clusterProfiler)
library(org.Hs.eg.db)
suppressMessages(library(ArchR))

source("/mnt/TRAINING/Zhongjiacheng/Customized Codes.R")

load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/Myeloid/Myeloid_BALF_reclustering.Rdata")

Markers.Myeloid.BALF.integrated$Delta <- Markers.Myeloid.BALF.integrated$pct.1 - Markers.Myeloid.BALF.integrated$pct.2
top10 <- Markers.Myeloid.BALF.integrated %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_log2FC)

Myeloid.BALF.integrated$Group <- ""
Myeloid.BALF.integrated$Group[grep(Myeloid.BALF.integrated$orig.ident, pattern = "Control")] <- "Control"
Myeloid.BALF.integrated$Group[Myeloid.BALF.integrated$orig.ident %in% c("PatientBALF")] <- "Patient"

png(filename = paste("BALFMyeloidIntegrated_Cell_Distribution.png", sep = ""), width = 2500, height = 2000, res = 400)
DimPlot(Myeloid.BALF.integrated, label = T)
dev.off()

png(filename = paste("BALFMyeloidIntegrated_Cell_Distribution_Split.png", sep = ""), width = 6000, height = 1000, res = 300)
DimPlot(Myeloid.BALF.integrated, split.by = "orig.ident")
dev.off()

Reanno.cluster.ids <- levels(Myeloid.BALF.integrated@active.ident)

Reanno.cluster.ids <- rep("φ", length(Reanno.cluster.ids))
Reanno.cluster.ids[c(0)+1] <- "Sterol Synthetic φ" # FDPS, FDFT1, IDI1 PMID: 23040070
Reanno.cluster.ids[c(25)+1] <- "CXCL8+ Proinflammatory φ" # common IL1B CCL4
Reanno.cluster.ids[c(30)+1] <- "CCL20+ Proinflammatory φ" # common IL1B CCL4
Reanno.cluster.ids[c(33)+1] <- "CXCL10+ φ" #
Reanno.cluster.ids[c(4, 11, 17, 26, 27, 34, 37)+1] <- "Profibrotic φ" # CCL18
Reanno.cluster.ids[c(36)+1] <- "Cycling φ" # MKI67
Reanno.cluster.ids[c(15, 39)+1] <- "PHKG1+ φ" # PHKG1 (glycogenolytic regulatory enzyme)
Reanno.cluster.ids[c(19)+1] <- "Metallothionein φ" #MT  MT1X, MT1E..
Reanno.cluster.ids[c(29)+1] <- "Dendritic" # CD1C
Reanno.cluster.ids[c(23)+1] <- "IFN-responsive φ" # IFIT PMID: 32725722
Reanno.cluster.ids[c(13)+1] <- "SOD2+ φ" # SOD2  #######
Reanno.cluster.ids[c(17)+1] <- "Regulatory φ" # A2M
Reanno.cluster.ids[c(28)+1] <- "mono-like φ" # common VCAN FCN1
Reanno.cluster.ids[c(24, 31)+1] <- "mono-like Regulatory φ" # A2M/CCL2 common VCAN FCN1 also CCL13
Reanno.cluster.ids[c(27, 35)+1] <- "mono-like Profibrotic φ" # CCL18 common VCAN FCN1

names(Reanno.cluster.ids) <- levels(Myeloid.BALF.integrated@active.ident)
Myeloid.BALF.integrated <- RenameIdents(Myeloid.BALF.integrated, Reanno.cluster.ids)

DimPlot(Myeloid.BALF.integrated)
Myeloid.BALF.integrated$FineGranularity <- Myeloid.BALF.integrated@active.ident
Myeloid.BALF.integrated$GroupCluster <- paste(Myeloid.BALF.integrated$Group, Myeloid.BALF.integrated$FineGranularity, sep = " ")

Cluster2Unicode <- readxl::read_xlsx(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Dec2UnicodeSymbol.xlsx")
BALFMyeloidcellsPalette <- c('φ' = "#6BC567", ####
                         'PHKG1+ φ' = "#489269",
                         'Profibrotic φ' = "#BEECC0",
                         'Regulatory φ' = "#D5FF5D",
                         'CXCL8+ Proinflammatory φ' = "#1EC048",
                         'CCL20+ Proinflammatory φ' = "#96AB2B",
                         'CXCL10+ φ' = "#608F21",
                         'SOD2+ φ' = "#8CE101",
                         'mono-like Regulatory φ' = "#C1F39F",
                         'mono-like Profibrotic φ' = "#00CC99",
                         'mono-like φ' = "#81F381",
                         'IFN-responsive φ' = "#B0D955",
                         'Dendritic' = "#93E9B6", ####
                         'Metallothionein φ' = "#9AF8C5",
                         'Sterol Synthetic φ' = "#009249",
                         'Cycling φ' = "#818E32"
)

svg(filename = "Cluster_distribution.svg", width = 9.5, height = 7.5)

tmp <- FetchData(Myeloid.BALF.integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident"))
tmp <- cbind(tmp, Myeloid.BALF.integrated@active.ident)
colnames(tmp)[ncol(tmp)] <- "FineGranularity"

tmp$Labels <- match(tmp$FineGranularity, names(BALFMyeloidcellsPalette))

LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
                             UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
LabelPositions$color <- "black"

ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = FineGranularity)) +
  geom_point(size = 0.2, alpha = 1) +
  scale_color_manual(values = BALFMyeloidcellsPalette[names(BALFMyeloidcellsPalette) %in% tmp$FineGranularity]) +
  annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rownames(LabelPositions),
           color = LabelPositions$color, size = 5,  fontface = "bold") +
  theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
                        legend.text = element_text(size = 15),
                                                axis.title = element_blank(),
                                                axis.text = element_blank(),
                                                axis.ticks = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 10,
                                                  shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$FineGranularity))]
  )))

dev.off()

#---------------------> Seurat-Style stacked violin plot <----------------------

svg(filename = paste("All_Cell_Violin.svg", sep = ""), width = 11.5, height = 7.5)
Features <- c("MARCO", "CCL18", "SOD2", "PHKG1", "CCL2", "A2M", "IL1B", "CCL4", "CXCL8", "CCL20", "CXCL10", "VCAN", "FCN1", "IFIT1", "IFIT2", "IFIT3",
              "MT1X", "MT1E", "CD1C", "FDPS", "FDFT1", "MKI67")

long <- FetchData(Myeloid.BALF.integrated, vars = Features, slot = "data")
long$MainType <- Myeloid.BALF.integrated@active.ident
long$Cell <- rownames(long)
long <- reshape2::melt(long, id.vars = c("Cell","MainType"), measure.vars = Features,
                       variable.name = "Feat", value.name = "Expr")

long$MainType <- factor(x = long$MainType, levels = rev(names(BALFMyeloidcellsPalette)))
long$Color <- BALFMyeloidcellsPalette[match(long$MainType, table = names(BALFMyeloidcellsPalette))]

ggplot(long, aes(Expr, factor(MainType))) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, aes(fill = MainType, color = MainType)) +
  scale_fill_manual("legend", values = BALFMyeloidcellsPalette) +
  scale_color_manual("legend", values = BALFMyeloidcellsPalette) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(Feat), scales = "free")  +
  cowplot::theme_cowplot(font_size = 15) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
dev.off()

svg(filename = "DimPlotSplit.svg", width = 9.5, height = 7.5)
tmp <- FetchData(Myeloid.BALF.integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident"))
tmp <- cbind(tmp, Myeloid.BALF.integrated@active.ident)
colnames(tmp)[ncol(tmp)] <- "MainCellType"

tmp$Labels <- match(tmp$MainCellType, names(BALFMyeloidcellsPalette))

LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
                             UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
LabelPositions$color <- "black"

ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = MainCellType)) +  geom_point(size = 0.2, alpha = 1) +

  scale_color_manual(values = BALFMyeloidcellsPalette[names(BALFMyeloidcellsPalette) %in% tmp$MainCellType]) +

        theme_ArchR() + facet_wrap(~orig.ident, scales = "free") + theme(legend.position = "none",
                                                                         strip.text.x = element_text(size = 12))
dev.off()

svg(filename = "RelativeFraction_BALF_Myeloid.svg", width = 7, height = 6)
RelativeFraction(Myeloid.BALF.integrated, ColData = "FineGranularity")
dev.off()

# DEG <- list()
#
# for (CellType in levels(Myeloid.BALF.integrated@active.ident)) {
#   FirstGroup <- paste("Patient", CellType, sep = " ")
#   SecondGroup <- paste("Control", CellType, sep = " ")
#
#   DEG[[CellType]] <- FindMarkers(object = Myeloid.BALF.integrated,
#                                  ident.1 = FirstGroup,
#                                  ident.2 = SecondGroup,
#                                  group.by = "GroupCluster",
#                                  logfc.threshold = 0.25,
#                                  test.use = "MAST")
#
#   Retain <- DEG[[CellType]]$p_val_adj < 0.05
#
#   DEG[[CellType]] <- DEG[[CellType]][Retain, ]
# }
#
# save(DEG, file = "MyeloidDEG.Rdata")
load("MyeloidDEG.Rdata")

Myeloid.BALF.integrated$NumDEG <- 0

for (CellType in names(DEG)) {
  Myeloid.BALF.integrated$NumDEG[Myeloid.BALF.integrated$FineGranularity == CellType] <- lapply(DEG, nrow)[[CellType]]
}

svg(filename = "NumberDiffDimPlot.svg", width = 9.5, height = 6.5)
tmp <- FetchData(Myeloid.BALF.integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident", "NumDEG"))
ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = NumDEG)) + geom_point(size = 0.3, alpha = 1) +
  theme_ArchR() + scale_color_gradient2(low = "grey",
                                        mid = "#6BC567",
                                        high = "#D5FF5D",
                                        midpoint = 2500) +
  facet_wrap(~orig.ident, scales = "free") + theme(legend.position = "right",
                                                   strip.text.x = element_text(size = 12),
                                                   legend.text = element_text(size = 12),
                                                   legend.title = element_text(size = 12),
                                                   axis.ticks.x = element_blank(),
                                                   axis.ticks.y = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   axis.title.y = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   axis.text.y = element_blank())
dev.off()


svg(filename = "NumberDiffDimPlotByGroup.svg",  width = 7.5, height = 4)
tmp <- FetchData(Myeloid.BALF.integrated, vars = c("UMAP_1", "UMAP_2", "Group", "NumDEG"))
ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = NumDEG)) + geom_point(size = 0.3, alpha = 1) +
  theme_ArchR() + scale_color_gradient2(low = "grey",
                                        mid = "#6BC567",
                                        high = "#D5FF5D",
                                        midpoint = 2500) +
  facet_wrap(~Group, scales = "free") + theme(legend.position = "right",
                                                   strip.text.x = element_text(size = 12),
                                                   legend.text = element_text(size = 12),
                                                   legend.title = element_text(size = 12),
                                                   axis.ticks.x = element_blank(),
                                                   axis.ticks.y = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   axis.title.y = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   axis.text.y = element_blank())
dev.off()


##### GESA and bubble plot focus on clusters with DEG > 1000, and φ

AnalyseSubTypes <- names(which(lapply(DEG, nrow) > 1000))

AnalyseSubTypes <- c('φ',
                     'PHKG1+ φ',
                     'mono-like Regulatory φ',
                     'CXCL8+ Proinflammatory φ',
                     'mono-like Profibrotic φ',
                     'Dendritic')

levels(AnalyseSubTypes) <- c('φ',
                             'PHKG1+ φ',
                             'mono-like Regulatory φ',
                             'CXCL8+ Proinflammatory φ',
                             'mono-like Profibrotic φ',
                             'Dendritic')

# GSEA.res <- List()
# for (CellType in AnalyseSubTypes) {
#   GSEA_input <- DEG[[CellType]]$avg_log2FC
#   names(GSEA_input) <- rownames(DEG[[CellType]])
#   GSEA_input <- sort(GSEA_input, decreasing = TRUE)
#   GSEA.res[[CellType]] <- gseGO(GSEA_input, OrgDb = org.Hs.eg.db,  keyType = "SYMBOL", ont = "BP",
#                                minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 0.1)
# }
#
# save(DEG, GSEA.res, file = "MyeloidGSEA.Rdata")
load("MyeloidGSEA.Rdata")

svg(filename = paste("GSEA/", "GSEACXCL8.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$`CXCL8+ Proinflammatory φ`, Descriptions = c("oxidative phosphorylation",
                                                                   "phagocytosis",
                                                                   "positive regulation of nitric-oxide synthase biosynthetic process",
                                                                   "defense response to bacterium"),
             base_size = 15, GSEAColor = c("#BEECC0", "#489269", "#8CE101", "#9AF8C5"),
             legendposition = c(0.43, 0.3), rel_heights = c(1.5, 0.5, 1), maxoverlaps = 55)
dev.off()

svg(filename = paste("GSEA/", "GSEADendritic.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$Dendritic, Descriptions = c("positive chemotaxis", "positive regulation of leukocyte chemotaxis", "lymphocyte chemotaxis"),
  base_size = 15, GSEAColor = as.character(BALFMyeloidcellsPalette)[5:7],
  legendposition = c(0.3, 0.3), rel_heights = c(1.5, 0.5, 1), maxoverlaps = 60)
dev.off()

# The following potential biomarkers were previously reported to be altered in hyper IGE syndrome

PotentialBioMarkers <- c("IL13",
                         "IL6",
                         "IFNG",

                         "CXCL3", # GROγ
                         "CCL5", # RANTES
                         "LEP", # Leptin

                         "TNF",
                         "IL4",
                         "IL10",
                         "IL12",
                         "IL7",
                         "IL18",
                         "IL17A",
                         "IL1B",
                         "IL8",
                         "CXCL10", # IP-10
                         "CXCL12" # DF-1alpha
)

ImmuneMarkers <- grep(paste(c("^CCL", "^CXCL", "^CXCR", "^CCR",
                              "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
                              "^TGF", "^IFN", "^STAT", "^DEFB",
                              "^LYZ", "S100A8", "S100A9", "^CAMP", "^LCN2", "^LTF", "CD86", "CD80",
                              paste("^IL", c(1:50), sep = "")), collapse = "|"),
                      rownames(Myeloid.BALF.integrated), value = T)

DiffImmuneMarkers <- list()

# For visualization, a log2FC curoff of 0.4 was set

for (CellType in AnalyseSubTypes) {
  DiffImmuneMarkers[[CellType]] <- rownames(DEG[[CellType]][rownames(DEG[[CellType]]) %in% c(PotentialBioMarkers, ImmuneMarkers) &
                                                              abs(DEG[[CellType]]$avg_log2FC) > 0.4,])
}

DiffImmuneMarkers <- Reduce(union, list(DiffImmuneMarkers))
DiffImmuneMarkers <- unique(unlist(DiffImmuneMarkers))

levels(DiffImmuneMarkers) <- rev(c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                                   "STAT1", "STAT2", "STAT3", "STAT5B",
                                   "CCL2", "CCL3",  "CCL3L1" , "CCL4", "CCL4L2", "CCL13", "CCL15", "CCL18", "CCL23",
                                   "CCR1", "CCR5", "CCRL2",
                                   "CXCL1", "CXCL2", "CXCL3", "CXCL10", "CXCL8", "CXCL16", "CXCR2", "CXCR4", "IL1B", "IL18",
                                   "IL1RAP", "IL1R2", "IL1RN", "IL2RG", "IL3RA", "IL4I1", "IL6R", "IL6ST", "IL7R", "IL10RA", "IL10RB", "IL17RA", "IL18BP",
                                   "IFNGR2", "IFNGR1", "IFNAR1", "TGFB1", "TGFBR2", "TGFBR1", "TGFBI", "TNF",
                                   "CAMP", "DEFB1", "S100A8", "S100A9", "LYZ", "CD86"))

svg(filename = "BALFMyeloidBubble.svg", width = 6.5, height = 12.5)
BubblePlot(Myeloid.BALF.integrated[,Myeloid.BALF.integrated$FineGranularity %in% AnalyseSubTypes], Features = DiffImmuneMarkers,
           CellTypes = AnalyseSubTypes, ColData = "FineGranularity", DiffRes = DEG, colorlow = "#71BC16", colorhigh = "#F6FFD9")
dev.off()

svg(filename = "CXCLSubGroupNoCD86.svg", width = 8, height = 4)
FeaturePlot(Myeloid.BALF.integrated, features = "CD86", split.by = "Group")
dev.off()

######################Combined BALF T cells Reclustering########################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/TNK/")
# 
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(emojifont)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# suppressMessages(library(ArchR))
# 
# source("/mnt/TRAINING/Zhongjiacheng/Customized Codes.R")
# 
# load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/MainCellTypes/TNK.BALF.Rdata")
# 
# TNK.BALF <- TNK.BALF[,TNK.BALF$orig.ident %in% names(which(table(TNK.BALF$orig.ident) > 200))] # Note, Control9 was filtered here
# 
# pc.use <- 30
# 
# TNK.BALF <- DietSeurat(TNK.BALF, assays = "RNA")
# 
# # Cell cycle scoring
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/cycle.rda")
# TNK.BALF <- NormalizeData(TNK.BALF) %>%
#   CellCycleScoring(g2m.features = g2m_genes, s.features = s_genes) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# TNK.BALF <- ScaleData(TNK.BALF)
# TNK.BALF <- RunPCA(TNK.BALF)
# 
# #Considerable differences due to cell cycle phase were observed. Based on this plot, we would regress out the variation due to cell cycle.
# 
# png(filename = paste("Cell_cycle_scoring.png", sep = ""), width = 6000, height = 2000, res = 300)
# DimPlot(TNK.BALF,
#         reduction = "pca",
#         group.by = "Phase",
#         split.by = "Phase")
# dev.off()
# 
# TNK.BALF <-RunUMAP(TNK.BALF, dims = 1:pc.use)
# 
# png(filename = paste("RawBALFPBMC_TNK_Before_batch_correction.png", sep = ""), width = 2500, height = 2000, res = 300)
# DimPlot(TNK.BALF, reduction = "umap", group.by = "orig.ident")
# dev.off()
# 
# DefaultAssay(TNK.BALF) <- "RNA"
# TNK.BALF <- DietSeurat(TNK.BALF, assays = "RNA")
# TNK.BALF.list <- SplitObject(TNK.BALF, split.by = "orig.ident")
# TNK.BALF.list <- lapply(TNK.BALF.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = T) # This step was intended to display the original normalized counts in subsequent visualization (SCTransform does not normalize the original RNA assay)
#   x <- SCTransform(x, vars.to.regress = c("percent.mt", "nCount_RNA", "S.Score", "G2M.Score"), verbose = T)
# })
# 
# integ_features <- SelectIntegrationFeatures(object.list = TNK.BALF.list,
#                                             nfeatures = 3000)
# 
# TNK.BALF.list <- PrepSCTIntegration(object.list = TNK.BALF.list,
#                                    anchor.features = integ_features)
# 
# integ_anchors <- FindIntegrationAnchors(object.list = TNK.BALF.list,
#                                         normalization.method = "SCT",
#                                         anchor.features = integ_features)
# 
# TNK.BALF.integrated <- IntegrateData(anchorset = integ_anchors,
#                                    normalization.method = "SCT")
# 
# TNK.BALF.integrated <- RunPCA(object = TNK.BALF.integrated)
# save(TNK.BALF.integrated, file = "TNK.BALF.integrated.Rdata")
# 
# 
# load("TNK.BALF.integrated.Rdata")
# 
# png(filename = paste("Elbow_Plot.png", sep = ""), width = 2500, height = 2000, res = 300)
# ElbowPlot(TNK.BALF.integrated, ndims = 50) # to determine number of dimensions for clustering
# dev.off()
# 
# pc.use <- 32
# TNK.BALF.integrated <- FindNeighbors(TNK.BALF.integrated, dims = 1:pc.use, verbose = T)
# TNK.BALF.integrated <- FindClusters(TNK.BALF.integrated, resolution = 3)
# TNK.BALF.integrated <- RunUMAP(TNK.BALF.integrated, dims = 1:pc.use)
# 
# DefaultAssay(TNK.BALF.integrated) <- "RNA"
# 
# Markers.TNK.BALF.integrated <- FindAllMarkers(TNK.BALF.integrated, only.pos = F)
# 
# save(TNK.BALF.integrated, Markers.TNK.BALF.integrated, file = "TNK_BALF_reclustering.Rdata")

################# Combined Analysis BALF TNK cells Annotation ##################
rm(list = ls())
setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/TNK/Analysis/")

library(Seurat)
library(dplyr)
library(ggplot2)
library(emojifont)
library(clusterProfiler)
library(org.Hs.eg.db)
suppressMessages(library(ArchR))

source("/mnt/TRAINING/Zhongjiacheng/Customized Codes.R")

load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFOnly/TNK/TNK_BALF_reclustering.Rdata")

Markers.TNK.BALF.integrated$Delta <- Markers.TNK.BALF.integrated$pct.1 - Markers.TNK.BALF.integrated$pct.2
top10 <- Markers.TNK.BALF.integrated %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_log2FC)

TNK.BALF.integrated$Group <- ""
TNK.BALF.integrated$Group[grep(TNK.BALF.integrated$orig.ident, pattern = "Control")] <- "Control"
TNK.BALF.integrated$Group[TNK.BALF.integrated$orig.ident %in% c("PatientBALF")] <- "Patient"

png(filename = paste("BALFTNKIntegrated_Cell_Distribution.png", sep = ""), width = 2500, height = 2000, res = 400)
DimPlot(TNK.BALF.integrated, label = T)
dev.off()

png(filename = paste("BALFTNKIntegrated_Cell_Distribution_Split.png", sep = ""), width = 6000, height = 1000, res = 300)
DimPlot(TNK.BALF.integrated, split.by = "orig.ident")
dev.off()

Reanno.cluster.ids <- levels(TNK.BALF.integrated@active.ident)

Reanno.cluster.ids <- rep("CD4", length(Reanno.cluster.ids))
Reanno.cluster.ids[c(4, 10)+1] <- "CD8"
Reanno.cluster.ids[c(3, 9)+1] <- "GZMK+ CD8"
Reanno.cluster.ids[c(0, 1)+1] <- "GZMH+ CD8"
Reanno.cluster.ids[c(5, 11)+1] <- "Naive CD4"
Reanno.cluster.ids[c(7)+1] <- "Naive CD8"
Reanno.cluster.ids[c(13)+1] <- "Treg" # FOXP3
Reanno.cluster.ids[c(12)+1] <- "NK" # KLRF1
Reanno.cluster.ids[c(14, 16)+1] <- "MARCO+ LC"
Reanno.cluster.ids[c(17)+1] <- "CPA3+ LC"

names(Reanno.cluster.ids) <- levels(TNK.BALF.integrated@active.ident)
TNK.BALF.integrated <- RenameIdents(TNK.BALF.integrated, Reanno.cluster.ids)

DimPlot(TNK.BALF.integrated)
TNK.BALF.integrated$FineGranularity <- TNK.BALF.integrated@active.ident
TNK.BALF.integrated$GroupCluster <- paste(TNK.BALF.integrated$Group, TNK.BALF.integrated$FineGranularity, sep = " ")

Cluster2Unicode <- readxl::read_xlsx(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Dec2UnicodeSymbol.xlsx")
BALFTNKPalette <- c('CD4' = "#77C6DD",
                    'Naive CD4' = "#B2E6EC",
                    'Naive CD8' = "#A8ABEA",
                    'Treg' = "#4AB4B1",
                    'CD8' = "#0D97D5",
                    'GZMK+ CD8' = "#7B93DF",
                    'GZMH+ CD8' = "#8F84B6",
                    'NK' = "#8672D0",
                    'CPA3+ LC' = "#B695BD",
                    'MARCO+ LC' = "#338C95"
)

svg(filename = "Cluster_distribution.svg", width = 6, height = 4.5)

tmp <- FetchData(TNK.BALF.integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident"))
tmp <- cbind(tmp, TNK.BALF.integrated@active.ident)
colnames(tmp)[ncol(tmp)] <- "FineGranularity"

tmp$Labels <- match(tmp$FineGranularity, names(BALFTNKPalette))

LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
                             UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
LabelPositions$color <- "black"

ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = FineGranularity)) +
  geom_point(size = 0.2, alpha = 1) +
  scale_color_manual(values = BALFTNKPalette[names(BALFTNKPalette) %in% tmp$FineGranularity]) +
  annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rownames(LabelPositions),
           color = LabelPositions$color, size = 5,  fontface = "bold") +
  theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
                        legend.text = element_text(size = 12),
                        legend.key.size = unit(2, 'lines'),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 7,
                                                  shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$FineGranularity))]
  )))

dev.off()

#---------------------> Seurat-Style stacked violin plot <----------------------

svg(filename = paste("All_Cell_Violin.svg", sep = ""), width = 6.5, height = 4.5)
Features <- c("CD4", "CD8A", "SELL", "LEF1", "FOXP3", "GZMH", "GZMK", "KLRF1", "MARCO", "CPA3")

long <- FetchData(TNK.BALF.integrated, vars = Features, slot = "data")
long$FineGranularity <- TNK.BALF.integrated$FineGranularity
long$Cell <- rownames(long)
long <- reshape2::melt(long, id.vars = c("Cell","FineGranularity"), measure.vars = Features,
                       variable.name = "Feat", value.name = "Expr")

long$FineGranularity <- factor(x = long$FineGranularity, levels = rev(names(BALFTNKPalette)))
long$Color <- BALFTNKPalette[match(long$FineGranularity, table = names(BALFTNKPalette))]

ggplot(long, aes(Expr, factor(FineGranularity))) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, aes(fill = FineGranularity, color = FineGranularity)) +
  scale_fill_manual("legend", values = BALFTNKPalette) +
  scale_color_manual("legend", values = BALFTNKPalette) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(Feat), scales = "free")  +
  cowplot::theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
dev.off()

svg(filename = "RelativeFraction_BALF_TNK.svg", width = 6, height = 6)
RelativeFraction(TNK.BALF.integrated, ColData = "FineGranularity")
dev.off()

# DEG <- list()
# 
# for (CellType in levels(TNK.BALF.integrated@active.ident)) {
#   FirstGroup <- paste("Patient", CellType, sep = " ")
#   SecondGroup <- paste("Control", CellType, sep = " ")
# 
#   DEG[[CellType]] <- FindMarkers(object = TNK.BALF.integrated,
#                                  ident.1 = FirstGroup,
#                                  ident.2 = SecondGroup,
#                                  group.by = "GroupCluster",
#                                  logfc.threshold = 0.25,
#                                  test.use = "MAST")
# 
#   Retain <- DEG[[CellType]]$p_val_adj < 0.05
# 
#   DEG[[CellType]] <- DEG[[CellType]][Retain, ]
# }
# 
# save(DEG, file = "TNKDEG.Rdata")
load("TNKDEG.Rdata")

TNK.BALF.integrated$NumDiff <- 0

for (CellType in names(DEG)) {
  TNK.BALF.integrated$NumDiff[TNK.BALF.integrated$FineGranularity == CellType] <- lapply(DEG, nrow)[[CellType]]
}

svg(filename = "NumberDiffDimPlot.svg", width = 9.5, height = 6.5)
tmp <- FetchData(TNK.BALF.integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident", "NumDiff"))
ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = NumDiff)) + geom_point(size = 0.3, alpha = 1) +
  theme_ArchR() + scale_color_gradient2(low = "grey",
                                        mid = "#2171B5",
                                        high = "#BDD7E7",
                                        midpoint = 100) +
  facet_wrap(~orig.ident, scales = "free") + theme(legend.position = "right",
                                                   strip.text.x = element_text(size = 12),
                                                   legend.text = element_text(size = 12),
                                                   legend.title = element_text(size = 12),
                                                   axis.ticks.x = element_blank(),
                                                   axis.ticks.y = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   axis.title.y = element_blank(),
                                                   axis.text.x = element_blank(),
                                                   axis.text.y = element_blank())
dev.off()

AnalyseSubTypes <- c('Naive CD4',
                     'CD4', 
                     'Naive CD8',
                     'CD8',
                     'NK',
                     'Treg')

levels(AnalyseSubTypes) <- c('Naive CD4',
                             'CD4', 
                             'Naive CD8',
                             'CD8',
                             'NK',
                             'Treg')

# The following potential biomarkers were previously reported to be altered in hyper IGE syndrome

PotentialBioMarkers <- c("IL13",
                         "IL6",
                         "IFNG",
                         
                         "CXCL3", # GROγ
                         "CCL5", # RANTES
                         "LEP", # Leptin
                         
                         "TNF",
                         "IL4",
                         "IL10",
                         "IL12",
                         "IL7",
                         "IL18",
                         "IL17A",
                         "IL1B",
                         "IL8",
                         "CXCL10", # IP-10
                         "CXCL12" # DF-1alpha
)

ImmuneInhibitors <- c("LAG3", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "SLAMF6", "BTLA", "TIM3")

ImmuneMarkers <- grep(paste(c(ImmuneInhibitors, "CD8A", "CD8B", "^CCL", "^CXCL", "^CXCR", "^CCR",
                              "^TNF", "^TGF", "^IFN", "^GZM", "^STAT", "^PRF", "^DEFB", 
                              "^CAMP", "^LCN2", "CD28", "CD40L",
                              paste("^IL", c(1:50), sep = "")), collapse = "|"), 
                      rownames(TNK.BALF.integrated), value = T)

DiffImmuneMarkers <- list()

# For visualization, a log2FC curoff of 0.4 was set

for (CellType in AnalyseSubTypes) {
  DiffImmuneMarkers[[CellType]] <- rownames(DEG[[CellType]][rownames(DEG[[CellType]]) %in% c(PotentialBioMarkers, ImmuneMarkers) &
                                                              abs(DEG[[CellType]]$avg_log2FC) > 0.4,])
}

DiffImmuneMarkers <- Reduce(union, list(DiffImmuneMarkers))
DiffImmuneMarkers <- unique(unlist(DiffImmuneMarkers))

levels(DiffImmuneMarkers) <- rev(c("LAG3", "TIGIT",
                                   "CXCL8",
                                   "IL4", "IL5", "IL13",
                                   "TNFRSF1B", "TNFRSF4",
                                   "IL2RG"))  

svg(filename = "BALFTNKBubble.svg", width = 5, height = 5)
BubblePlot(TNK.BALF.integrated[,TNK.BALF.integrated$FineGranularity %in% AnalyseSubTypes], Features = DiffImmuneMarkers,
           CellTypes = AnalyseSubTypes, ColData = "FineGranularity", DiffRes = DEG, colorlow = "#3587A4", colorhigh = "#DAE3FE", colorbarwidth = 6)
dev.off()

svg(filename = "CD28.svg", width = 8, height = 4)
FeaturePlot(TNK.BALF.integrated, features = "CD28", split.by = "Group")
dev.off()
