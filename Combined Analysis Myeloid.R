################# Combined Analysis M cells Preprocessing 1#####################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Preprocessing/")
# library(Seurat)
# library(dplyr)
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/ReannotationAndElimination/SeuratPBMC_Annotated_Myeloid.Rdata")
# 
# M.RNA <- DietSeurat(Myeloid.RNA, assays = "RNA")
# 
# M.list <- SplitObject(M.RNA, split.by = "orig.ident")
# M.list <- lapply(M.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = T) # This step was intended to display the original normalized counts in subsequent visualization (SCTransform does not normalize the original RNA assay)
#   x <- SCTransform(x, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = T)
# })
# 
# integ_features <- SelectIntegrationFeatures(object.list = M.list,
#                                             nfeatures = 3000)
# 
# M.list <- PrepSCTIntegration(object.list = M.list,
#                                       anchor.features = integ_features)
# 
# integ_anchors <- FindIntegrationAnchors(object.list = M.list,
#                                         normalization.method = "SCT",
#                                         anchor.features = integ_features)
# 
# M.RNA <- IntegrateData(anchorset = integ_anchors,
#                          normalization.method = "SCT")
# 
# M.RNA <- RunPCA(M.RNA)
# 
# save(M.RNA, file = "M.RNA.Rdata")
# 
################### Combined Analysis M cells Preprocessing2 ###################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Preprocessing/")
# library(Seurat)
# library(dplyr)
# 
# load("M.RNA.Rdata")
# 
# svg(filename = paste("Elbow_Plot.svg", sep = ""), width = 7, height = 5)
# ElbowPlot(M.RNA, ndims = 50) # to determine number of dimensions for clustering
# dev.off()
# 
# pc.use <- 30
# M.RNA <- FindNeighbors(M.RNA, dims = 1:pc.use, verbose = T)
# M.RNA <- FindClusters(M.RNA, resolution = 1.5)
# M.RNA <- RunUMAP(M.RNA, dims = 1:pc.use)
# 
# DimPlot(M.RNA, label = T)
# 
# DefaultAssay(M.RNA) <- "RNA"
# Cluster.markers.m <- FindAllMarkers(M.RNA, only.pos = F, test.use = "negbinom", latent.vars = "orig.ident")
# 
# save(M.RNA, Cluster.markers.m, file = "Myeloid_reclustering.Rdata")

################## Combined Analysis Myeloid cells Annotation ##################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/")
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(emojifont)
# library(ComplexHeatmap)
# library(VennDiagram)
# library(org.Hs.eg.db)
# suppressMessages(library(ArchR))
# addArchRThreads(threads = 16)
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Preprocessing/Myeloid_reclustering.Rdata")
# 
# Cluster.markers.m$Delta <- Cluster.markers.m$pct.1 - Cluster.markers.m$pct.2
# top10 <- Cluster.markers.m %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_log2FC)
# 
# M.RNA$FineGranularity <- "CD14+ mono"
# M.RNA$IntermediateGranularity <- ""
# 
# M.RNA$FineGranularity[M.RNA$seurat_clusters %in% c("10", "11")] <- "CD16+ mono" #
# M.RNA$FineGranularity[M.RNA$seurat_clusters %in% c("16")] <- "CD3E+ mono" #
# M.RNA$FineGranularity[M.RNA$seurat_clusters %in% c("17")] <- "pDC" # LILRA4
# M.RNA$FineGranularity[M.RNA$seurat_clusters %in% c("7")] <- "IFN-responsive mono" # IFIT1, IFIT2, IFIT3
# M.RNA$FineGranularity[M.RNA$seurat_clusters %in% c("14")] <- "Dendritic" # CD1C
# M.RNA$FineGranularity[M.RNA$seurat_clusters %in% c("6")] <- "Intermediate mono" # CD14, FCGR3A
# M.RNA$FineGranularity[M.RNA$seurat_clusters %in% c("5")] <- "CXCL8+ mono"
# M.RNA$FineGranularity[M.RNA$seurat_clusters %in% c("12")] <- "IL1B+ mono"
# 
# M.RNA$IntermediateGranularity <- M.RNA$FineGranularity
# 
# Cluster2Unicode <- readxl::read_xlsx(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Dec2UnicodeSymbol.xlsx")
# McellsPalette <- c('CD14+ mono' = "#C9EF75",
#                    'CD16+ mono' = "#73CF9F",
#                    'CD3E+ mono' = "#82CC94",
# 
#                    'pDC' = "#BBF785",
#                    'IFN-responsive mono' = "#B0D955",
#                    'Dendritic' = "#93E9B6",
# 
#                    'Intermediate mono' = "#BCCC82",
#                    'CXCL8+ mono' = "#1EC048",
#                    'IL1B+ mono' = "#72B701"
# )
# 
# NumCluster <- length(table(M.RNA$seurat_clusters)) - 1
# 
# M.RNA$IntermediateGranularity[M.RNA$FineGranularity %in% c("CD14+ mono", "IFN-responsive mono", "Intermediate mono",
#                                                            "CXCL8+ mono", "IL1B+ mono")] <- "CD14+ mono"
# 
# Reanno.cluster.ids <- M.RNA$FineGranularity[match(as.character(0:NumCluster), table = M.RNA$seurat_clusters)]
# names(Reanno.cluster.ids) <- levels(M.RNA@active.ident)
# M.RNA <- RenameIdents(M.RNA, Reanno.cluster.ids)
# 
# DimPlot(M.RNA)
# 
# svg(filename = "Myeloid_IntermediateGranularity.svg", width = 9.5, height = 7.5)
# tmp <- FetchData(M.RNA, vars = c("UMAP_1", "UMAP_2", "orig.ident", "IntermediateGranularity"))
# 
# tmp$Labels <- match(tmp$IntermediateGranularity, names(McellsPalette))
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = IntermediateGranularity)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = McellsPalette[names(McellsPalette) %in% tmp$IntermediateGranularity]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rank(match(rownames(LabelPositions), names(McellsPalette))),
#            color = LabelPositions$color, size = 7,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
#                         legend.key.size = unit(2.5, 'cm'),
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(override.aes = list(size = 9,
#                                                   shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$IntermediateGranularity))]
#   )))
# dev.off()
# 
# svg(filename = "MyeloidFineGranularity.svg", width = 9.5, height = 7.5)
# tmp <- FetchData(M.RNA, vars = c("UMAP_1", "UMAP_2", "orig.ident", "FineGranularity"))
# 
# tmp$Labels <- match(tmp$FineGranularity, names(McellsPalette))
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = FineGranularity)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = McellsPalette[names(McellsPalette) %in% tmp$FineGranularity]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rownames(LabelPositions),
#            color = LabelPositions$color, size = 7,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
#                         legend.key.size = unit(1.8, 'cm'),
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(override.aes = list(size = 9,
#                                                   shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$FineGranularity))]
#   )))
# dev.off()
# 
# #---------------------> Seurat-Style stacked violin plot <----------------------
# svg(filename = paste("Myeloid_FineGranularity_Violin.svg", sep = ""), width = 12, height = 9)
# Features <- c("CD14", "FCGR3A", "CD3E", "LILRA4", "IFIT1", "IFIT2", "IFIT3", "CD1C", "CXCL8", "IL1B", "FCGR3A")
# 
# long <- FetchData(M.RNA, vars = Features, slot = "data")
# long$FineGranularity <- M.RNA$FineGranularity
# long$Cell <- rownames(long)
# long <- reshape2::melt(long, id.vars = c("Cell","FineGranularity"), measure.vars = Features,
#                        variable.name = "Feat", value.name = "Expr")
# 
# long$FineGranularity <- factor(x = long$FineGranularity, levels = rev(names(McellsPalette)))
# long$Color <- McellsPalette[match(long$FineGranularity, table = names(McellsPalette))]
# 
# ggplot(long, aes(Expr, factor(FineGranularity))) +
#   geom_violin(scale = "width", adjust = 1, trim = TRUE, aes(fill = FineGranularity, color = FineGranularity)) +
#   scale_fill_manual("legend", values = McellsPalette) +
#   scale_color_manual("legend", values = McellsPalette) +
#   scale_x_continuous(expand = c(0, 0), labels = function(x)
#     c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
#   facet_grid(cols = vars(Feat), scales = "free")  +
#   cowplot::theme_cowplot(font_size = 20) +
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
# # ATAC
# 
# M.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/AnnotatedATAC/Myeloid/"))
# 
# M.ATAC <- addIterativeLSI(ArchRProj = M.ATAC, force = T)
# 
# M.ATAC <- addClusters(
#   input = M.ATAC,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "Clusters",
#   resolution = 0.8
# )
# 
# # The following tmp object was only used for peak calling,
# # wherein each cluster contains at least 200 cells
# tmp <- M.ATAC[M.ATAC$Clusters %in% names(which(table(M.ATAC$Clusters) > 200)),]
# tmp <- addGroupCoverages(ArchRProj = tmp, groupBy = "Clusters")
# tmp <- addReproduciblePeakSet(
#   ArchRProj = tmp, reproducibility = "1",
#   groupBy = "Clusters", pathToMacs2 = "/mnt/TRAINING/Zhongjiacheng/anaconda2/bin/macs2"
# )
# 
# M.ATAC <- addPeakSet(M.ATAC, peakSet = tmp@peakSet, force = T)
# M.ATAC <- addPeakMatrix(M.ATAC)
# 
# M.ATAC <- addIterativeLSI(
#   ArchRProj = M.ATAC,
#   useMatrix = "PeakMatrix",
#   name = "PeakIterativeLSI",
#   force = T
# )
# 
# M.ATAC <- addHarmony(
#   ArchRProj = M.ATAC,
#   reducedDims = "PeakIterativeLSI",
#   name = "Harmony",
#   groupBy = "Sample", force = T
# )
# 
# M.ATAC <- addUMAP(
#   ArchRProj = M.ATAC,
#   reducedDims = "Harmony",
#   name = "HarmonyUMAP",
#   nNeighbors = 30,
#   minDist = 0.5,
#   metric = "cosine",
#   seed = 2, force = T
# )
# 
# groupList <- SimpleList()
# for (i in c("Control1", "Control2", "Patient")) {
#   RNA.cells <- colnames(M.RNA[,grep(M.RNA$orig.ident, pattern = i)])
#   ATAC.cells <- grep(M.ATAC$cellNames, pattern = i, value = T)
# 
#   groupList[[i]] <- SimpleList(
#     ATAC = ATAC.cells,
#     RNA = RNA.cells
#   )
# }
# 
# M.ATAC <- addGeneIntegrationMatrix(
#   ArchRProj = M.ATAC,
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "PeakIterativeLSI",
#   seRNA = M.RNA,
#   addToArrow = T,
#   groupList = groupList,
#   groupRNA = "IntermediateGranularity",
#   nameCell = "predictedCell_Co",
#   nameGroup = "predictedGroup_Co",
#   nameScore = "predictedScore_Co",
#   force = T, useImputation = F
# )
# 
# png(filename =  "Constrained_Myeloid_Cell_Type_Transfer_Score.png", width = 2000, height = 2000, res = 300)
# print(plotEmbedding(
#   M.ATAC,
#   colorBy = "cellColData",
#   name = "predictedScore_Co",
#   embedding = "HarmonyUMAP"
# ))
# dev.off()
# 
# M.ATAC <- M.ATAC[M.ATAC$predictedScore_Co > 0.5,]
# 
# svg(filename =  "Constrained_Myeloid_Cell_Type_Transfer_Score_Filtered.svg", width = 6, height = 6)
# print(plotEmbedding(
#   M.ATAC,
#   colorBy = "cellColData",
#   name = "predictedScore_Co",
#   embedding = "HarmonyUMAP"
# ))
# dev.off()
# 
# # The following step filters out cell types that were not shared across all samples
# 
# ValidPseudoBulk <- Reduce(intersect, list(names(table(M.ATAC$predictedGroup_Co[M.ATAC$Sample == "Control1-scATAC"])),
#                                           names(table(M.ATAC$predictedGroup_Co[M.ATAC$Sample == "Control2-scATAC"])),
#                                           names(table(M.ATAC$predictedGroup_Co[M.ATAC$Sample == "Patient-scATAC"]))))
# 
# # The following step filters out cell types with number less than 40 in at least one sample
# 
# ValidPseudoBulk <- Reduce(intersect, list(names(which(table(M.ATAC$predictedGroup_Co[M.ATAC$Sample == "Control1-scATAC" & M.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40)),
#                                           names(which(table(M.ATAC$predictedGroup_Co[M.ATAC$Sample == "Control2-scATAC" & M.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40)),
#                                           names(which(table(M.ATAC$predictedGroup_Co[M.ATAC$Sample == "Patient-scATAC" & M.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40))))
# 
# M.ATAC <- M.ATAC[M.ATAC$predictedGroup_Co %in% ValidPseudoBulk,]
# 
# svg(filename = "Constrained_Myeloid_Cell_Type_Transfer_Filtered_Group.svg", width = 8, height = 9)
# tmp <- getEmbedding(ArchRProj = M.ATAC, embedding = "HarmonyUMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$predictedGroup_Co <- M.ATAC$predictedGroup_Co
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[grep(rownames(tmp), pattern = "Control"),"UMAP_1"]), INDEX = tmp[grep(rownames(tmp), pattern = "Control"),"predictedGroup_Co"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[grep(rownames(tmp), pattern = "Control"),"UMAP_2"]), INDEX = tmp[grep(rownames(tmp), pattern = "Control"),"predictedGroup_Co"], FUN = mean)))
# 
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = predictedGroup_Co)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = McellsPalette[names(McellsPalette) %in% tmp$predictedGroup_Co]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rank(match(rownames(LabelPositions), names(McellsPalette))),
#            color = LabelPositions$color, size = 7,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "bottom", legend.title = element_blank(),
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(nrow = 2, override.aes = list(size = 9,
#                                                   shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$predictedGroup_Co))]
#   )))
# dev.off()
# 
# svg(filename = "MyeloidATAC_Reclustering_Colored_Samples.svg", width = 8, height = 8.5)
# tmp <- getEmbedding(ArchRProj = M.ATAC, embedding = "HarmonyUMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$Sample <- M.ATAC$Sample
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
# # # Double-check canonical cell marker
# #
# markerGenes  <- c(
#   "CD14",
#   "FCGR3A",
#   "LYZ",
#   "CD3E"
# )
# 
# M.ATAC <- addImputeWeights(M.ATAC)
# 
# p <- plotEmbedding(
#   ArchRProj = M.ATAC,
#   colorBy = "GeneScoreMatrix",
#   name = markerGenes,
#   continuousSet = "horizonExtra",
#   embedding = "HarmonyUMAP",
#   imputeWeights = getImputeWeights(M.ATAC)
# )
# 
# for (i in markerGenes) {
#   png(filename = paste("CheckMarkerLabelTransfer/", i, "embedding.png", sep = ""),  width = 2000, height = 2000, res = 300)
#   print(p[i])
#   dev.off()
# }
# 
# M.ATAC <- addMotifAnnotations(ArchRProj = M.ATAC, motifSet = "cisbp", name = "Motif")
# 
# M.ATAC <- addBgdPeaks(M.ATAC)
# M.ATAC <- addDeviationsMatrix(
#   ArchRProj = M.ATAC,
#   peakAnnotation = "Motif",
#   force = T
# )
# 
# M.ATAC$Group <- ""
# M.ATAC$Group[grep(M.ATAC$Sample, pattern = "Control")] <- "Control"
# M.ATAC$Group[grep(M.ATAC$Sample, pattern = "Patient")] <- "Patient"
# 
# M.ATAC$GroupSpecificPseudoBulk <- paste(M.ATAC$Group, M.ATAC$predictedGroup_Co)
# M.ATAC <- addGroupCoverages(ArchRProj = M.ATAC, groupBy = "GroupSpecificPseudoBulk")
# 
# M.RNA$Group <- ""
# M.RNA$Group[grep(M.RNA$orig.ident, pattern = "Control")] <- "Control"
# M.RNA$Group[M.RNA$orig.ident %in% c("PatientPBMC")] <- "Patient"
# 
# M.RNA$GroupCluster <- paste(M.RNA$Group, M.RNA$IntermediateGranularity, sep = " ")
# 
# saveArchRProject(ArchRProj = M.ATAC, outputDirectory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/ATAC/", load = F, overwrite = T)
# save(M.RNA, file = "AnnotatedMyeloid.Rdata")

####################### Differential Analyses Myeloid ##########################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/")
# library(clusterProfiler)
# library(Seurat)
# library(org.Hs.eg.db)
# suppressMessages(library(ArchR))
# 
# M.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/ATAC/"))
# load(file = "AnnotatedMyeloid.Rdata")
# 
# ##### DEG Differential Expressed Genes
# 
# DEG <- list()
# 
# for (CellType in c(names(table(M.ATAC$predictedGroup_Co)), "Dendritic")) {
#   FirstGroup <- paste("Patient", CellType, sep = " ")
#   SecondGroup <- paste("Control", CellType, sep = " ")
# 
#   DEGPv1To2 <- FindMarkers(object = M.RNA[,M.RNA$orig.ident %in% c("PatientPBMC",
#                                                                        "Control1PBMC",
#                                                                        "Control2PBMC")],
#                            ident.1 = FirstGroup,
#                            ident.2 = SecondGroup,
#                            group.by = "GroupCluster",
#                            logfc.threshold = 0.25,
#                            test.use = "MAST")
# 
#   DEGPv3To5 <- FindMarkers(object = M.RNA[,M.RNA$orig.ident %in% c("PatientPBMC",
#                                                                        "Control3PBMC",
#                                                                        "Control4PBMC",
#                                                                        "Control5PBMC")],
#                            ident.1 = FirstGroup,
#                            ident.2 = SecondGroup,
#                            group.by = "GroupCluster",
#                            logfc.threshold = 0.25,
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
# ##### DAP Differential Accessibility Regions/DET Differential Enrichment of TF Motif
# # Fix the mColSums issue in ArchR 1.0.1 (https://github.com/GreenleafLab/ArchR/issues/562)
# source("/mnt/TRAINING/Zhongjiacheng/MarkerSCFix.R")
# 
# DAP <- list()
# DET <- list()
# 
# for (CellType in names(table(M.ATAC$predictedGroup_Co))) {
#   FirstGroup <- paste("Patient", CellType)
#   SecondGroup <- paste("Control", CellType)
# 
#   CutOffUp <- "FDR < 0.05 & Log2FC > 0.25"
#   CutOffDown <- "FDR < 0.05 & Log2FC < -0.25"
# 
#   AllPeak <- getMarkerFeatures(
#     ArchRProj = M.ATAC,
#     useMatrix = "PeakMatrix",
#     groupBy = "GroupSpecificPseudoBulk",
#     testMethod = "wilcoxon",
#     bias = c("TSSEnrichment", "log10(nFrags)"),
#     useGroups = FirstGroup,
#     bgdGroups = SecondGroup)
# 
#   DAP[[CellType]] <- GRangesList(Up = getMarkers(AllPeak, cutOff = CutOffUp, returnGR = TRUE)@listData[[1]],
#                                             Down = getMarkers(AllPeak, cutOff = CutOffDown, returnGR = TRUE)@listData[[1]])
# 
#   DET[[CellType]]$motifsUp <- peakAnnoEnrichment(
#     seMarker = AllPeak,
#     ArchRProj = M.ATAC,
#     peakAnnotation = "Motif",
#     cutOff = CutOffUp
#   )
# 
#   DET[[CellType]]$motifsDown <- peakAnnoEnrichment(
#     seMarker = AllPeak,
#     ArchRProj = M.ATAC,
#     peakAnnotation = "Motif",
#     cutOff = CutOffDown
#   )
# }
# 
# GSEA.res <- List()
# for (CellType in c(names(table(M.ATAC$predictedGroup_Co)), "Dendritic")) {
#   GSEA_input <- DEG[[CellType]]$avg_log2FC
#   names(GSEA_input) <- rownames(DEG[[CellType]])
#   GSEA_input <- sort(GSEA_input, decreasing = TRUE)
#   GSEA.res[[CellType]] <- gseGO(GSEA_input, OrgDb = org.Hs.eg.db,  keyType = "SYMBOL", ont = "BP",
#                                minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 0.1)
# }
# 
# save(DEG, DAP, DET, GSEA.res, file = "Differential_Test_Results.Rdata")
# 
# # Get necessary genomic information: GeneBody, TSS, Promoter and Motif Position
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/")
# 
# library(org.Hs.eg.db)
# library(annotate)
# library(GenomicFeatures)
# library(EnsDb.Hsapiens.v86) # EnsDb.Hsapiens.v86 for human gene definitions of the Ensembl code database version 86 that based on the GRCh38 genome
# suppressMessages(library(ArchR))
# 
# M.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/ATAC/"))
# 
# TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# motifPositions <- getPositions(M.ATAC)
# GeneBody <- getGeneAnnotation(M.ATAC)$genes
# TSS <- compEpiTools::TSS(TxDb)
# 
# # The result uses "." followed by a number to indicate different transcript,
# # We refer these transcripts back to the original transcript ID by eliminating
# # information after "."
# 
# TSS@elementMetadata$TranscriptID <- unlist(lapply(strsplit(TSS$tx_name, split = "\\."), '[[', 1))
# 
# # Covert transcript ID to symbols, and eliminate transcripts without matched symbols
# 
# symbols <- select(EnsDb.Hsapiens.v86, columns = "GENENAME", keytype = "TXID", keys = TSS@elementMetadata$TranscriptID)
# TSS@elementMetadata$symbol <- symbols$GENENAME[match(TSS@elementMetadata$TranscriptID, table = symbols$TXID)]
# TSS <- TSS[!is.na(TSS@elementMetadata$symbol)]
# 
# Promoter <- promoters(genes(TxDb), upstream = 1500, downstream = 500)
# Promoter@elementMetadata$symbol <- unlist(lookUp(unlist(Promoter@elementMetadata$gene_id),'org.Hs.eg','SYMBOL'))
# Promoter@elementMetadata <- Promoter@elementMetadata[-1] # Remove redundant information such as gene_id etc
# Promoter <- Promoter[!is.na(Promoter$symbol)]
# 
# Promoter <- Promoter[Promoter@seqnames %in% GeneBody@seqnames]
# TSS <- TSS[TSS@seqnames %in% GeneBody@seqnames]
# 
# save(TSS, Promoter, GeneBody, motifPositions, file = "GenomicInformation.Rdata")

######################### DownStream Analysis Myeloid ##########################
rm(list = ls())
setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/")

library(dplyr)
library(ggrepel)
library(Seurat)
library(corrplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)
library(ggalluvial)
library(ComplexHeatmap)
library(cowplot)
library(monocle3)
library(SeuratWrappers)
library(colorspace)
suppressMessages(library(ArchR))
addArchRThreads(threads = 16)
source("/mnt/TRAINING/Zhongjiacheng/Customized Codes.R")
source("/mnt/TRAINING/Zhongjiacheng/Step4_functions.R")

M.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/ATAC/"))
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/AnnotatedMyeloid.Rdata")
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/Differential_Test_Results.Rdata")
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/Annotation/GenomicInformation.Rdata")

McellsPalette <- c('CD14+ mono' = "#C9EF75",
                   'CD16+ mono' = "#73CF9F",
                   'CD3E+ mono' = "#82CC94",

                   'pDC' = "#BBF785",
                   'IFN-responsive mono' = "#B0D955",
                   'Dendritic' = "#93E9B6",

                   'Intermediate mono' = "#BCCC82",
                   'CXCL8+ mono' = "#1EC048",
                   'IL1B+ mono' = "#72B701"
)

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


ImmuneMarkers <- grep(paste(c("^CCL", "^CXCL", "^CXCR", "^CCR", "TNF",
                              "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
                              "^TGF", "^IFN", "^STAT", "^DEFB", 
                              "^LYZ", "S100A8", "S100A9", "^CAMP", "^LCN2", "^LTF", "CD86", "CD80",
                              paste("^IL", c(1:50), sep = "")), collapse = "|"), 
                      rownames(M.RNA), value = T)

DiffImmuneMarkers <- list()

# # For visualization, a log2FC curoff of 0.4 was set

for (CellType in c(names(table(M.ATAC$predictedGroup_Co)), "Dendritic")) {
  DiffImmuneMarkers[[CellType]] <- rownames(DEG[[CellType]][rownames(DEG[[CellType]]) %in% c(PotentialBioMarkers, ImmuneMarkers) &
                                                              abs(DEG[[CellType]]$avg_log2FC) > 0.4,])
}

DiffImmuneMarkers <- Reduce(union, list(DiffImmuneMarkers))
DiffImmuneMarkers <- unique(unlist(DiffImmuneMarkers))

levels(DiffImmuneMarkers) <- rev(c("HLA-A", "HLA-B", "HLA-C", "HLA-E",
                                   "STAT1", "STAT2", "STAT3", "STAT5B", "STAT6",
                                   "IFNAR2", "IFNGR1", "CXCL16", "CXCR4",
                                   "IL1B", "IL15", "IL6R", "IL17RA", "IL13RA1", "IL6ST",
                                   "TGFBR2", "TNF", "TNFAIP2","TNFAIP3", "TNFAIP8",
                                   "TNFRSF1B", "TNFSF10", "TNFSF13", "TNFRSF14",
                                   "S100A8"))

AnalyseSubTypes <- c(names(table(M.ATAC$predictedGroup_Co)), "Dendritic")
levels(AnalyseSubTypes) <- c("CD14+ mono", "CD16+ mono", "Dendritic")

svg(filename = "MyeloidBubble.svg", width = 3.5, height = 6.5)
BubblePlot(M.RNA[,M.RNA$orig.ident %in% c("Control1PBMC", "Control2PBMC", "PatientPBMC")], Features = DiffImmuneMarkers,
           CellTypes = AnalyseSubTypes, ColData = "IntermediateGranularity", DiffRes = DEG, colorlow = "#71BC16", colorhigh = "#F6FFD9")
dev.off()

svg(filename = "RelativeFraction_Myeloid.svg", width = 6, height = 6)
RelativeFraction(M.RNA, ColData = "FineGranularity")
dev.off()

svg(filename = paste("GSEA/", "GSEACD14.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$`CD14+ mono`, Descriptions = c("oxidative phosphorylation", 
                                                     "positive regulation of nitric oxide biosynthetic process",
                                                     "Fc-gamma receptor signaling pathway involved in phagocytosis"),
             base_size = 15,
             legendposition = c(0.45, 0.3), 
             rel_heights = c(1.5, 0.5, 1),
             GSEAColor = c("#8CE101", "#489269", "#9AF8C5"), maxoverlaps = 60
)
dev.off()

svg(filename = paste("GSEA/", "GSEACD16.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$`CD16+ mono`, Descriptions = c("oxidative phosphorylation", 
                                                     "positive regulation of nitric-oxide synthase biosynthetic process",
                                                     "Fc-gamma receptor signaling pathway involved in phagocytosis"),
             base_size = 15,
             legendposition = c(0.45, 0.3), 
             rel_heights = c(1.5, 0.5, 1),
             GSEAColor = c("#8CE101", "#489269", "#9AF8C5"), maxoverlaps = 60
)
dev.off()

GSEACD14nit <- CoreEnrichment(GSEA.res$`CD14+ mono`@result$core_enrichment[GSEA.res$`CD14+ mono`@result$Description %in% "positive regulation of nitric oxide biosynthetic process"])

GSEACD16nit <- CoreEnrichment(GSEA.res$`CD16+ mono`@result$core_enrichment[GSEA.res$`CD16+ mono`@result$Description %in%
                                                                             "positive regulation of nitric-oxide synthase biosynthetic process"])

GSEACD14oxi <- CoreEnrichment(GSEA.res$`CD14+ mono`@result$core_enrichment[GSEA.res$`CD14+ mono`@result$Description %in% "oxidative phosphorylation"])

GSEACD16oxi <- CoreEnrichment(GSEA.res$`CD16+ mono`@result$core_enrichment[GSEA.res$`CD16+ mono`@result$Description %in% "oxidative phosphorylation"])


GSEACD14pha <- CoreEnrichment(GSEA.res$`CD14+ mono`@result$core_enrichment[GSEA.res$`CD14+ mono`@result$Description %in% "immune response-regulating cell surface receptor signaling pathway involved in phagocytosis"])

GSEACD16pha <- CoreEnrichment(GSEA.res$`CD16+ mono`@result$core_enrichment[GSEA.res$`CD16+ mono`@result$Description %in% c("Fc-gamma receptor signaling pathway involved in phagocytosis")])


GSEACD14MHC <- CoreEnrichment(GSEA.res$`CD14+ mono`@result$core_enrichment[GSEA.res$`CD14+ mono`@result$Description %in% c("antigen processing and presentation of peptide antigen via MHC class I",
                                                                                                                           "antigen processing and presentation of peptide antigen via MHC class II")])
GSEACD16MHC <- CoreEnrichment(GSEA.res$`CD16+ mono`@result$core_enrichment[GSEA.res$`CD16+ mono`@result$Description %in% c("antigen processing and presentation of peptide antigen via MHC class I",
                                                                                                                           "antigen processing and presentation of peptide antigen via MHC class II")])
################################## ATAC ########################################
ControlPalette <- McellsPalette[names(McellsPalette) %in% M.ATAC$predictedGroup_Co]
cols1 <- readhex(file = textConnection(paste(ControlPalette, collapse = "\n")),
                 class = "RGB")

#transform to hue/lightness/saturation colorspace

cols1 <- as(cols1, "HLS")
PatientPalette <- cols1 -> ControlPalette

#multiplicative decrease of lightness

PatientPalette@coords[, "L"] <- PatientPalette@coords[, "L"] * 0.75
ControlPalette@coords[, "L"] <- ControlPalette@coords[, "L"] * 1.25

#going via rgb seems to work better

cols1 <- as(cols1, "RGB") %>% hex()
PatientPalette <- as(PatientPalette, "RGB") %>% hex()
ControlPalette <- as(ControlPalette, "RGB") %>% hex()

svg(filename = "ATACGroupLegendKeys.svg", width = 7, height = 5)
plot.new()
plot(x = seq_along(ControlPalette), y = rep(1, length(ControlPalette)),
     col = ControlPalette, pch = 16, ylim = c(0, 4.5), cex = 5,
     xlab = "", ylab = "")
points(x = seq_along(PatientPalette), y = rep(2, length(PatientPalette)),
       col = PatientPalette, pch = 16, cex = 5)
dev.off()

names(PatientPalette) <- names(McellsPalette[names(McellsPalette) %in% M.ATAC$predictedGroup_Co]) -> names(ControlPalette)
names(PatientPalette) <- paste("Patient", names(PatientPalette))
names(ControlPalette) <- paste("Control", names(ControlPalette))
ATACPalette <- c(PatientPalette, ControlPalette)

orderC <- c("Control CD14+ mono",
            "Control CD16+ mono")

orderP <- c("Patient CD14+ mono",
            "Patient CD16+ mono")

################################################################################
Patient.M.ATAC <- M.ATAC[M.ATAC$Sample == "Patient-scATAC",]
Control.M.ATAC <- M.ATAC[!M.ATAC$Sample == "Patient-scATAC",]

# Create a MotifMatrix of chromVAR deviations and deviation z-scores for all motifs. 
# This data, averaged by clusters, was obtained by using the getGroupSE() function 
# which returns a SummarizedExperiment. 

seGroupMotif <- getGroupSE(ArchRProj = M.ATAC, useMatrix = "MotifMatrix", groupBy = "GroupSpecificPseudoBulk")

# This SummarizedExperiment was further subset to just the deviation z-scores.

DeviationZ <- seGroupMotif[rowData(seGroupMotif)$seqnames == "z",]

for (CellType in names(table(M.ATAC$predictedGroup_Co))) {
  svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/EnrichedMotif/Up-", CellType, ".svg", sep = ""), width = 6, height = 7)
  print(PlotEnrichedTFMotif(DET[[CellType]]$motifsUp, dotsize = 3, labelsize = 4, repulsion = 25, DeviationZ = DeviationZ, UpIn = "Patient"))
  dev.off()
  
  svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/EnrichedMotif/Down-", CellType, ".svg", sep = ""), width = 6, height = 7)
  print(PlotEnrichedTFMotif(DET[[CellType]]$motifsDown, dotsize = 3, labelsize = 4, repulsion = 25, DeviationZ = DeviationZ, UpIn = "Control"))
  dev.off()
  
  DET[[CellType]]$motifsUpTop <- GetRepresentativeTF(DET[[CellType]]$motifsUp, topTF = 20, DeviationZ = DeviationZ, UpIn = "Patient")
  DET[[CellType]]$motifsDownTop <- GetRepresentativeTF(DET[[CellType]]$motifsDown, topTF = 20, DeviationZ = DeviationZ, UpIn = "Control")
}

AllSigTFMotif <- list()
for (CellType in names(table(M.ATAC$predictedGroup_Co))) {
  AllSigTFMotif[[paste(CellType, "Up")]] <- names(DET[[CellType]]$motifsUpTop)
  AllSigTFMotif[[paste(CellType, "Down")]] <- names(DET[[CellType]]$motifsDownTop)
}

UnionSigTFMotif <- Reduce(union, AllSigTFMotif)

IntersectSigTFMotifUp <- Reduce(intersect, AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Up")])
IntersectSigTFMotifDown <- Reduce(intersect, AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Down")])

svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/EnrichedMotif/", "DownUpset.svg", sep = ""), width = 5, height = 3)
UpSet(make_comb_mat(AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Down")]), 
      comb_col = "#C0E2FC", bg_col = "#F0F0FF", bg_pt_col = "#1C97F4")
dev.off()

svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/EnrichedMotif/", "UUppset.svg", sep = ""), width = 5, height = 3)
UpSet(make_comb_mat(AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Up")]), 
      comb_col = "#FCC8C4", bg_col = "#FFF1F3", bg_pt_col = "#F52B1B")
dev.off()

################# Associate CAR and TF and conversion to bam ###################
# Fragment_By_CellType <- list()
# for (CellType in names(table(M.ATAC$GroupSpecificPseudoBulk))) {
#   SubSlot <- paste(unlist(strsplit(CellType, split = " ")), collapse = "_")
# 
#   Fragment_By_CellType[[SubSlot]] <- unlist(getFragmentsFromProject(M.ATAC, cellNames = rownames(M.ATAC[M.ATAC$GroupSpecificPseudoBulk == CellType,])))
#   Fragment_By_CellType[[SubSlot]]$index <- as.character(Fragment_By_CellType[[SubSlot]]$RG)
# }
# 
# library(GenomicRanges)
# 
# Output_to_bed_files(Fragment_list_cl = Fragment_By_CellType, folder = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/FragmentsByCellTypes/")
# 
# # Execute bedpetobam in /FragmentsByCellTypes, output to /BAM
# 
# bedpetobam <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   outdir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/BAM/"
#   bedpetobam <- c(bedpetobam, paste("bedtools bedpetobam -i ", i, "_fragments_cl_bamGR_pe.bed -g hg38.chrom.sizes.txt > ", outdir, i, "_fragments_cl_bamGR_pe.bam", sep = ""))
# }
# 
# bedpetobam <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/FragmentsByCellTypes/"), bedpetobam)
# 
# # Execute samtools_sort and samtools_index in /BAM, output to /sBAM
# 
# samtools_sort <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   outdir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/sBAM/"
#   samtools_sort <- c(samtools_sort, (paste("samtools sort -o ", outdir, i, "_fragments_cl_bamGR_pe_s.bam ", i, "_fragments_cl_bamGR_pe.bam", sep = "")))
# }
# 
# samtools_sort <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/BAM/"), samtools_sort)
# 
# # Execute samtools_sort and samtools_index in /sBAM, output to /sBAM
# 
# samtools_index <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   samtools_index <- c(samtools_index, paste("samtools index ", i, "_fragments_cl_bamGR_pe_s.bam", sep = ""))
# }
# 
# samtools_index <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/sBAM/"), samtools_index)
# 
# write.table(c(bedpetobam, samtools_sort, samtools_index), file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/fragments2bams.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

########################## ATAC signal correction ##############################
# # The precompiled version of the hg38 genome in ArchR uses BSgenome.Hsapiens.UCSC.hg38,
# # TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db,
# # and a blacklist that was merged using ArchR::mergeGR() from the hg38 v2 blacklist regions
# # and from mitochondrial regions that show high mappability to the hg38 nuclear genome from Caleb Lareau and Jason Buenrostro.
# # To set a global genome default to the precompiled hg38 genome:
# 
# library(rtracklayer)
# Hg38blacklist <- import("/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/BlackList/hg38-blacklist.v2.bed.gz", format = "bed")
# Hg38MitochondrialRegions <- read.table("/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/BlackList/hg38_peaks.narrowPeak.txt")
# Hg38MitochondrialRegions <- Hg38MitochondrialRegions[,c(1, 2, 3, 6)]
# colnames(Hg38MitochondrialRegions) <- c("chr", "start", "end", "strand")
# Hg38MitochondrialRegions <- makeGRangesFromDataFrame(Hg38MitochondrialRegions)
# 
# # reduce method will align the ranges and merge overlapping ranges to produce a simplified set.
# Hg38TotalBlackList <- reduce(c(Hg38blacklist, Hg38MitochondrialRegions))
# 
# export.bed(Hg38TotalBlackList, con = "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/BlackList/Hg38TotalBlackList.bed")
# export.bed(getPeakSet(M.ATAC), con = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ATACorrect/AllPeaks.bed")
# 
# ATACorrect <- c()
# for(i in names(Fragment_By_CellType)){
#   Input <- paste(i, "_fragments_cl_bamGR_pe_s.bam", sep = "")
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   Genome <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/hg38.fa"
#   BlackList <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/BlackList/Hg38TotalBlackList.bed"
#   OutPut <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ATACorrect/"
# 
#   if (!dir.exists(paste(OutPut, i, sep = ""))) {
#     dir.create(paste(OutPut, i, sep = ""))
#   }
# 
#   outdir <- paste(OutPut, i, sep = "")
# 
#   ATACorrect <- c(ATACorrect, paste("nohup", "TOBIAS ATACorrect --read_shift 0 0 --bam", Input, "--genome", Genome, "--peaks", AllPeaks, "--blacklist", BlackList, "--outdir", outdir, "--cores 16 &"))
# }
# 
# ATACorrect <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/sBAM/"), ATACorrect)
# 
# write.table(ATACorrect, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ATACorrectScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

############################### ScoreBigwig ####################################
# ScoreBigwig <- c()
# for(i in names(Fragment_By_CellType)){
#   Input <- paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ATACorrect/", i, "/", i, "_fragments_cl_bamGR_pe_s_corrected.bw", sep = "")
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   OutPut <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ScoreBigwig/"
#   OutFileName <- paste(i, "_footprints.bw", sep = "")
# 
#   ScoreBigwig <- c(ScoreBigwig, paste("nohup", "TOBIAS FootprintScores --signal", Input, "--regions", AllPeaks, "--output", OutFileName, "--cores 16 &"))
# }
# 
# ScoreBigwig <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ScoreBigwig/"), ScoreBigwig)
# 
# write.table(ScoreBigwig, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ScoreBigwigScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

##################### BINDetect differential TF occupancy ######################
# # If error occurs, enter ulimit -n 50000 in console
# BINDetect <- c()
# for(CellType in names(DEG)){
#   FirstGroup <- paste(c("Patient", paste(unlist(strsplit(CellType, split = " ")), collapse = "_")), collapse = "_")
#   SecondGroup <- paste(c("Control", paste(unlist(strsplit(CellType, split = " ")), collapse = "_")), collapse = "_")
# 
#   motifsMEME <- "/mnt/TRAINING/Zhongjiacheng/pwm_to_meme/Homo.cisbp.meme"
#   Genome <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/hg38.fa"
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   BigwigDir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/ScoreBigwig/"
#   OutDir <- paste(FirstGroup, SecondGroup, sep = "_vs_")
# 
#   BINDetect <- c(BINDetect, paste("TOBIAS BINDetect", "--motifs", motifsMEME, "--signals", paste(BigwigDir, FirstGroup, "_footprints.bw", sep = ""),
#                                   paste(BigwigDir, SecondGroup, "_footprints.bw", sep = ""), "--genome", Genome, "--peaks", AllPeaks, "--outdir"
#                                   , OutDir, "--cond_names", paste(FirstGroup, SecondGroup) , "--cores 32"))
# }
# 
# BINDetect <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/DiffTFOccupancy/"), BINDetect)
# 
# write.table(BINDetect, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/BINDetectScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

################ Visualization of differential TF binding ######################

ValidTFB <- names(table(Patient.M.ATAC$predictedGroup_Co))
ValidTFB <- unlist(lapply(strsplit(ValidTFB, split = " "), paste, collapse = "_"))

################ Subset the binding events to TF binding only ##################

DiffTFBind <- list()

for (CellType in ValidTFB) {
  tmp <- readxl::read_xlsx(paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/Myeloid/DownStream/TOBIAS/DiffTFOccupancy/", 
                                 "Patient_", CellType, "_vs_Control_", CellType, "/bindetect_results.xlsx", sep = ""), sheet = 2) %>% as.data.frame()
  tmp$clustername <- unlist(lapply(strsplit(tmp$cluster, split = "_"), '[[', 2))
  DiffTFBind[[CellType]] <- tmp[tmp$clustername %in% unlist(lapply(strsplit(names(motifPositions), split = "_"), '[[', 1)),]
}

svg(filename = paste("Volcano/", "CD14Vol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["CD14+_mono"]], 
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control CD14+ mono"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient CD14+ mono"])))
dev.off()

svg(filename = paste("Volcano/", "CD16Vol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["CD16+_mono"]], 
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control CD16+ mono"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient CD16+ mono"])))
dev.off()


suppressMessages(InspectEmbedding(M.ATAC, QueryGene = "ZNF524", 
                                  UseMatrix = getAvailableMatrices(M.ATAC)[3],
                                  Embed = "HarmonyUMAP"))


