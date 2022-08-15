################# Combined Analysis B cells Preprocessing 1#####################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Preprocessing/")
# library(Seurat)
# library(dplyr)
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/ReannotationAndElimination/SeuratPBMC_Annotated_B.Rdata")
# 
# B.RNA <- DietSeurat(B.RNA, assays = "RNA")
# 
# B.list <- SplitObject(B.RNA, split.by = "orig.ident")
# B.list <- lapply(B.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = T) # This step was intended to display the original normalized counts in subsequent visualization (SCTransform does not normalize the original RNA assay)
#   x <- SCTransform(x, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = T)
# })
# 
# integ_features <- SelectIntegrationFeatures(object.list = B.list,
#                                             nfeatures = 3000)
# 
# B.list <- PrepSCTIntegration(object.list = B.list,
#                                       anchor.features = integ_features)
# 
# integ_anchors <- FindIntegrationAnchors(object.list = B.list,
#                                         normalization.method = "SCT",
#                                         anchor.features = integ_features)
# 
# B.RNA <- IntegrateData(anchorset = integ_anchors,
#                          normalization.method = "SCT")
# 
# B.RNA <- RunPCA(B.RNA)
# 
# save(B.RNA, file = "B.RNA.Rdata")

################### Combined Analysis B cells Preprocessing2 ###################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Preprocessing/")
# library(Seurat)
# library(dplyr)
# 
# load("B.RNA.Rdata")
# 
# svg(filename = paste("Elbow_Plot.svg", sep = ""), width = 7, height = 5)
# ElbowPlot(B.RNA, ndims = 50) # to determine number of dimensions for clustering
# dev.off()
# 
# pc.use <- 32
# B.RNA <- FindNeighbors(B.RNA, dims = 1:pc.use, verbose = T)
# B.RNA <- FindClusters(B.RNA, resolution = 1.5)
# B.RNA <- RunUMAP(B.RNA, dims = 1:pc.use)
# 
# DimPlot(B.RNA, label = T)
# 
# DefaultAssay(B.RNA) <- "RNA"
# Cluster.markers.b <- FindAllMarkers(B.RNA, only.pos = F, test.use = "negbinom", latent.vars = "orig.ident")
# 
# save(B.RNA, Cluster.markers.b, file = "B_reclustering.Rdata")

#################### Combined Analysis B cells Annotation ######################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/")
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
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Preprocessing/B_reclustering.Rdata")
# 
# Cluster.markers.b$Delta <- Cluster.markers.b$pct.1 - Cluster.markers.b$pct.2
# top10 <- Cluster.markers.b %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_log2FC)
# 
# B.RNA$FineGranularity <- ""
# B.RNA$IntermediateGranularity <- ""
# 
# B.RNA$FineGranularity[B.RNA$seurat_clusters %in% c("10")] <- "Plasma" # MZB1
# B.RNA$FineGranularity[B.RNA$seurat_clusters %in% c("9")] <- "Atypical memory B" # FCRL5
# B.RNA$FineGranularity[B.RNA$seurat_clusters %in% c("1", "3", "4")] <- "Naive B" # IL4R
# B.RNA$FineGranularity[B.RNA$seurat_clusters %in% c("2")] <- "Memory B" # COCH, SSPN PMID: 32351704
# B.RNA$FineGranularity[B.RNA$seurat_clusters %in% c("8")] <- "CD3E+ B" # CD3E
# B.RNA$FineGranularity[B.RNA$seurat_clusters %in% c("11")] <- "LYZ+ B"
# B.RNA$FineGranularity[B.RNA$seurat_clusters %in% c("5", "0", "7", "6")] <- "B"
# 
# B.RNA$IntermediateGranularity <- B.RNA$FineGranularity
# 
# Cluster2Unicode <- readxl::read_xlsx(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Dec2UnicodeSymbol.xlsx")
# BcellsPalette <- c('Plasma' = "#FF9F9F",
#                    'Atypical memory B' = "#FFCB6D",
#                    'Naive B' = "#FCCFCC",
#                    'Memory B' = "#FF7979",
# 
#                    'CD3E+ B' = "#FDA9CD",
#                    'LYZ+ B' = "#BF6623",
#                    'B' = "#FFA771"
# )
# 
# NumCluster <- length(table(B.RNA$seurat_clusters)) - 1
# 
# B.RNA$IntermediateGranularity[B.RNA$FineGranularity %in% c("Memory B", "B", "LYZ+ B")] <- "B"
# 
# Reanno.cluster.ids <- B.RNA$FineGranularity[match(as.character(0:NumCluster), table = B.RNA$seurat_clusters)]
# names(Reanno.cluster.ids) <- levels(B.RNA@active.ident)
# B.RNA <- RenameIdents(B.RNA, Reanno.cluster.ids)
# 
# DimPlot(B.RNA)
# 
# svg(filename = "B_IntermediateGranularity.svg", width = 6.5, height = 4)
# tmp <- FetchData(B.RNA, vars = c("UMAP_1", "UMAP_2", "orig.ident", "IntermediateGranularity"))
# 
# tmp$Labels <- match(tmp$IntermediateGranularity, names(BcellsPalette))
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = IntermediateGranularity)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = BcellsPalette[names(BcellsPalette) %in% tmp$IntermediateGranularity]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rank(match(rownames(LabelPositions), names(BcellsPalette))),
#            color = LabelPositions$color, size = 7,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
#                         legend.key.size = unit(1.5, 'cm'),
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(override.aes = list(size = 9,
#                                                   shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$IntermediateGranularity))]
#   )))
# dev.off()
# 
# svg(filename = "B_FineGranularity.svg", width = 6.5, height = 4)
# tmp <- FetchData(B.RNA, vars = c("UMAP_1", "UMAP_2", "orig.ident", "FineGranularity"))
# 
# tmp$Labels <- match(tmp$FineGranularity, names(BcellsPalette))
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = FineGranularity)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = BcellsPalette[names(BcellsPalette) %in% tmp$FineGranularity]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rownames(LabelPositions),
#            color = LabelPositions$color, size = 7,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
#                         legend.text = element_text(size = 17),
#                         legend.key.size = unit(2.5, 'lines'),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(override.aes = list(size = 9,
#                                                   shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$FineGranularity))]
#   )))
# dev.off()
# 
# #---------------------> Seurat-Style stacked violin plot <----------------------
# svg(filename = paste("B_FineGranularity_Violin.svg", sep = ""), width = 9, height = 6)
# Features <- c("MZB1", "FCRL5", "IL4R", "IGHD", "COCH", "SSPN", "CD3E", "LYZ")
# 
# long <- FetchData(B.RNA, vars = Features, slot = "data")
# long$FineGranularity <- B.RNA$FineGranularity
# long$Cell <- rownames(long)
# long <- reshape2::melt(long, id.vars = c("Cell","FineGranularity"), measure.vars = Features,
#                        variable.name = "Feat", value.name = "Expr")
# 
# long$FineGranularity <- factor(x = long$FineGranularity, levels = rev(names(BcellsPalette)))
# long$Color <- BcellsPalette[match(long$FineGranularity, table = names(BcellsPalette))]
# 
# ggplot(long, aes(Expr, factor(FineGranularity))) +
#   geom_violin(scale = "width", adjust = 1, trim = TRUE, aes(fill = FineGranularity, color = FineGranularity)) +
#   scale_fill_manual("legend", values = BcellsPalette) +
#   scale_color_manual("legend", values = BcellsPalette) +
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
# B.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/AnnotatedATAC/B/"))
# 
# B.ATAC <- addIterativeLSI(ArchRProj = B.ATAC, force = T)
# 
# B.ATAC <- addClusters(
#   input = B.ATAC,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "Clusters",
#   resolution = 0.8
# )
# 
# # The following tmp object was only used for peak calling,
# # wherein each cluster contains at least 200 cells
# tmp <- B.ATAC[B.ATAC$Clusters %in% names(which(table(B.ATAC$Clusters) > 200)),]
# tmp <- addGroupCoverages(ArchRProj = tmp, groupBy = "Clusters")
# tmp <- addReproduciblePeakSet(
#   ArchRProj = tmp, reproducibility = "1",
#   groupBy = "Clusters", pathToMacs2 = "/mnt/TRAINING/Zhongjiacheng/anaconda2/bin/macs2"
# )
# 
# B.ATAC <- addPeakSet(B.ATAC, peakSet = tmp@peakSet, force = T)
# B.ATAC <- addPeakMatrix(B.ATAC)
# 
# B.ATAC <- addIterativeLSI(
#   ArchRProj = B.ATAC,
#   useMatrix = "PeakMatrix",
#   name = "PeakIterativeLSI",
#   force = T
# )
# 
# B.ATAC <- addHarmony(
#   ArchRProj = B.ATAC,
#   reducedDims = "PeakIterativeLSI",
#   name = "Harmony",
#   groupBy = "Sample", force = T
# )
# 
# B.ATAC <- addUMAP(
#   ArchRProj = B.ATAC,
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
#   RNA.cells <- colnames(B.RNA[,grep(B.RNA$orig.ident, pattern = i)])
#   ATAC.cells <- grep(B.ATAC$cellNames, pattern = i, value = T)
# 
#   groupList[[i]] <- SimpleList(
#     ATAC = ATAC.cells,
#     RNA = RNA.cells
#   )
# }
# 
# B.ATAC <- addGeneIntegrationMatrix(
#   ArchRProj = B.ATAC,
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "PeakIterativeLSI",
#   seRNA = B.RNA,
#   addToArrow = T,
#   groupList = groupList,
#   groupRNA = "IntermediateGranularity",
#   nameCell = "predictedCell_Co",
#   nameGroup = "predictedGroup_Co",
#   nameScore = "predictedScore_Co",
#   force = T, useImputation = F
# )
# 
# png(filename =  "Constrained_B_Cell_Type_Transfer_Score.png", width = 2000, height = 2000, res = 300)
# print(plotEmbedding(
#   B.ATAC,
#   colorBy = "cellColData",
#   name = "predictedScore_Co",
#   embedding = "HarmonyUMAP"
# ))
# dev.off()
# 
# B.ATAC <- B.ATAC[B.ATAC$predictedScore_Co > 0.5,]
# 
# svg(filename =  "Constrained_B_Cell_Type_Transfer_Score_Filtered.svg", width = 6, height = 6)
# print(plotEmbedding(
#   B.ATAC,
#   colorBy = "cellColData",
#   name = "predictedScore_Co",
#   embedding = "HarmonyUMAP"
# ))
# dev.off()
# 
# # The following step filters out cell types that were not shared across all samples
# 
# ValidPseudoBulk <- Reduce(intersect, list(names(table(B.ATAC$predictedGroup_Co[B.ATAC$Sample == "Control1-scATAC"])),
#                                           names(table(B.ATAC$predictedGroup_Co[B.ATAC$Sample == "Control2-scATAC"])),
#                                           names(table(B.ATAC$predictedGroup_Co[B.ATAC$Sample == "Patient-scATAC"]))))
# 
# # The following step filters out cell types with number less than 40 in at least one sample
# 
# ValidPseudoBulk <- Reduce(intersect, list(names(which(table(B.ATAC$predictedGroup_Co[B.ATAC$Sample == "Control1-scATAC" & B.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40)),
#                                           names(which(table(B.ATAC$predictedGroup_Co[B.ATAC$Sample == "Control2-scATAC" & B.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40)),
#                                           names(which(table(B.ATAC$predictedGroup_Co[B.ATAC$Sample == "Patient-scATAC" & B.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40))))
# 
# B.ATAC <- B.ATAC[B.ATAC$predictedGroup_Co %in% ValidPseudoBulk,]
# 
# svg(filename = "Constrained_B_Cell_Type_Transfer_Filtered_Group.svg", width = 4, height = 5)
# tmp <- getEmbedding(ArchRProj = B.ATAC, embedding = "HarmonyUMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$predictedGroup_Co <- B.ATAC$predictedGroup_Co
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[grep(rownames(tmp), pattern = "Control"),"UMAP_1"]), INDEX = tmp[grep(rownames(tmp), pattern = "Control"),"predictedGroup_Co"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[grep(rownames(tmp), pattern = "Control"),"UMAP_2"]), INDEX = tmp[grep(rownames(tmp), pattern = "Control"),"predictedGroup_Co"], FUN = mean)))
# 
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = predictedGroup_Co)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = BcellsPalette[names(BcellsPalette) %in% tmp$predictedGroup_Co]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rank(match(rownames(LabelPositions), names(BcellsPalette))),
#            color = LabelPositions$color, size = 7,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "bottom", legend.title = element_blank(),
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(nrow = 1, override.aes = list(size = 9,
#                                                   shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$predictedGroup_Co))]
#   )))
# dev.off()
# 
# svg(filename = "BATAC_Reclustering_Colored_Samples.svg", width = 8, height = 8.5)
# tmp <- getEmbedding(ArchRProj = B.ATAC, embedding = "HarmonyUMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$Sample <- B.ATAC$Sample
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
# # Double-check canonical cell marker
# 
# markerGenes  <- c(
#   "IL4R",
#   "MS4A1"
# )
# 
# B.ATAC <- addImputeWeights(B.ATAC)
# 
# p <- plotEmbedding(
#   ArchRProj = B.ATAC,
#   colorBy = "GeneScoreMatrix",
#   name = markerGenes,
#   continuousSet = "horizonExtra",
#   embedding = "HarmonyUMAP",
#   imputeWeights = getImputeWeights(B.ATAC)
# )
# 
# for (i in markerGenes) {
#   png(filename = paste("CheckMarkerLabelTransfer/", i, "embedding.png", sep = ""),  width = 2000, height = 2000, res = 300)
#   print(p[i])
#   dev.off()
# }
# 
# B.ATAC <- addMotifAnnotations(ArchRProj = B.ATAC, motifSet = "cisbp", name = "Motif")
# 
# B.ATAC <- addBgdPeaks(B.ATAC)
# B.ATAC <- addDeviationsMatrix(
#   ArchRProj = B.ATAC,
#   peakAnnotation = "Motif",
#   force = T
# )
# 
# B.ATAC$Group <- ""
# B.ATAC$Group[grep(B.ATAC$Sample, pattern = "Control")] <- "Control"
# B.ATAC$Group[grep(B.ATAC$Sample, pattern = "Patient")] <- "Patient"
# 
# B.ATAC$GroupSpecificPseudoBulk <- paste(B.ATAC$Group, B.ATAC$predictedGroup_Co)
# B.ATAC <- addGroupCoverages(ArchRProj = B.ATAC, groupBy = "GroupSpecificPseudoBulk")
# 
# B.RNA$Group <- ""
# B.RNA$Group[grep(B.RNA$orig.ident, pattern = "Control")] <- "Control"
# B.RNA$Group[B.RNA$orig.ident %in% c("PatientPBMC")] <- "Patient"
# 
# B.RNA$GroupCluster <- paste(B.RNA$Group, B.RNA$IntermediateGranularity, sep = " ")
# 
# saveArchRProject(ArchRProj = B.ATAC, outputDirectory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/ATAC/", load = F, overwrite = T)
# save(B.RNA, file = "AnnotatedB.Rdata")

######################### Differential Analyses B ##############################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/")
# library(clusterProfiler)
# library(Seurat)
# library(org.Hs.eg.db)
# suppressMessages(library(ArchR))
# 
# B.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/ATAC/"))
# load(file = "AnnotatedB.Rdata")
# 
# ##### DEG Differential Expressed Genes
# 
# DEG <- list()
# 
# for (CellType in names(table(B.ATAC$predictedGroup_Co))) {
#   FirstGroup <- paste("Patient", CellType, sep = " ")
#   SecondGroup <- paste("Control", CellType, sep = " ")
# 
#   DEGPv1To2 <- FindMarkers(object = B.RNA[,B.RNA$orig.ident %in% c("PatientPBMC",
#                                                                        "Control1PBMC",
#                                                                        "Control2PBMC")],
#                            ident.1 = FirstGroup,
#                            ident.2 = SecondGroup,
#                            group.by = "GroupCluster",
#                            logfc.threshold = 0.25,
#                            test.use = "MAST")
# 
#   DEGPv3To5 <- FindMarkers(object = B.RNA[,B.RNA$orig.ident %in% c("PatientPBMC",
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
# for (CellType in names(table(B.ATAC$predictedGroup_Co))) {
#   FirstGroup <- paste("Patient", CellType)
#   SecondGroup <- paste("Control", CellType)
# 
#   CutOffUp <- "FDR < 0.05 & Log2FC > 0.25"
#   CutOffDown <- "FDR < 0.05 & Log2FC < -0.25"
# 
#   AllPeak <- getMarkerFeatures(
#     ArchRProj = B.ATAC,
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
#     ArchRProj = B.ATAC,
#     peakAnnotation = "Motif",
#     cutOff = CutOffUp
#   )
# 
#   DET[[CellType]]$motifsDown <- peakAnnoEnrichment(
#     seMarker = AllPeak,
#     ArchRProj = B.ATAC,
#     peakAnnotation = "Motif",
#     cutOff = CutOffDown
#   )
# }
# 
# GSEA.res <- List()
# for (CellType in names(table(B.ATAC$predictedGroup_Co))) {
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
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/")
# 
# library(org.Hs.eg.db)
# library(annotate)
# library(GenomicFeatures)
# library(EnsDb.Hsapiens.v86) # EnsDb.Hsapiens.v86 for human gene definitions of the Ensembl code database version 86 that based on the GRCh38 genome
# suppressMessages(library(ArchR))
# 
# B.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/ATAC/"))
# 
# TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# motifPositions <- getPositions(B.ATAC)
# GeneBody <- getGeneAnnotation(B.ATAC)$genes
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

########################## DownStream Analysis B #############################
rm(list = ls())
setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/")

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

BCR10Xp <- immunarch::repLoad(.path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ImmuneRepertoire/ProcessedData/BCR/")

B.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/ATAC/"))
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/AnnotatedB.Rdata")
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/Differential_Test_Results.Rdata")
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/Annotation/GenomicInformation.Rdata")

BcellsPalette <- c('Plasma' = "#FF9F9F",
                   'Atypical memory B' = "#FFCB6D",
                   'Naive B' = "#FCCFCC",
                   'Memory B' = "#FF7979",

                   'CD3E+ B' = "#FDA9CD",
                   'LYZ+ B' = "#BF6623",
                   'B' = "#FFA771"
)

# Integrate BCR and scRNA
Immune <-B.RNA[,B.RNA$orig.ident %in% c("Control3PBMC", "Control4PBMC", "PatientPBMC")]

Immune$CloneType <- NA
Immune$Combination <- NA
Immune$V <- NA
Immune$D <- NA
Immune$J <- NA
Immune$CloneSize <- NA
Immune$CD3AA <- NA

# Prevalent  Clones (Number of clone > 1) were assigned unique clonetype names
for (i in names(BCR10Xp$data)) {
  BCR10Xp$data[[i]]$CloneType <- paste("CloneType", c(1:nrow(BCR10Xp$data[[i]])))
  for (j in which(BCR10Xp$data[[i]]$Clones > 1)) {
    SampleBarcode <- colnames(Immune)[grep(Immune$orig.ident, pattern = i)]
    mergeindex <- match(unlist(strsplit(BCR10Xp$data[[i]]$Barcode[j], ";")), 
                        substr(SampleBarcode, start = 1, stop = 18)) %>% 
      na.omit() %>% as.numeric()
    Immune$CloneType[colnames(Immune) %in% SampleBarcode][mergeindex] <- paste(i, BCR10Xp$data[[i]]$CloneType[j])
    Immune$Combination[colnames(Immune) %in% SampleBarcode][mergeindex] <- paste(BCR10Xp$data[[i]]$V.name[j], BCR10Xp$data[[i]]$J.name[j])
    Immune$V[colnames(Immune) %in% SampleBarcode][mergeindex] <- BCR10Xp$data[[i]]$V.name[j]
    Immune$D[colnames(Immune) %in% SampleBarcode][mergeindex] <- BCR10Xp$data[[i]]$D.name[j]
    Immune$J[colnames(Immune) %in% SampleBarcode][mergeindex] <- BCR10Xp$data[[i]]$J.name[j]
    Immune$CD3AA[colnames(Immune) %in% SampleBarcode][mergeindex] <- BCR10Xp$data[[i]]$CDR3.aa[j]
  }
}

# Clonal space homeostasis
svg(filename = "Clonal_space_homeostasis.svg", width = 10.5, height = 4)
print(ClonalSpaceHomeostasis(BCR10Xp$data, 
                             ColSmall = "#F51717", ColMedium = "#F08484", 
                             ColLarge = "#FCC4C4", ColHyperexpanded = "#FDEDED"))
dev.off()

for (i in names(table(Immune$CloneType))) {
  Immune$CloneSize[Immune$CloneType == i] <- table(Immune$CloneType)[i]
}

svg(filename = "CloneSizeDimPlot.svg", width = 10.5, height = 4)
tmp <- FetchData(Immune, vars = c("UMAP_1", "UMAP_2", "orig.ident", "CloneType", "CloneSize"))
ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = CloneSize)) + geom_point(size = 0.3, alpha = 1) + 
  theme_ArchR() + scale_color_gradient(low = "#FF7979", high = "#FFD9D9") + 
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
                              "^TNF", "^TGF", "^IFN", "^GZM", "^STAT", 
                              "CD80", "CD86", "ICOSLG", "CD37", "CD22", "CD40",
                              paste("^IL", c(1:50), sep = "")), collapse = "|"), 
                      rownames(B.RNA), value = T)

DiffImmuneMarkers <- list()

# For visualization, a log2FC curoff of 0.4 was set

for (CellType in names(table(B.ATAC$predictedGroup_Co))) {
  DiffImmuneMarkers[[CellType]] <- rownames(DEG[[CellType]][rownames(DEG[[CellType]]) %in% c(PotentialBioMarkers, ImmuneMarkers) &
                                                              abs(DEG[[CellType]]$avg_log2FC) > 0.4,])
}

DiffImmuneMarkers <- Reduce(union, list(DiffImmuneMarkers))
DiffImmuneMarkers <- unique(unlist(DiffImmuneMarkers))

levels(DiffImmuneMarkers) <- rev(c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                                   "STAT1", "STAT2", "STAT3", 
                                   "CXCR4", 
                                   "IFNG-AS1", 
                                   "IL16",
                                   "IL2RA", "IL2RG", "IL4R", 
                                   "TGFB1", "TGFBR2",
                                   "TNFRSF13C", "TNFRSF13B", "TNFRSF14",
                                   "CD37", "CD40"))  

AnalyseSubTypes <- names(table(B.ATAC$predictedGroup_Co))
levels(AnalyseSubTypes) <- c("Naive B", "B")

svg(filename = "BBubble.svg", width = 5, height = 7)
BubblePlot(B.RNA[,B.RNA$orig.ident %in% c("Control1PBMC", "Control2PBMC", "PatientPBMC")], Features = DiffImmuneMarkers,
           CellTypes = AnalyseSubTypes, ColData = "IntermediateGranularity", DiffRes = DEG, colorlow = "#F65E00", colorhigh = "#FEE2E2")
dev.off()

svg(filename = "RelativeFraction_B.svg", width = 6.5, height = 5.5)
RelativeFraction(B.RNA, ColData = "FineGranularity")
dev.off()

svg(filename = "ICOSLG.svg", width = 8, height = 4)
FeaturePlot(B.RNA, features = "ICOSLG", split.by = "Group")
dev.off()

svg(filename = "CD80.svg", width = 8, height = 4)
FeaturePlot(B.RNA, features = "CD80", split.by = "Group")
dev.off()

svg(filename = "CD86.svg", width = 8, height = 4)
FeaturePlot(B.RNA, features = "CD86", split.by = "Group")
dev.off()

svg(filename = "CD40.svg", width = 8, height = 4)
FeaturePlot(B.RNA, features = "CD40", split.by = "Group")
dev.off()

svg(filename = paste("GSEA/", "GSEAB.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$`B`, Descriptions = c("immunoglobulin production", 
                                            "immunoglobulin production involved in immunoglobulin mediated immune response", 
                                            "B cell activation", "B cell receptor signaling pathway"), 
             base_size = 14, legendposition = c(0.5, 0.3), rel_heights = c(1.5, 0.5, 1.0), GSEAColor = c("#FDA9CD", "#FCCFCC", 
                                                                                                         "#FF7979", "#BF6623"))
dev.off()

svg(filename = paste("GSEA/", "GSEANaiveB.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$`Naive B`, Descriptions = c("B cell receptor signaling pathway"), 
             base_size = 15, legendposition = c(0.25, 0.2), rel_heights = c(1.5, 0.5, 1.0), GSEAColor = c("#FCCFCC"))
dev.off()


GSEANBrec <- CoreEnrichment(GSEA.res$`Naive B`@result$core_enrichment[GSEA.res$`Naive B`@result$Description %in%
                                                                        c("B cell receptor signaling pathway")])

GSEANBbac <- CoreEnrichment(GSEA.res$`Naive B`@result$core_enrichment[GSEA.res$`Naive B`@result$Description %in%
                                                                        c("response to bacterium")])

GSEABglo <- CoreEnrichment(GSEA.res$`B`@result$core_enrichment[GSEA.res$`B`@result$Description %in%
                                                                        c("immunoglobulin production", 
                                                                          "immunoglobulin production involved in immunoglobulin mediated immune response")])

GSEABact <- CoreEnrichment(GSEA.res$`B`@result$core_enrichment[GSEA.res$`B`@result$Description %in%
                                                                  c("B cell activation")])

GSEABrec <- CoreEnrichment(GSEA.res$`B`@result$core_enrichment[GSEA.res$`B`@result$Description %in%
                                                                  c("B cell receptor signaling pathway")])

GSEABbac <- CoreEnrichment(GSEA.res$`B`@result$core_enrichment[GSEA.res$`B`@result$Description %in%
                                                                  c("response to bacterium")])

################################## ATAC ########################################
ControlPalette <- BcellsPalette[names(BcellsPalette) %in% B.ATAC$predictedGroup_Co]
cols1 <- readhex(file = textConnection(paste(ControlPalette, collapse = "\n")),
                 class = "RGB")

#transform to hue/lightness/saturation colorspace

cols1 <- as(cols1, "HLS")
PatientPalette <- cols1 -> ControlPalette

#multiplicative decrease of lightness

PatientPalette@coords[, "L"] <- PatientPalette@coords[, "L"] * 0.75
ControlPalette@coords[, "L"] <- ControlPalette@coords[, "L"] * 1.05

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

names(PatientPalette) <- names(BcellsPalette[names(BcellsPalette) %in% B.ATAC$predictedGroup_Co]) -> names(ControlPalette)
names(PatientPalette) <- paste("Patient", names(PatientPalette))
names(ControlPalette) <- paste("Control", names(ControlPalette))
ATACPalette <- c(PatientPalette, ControlPalette)

orderC <- c("Control Naive B",
            "Control B")

orderP <- c("Patient Naive B",
            "Patient B")

################################################################################
Patient.B.ATAC <- B.ATAC[B.ATAC$Sample == "Patient-scATAC",]
Control.B.ATAC <- B.ATAC[!B.ATAC$Sample == "Patient-scATAC",]

# Create a MotifMatrix of chromVAR deviations and deviation z-scores for all motifs. 
# This data, averaged by clusters, was obtained by using the getGroupSE() function 
# which returns a SummarizedExperiment. 

seGroupMotif <- getGroupSE(ArchRProj = B.ATAC, useMatrix = "MotifMatrix", groupBy = "GroupSpecificPseudoBulk")

# This SummarizedExperiment was further subset to just the deviation z-scores.

DeviationZ <- seGroupMotif[rowData(seGroupMotif)$seqnames == "z",]

for (CellType in names(table(B.ATAC$predictedGroup_Co))) {
  svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/EnrichedMotif/Up-", CellType, ".svg", sep = ""), width = 6, height = 7)
  print(PlotEnrichedTFMotif(DET[[CellType]]$motifsUp, dotsize = 3, labelsize = 4, repulsion = 25, DeviationZ = DeviationZ, UpIn = "Patient"))
  dev.off()
  
  svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/EnrichedMotif/Down-", CellType, ".svg", sep = ""), width = 6, height = 7)
  print(PlotEnrichedTFMotif(DET[[CellType]]$motifsDown, dotsize = 3, labelsize = 4, repulsion = 25, DeviationZ = DeviationZ, UpIn = "Control"))
  dev.off()
  
  DET[[CellType]]$motifsUpTop <- GetRepresentativeTF(DET[[CellType]]$motifsUp, topTF = 20, DeviationZ = DeviationZ, UpIn = "Patient")
  DET[[CellType]]$motifsDownTop <- GetRepresentativeTF(DET[[CellType]]$motifsDown, topTF = 20, DeviationZ = DeviationZ, UpIn = "Control")
}

AllSigTFMotif <- list()
for (CellType in names(table(B.ATAC$predictedGroup_Co))) {
  AllSigTFMotif[[paste(CellType, "Up")]] <- names(DET[[CellType]]$motifsUpTop)
  AllSigTFMotif[[paste(CellType, "Down")]] <- names(DET[[CellType]]$motifsDownTop)
}

UnionSigTFMotif <- Reduce(union, AllSigTFMotif)

IntersectSigTFMotifUp <- Reduce(intersect, AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Up")])
IntersectSigTFMotifDown <- Reduce(intersect, AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Down")])

################# Associate CAR and TF and conversion to bam ###################
# Fragment_By_CellType <- list()
# for (CellType in names(table(B.ATAC$GroupSpecificPseudoBulk))) {
#   SubSlot <- paste(unlist(strsplit(CellType, split = " ")), collapse = "_")
# 
#   Fragment_By_CellType[[SubSlot]] <- unlist(getFragmentsFromProject(B.ATAC, cellNames = rownames(B.ATAC[B.ATAC$GroupSpecificPseudoBulk == CellType,])))
#   Fragment_By_CellType[[SubSlot]]$index <- as.character(Fragment_By_CellType[[SubSlot]]$RG)
# }
# 
# library(GenomicRanges)
# 
# Output_to_bed_files(Fragment_list_cl = Fragment_By_CellType, folder = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/FragmentsByCellTypes/")
# 
# # Execute bedpetobam in /FragmentsByCellTypes, output to /BAM
# 
# bedpetobam <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   outdir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/BAM/"
#   bedpetobam <- c(bedpetobam, paste("bedtools bedpetobam -i ", i, "_fragments_cl_bamGR_pe.bed -g hg38.chrom.sizes.txt > ", outdir, i, "_fragments_cl_bamGR_pe.bam", sep = ""))
# }
# 
# bedpetobam <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/FragmentsByCellTypes/"), bedpetobam)
# 
# # Execute samtools_sort and samtools_index in /BAM, output to /sBAM
# 
# samtools_sort <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   outdir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/sBAM/"
#   samtools_sort <- c(samtools_sort, (paste("samtools sort -o ", outdir, i, "_fragments_cl_bamGR_pe_s.bam ", i, "_fragments_cl_bamGR_pe.bam", sep = "")))
# }
# 
# samtools_sort <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/BAM/"), samtools_sort)
# 
# # Execute samtools_sort and samtools_index in /sBAM, output to /sBAM
# 
# samtools_index <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   samtools_index <- c(samtools_index, paste("samtools index ", i, "_fragments_cl_bamGR_pe_s.bam", sep = ""))
# }
# 
# samtools_index <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/sBAM/"), samtools_index)
# 
# write.table(c(bedpetobam, samtools_sort, samtools_index), file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/fragments2bams.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

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
# export.bed(getPeakSet(B.ATAC), con = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ATACorrect/AllPeaks.bed")
# 
# ATACorrect <- c()
# for(i in names(Fragment_By_CellType)){
#   Input <- paste(i, "_fragments_cl_bamGR_pe_s.bam", sep = "")
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   Genome <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/hg38.fa"
#   BlackList <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/BlackList/Hg38TotalBlackList.bed"
#   OutPut <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ATACorrect/"
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
# ATACorrect <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/sBAM/"), ATACorrect)
# 
# write.table(ATACorrect, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ATACorrectScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

############################### ScoreBigwig ####################################
# ScoreBigwig <- c()
# for(i in names(Fragment_By_CellType)){
#   Input <- paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ATACorrect/", i, "/", i, "_fragments_cl_bamGR_pe_s_corrected.bw", sep = "")
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   OutPut <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ScoreBigwig/"
#   OutFileName <- paste(i, "_footprints.bw", sep = "")
# 
#   ScoreBigwig <- c(ScoreBigwig, paste("nohup", "TOBIAS FootprintScores --signal", Input, "--regions", AllPeaks, "--output", OutFileName, "--cores 16 &"))
# }
# 
# ScoreBigwig <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ScoreBigwig/"), ScoreBigwig)
# 
# write.table(ScoreBigwig, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ScoreBigwigScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

##################### BINDetect differential TF occupancy ######################
# If error occurs, enter ulimit -n 50000 in console
# BINDetect <- c()
# for(CellType in names(DEG)){
#   FirstGroup <- paste(c("Patient", paste(unlist(strsplit(CellType, split = " ")), collapse = "_")), collapse = "_")
#   SecondGroup <- paste(c("Control", paste(unlist(strsplit(CellType, split = " ")), collapse = "_")), collapse = "_")
# 
#   motifsMEME <- "/mnt/TRAINING/Zhongjiacheng/pwm_to_meme/Homo.cisbp.meme"
#   Genome <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/hg38.fa"
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   BigwigDir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/ScoreBigwig/"
#   OutDir <- paste(FirstGroup, SecondGroup, sep = "_vs_")
# 
#   BINDetect <- c(BINDetect, paste("TOBIAS BINDetect", "--motifs", motifsMEME, "--signals", paste(BigwigDir, FirstGroup, "_footprints.bw", sep = ""),
#                                   paste(BigwigDir, SecondGroup, "_footprints.bw", sep = ""), "--genome", Genome, "--peaks", AllPeaks, "--outdir"
#                                   , OutDir, "--cond_names", paste(FirstGroup, SecondGroup) , "--cores 16"))
# }
# 
# BINDetect <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/DiffTFOccupancy/"), BINDetect)
# 
# write.table(BINDetect, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/BINDetectScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

################ Visualization of differential TF binding ######################

ValidTFB <- names(table(Patient.B.ATAC$predictedGroup_Co))
ValidTFB <- unlist(lapply(strsplit(ValidTFB, split = " "), paste, collapse = "_"))

################ Subset the binding events to TF binding only ##################

DiffTFBind <- list()

for (CellType in ValidTFB) {
  tmp <- readxl::read_xlsx(paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/B/DownStream/TOBIAS/DiffTFOccupancy/", 
                                 "Patient_", CellType, "_vs_Control_", CellType, "/bindetect_results.xlsx", sep = ""), sheet = 2) %>% as.data.frame()
  tmp$clustername <- unlist(lapply(strsplit(tmp$cluster, split = "_"), '[[', 2))
  DiffTFBind[[CellType]] <- tmp[tmp$clustername %in% unlist(lapply(strsplit(names(motifPositions), split = "_"), '[[', 1)),]
}

svg(filename = paste("Volcano/", "BVol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["B"]], 
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control B"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient B"])))
dev.off()

svg(filename = paste("Volcano/", "NBVol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["Naive_B"]], 
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control Naive B"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient Naive B"])))
dev.off()

suppressMessages(InspectEmbedding(B.ATAC, QueryGene = "ZNF524", 
                                  UseMatrix = getAvailableMatrices(B.ATAC)[3],
                                  Embed = "HarmonyUMAP"))



