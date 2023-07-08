################# Combined Analysis T cells Preprocessing 1#####################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Preprocessing/")
# library(Seurat)
# library(dplyr)
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/SeuratPBMC/ReannotationAndElimination/SeuratPBMC_Annotated_TNK.Rdata")
# 
# TNK.RNA <- DietSeurat(TNK.RNA, assays = "RNA")
# 
# TNK.list <- SplitObject(TNK.RNA, split.by = "orig.ident")
# TNK.list <- lapply(TNK.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = T) # This step was intended to display the original normalized counts in subsequent visualization (SCTransform does not normalize the original RNA assay)
#   x <- SCTransform(x, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = T)
# })
# 
# integ_features <- SelectIntegrationFeatures(object.list = TNK.list,
#                                             nfeatures = 3000)
# 
# TNK.list <- PrepSCTIntegration(object.list = TNK.list,
#                                       anchor.features = integ_features)
# 
# integ_anchors <- FindIntegrationAnchors(object.list = TNK.list,
#                                         normalization.method = "SCT",
#                                         anchor.features = integ_features)
# 
# TNK.RNA <- IntegrateData(anchorset = integ_anchors,
#                          normalization.method = "SCT")
# 
# TNK.RNA <- RunPCA(TNK.RNA)
# 
# save(TNK.RNA, file = "TNK.RNA.Rdata")

################### Combined Analysis T cells Preprocessing2 ###################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Preprocessing/")
# library(Seurat)
# library(dplyr)
#
# load("TNK.RNA.Rdata")
#
# svg(filename = paste("Elbow_Plot.svg", sep = ""), width = 7, height = 5)
# ElbowPlot(TNK.RNA, ndims = 50) # to determine number of dimensions for clustering
# dev.off()
#
# pc.use <- 27
# TNK.RNA <- FindNeighbors(TNK.RNA, dims = 1:pc.use, verbose = T)
# TNK.RNA <- FindClusters(TNK.RNA, resolution = 2.8)
# TNK.RNA <- RunUMAP(TNK.RNA, dims = 1:pc.use)
#
# DimPlot(TNK.RNA, label = T)
#
# DefaultAssay(TNK.RNA) <- "RNA"
# Cluster.markers.t <- FindAllMarkers(TNK.RNA, only.pos = F, test.use = "negbinom", latent.vars = "orig.ident")
#
# save(TNK.RNA, Cluster.markers.t, file = "TNK_reclustering.Rdata")

#####################Combined Analysis T cells Annotation#######################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/")
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
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Preprocessing/TNK_reclustering.Rdata")
# 
# Cluster.markers.t$Delta <- Cluster.markers.t$pct.1 - Cluster.markers.t$pct.2
# top10 <- Cluster.markers.t %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_log2FC)
# 
# TNK.RNA$FineGranularity <- ""
# TNK.RNA$IntermediateGranularity <- ""
# 
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("3", "4", "6", "8", "12", "14", "17", "33", "35")] <- "Naive CD4" #LEF1, SELL
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("11", "13", "26")] <- "CD4"
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("21", "22")] <- "Memory CD4" #S100A4 PMID: 22421787
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("19")] <- "GZMK+ CD4"
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("5", "31")] <- "Th2-like" #GATA3
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("2")] <- "Tfh-like" #MAF, ICOS
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("10")] <- "Th17-like" #CCR6 (S100A4)
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("29", "30")] <- "Treg" #FOXP3
# 
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("0", "7", "18", "37")] <- "Naive CD8" #LEF1, SELL
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("28")] <- "CD8"
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("16", "20")] <- "GZMK+ CD8"
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("1", "15")] <- "GZMH+ CD8" # GZMH S100A4
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("9", "25")] <- "NK" #
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("23", "27")] <- "NKT" #
# 
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("24")] <- "MAIT" #SLC4A10
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("32")] <- "γδ"  # TRDC, TRGC1, TRDV2 gamma/delta T cells are usually CD4-CD8-
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("39")] <- "Cycling LC" #MKI67
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("34", "36")] <- "Other LC"
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("40")] <- "LYZ+ LC"
# TNK.RNA$FineGranularity[TNK.RNA$seurat_clusters %in% c("38")] <- "Unstimulated NK" #XCL1 and XCL2 are constitutively expressed by unstimulated natural killer (NK) cells
# 
# TNK.RNA$IntermediateGranularity <- TNK.RNA$FineGranularity
# 
# Cluster2Unicode <- readxl::read_xlsx(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Dec2UnicodeSymbol.xlsx")
# TcellsPalette <- c('Naive CD8' = "#A8ABEA",
#                    'Naive CD4' = "#B2E6EC",
#                    'CD4' = "#77C6DD",
# 
#                    'Tfh-like' = "#B4BDDE",
#                    'Th2-like' = "#C1EFFF",
#                    'Th17-like' = "#9ECFF8",
# 
#                    'Treg' = "#4AB4B1",
#                    'Memory CD4' = "#94B0F0",
#                    'GZMK+ CD4' = "#55C0F5",
# 
#                    'CD8' = "#0D97D5",
#                    'GZMK+ CD8' = "#7B93DF",
#                    'GZMH+ CD8' = "#8F84B6",
# 
#                    'γδ' = "#CCBFFD",
#                    'MAIT' = "#B494E4",
#                    'NKT' = "#6F9DB9",
# 
#                    'NK' = "#8672D0",
#                    'Unstimulated NK' = "#8B41ED",
#                    'Cycling LC' = "#6574FF",
#                    'Other LC' = "#93A0BF",
#                    'LYZ+ LC' = "#6FABB5"
# )
# 
# NumCluster <- length(table(TNK.RNA$seurat_clusters)) - 1
# 
# TNK.RNA$IntermediateGranularity[TNK.RNA$FineGranularity %in% c("CD8", "GZMH+ CD8", "GZMK+ CD8")] <- "CD8"
# TNK.RNA$IntermediateGranularity[TNK.RNA$FineGranularity %in% c("Tfh-like", "Th17-like", "Th2-like", "CD4", "Memory CD4", "GZMK+ CD4")] <- "CD4"
# 
# Reanno.cluster.ids <- TNK.RNA$FineGranularity[match(as.character(0:NumCluster), table = TNK.RNA$seurat_clusters)]
# names(Reanno.cluster.ids) <- levels(TNK.RNA@active.ident)
# TNK.RNA <- RenameIdents(TNK.RNA, Reanno.cluster.ids)
# 
# DimPlot(TNK.RNA)
# 
# svg(filename = "TNK_IntermediateGranularity.svg", width = 9.5, height = 7.5)
# tmp <- FetchData(TNK.RNA, vars = c("UMAP_1", "UMAP_2", "orig.ident", "IntermediateGranularity"))
# 
# tmp$Labels <- match(tmp$IntermediateGranularity, names(TcellsPalette))
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = IntermediateGranularity)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = TcellsPalette[names(TcellsPalette) %in% tmp$IntermediateGranularity]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rank(match(rownames(LabelPositions), names(TcellsPalette))),
#            color = LabelPositions$color, size = 7,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
#                         legend.key.size = unit(1.3, 'cm'),
#                         legend.text = element_text(size = 17),
#                         axis.title = element_blank(),
#                         axis.text = element_blank(),
#                         axis.ticks = element_blank()) +
#   guides(color = guide_legend(override.aes = list(size = 9,
#                                                   shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$IntermediateGranularity))]
#   )))
# dev.off()
# 
# svg(filename = "TNKFineGranularity.svg", width = 9.5, height = 7.5)
# tmp <- FetchData(TNK.RNA, vars = c("UMAP_1", "UMAP_2", "orig.ident", "FineGranularity"))
# 
# tmp$Labels <- match(tmp$FineGranularity, names(TcellsPalette))
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = FineGranularity)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = TcellsPalette[names(TcellsPalette) %in% tmp$FineGranularity]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rownames(LabelPositions),
#            color = LabelPositions$color, size = 7,  fontface = "bold") +
#   theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
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
# svg(filename = paste("TNK_FineGranularity_Violin.svg", sep = ""), width = 12, height = 9)
# Features <- c("CD4", "CD8A", "SELL", "LEF1", "S100A4", "CCR6", "RORC", "GATA3", "ICOS", "MAF", "FOXP3", "GZMH", "GZMK", "ZEB2",
#               "KLRF1",  "XCL1", "XCL2", "SLC4A10", "TRDC", "TRGC1", "TRDV2", "MKI67", "LYZ")
# 
# long <- FetchData(TNK.RNA, vars = Features, slot = "data")
# long$FineGranularity <- TNK.RNA$FineGranularity
# long$Cell <- rownames(long)
# long <- reshape2::melt(long, id.vars = c("Cell","FineGranularity"), measure.vars = Features,
#                        variable.name = "Feat", value.name = "Expr")
# 
# long$FineGranularity <- factor(x = long$FineGranularity, levels = rev(names(TcellsPalette)))
# long$Color <- TcellsPalette[match(long$FineGranularity, table = names(TcellsPalette))]
# 
# ggplot(long, aes(Expr, factor(FineGranularity))) +
#   geom_violin(scale = "width", adjust = 1, trim = TRUE, aes(fill = FineGranularity, color = FineGranularity)) +
#   scale_fill_manual("legend", values = TcellsPalette) +
#   scale_color_manual("legend", values = TcellsPalette) +
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
# TNK.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ATAC/AnnotatedATAC/TNK/"))
# 
# TNK.ATAC <- addIterativeLSI(ArchRProj = TNK.ATAC, force = T)
# 
# TNK.ATAC <- addClusters(
#   input = TNK.ATAC,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "Clusters",
#   resolution = 0.8
# )
# 
# # The following tmp object was only used for peak calling,
# # wherein each cluster contains at least 200 cells
# tmp <- TNK.ATAC[TNK.ATAC$Clusters %in% names(which(table(TNK.ATAC$Clusters) > 200)),]
# tmp <- addGroupCoverages(ArchRProj = tmp, groupBy = "Clusters")
# tmp <- addReproduciblePeakSet(
#   ArchRProj = tmp, reproducibility = "1",
#   groupBy = "Clusters", pathToMacs2 = "/mnt/TRAINING/Zhongjiacheng/anaconda2/bin/macs2"
# )
# 
# TNK.ATAC <- addPeakSet(TNK.ATAC, peakSet = tmp@peakSet, force = T)
# TNK.ATAC <- addPeakMatrix(TNK.ATAC)
# 
# TNK.ATAC <- addIterativeLSI(
#   ArchRProj = TNK.ATAC,
#   useMatrix = "PeakMatrix",
#   name = "PeakIterativeLSI",
#   force = T
# )
# 
# TNK.ATAC <- addHarmony(
#   ArchRProj = TNK.ATAC,
#   reducedDims = "PeakIterativeLSI",
#   name = "Harmony",
#   groupBy = "Sample", force = T
# )
# 
# TNK.ATAC <- addUMAP(
#   ArchRProj = TNK.ATAC,
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
#   RNA.cells <- colnames(TNK.RNA[,grep(TNK.RNA$orig.ident, pattern = i)])
#   ATAC.cells <- grep(TNK.ATAC$cellNames, pattern = i, value = T)
# 
#   groupList[[i]] <- SimpleList(
#     ATAC = ATAC.cells,
#     RNA = RNA.cells
#   )
# }
# 
# TNK.ATAC <- addGeneIntegrationMatrix(
#   ArchRProj = TNK.ATAC,
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "PeakIterativeLSI",
#   seRNA = TNK.RNA,
#   addToArrow = T,
#   groupList = groupList,
#   groupRNA = "IntermediateGranularity",
#   nameCell = "predictedCell_Co",
#   nameGroup = "predictedGroup_Co",
#   nameScore = "predictedScore_Co",
#   force = T, useImputation = F
# )
# 
# png(filename =  "Constrained_T_Cell_Type_Transfer_Score.png", width = 2000, height = 2000, res = 300)
# print(plotEmbedding(
#   TNK.ATAC,
#   colorBy = "cellColData",
#   name = "predictedScore_Co",
#   embedding = "HarmonyUMAP"
# ))
# dev.off()
# 
# TNK.ATAC <- TNK.ATAC[TNK.ATAC$predictedScore_Co > 0.5,]
# 
# svg(filename =  "Constrained_T_Cell_Type_Transfer_Score_Filtered.svg", width = 6, height = 6)
# print(plotEmbedding(
#   TNK.ATAC,
#   colorBy = "cellColData",
#   name = "predictedScore_Co",
#   embedding = "HarmonyUMAP"
# ))
# dev.off()
# 
# # The following step filters out cell types that were not shared across all samples
# 
# ValidPseudoBulk <- Reduce(intersect, list(names(table(TNK.ATAC$predictedGroup_Co[TNK.ATAC$Sample == "Control1-scATAC"])),
#                                           names(table(TNK.ATAC$predictedGroup_Co[TNK.ATAC$Sample == "Control2-scATAC"])),
#                                           names(table(TNK.ATAC$predictedGroup_Co[TNK.ATAC$Sample == "Patient-scATAC"]))))
# 
# # The following step filters out cell types with number less than 40 in at least one sample
# 
# ValidPseudoBulk <- Reduce(intersect, list(names(which(table(TNK.ATAC$predictedGroup_Co[TNK.ATAC$Sample == "Control1-scATAC" & TNK.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40)),
#                                           names(which(table(TNK.ATAC$predictedGroup_Co[TNK.ATAC$Sample == "Control2-scATAC" & TNK.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40)),
#                                           names(which(table(TNK.ATAC$predictedGroup_Co[TNK.ATAC$Sample == "Patient-scATAC" & TNK.ATAC$predictedGroup_Co %in% ValidPseudoBulk]) > 40))))
# 
# TNK.ATAC <- TNK.ATAC[TNK.ATAC$predictedGroup_Co %in% ValidPseudoBulk,]
# 
# svg(filename = "Constrained_T_Cell_Type_Transfer_Filtered_Group.svg", width = 8, height = 9)
# tmp <- getEmbedding(ArchRProj = TNK.ATAC, embedding = "HarmonyUMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$predictedGroup_Co <- TNK.ATAC$predictedGroup_Co
# 
# LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[grep(rownames(tmp), pattern = "Control"),"UMAP_1"]), INDEX = tmp[grep(rownames(tmp), pattern = "Control"),"predictedGroup_Co"], FUN = mean)),
#                              UMAP2 = data.frame(UMAP2 = tapply((tmp[grep(rownames(tmp), pattern = "Control"),"UMAP_2"]), INDEX = tmp[grep(rownames(tmp), pattern = "Control"),"predictedGroup_Co"], FUN = mean)))
# 
# LabelPositions$color <- "black"
# 
# ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = predictedGroup_Co)) +
#   geom_point(size = 0.2, alpha = 1) +
#   scale_color_manual(values = TcellsPalette[names(TcellsPalette) %in% tmp$predictedGroup_Co]) +
#   annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rank(match(rownames(LabelPositions), names(TcellsPalette))),
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
# svg(filename = "TNKATAC_Reclustering_Colored_Samples.svg", width = 8, height = 8.5)
# tmp <- getEmbedding(ArchRProj = TNK.ATAC, embedding = "HarmonyUMAP", returnDF = TRUE)
# colnames(tmp) <- c("UMAP_1", "UMAP_2")
# tmp$Sample <- TNK.ATAC$Sample
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
#   "CD3E",
#   "CD4",
#   "CD8A",
#   "CD8B",
#   "KLRF1",
#   "SELL",
#   "LEF1",
#   "FOXP3"
# )
# 
# TNK.ATAC <- addImputeWeights(TNK.ATAC)
# 
# p <- plotEmbedding(
#   ArchRProj = TNK.ATAC,
#   colorBy = "GeneScoreMatrix",
#   name = markerGenes,
#   continuousSet = "horizonExtra",
#   embedding = "HarmonyUMAP",
#   imputeWeights = getImputeWeights(TNK.ATAC)
# )
# 
# for (i in markerGenes) {
#   png(filename = paste("CheckMarkerLabelTransfer/", i, "embedding.png", sep = ""),  width = 2000, height = 2000, res = 300)
#   print(p[i])
#   dev.off()
# }
# 
# TNK.ATAC <- addMotifAnnotations(ArchRProj = TNK.ATAC, motifSet = "cisbp", name = "Motif")
# 
# TNK.ATAC <- addBgdPeaks(TNK.ATAC)
# TNK.ATAC <- addDeviationsMatrix(
#   ArchRProj = TNK.ATAC,
#   peakAnnotation = "Motif",
#   force = T
# )
# 
# TNK.ATAC$Group <- ""
# TNK.ATAC$Group[grep(TNK.ATAC$Sample, pattern = "Control")] <- "Control"
# TNK.ATAC$Group[grep(TNK.ATAC$Sample, pattern = "Patient")] <- "Patient"
# 
# TNK.ATAC$GroupSpecificPseudoBulk <- paste(TNK.ATAC$Group, TNK.ATAC$predictedGroup_Co)
# TNK.ATAC <- addGroupCoverages(ArchRProj = TNK.ATAC, groupBy = "GroupSpecificPseudoBulk")
# 
# TNK.RNA$Group <- ""
# TNK.RNA$Group[grep(TNK.RNA$orig.ident, pattern = "Control")] <- "Control"
# TNK.RNA$Group[TNK.RNA$orig.ident %in% c("PatientPBMC")] <- "Patient"
# 
# TNK.RNA$GroupCluster <- paste(TNK.RNA$Group, TNK.RNA$IntermediateGranularity, sep = " ")
# 
# saveArchRProject(ArchRProj = TNK.ATAC, outputDirectory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/ATAC/", load = F, overwrite = T)
# save(TNK.RNA, file = "AnnotatedTNK.Rdata")

######################## Differential Analyses TNK##############################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/")
# library(clusterProfiler)
# library(Seurat)
# library(org.Hs.eg.db)
# suppressMessages(library(ArchR))
# 
# TNK.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/ATAC/"))
# load(file = "AnnotatedTNK.Rdata")
# 
# ##### DEG Differential Expressed Genes
# 
# DEG <- list()
# 
# for (CellType in c(names(table(TNK.ATAC$predictedGroup_Co)), "NK")) {
#   FirstGroup <- paste("Patient", CellType, sep = " ")
#   SecondGroup <- paste("Control", CellType, sep = " ")
# 
#   DEGPv1To2 <- FindMarkers(object = TNK.RNA[,TNK.RNA$orig.ident %in% c("PatientPBMC",
#                                                                        "Control1PBMC",
#                                                                        "Control2PBMC")],
#                            ident.1 = FirstGroup,
#                            ident.2 = SecondGroup,
#                            group.by = "GroupCluster",
#                            logfc.threshold = 0.25,
#                            test.use = "MAST")
# 
#   DEGPv3To5 <- FindMarkers(object = TNK.RNA[,TNK.RNA$orig.ident %in% c("PatientPBMC",
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
# for (CellType in names(table(TNK.ATAC$predictedGroup_Co))) {
#   FirstGroup <- paste("Patient", CellType)
#   SecondGroup <- paste("Control", CellType)
# 
#   CutOffUp <- "FDR < 0.05 & Log2FC > 0.25"
#   CutOffDown <- "FDR < 0.05 & Log2FC < -0.25"
# 
#   AllPeak <- getMarkerFeatures(
#     ArchRProj = TNK.ATAC,
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
#     ArchRProj = TNK.ATAC,
#     peakAnnotation = "Motif",
#     cutOff = CutOffUp
#   )
# 
#   DET[[CellType]]$motifsDown <- peakAnnoEnrichment(
#     seMarker = AllPeak,
#     ArchRProj = TNK.ATAC,
#     peakAnnotation = "Motif",
#     cutOff = CutOffDown
#   )
# }
# 
# GSEA.res <- List()
# for (CellType in c(names(table(TNK.ATAC$predictedGroup_Co)), "NK")) {
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
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/")
# 
# library(org.Hs.eg.db)
# library(annotate)
# library(GenomicFeatures)
# library(EnsDb.Hsapiens.v86) # EnsDb.Hsapiens.v86 for human gene definitions of the Ensembl code database version 86 that based on the GRCh38 genome
# suppressMessages(library(ArchR))
# 
# TNK.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/ATAC/"))
# 
# TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# motifPositions <- getPositions(TNK.ATAC)
# GeneBody <- getGeneAnnotation(TNK.ATAC)$genes
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

########################## DownStream Analysis TNK #############################
rm(list = ls())
setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/")

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

TCR10Xp <- immunarch::repLoad(.path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/ImmuneRepertoire/ProcessedData/TCR/")

TNK.ATAC <- suppressMessages(loadArchRProject(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/ATAC/"))
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/AnnotatedTNK.Rdata")
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/Differential_Test_Results.Rdata")
load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/Annotation/GenomicInformation.Rdata")

TcellsPalette <- c('Naive CD8' = "#A8ABEA",
                   'Naive CD4' = "#B2E6EC",
                   'CD4' = "#77C6DD",

                   'Tfh-like' = "#B4BDDE",
                   'Th2-like' = "#C1EFFF",
                   'Th17-like' = "#9ECFF8",

                   'Treg' = "#4AB4B1",
                   'Memory CD4' = "#94B0F0",
                   'GZMK+ CD4' = "#55C0F5",

                   'CD8' = "#0D97D5",
                   'GZMK+ CD8' = "#7B93DF",
                   'GZMH+ CD8' = "#8F84B6",

                   'γδ' = "#CCBFFD",
                   'MAIT' = "#B494E4",
                   'NKT' = "#6F9DB9",

                   'NK' = "#8672D0",
                   'Unstimulated NK' = "#8B41ED",
                   'Cycling LC' = "#6574FF",
                   'Other LC' = "#93A0BF",
                   'LYZ+ LC' = "#6FABB5"
)

# Integrate TCR and scRNA
Immune <-TNK.RNA[,TNK.RNA$orig.ident %in% c("Control3PBMC", "Control4PBMC", "PatientPBMC")]

Immune$CloneType <- NA
Immune$Combination <- NA
Immune$V <- NA
Immune$D <- NA
Immune$J <- NA
Immune$CloneSize <- NA
Immune$CD3AA <- NA

# Clones (Number of clone > 1) were assigned unique clonetype names
for (i in names(TCR10Xp$data)) {
  TCR10Xp$data[[i]]$CloneType <- paste("CloneType", c(1:nrow(TCR10Xp$data[[i]])))
  for (j in which(TCR10Xp$data[[i]]$Clones > 1)) {
    SampleBarcode <- colnames(Immune)[grep(Immune$orig.ident, pattern = i)]
    mergeindex <- match(unlist(strsplit(TCR10Xp$data[[i]]$Barcode[j], ";")), 
                        substr(SampleBarcode, start = 1, stop = 18)) %>% 
      na.omit() %>% as.numeric()
    Immune$CloneType[colnames(Immune) %in% SampleBarcode][mergeindex] <- paste(i, TCR10Xp$data[[i]]$CloneType[j])
    Immune$Combination[colnames(Immune) %in% SampleBarcode][mergeindex] <- paste(TCR10Xp$data[[i]]$V.name[j], TCR10Xp$data[[i]]$J.name[j])
    Immune$V[colnames(Immune) %in% SampleBarcode][mergeindex] <- TCR10Xp$data[[i]]$V.name[j]
    Immune$D[colnames(Immune) %in% SampleBarcode][mergeindex] <- TCR10Xp$data[[i]]$D.name[j]
    Immune$J[colnames(Immune) %in% SampleBarcode][mergeindex] <- TCR10Xp$data[[i]]$J.name[j]
    Immune$CD3AA[colnames(Immune) %in% SampleBarcode][mergeindex] <- TCR10Xp$data[[i]]$CDR3.aa[j]
  }
}

# Clonal space homeostasis
svg(filename = "Clonal_space_homeostasis.svg", width = 10.5, height = 4)
print(ClonalSpaceHomeostasis(TCR10Xp$data, 
                             ColSmall = "#2171B5", ColMedium = "#6BAED6", 
                             ColLarge = "#BDD7E7", ColHyperexpanded = "#EFF3FF"))
dev.off()

for (i in names(table(Immune$CloneType))) {
  Immune$CloneSize[Immune$CloneType == i] <- table(Immune$CloneType)[i]
}

svg(filename = "CloneSizeDimPlot.svg", width = 10.5, height = 4)
tmp <- FetchData(Immune, vars = c("UMAP_1", "UMAP_2", "orig.ident", "CloneType", "CloneSize"))
ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = CloneSize)) + geom_point(size = 0.3, alpha = 1) + 
  theme_ArchR() + scale_color_gradient(low = "#2171B5", high = "#BDD7E7") + 
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

# Cells with prevalent TCR/BCR clonotypes (number of clones > 3) were deemed as expanded.

Immune$Expanded <- FALSE
Immune$Expanded[Immune$CloneSize > 3] <- TRUE

# monocle3 Inference on trajectory

Trajectory <- function(seu, rootcluster, directory, filename, palette, mainplottitle, insertplotfontsize = 3, insertplotheight, insertplotwidth) {
  df <- c()
  for (i in names(table(seu$FineGranularity))) {
    Expanded <- seu$Expanded[seu$FineGranularity == i]
    df <- rbind(df, c(i, length(which(Expanded))/length(Expanded)))
  }
  
  df <- as.data.frame(df)
  colnames(df) <- c("CellType", "ExpanPercent")
  df <- arrange(df, ExpanPercent)
  df$CellType <- factor(df$CellType, levels = as.character(df$CellType))
  df$ExpanPercent <- as.numeric(as.character(df$ExpanPercent))
  df$Color <- palette[match(df$CellType, names(palette))]
  
  svg(filename = paste(directory, "/", filename, "insertplot", ".svg", sep = ""), width = insertplotwidth, height = insertplotheight)
  print(ggplot(df, aes(x = CellType, y = ExpanPercent, fill = CellType)) + ylim(0, 1) + 
          geom_bar(stat = "identity", fill = df$Color) + theme_ArchR() + geom_text(aes(label = CellType), 
                                                                                   position = position_dodge(width = 0.9), hjust = "inward",
                                                                                   size = insertplotfontsize) +
          coord_flip() + theme(legend.position = "none", axis.title.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.text.y = element_blank(), 
                               axis.title.x = element_blank()))
  dev.off()
  
  cds <- as.cell_data_set(seu)
  cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
  cds <- learn_graph(cds, verbose = FALSE, use_partition = FALSE, close_loop = FALSE,
                     learn_graph_control = list(minimal_branch_len = 1))
  cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds)[cds@colData$FineGranularity == rootcluster])
  
  svg(filename = paste(directory, "/", filename, "trajectory", ".svg", sep = ""), width = 6, height = 6)

  print(plot_cells(cds, cell_size = 1, 
                   group_label_size = 0,
                   color_cells_by = "FineGranularity",
                   label_groups_by_cluster = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   label_roots = FALSE,
                   trajectory_graph_color = "#C0C0C0",
                   trajectory_graph_segment_size = 1.25) +  ggtitle(mainplottitle) +
          scale_color_manual(values = palette) + theme(plot.title = element_blank())) 
  dev.off()
  
  # return(pseudotime(cds, reduction_method = "UMAP"))
}

Trajectory(Immune[,Immune$orig.ident == "PatientPBMC" & Immune$IntermediateGranularity %in% c("CD4", "Naive CD4")], rootcluster = "Naive CD4", 
           directory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/RNATrajectory/", palette = TcellsPalette,
           filename = "CD4", mainplottitle = "Patient CD4", insertplotfontsize = 3, insertplotheight = 3, insertplotwidth = 3)

Trajectory(Immune[,Immune$orig.ident == "PatientPBMC" & Immune$IntermediateGranularity %in% c("CD8", "Naive CD8")], rootcluster = "Naive CD8", 
           directory = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/RNATrajectory/", palette = TcellsPalette,
           filename = "CD8", mainplottitle = "Patient CD8", insertplotfontsize = 3, insertplotheight = 2, insertplotwidth = 3)

# Clone type flows between compartments

PatientCD8Clone <- sort(table(grep(Immune$CloneType[Immune$IntermediateGranularity == "CD8"], pattern = "Patient", value = T)), decreasing = T)/
  length(grep(Immune$CloneType[Immune$IntermediateGranularity == "CD8"], pattern = "Patient", value = T))
PatientCD4Clone <- sort(table(grep(Immune$CloneType[Immune$IntermediateGranularity == "CD4"], pattern = "Patient", value = T)), decreasing = T)/
  length(grep(Immune$CloneType[Immune$IntermediateGranularity == "CD4"], pattern = "Patient", value = T))
PatientNKTClone <- sort(table(grep(Immune$CloneType[Immune$IntermediateGranularity == "NKT"], pattern = "Patient", value = T)), decreasing = T)/
  length(grep(Immune$CloneType[Immune$IntermediateGranularity == "NKT"], pattern = "Patient", value = T))

CombineCloneCD8 <- c()
for (i in names(PatientCD8Clone)) {
  CombineCloneCD8 <- rbind(CombineCloneCD8, data.frame(CellType = "PatientCD8",
                                                     CloneType = i, 
                                                     Proportion = PatientCD8Clone[names(PatientCD8Clone) == i]))
}

CombineCloneCD4 <- c()
for (i in names(PatientCD4Clone)) {
  CombineCloneCD4 <- rbind(CombineCloneCD4, data.frame(CellType = "PatientCD4",
                                                     CloneType = i, 
                                                     Proportion = PatientCD4Clone[names(PatientCD4Clone) == i]))
}

CombineCloneNKT <- c()
for (i in names(PatientNKTClone)) {
  CombineCloneCD8 <- rbind(CombineCloneCD8, data.frame(CellType = "PatientNKT",
                                                       CloneType = i, 
                                                       Proportion = PatientNKTClone[names(PatientNKTClone) == i]))
}

CombineCloneAll <- rbind(CombineCloneCD8, CombineCloneCD4, CombineCloneNKT)

CombineCloneAll$NumericOrder <- as.numeric(sapply(strsplit(as.character(CombineCloneAll$CloneType), split = "CloneType"), '[[', 2))

CombineCloneAll <- arrange(CombineCloneAll, CellType, desc(NumericOrder))

CombineCloneAll$Size <- CombineCloneAll$Proportion*10

CombineCloneAll$CloneTypeIdx <- as.numeric(unlist(lapply(strsplit(CombineCloneAll$CloneType, split = "Patient CloneType "), '[[', 2)))

CombineClonePalette <- colorRampPalette(colors = c("#EFF3FF", "#6BAED6", "#2171B5", "#BDD7E7"), bias = 0.9)(length(unique(CombineCloneAll$CloneTypeIdx)))

Divide <- round(length(unique(CombineCloneAll$CloneTypeIdx))/10)
SortedIdx <- sort(unique(CombineCloneAll$CloneTypeIdx))

RearrangeIdx <- SortedIdx[-c(1:10)]

for (i in 1:10) {
  AfterIndex <- i*Divide
  RearrangeIdx <- append(RearrangeIdx, SortedIdx[i], after = AfterIndex)
}

RearrangeIdx <- paste("Patient CloneType", RearrangeIdx)
names(CombineClonePalette) <- RearrangeIdx

DimPlot(Immune[,Immune$orig.ident == "PatientPBMC" & 
                 Immune$CloneType %in% paste(rep("Patient CloneType", 10), c(1:10))], group.by = "CloneType") + DimPlot(Immune[,Immune$orig.ident == "PatientPBMC"])

svg(filename = "SharedCloneTypes.svg", width = 6, height = 6)
ggplot(CombineCloneAll,
       aes(x = CellType, stratum = CloneType, alluvium = CloneType,
           y = Proportion,
           fill = CloneType, label = CloneType)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow(alpha = .6, width = 1/2, decreasing = FALSE) +
  geom_stratum(alpha = .6, width = 1/2, size = 0.2, decreasing = FALSE) +
  ggfittext::geom_fit_text(stat = "stratum", min.size = 0, width = 1/2, grow = TRUE, decreasing = FALSE) +
  ggtitle("CloneType flows between compartments") + theme_ArchR() + theme_minimal() + scale_fill_manual(values = CombineClonePalette) +
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(vjust = 5, size = 12),
        plot.title = element_text(hjust = 0.5, vjust = -3),
        axis.title.x = element_blank())
dev.off()

# Export top clone types

tmp <- TCR10Xp$data$Patient[TCR10Xp$data$Patient$CloneType %in% paste("CloneType", c(1:10)), c(1,3,4,5,6,7,20)]

BetaChain <- cbind(tmp$Clones,
                   sapply(strsplit(tmp$CDR3.nt, split = ";"), '[[', 2),
                   sapply(strsplit(tmp$CDR3.aa, split = ";"), '[[', 2),
                   sapply(strsplit(tmp$V.name, split = ";"), '[[', 2),
                   sapply(strsplit(tmp$D.name, split = ";"), '[[', 2),
                   sapply(strsplit(tmp$J.name, split = ";"), '[[', 2),
                   tmp$CloneType)
AlphaChain <- cbind(tmp$Clones,
                    sapply(strsplit(tmp$CDR3.nt, split = ";"), '[[', 1),
                    sapply(strsplit(tmp$CDR3.aa, split = ";"), '[[', 1),
                    sapply(strsplit(tmp$V.name, split = ";"), '[[', 1),
                    sapply(strsplit(tmp$D.name, split = ";"), '[[', 1),
                    sapply(strsplit(tmp$J.name, split = ";"), '[[', 1),
                    tmp$CloneType)

tmp <- cbind(AlphaChain, BetaChain)

write.csv(tmp, "Top10CloneType.csv", row.names = F)

################################ scVDJ annotation ##############################
# # Betavdjdb <- dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb",
# #                .species = "HomoSapiens", .chain = "TRB", .pathology = NA)
# # Alphavdjdb <- dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb",
# #                 .species = "HomoSapiens", .chain = "TRA", .pathology = NA)
# #
# # save(Betavdjdb, Alphavdjdb, file = "vdjdb.Rdata")
# load(file = "vdjdb.Rdata")
# 
# PatientTCRdf <- TCR10Xp$data$Patient
# 
# for (i in 1:nrow(PatientTCRdf)) {
#   PatientTCRdf$VGeneA[i] <- unlist(lapply(strsplit(PatientTCRdf$V.name, split = ";"), '[[', 1))[i]
#   PatientTCRdf$DGeneA[i] <- unlist(lapply(strsplit(PatientTCRdf$D.name, split = ";"), '[[', 1))[i]
#   PatientTCRdf$JGeneA[i] <- unlist(lapply(strsplit(PatientTCRdf$J.name, split = ";"), '[[', 1))[i]
#   PatientTCRdf$CDR3AA_Alpha[i] <- unlist(lapply(strsplit(PatientTCRdf$CDR3.aa, split = ";"), '[[', 1))[i]
# 
#     if (PatientTCRdf$DGeneA[i] == "NA") {
#     PatientTCRdf$AlphaChain[i] <- paste(c(PatientTCRdf$VGeneA[i], PatientTCRdf$JGeneA[i]), collapse = "/")
#   } else {
#     PatientTCRdf$AlphaChain[i] <- paste(c(PatientTCRdf$VGeneA[i], PatientTCRdf$DGeneA[i], PatientTCRdf$JGeneA[i]), collapse = "/")
#   }
# }
# 
# for (i in 1:nrow(PatientTCRdf)) {
#   PatientTCRdf$VGeneB[i] <- unlist(lapply(strsplit(PatientTCRdf$V.name, split = ";"), '[[', 2))[i]
#   PatientTCRdf$DGeneB[i] <- unlist(lapply(strsplit(PatientTCRdf$D.name, split = ";"), '[[', 2))[i]
#   PatientTCRdf$JGeneB[i] <- unlist(lapply(strsplit(PatientTCRdf$J.name, split = ";"), '[[', 2))[i]
#   PatientTCRdf$CDR3AA_Beta[i] <- unlist(lapply(strsplit(PatientTCRdf$CDR3.aa, split = ";"), '[[', 2))[i]
# 
#   if (PatientTCRdf$DGeneB[i] == "NA") {
#     PatientTCRdf$BetaChain[i] <- paste(c(PatientTCRdf$VGeneB[i], PatientTCRdf$JGeneB[i]), collapse = "/")
#   } else {
#     PatientTCRdf$BetaChain[i] <- paste(c(PatientTCRdf$VGeneB[i], PatientTCRdf$DGeneB[i], PatientTCRdf$JGeneB[i]), collapse = "/")
#   }
# }
# 
# PatientTCRdf[,c("Alpha.antigen.epitope", "Alpha.antigen.gene", "Alpha.antigen.species",
#                 "Beta.antigen.epitope", "Beta.antigen.gene", "Beta.antigen.species")] <- ""
# 
# # Nomenclature for T-cell receptor (TCR) gene segments of the immune system
# # The distinct NLV- and GIL- specific TCRs (defined by unique combinations of the V gene,
# # CDR3 amino acid sequence, and J gene) ranged from 21–784 for TCRα and from 13–1,030 for
# # TCRβ per person of freshly isolated cells
# 
# for (i in 1:nrow(PatientTCRdf)) {
#   tmpRowA <- c()
#   for (j in 1:nrow(Alphavdjdb)) {
#     if (PatientTCRdf$VGeneA[i] %in% unlist(strsplit(Alphavdjdb$v.segm[j], split = ",")) &
#         PatientTCRdf$JGeneA[i] %in% unlist(strsplit(Alphavdjdb$j.segm[j], split = ",")) &
#         PatientTCRdf$CDR3AA_Alpha[i] %in% unlist(strsplit(Alphavdjdb$cdr3[j], split = ","))) {
#       tmpRowA <- rbind(tmpRowA, Alphavdjdb[j,c("antigen.epitope", "antigen.gene", "antigen.species", "mhc.class")])
#       tmpRowA <- arrange(tmpRowA, antigen.species, antigen.gene)
#       PatientTCRdf[i,c("Alpha.antigen.epitope", "Alpha.antigen.gene", "Alpha.antigen.species", "Alpha.mhc.class")] <- unlist(lapply(tmpRowA, paste, collapse = "/"))
#     }
#   }
# 
#   tmpRowB <- c()
#   for (j in 1:nrow(Betavdjdb)) {
#     if (PatientTCRdf$VGeneB[i] %in% unlist(strsplit(Betavdjdb$v.segm[j], split = ",")) &
#         PatientTCRdf$JGeneB[i] %in% unlist(strsplit(Betavdjdb$j.segm[j], split = ",")) &
#         PatientTCRdf$CDR3AA_Beta[i] %in% unlist(strsplit(Betavdjdb$cdr3[j], split = ","))) {
#       tmpRowB <- rbind(tmpRowB, Betavdjdb[j,c("antigen.epitope", "antigen.gene", "antigen.species", "mhc.class")])
#       tmpRowB <- arrange(tmpRowB, antigen.species, antigen.gene)
#       PatientTCRdf[i,c("Beta.antigen.epitope", "Beta.antigen.gene", "Beta.antigen.species", "Beta.mhc.class")] <- unlist(lapply(tmpRowB, paste, collapse = "/"))
#     }
#   }
# }
# 
# save(PatientTCRdf, file = "PatientTCRdf.Rdata")

load(file = "PatientTCRdf.Rdata")

# Check whether an entire TRA-TRB was annotated
PatientTCRdf[which((!PatientTCRdf$Alpha.antigen.species == "") & (!PatientTCRdf$Beta.antigen.species == "")),]

AlluvialPlotA <- PatientTCRdf[which(!PatientTCRdf$Alpha.antigen.epitope == ""), c("CloneType", "AlphaChain", "Clones", "Alpha.antigen.gene", "Alpha.antigen.species", "Alpha.mhc.class")]
AlluvialPlotB <- PatientTCRdf[which(!PatientTCRdf$Beta.antigen.epitope == ""), c("CloneType", "BetaChain", "Clones", "Beta.antigen.gene", "Beta.antigen.species", "Beta.mhc.class")]

write.csv(AlluvialPlotA, file = "AlluvialPlotA.csv", row.names = F)
write.csv(AlluvialPlotB, file = "AlluvialPlotB.csv", row.names = F)

colnames(AlluvialPlotA) <- c("CloneType", "Gene usage", "Clones", "Antigen.gene", "Antigen.species", "MHC.class") -> colnames(AlluvialPlotB)

# It is unnecessary to visualize MHC class (almost all of them are MHCI)

AlluvialPlotA <- AlluvialPlotA[,!colnames(AlluvialPlotA) == "MHC.class"]
AlluvialPlotB <- AlluvialPlotB[,!colnames(AlluvialPlotB) == "MHC.class"]

AlluvialPlotDf <- rbind(as.data.frame(lapply(AlluvialPlotA, rep, AlluvialPlotA$Clones)),
                        as.data.frame(lapply(AlluvialPlotB, rep, AlluvialPlotB$Clones)))

AlluvialPlotDf <- AlluvialPlotDf[-which(colnames(AlluvialPlotDf) == "Clones")]
AlluvialPlotDf <- to_lodes_form(data.frame(AlluvialPlotDf), key = "Xaxes", axes = 1:ncol(AlluvialPlotDf)) # Here axes defined which column(s) to be used as x axis labels
AlluvialPlotDf$Colors <- rep(AlluvialPlotDf$stratum[AlluvialPlotDf$Xaxes == "CloneType"], length(table(AlluvialPlotDf$Xaxes)))

VDJAnnoPalette <- colorRampPalette(colors = c("#EFF3FF", "#6BAED6", "#2171B5", "#BDD7E7"), bias = 0.9)(length(unique(AlluvialPlotDf$Colors)))

Divide <- round(length(unique(AlluvialPlotA$CloneType))/6)
SortedIdx <- arrange(AlluvialPlotA, desc(Clones))$CloneType

RearrangeIdx <- SortedIdx[-c(1:10)]

for (i in 1:10) {
  AfterIndex <- i*Divide
  RearrangeIdx <- append(RearrangeIdx, SortedIdx[i], after = AfterIndex)
}

names(VDJAnnoPalette) <- RearrangeIdx

svg(filename = "VDJAnnotation.svg", width = 6, height = 3)
ggplot(data = AlluvialPlotDf,
       aes(x = Xaxes, stratum = stratum, alluvium = alluvium, fill = Colors,
           label = stratum)) +
  geom_flow(alpha = .7, width = 1/2, decreasing = FALSE) +
  geom_stratum(fill = 'white', color = 'skyblue', alpha= .7, width = 1/2, size = 0.2, decreasing = FALSE)  + theme_minimal() +
  ggfittext::geom_fit_text(stat = "stratum", min.size = 0, width = 1/2, grow = TRUE, decreasing = FALSE) +
  scale_fill_manual(values = VDJAnnoPalette) +
  theme(legend.position = "None", axis.title.x = element_blank(),
        axis.text.y =  element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
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

ImmuneInhibitors <- c("LAG3", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "SLAMF6", "BTLA", "TIM3")

ImmuneMarkers <- grep(paste(c(ImmuneInhibitors, "CD8A", "CD8B", "^CCL", "^CXCL", "^CXCR", "^CCR",
                              "^TNF", "^TGF", "^IFN", "^GZM", "^STAT", "^PRF", "^DEFB", 
                              "^CAMP", "^LCN2", "CD28", "CD40L", ExactMatch("PTPRC"),
                              paste("^IL", c(1:50), sep = "")), collapse = "|"), 
                      rownames(TNK.RNA), value = T)

DiffImmuneMarkers <- list()

# For visualization, a log2FC curoff of 0.4 was set

for (CellType in c(names(table(TNK.ATAC$predictedGroup_Co)), "NK")) {
  DiffImmuneMarkers[[CellType]] <- rownames(DEG[[CellType]][rownames(DEG[[CellType]]) %in% c(PotentialBioMarkers, ImmuneMarkers) &
                                                              abs(DEG[[CellType]]$avg_log2FC) > 0.4,])
}

DiffImmuneMarkers <- Reduce(union, list(DiffImmuneMarkers))
DiffImmuneMarkers <- unique(unlist(DiffImmuneMarkers))

levels(DiffImmuneMarkers) <- rev(c("CD8A", "CD8B", 
                                   "LAG3", "TIGIT", "CTLA4", 
                                   "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1",
                                   "STAT2", "STAT3", "STAT4", "STAT5B", "STAT6", 
                                   "CCL3", "CCL4", "CCL4L2", "CCL5", "CCR6","CCR7", "CXCR3",
                                   "IFNG", "IFNG-AS1", "IFNAR2", "IFNGR1",
                                   "IL16", "IL32", "IL2RA", "IL2RB", "IL2RG", "IL18R1", "IL18RAP",
                                   "IL6R", "IL7R", "IL6ST",
                                   "TGFB1", "TGFBR1", "TGFBR2", "TGFBR3",
                                   "TNFAIP3", "TNFRSF1A", "TNFRSF14","TNFRSF25","TNFAIP8",   
                                   "CD28", "PTPRC"))  

AnalyseSubTypes <- c(names(table(TNK.ATAC$predictedGroup_Co)), "NK")
levels(AnalyseSubTypes) <- c("Naive CD4", "CD4", "Treg", 
                             "Naive CD8", "CD8", "NKT", "NK")

svg(filename = "TNKBubble.svg", width = 5.5, height = 9)
BubblePlot(TNK.RNA[,TNK.RNA$orig.ident %in% c("Control1PBMC", "Control2PBMC", "PatientPBMC")], Features = DiffImmuneMarkers,
           CellTypes = AnalyseSubTypes, ColData = "IntermediateGranularity", DiffRes = DEG, colorlow = "#3587A4", colorhigh = "#DAE3FE")
dev.off()

svg(filename = "RelativeFraction_TNK.svg", width = 6, height = 6)
RelativeFraction(TNK.RNA, ColData = "FineGranularity")
dev.off()

svg(filename = "CD28.svg", width = 8, height = 4)
FeaturePlot(TNK.RNA, features = "CD28", split.by = "Group")
dev.off()

svg(filename = "ICOS.svg", width = 8, height = 4)
FeaturePlot(TNK.RNA, features = "ICOS", split.by = "Group")
dev.off()

svg(filename = "CD40LG.svg", width = 8, height = 4)
FeaturePlot(TNK.RNA, features = "CD40LG", split.by = "Group")
dev.off()

svg(filename = paste("GSEA/", "GSEACD4.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$CD4, Descriptions = c("negative regulation of cell cycle phase transition",
                                           "T cell mediated cytotoxicity", 
                                           "cellular defense response"), 
             base_size = 15, legendposition = c(0.35, 0.2), 
             rel_heights = c(1.5, 0.5, 1.0), GSEAColor = c("#77C6DD", "#B4BDDE", "#C1EFFF"))
dev.off()

svg(filename = paste("GSEA/", "GSEACD8.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$CD8, Descriptions = c("negative regulation of cell cycle phase transition",
                                            "positive regulation of T cell mediated cytotoxicity",
                                            "cellular defense response"),
             base_size = 15,
             legendposition = c(0.38, 0.2), rel_heights = c(1.5, 0.5, 1.0), GSEAColor = c("#A8ABEA", "#6574FF", "#B2E6EC"))
dev.off()

svg(filename = paste("GSEA/", "GSEANaiveCD4.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$`Naive CD4`, Descriptions = c("alpha-beta T cell differentiation",
                                                    "lymphocyte differentiation"), 
             base_size = 15,
             legendposition = c(0.25, 0.2), 
             rel_heights = c(1.5, 0.5, 1.0), GSEAColor = c("#B2E6EC", "#9ECFF8"))
dev.off()

svg(filename = paste("GSEA/", "GSEANaiveCD8.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$`Naive CD8`, Descriptions = c("alpha-beta T cell differentiation",
                                                    "lymphocyte differentiation"), 
             base_size = 15,
             legendposition = c(0.25, 0.2), 
             rel_heights = c(1.5, 0.5, 1.0), GSEAColor = c("#A8ABEA", "#6F9DB9"))
dev.off()

svg(filename = paste("GSEA/", "GSEATreg.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$Treg, Descriptions = c("response to transforming growth factor beta"), 
             base_size = 15, legendposition = c(0.35, 0.2), 
             rel_heights = c(1.5, 0.5, 1.0), GSEAColor = c("#4AB4B1"))
dev.off()

svg(filename = paste("GSEA/", "GSEANKT.svg", sep = ""), width = 7, height = 5)
ImprovedGSEA(GSEA.res$NKT, Descriptions = c("negative regulation of cell cycle phase transition",
                                            "positive regulation of T cell mediated cytotoxicity",
                                            "Fc receptor signaling pathway",
                                            "Fc receptor mediated stimulatory signaling pathway"), 
             base_size = 15, legendposition = c(0.37, 0.22), 
             rel_heights = c(1.5, 0.5, 1.0), GSEAColor = c("#6F9DB9", "#93A0BF", "#8672D0", "#CCBFFD"))
dev.off()


GSEACD4Tox <- CoreEnrichment(GSEA.res$CD4@result$core_enrichment[GSEA.res$CD4@result$Description %in% c("T cell mediated cytotoxicity", 
                                                                                                        "cellular defense response")])

GSEACD4Cyc <- CoreEnrichment(GSEA.res$CD4@result$core_enrichment[GSEA.res$CD4@result$Description %in% c("negative regulation of cell cycle phase transition")])

GSEACD8Tox <- CoreEnrichment(GSEA.res$CD8@result$core_enrichment[GSEA.res$CD8@result$Description %in% c("positive regulation of T cell mediated cytotoxicity",
                                                                                                        "cellular defense response")])

GSEACD8Cyc <- CoreEnrichment(GSEA.res$CD8@result$core_enrichment[GSEA.res$CD8@result$Description %in% c("negative regulation of cell cycle phase transition")])

GSEANCD4Dif <- CoreEnrichment(GSEA.res$`Naive CD4`@result$core_enrichment[GSEA.res$`Naive CD4`@result$Description %in% c("alpha-beta T cell differentiation", "lymphocyte differentiation")])

GSEANCD8Dif <- CoreEnrichment(GSEA.res$`Naive CD8`@result$core_enrichment[GSEA.res$`Naive CD8`@result$Description %in% c("alpha-beta T cell differentiation", "lymphocyte differentiation")])

GSEATrg <- CoreEnrichment(GSEA.res$Treg@result$core_enrichment[GSEA.res$Treg@result$Description %in% c("response to transforming growth factor beta")])

GSEANKTFc <- CoreEnrichment(GSEA.res$NK@result$core_enrichment[GSEA.res$NK@result$Description %in% c("Fc receptor signaling pathway", "Fc receptor mediated stimulatory signaling pathway")]) 

GSEANKTox <- CoreEnrichment(GSEA.res$NK@result$core_enrichment[GSEA.res$NK@result$Description %in% c("positive regulation of T cell mediated cytotoxicity")]) 

GSEANKTCyc <- CoreEnrichment(GSEA.res$NK@result$core_enrichment[GSEA.res$NKT@result$Description %in% c("negative regulation of cell cycle phase transition")]) 

################################## ATAC ########################################
ControlPalette <- TcellsPalette[names(TcellsPalette) %in% TNK.ATAC$predictedGroup_Co]
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

names(PatientPalette) <- names(TcellsPalette[names(TcellsPalette) %in% TNK.ATAC$predictedGroup_Co]) -> names(ControlPalette)
names(PatientPalette) <- paste("Patient", names(PatientPalette))
names(ControlPalette) <- paste("Control", names(ControlPalette))
ATACPalette <- c(PatientPalette, ControlPalette)

orderC <- c("Control Naive CD4",
            "Control Naive CD8",
            "Control CD4",
            "Control CD8",
            "Control NKT",
            "Control Treg")

orderP <- c("Patient Naive CD4",
            "Patient Naive CD8",
            "Patient CD4",
            "Patient CD8",
            "Patient NKT",
            "Patient Treg")

################## Identification of Positive TF-Regulators ####################
Patient.TNK.ATAC <- TNK.ATAC[TNK.ATAC$Sample == "Patient-scATAC",]
Control.TNK.ATAC <- TNK.ATAC[!TNK.ATAC$Sample == "Patient-scATAC",]

PatientPTRMatrix <- PositiveTFRegulator(Patient.TNK.ATAC)
ControlPTRMatrix <- PositiveTFRegulator(Control.TNK.ATAC)
UnionPositiveTFRegulator <- union(PatientPTRMatrix$PositiveTFs,
                                  ControlPTRMatrix$PositiveTFs)

PatientSpecificPositiveTFRegulator <- PatientPTRMatrix$PositiveTFs[!PatientPTRMatrix$PositiveTFs %in% ControlPTRMatrix$PositiveTFs]
ControlSpecificPositiveTFRegulator <- ControlPTRMatrix$PositiveTFs[!ControlPTRMatrix$PositiveTFs %in% PatientPTRMatrix$PositiveTFs]
PatientTF <- UnionPositiveTFRegulator[!UnionPositiveTFRegulator %in% ControlSpecificPositiveTFRegulator]

# Select maximum scale gene expression/motif z score across all patient and control cell sub-clusters to draw common color scale bar
GeneMatrixRescale <- cbind(PatientPTRMatrix$PositiveTFsGeneMatrix[UnionPositiveTFRegulator,],
                           ControlPTRMatrix$PositiveTFsGeneMatrix[UnionPositiveTFRegulator,])
GeneMotifRescale <- cbind(PatientPTRMatrix$PositiveTFsMotifMatrix[UnionPositiveTFRegulator,],
                          ControlPTRMatrix$PositiveTFsMotifMatrix[UnionPositiveTFRegulator,])
GeneMatrixRescale <- t(scale(t(GeneMatrixRescale)))
GeneMotifRescale <- t(scale(t(GeneMotifRescale)))
ExpLim <- max(abs(ceiling(max(apply(GeneMatrixRescale, 2, max)))),
              abs(floor(min(apply(GeneMatrixRescale, 2, min)))))
zLim <- max(abs(ceiling(max(apply(GeneMotifRescale, 2, max)))),
            abs(floor(min(apply(GeneMotifRescale, 2, min)))))

# Main function for visualization
DiagPlot <- function(GeneMatrix, MotifMatrix, Order, Color1, Color2) {
  GM <- GeneMatrix
  MM <-  MotifMatrix
  GM <- GM[, Order]
  MM <- MM[, Order]
  MM <- MM[rownames(GM),]
  UpColor <- colorRamp2(breaks = c(-ExpLim, ExpLim), colors = c("#FFFADD","#AB221F")) # Upper diagonal Red
  DnColor <- colorRamp2(breaks = c(-zLim, zLim), colors = c("#FFFADD","#3878C1")) # Lower diagonal Blue
  DiagFunc <- function(up, down){
    function(j, i, x, y, width, height, fill){
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width),
                   unit.c(y - 0.5*height, y + 0.5*height, y + 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "white"))
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width),
                   unit.c(y + 0.5*height, y - 0.5*height, y - 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "white"))
    }
  }
  rowcolors <- rep("black", length(GM))
  rowcolors[rownames(GM) %in% Color1] <- "#AB221F"
  rowcolors[rownames(GM) %in% Color2] <- "#3878C1"
  rowfonts <- rep("plain", length(GM))
  rowfonts[!rowcolors == "black"] <- "bold.italic"
  Heatmap(GM, rect_gp = gpar(type = "none"), show_column_dend = F,
          show_row_dend = F,
          show_heatmap_legend = F, cluster_rows = F,
          cluster_columns = F, row_names_gp = gpar(col = rowcolors, fontface = rowfonts),
          cell_fun = DiagFunc(up = GM, down = MM))
}

PDiag <- DiagPlot(GeneMatrix = GeneMatrixRescale[,grep(colnames(GeneMatrixRescale), pattern = "Patient")],
                  MotifMatrix = GeneMotifRescale[,grep(colnames(GeneMotifRescale), pattern = "Patient")],
                  Order = orderP, Color1 = PatientSpecificPositiveTFRegulator,
                  Color2 = ControlSpecificPositiveTFRegulator)

CDiag <- DiagPlot(GeneMatrix = GeneMatrixRescale[,grep(colnames(GeneMatrixRescale), pattern = "Control")],
                  MotifMatrix = GeneMotifRescale[,grep(colnames(GeneMotifRescale), pattern = "Control")],
                  Order = orderC, Color1 = PatientSpecificPositiveTFRegulator,
                  Color2 = ControlSpecificPositiveTFRegulator)

lgd <- list(Legend(title = "Scaled Expression",
                   col_fun = colorRamp2(breaks = c(-ExpLim, ExpLim), colors = c("#FFFADD","#AB221F")),
                   at = c(-ExpLim, 0, ExpLim),
                   direction = "horizontal"),
            Legend(title = "Scaled ChromVAR z-score",
                   col_fun = colorRamp2(breaks = c(-zLim, zLim), colors = c("#FFFADD","#3878C1")),
                   at = c(-zLim, 0, zLim),
                   direction = "horizontal"))

tmp <- PDiag + CDiag

svg(filename = "PositiveTFRegulatorDiagonal.svg", width = 6, height = 7.5)
draw(tmp, heatmap_legend_side = "top", cluster_rows = T,
     column_title = "Union of positive TF regulators (PTFR)",
     annotation_legend_list = lgd,
     annotation_legend_side = "bottom")
dev.off()

svg("Venn.svg", width = 6, height = 4)
display_venn(x = list(PatientPTR = PatientPTRMatrix$PositiveTFs,
                      ControlPTR = ControlPTRMatrix$PositiveTFs),
             col ="transparent", fill = c("#F98077","#6EBDF8"), alpha = 0.5,
             cex = 2, fontfamily = "serif", fontface = "bold", cat.cex = 0, cat.pos = 0,
             cat.dist = 0.04, cat.fontfamily = "serif", rotation.degree = 360)
dev.off()

# Create a MotifMatrix of chromVAR deviations and deviation z-scores for all motifs. 
# This data, averaged by clusters, was obtained by using the getGroupSE() function 
# which returns a SummarizedExperiment. 

seGroupMotif <- getGroupSE(ArchRProj = TNK.ATAC, useMatrix = "MotifMatrix", groupBy = "GroupSpecificPseudoBulk")

# This SummarizedExperiment was further subset to just the deviation z-scores.

DeviationZ <- seGroupMotif[rowData(seGroupMotif)$seqnames == "z",]

for (CellType in names(table(TNK.ATAC$predictedGroup_Co))) {
  svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/EnrichedMotif/Up-", CellType, ".svg", sep = ""), width = 6, height = 7)
  print(PlotEnrichedTFMotif(DET[[CellType]]$motifsUp, dotsize = 3, labelsize = 4, repulsion = 25, DeviationZ = DeviationZ, UpIn = "Patient"))
  dev.off()
  
  svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/EnrichedMotif/Down-", CellType, ".svg", sep = ""), width = 6, height = 7)
  print(PlotEnrichedTFMotif(DET[[CellType]]$motifsDown, dotsize = 3, labelsize = 4, repulsion = 25, DeviationZ = DeviationZ, UpIn = "Control"))
  dev.off()
  
  DET[[CellType]]$motifsUpTop <- GetRepresentativeTF(DET[[CellType]]$motifsUp, topTF = 20, DeviationZ = DeviationZ, UpIn = "Patient")
  DET[[CellType]]$motifsDownTop <- GetRepresentativeTF(DET[[CellType]]$motifsDown, topTF = 20, DeviationZ = DeviationZ, UpIn = "Control")
}

suppressMessages(InspectEmbedding(TNK.ATAC, QueryGene = "HOXB7",
                                  UseMatrix = getAvailableMatrices(TNK.ATAC)[3],
                                  Embed = "HarmonyUMAP"))

AllSigTFMotif <- list()
for (CellType in names(table(TNK.ATAC$predictedGroup_Co))) {
  AllSigTFMotif[[paste(CellType, "Up")]] <- names(DET[[CellType]]$motifsUpTop)
  AllSigTFMotif[[paste(CellType, "Down")]] <- names(DET[[CellType]]$motifsDownTop)
}

UnionSigTFMotif <- Reduce(union, AllSigTFMotif)

IntersectSigTFMotifUp <- Reduce(intersect, AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Up")])
IntersectSigTFMotifDown <- Reduce(intersect, AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Down")])

svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/EnrichedMotif/", "DownUpset.svg", sep = ""), width = 5, height = 3)
UpSet(make_comb_mat(AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Down")]), 
      comb_col = "#C0E2FC", bg_col = "#F0F0FF", bg_pt_col = "#1C97F4")
dev.off()

svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/EnrichedMotif/", "UUppset.svg", sep = ""), width = 5, height = 3)
UpSet(make_comb_mat(AllSigTFMotif[grep(names(AllSigTFMotif), pattern = "Up")]), 
      comb_col = "#FCC8C4", bg_col = "#FFF1F3", bg_pt_col = "#F52B1B")
dev.off()

############ Identification of Correlated Accessible Regions (CAR) #############
# FilteredPatientP2G <- list()
# 
# # nabor::knn(data = data, query = query, k = k, ...) :  Requesting more points (100) than available in cloud (94)
# # addPeak2GeneLinks require at least 100 cells for each subcluster, therefore, Treg was filtered
# 
# ValidP2G <- names(which(table(Patient.TNK.ATAC$predictedGroup_Co) > 100))
# 
# for (CellType in ValidP2G) {
#   FilteredPatientP2G[[CellType]] <- Get_p2g_fun(addPeak2GeneLinks(Patient.TNK.ATAC[Patient.TNK.ATAC$predictedGroup_Co == CellType,],
#                                                                               reducedDims = "PeakIterativeLSI")) %>% FilterP2G(TSS = TSS, Promoter = Promoter)
# }
# 
# save(FilteredPatientP2G, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/FilteredPatientP2G.Rdata")

load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/FilteredPatientP2G.Rdata")

################# Associate CAR and TF and conversion to bam ###################
# Fragment_By_CellType <- list()
# for (CellType in names(table(TNK.ATAC$GroupSpecificPseudoBulk))) {
#   SubSlot <- paste(unlist(strsplit(CellType, split = " ")), collapse = "_")
# 
#   Fragment_By_CellType[[SubSlot]] <- unlist(getFragmentsFromProject(TNK.ATAC, cellNames = rownames(TNK.ATAC[TNK.ATAC$GroupSpecificPseudoBulk == CellType,])))
#   Fragment_By_CellType[[SubSlot]]$index <- as.character(Fragment_By_CellType[[SubSlot]]$RG)
# }
# 
# library(GenomicRanges)
# 
# Output_to_bed_files(Fragment_list_cl = Fragment_By_CellType, folder = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/FragmentsByCellTypes/")
# 
# # Execute bedpetobam in /FragmentsByCellTypes, output to /BAM
# 
# bedpetobam <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   outdir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/BAM/"
#   bedpetobam <- c(bedpetobam, paste("bedtools bedpetobam -i ", i, "_fragments_cl_bamGR_pe.bed -g hg38.chrom.sizes.txt > ", outdir, i, "_fragments_cl_bamGR_pe.bam", sep = ""))
# }
# 
# bedpetobam <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/FragmentsByCellTypes/"), bedpetobam)
# 
# # Execute samtools_sort and samtools_index in /BAM, output to /sBAM
# 
# samtools_sort <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   outdir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/sBAM/"
#   samtools_sort <- c(samtools_sort, (paste("samtools sort -o ", outdir, i, "_fragments_cl_bamGR_pe_s.bam ", i, "_fragments_cl_bamGR_pe.bam", sep = "")))
# }
# 
# samtools_sort <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/BAM/"), samtools_sort)
# 
# # Execute samtools_sort and samtools_index in /sBAM, output to /sBAM
# 
# samtools_index <- c()
# for(i in names(Fragment_By_CellType)){
#   print(i)
#   samtools_index <- c(samtools_index, paste("samtools index ", i, "_fragments_cl_bamGR_pe_s.bam", sep = ""))
# }
# 
# samtools_index <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/sBAM/"), samtools_index)
# 
# write.table(c(bedpetobam, samtools_sort, samtools_index), file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/fragments2bams.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

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
# export.bed(getPeakSet(TNK.ATAC), con = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ATACorrect/AllPeaks.bed")
# 
# ATACorrect <- c()
# for(i in names(Fragment_By_CellType)){
#   Input <- paste(i, "_fragments_cl_bamGR_pe_s.bam", sep = "")
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   Genome <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/hg38.fa"
#   BlackList <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/BlackList/Hg38TotalBlackList.bed"
#   OutPut <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ATACorrect/"
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
# ATACorrect <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/sBAM/"), ATACorrect)
# 
# write.table(ATACorrect, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ATACorrectScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

############################### ScoreBigwig ####################################
# ScoreBigwig <- c()
# for(i in names(Fragment_By_CellType)){
#   Input <- paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ATACorrect/", i, "/", i, "_fragments_cl_bamGR_pe_s_corrected.bw", sep = "")
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   OutPut <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ScoreBigwig/"
#   OutFileName <- paste(i, "_footprints.bw", sep = "")
# 
#   ScoreBigwig <- c(ScoreBigwig, paste("nohup", "TOBIAS FootprintScores --signal", Input, "--regions", AllPeaks, "--output", OutFileName, "--cores 16 &"))
# }
# 
# ScoreBigwig <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ScoreBigwig/"), ScoreBigwig)
# 
# write.table(ScoreBigwig, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ScoreBigwigScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

##################### BINDetect differential TF occupancy ######################
# # If error occurs, enter ulimit -n 50000 in console
# BINDetect <- c()
# for(CellType in names(DEG)){
#   FirstGroup <- paste(c("Patient", paste(unlist(strsplit(CellType, split = " ")), collapse = "_")), collapse = "_")
#   SecondGroup <- paste(c("Control", paste(unlist(strsplit(CellType, split = " ")), collapse = "_")), collapse = "_")
# 
#   motifsMEME <- "/mnt/TRAINING/Zhongjiacheng/pwm_to_meme/Homo.cisbp.meme"
#   Genome <- "/mnt/TRAINING/Zhongjiacheng/ReferenceGenome/hg38.fa"
#   AllPeaks <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ATACorrect/AllPeaks.bed"
#   BigwigDir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ScoreBigwig/"
#   OutDir <- paste(FirstGroup, SecondGroup, sep = "_vs_")
# 
#   BINDetect <- c(BINDetect, paste("TOBIAS BINDetect", "--motifs", motifsMEME, "--signals", paste(BigwigDir, FirstGroup, "_footprints.bw", sep = ""),
#                                   paste(BigwigDir, SecondGroup, "_footprints.bw", sep = ""), "--genome", Genome, "--peaks", AllPeaks, "--outdir"
#                                   , OutDir, "--cond_names", paste(FirstGroup, SecondGroup) , "--cores 16"))
# }
# 
# BINDetect <- c(paste("cd", "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/DiffTFOccupancy/"), BINDetect)
# 
# write.table(BINDetect, file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/BINDetectScript.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

################ Visualization of differential TF binding ######################

ValidTFB <- names(table(Patient.TNK.ATAC$predictedGroup_Co))
ValidTFB <- unlist(lapply(strsplit(ValidTFB, split = " "), paste, collapse = "_"))

################ Subset the binding events to TF binding only ##################

DiffTFBind <- list()

for (CellType in ValidTFB) {
  tmp <- readxl::read_xlsx(paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/DiffTFOccupancy/", 
                                 "Patient_", CellType, "_vs_Control_", CellType, "/bindetect_results.xlsx", sep = ""), sheet = 2) %>% as.data.frame()
  tmp$clustername <- unlist(lapply(strsplit(tmp$cluster, split = "_"), '[[', 2))
  DiffTFBind[[CellType]] <- tmp[tmp$clustername %in% unlist(lapply(strsplit(names(motifPositions), split = "_"), '[[', 1)),]
}

svg(filename = paste("Volcano/", "CD4Vol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["CD4"]], 
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control CD4"]),"grey",
                                                  as.character(ATACPalette[names(ATACPalette) == "Patient CD4"])))
dev.off()

svg(filename = paste("Volcano/", "CD8Vol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["CD8"]], 
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control CD8"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient CD8"])))
dev.off()

svg(filename = paste("Volcano/", "CD8VolEOMES.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["CD8"]],showgenes = "EOMES",
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control CD8"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient CD8"])))
dev.off()

svg(filename = paste("Volcano/", "NaiveCD4Vol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["Naive_CD4"]],
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control Naive CD4"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient Naive CD4"])))
dev.off()

svg(filename = paste("Volcano/", "NaiveCD8Vol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["Naive_CD8"]],
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control Naive CD8"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient Naive CD8"])))
dev.off()


svg(filename = paste("Volcano/", "TregVol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["Treg"]],
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control Treg"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient Treg"])))
dev.off()

svg(filename = paste("Volcano/", "NKTVol.svg", sep = ""), width = 6, height = 6)
Volcanoplot(DiffTFBind[["NKT"]], 
            Sig = 10^-2, Change = 0.025, Color = c(as.character(ATACPalette[names(ATACPalette) == "Control NKT"]),"grey",
                                                   as.character(ATACPalette[names(ATACPalette) == "Patient NKT"])))
dev.off()


svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/Volcano/", "UUppset.svg", sep = ""), width = 5, height = 3)
UpSet(make_comb_mat(AllMotifUpList), comb_col = "#FCC8C4", bg_col = "#FFF1F3", bg_pt_col = "#F52B1B")
dev.off()

svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/Volcano/", "DownUpset.svg", sep = ""), width = 5, height = 3)
UpSet(make_comb_mat(AllMotifDownList), comb_col = "#C0E2FC", bg_col = "#F0F0FF", bg_pt_col = "#1C97F4")
dev.off()


AllSigTFBindUp <- unique(c(GetTopDiffBindTF(DiffTFBind[["CD4"]], Up = T),
                    GetTopDiffBindTF(DiffTFBind[["CD8"]], Up = T),
                    GetTopDiffBindTF(DiffTFBind[["Naive_CD4"]], Up = T),
                    GetTopDiffBindTF(DiffTFBind[["Naive_CD8"]], Up = T),
                    GetTopDiffBindTF(DiffTFBind[["NKT"]], Up = T),
                    GetTopDiffBindTF(DiffTFBind[["Treg"]], Up = T)))

AllSigTFBindDown <- unique(c(GetTopDiffBindTF(DiffTFBind[["CD4"]], Up = F),
                           GetTopDiffBindTF(DiffTFBind[["CD8"]], Up = F),
                           GetTopDiffBindTF(DiffTFBind[["Naive_CD4"]], Up = F),
                           GetTopDiffBindTF(DiffTFBind[["Naive_CD8"]], Up = F),
                           GetTopDiffBindTF(DiffTFBind[["NKT"]], Up = F),
                           GetTopDiffBindTF(DiffTFBind[["Treg"]], Up = F)))  

AllMotifUpList <-  list(CD4Up = GetTopDiffBindTF(DiffTFBind[["CD4"]], Up = T),
                        CD8Up = GetTopDiffBindTF(DiffTFBind[["CD8"]], Up = T),
                        NaiveCD4Up = GetTopDiffBindTF(DiffTFBind[["Naive_CD4"]], Up = T),
                        NaiveCD8Up = GetTopDiffBindTF(DiffTFBind[["Naive_CD8"]], Up = T),
                        NKTUp = GetTopDiffBindTF(DiffTFBind[["NKT"]], Up = T),
                        TregUp = GetTopDiffBindTF(DiffTFBind[["Treg"]], Up = T))

AllMotifDownList <-  list(CD4Down = GetTopDiffBindTF(DiffTFBind[["CD4"]], Up = F),
                        CD8Down = GetTopDiffBindTF(DiffTFBind[["CD8"]], Up = F),
                        NaiveCD4Down = GetTopDiffBindTF(DiffTFBind[["Naive_CD4"]], Up = F),
                        NaiveCD8Down = GetTopDiffBindTF(DiffTFBind[["Naive_CD8"]], Up = F),
                        NKTDown = GetTopDiffBindTF(DiffTFBind[["NKT"]], Up = F),
                        TregDown = GetTopDiffBindTF(DiffTFBind[["Treg"]], Up = F))

IntersectTFBindUp <- unique(Reduce(intersect, AllMotifUpList))

IntersectTFBindDown <- unique(Reduce(intersect, AllMotifDownList))

CD4SpecificUp <- AllMotifUpList$CD4Up[!AllMotifUpList$CD4Up %in% unlist(AllMotifUpList[which(!names(AllMotifUpList) == "CD4Up")])]

TregSpecificUp <- AllMotifUpList$TregUp[!AllMotifUpList$TregUp %in% unlist(AllMotifUpList[which(!names(AllMotifUpList) == "TregUp")])]

NKTSpecificUp <- AllMotifUpList$NKTUp[!AllMotifUpList$NKTUp %in% unlist(AllMotifUpList[which(!names(AllMotifUpList) == "NKTUp")])]

SpecificUpExceptNKT <- Reduce(intersect, list(AllMotifUpList$CD4Up,
                                              AllMotifUpList$CD8Up,
                                              AllMotifUpList$NaiveCD4Up,
                                              AllMotifUpList$NaiveCD8Up,
                                              AllMotifUpList$TregUp))[which(!Reduce(intersect, list(AllMotifUpList$CD4Up,
                                                                                                    AllMotifUpList$CD8Up,
                                                                                                    AllMotifUpList$NaiveCD4Up,
                                                                                                    AllMotifUpList$NaiveCD8Up,
                                                                                                    AllMotifUpList$TregUp)) %in% unlist(AllMotifUpList[which(!names(AllMotifUpList) %in% 
                                                                                                                                                               c("CD4Up", "CD8Up", "NaiveCD4Up", "NaiveCD8Up", "TregUp"))]))]



NKTSpecificDown <- AllMotifDownList$NKTDown[!AllMotifDownList$NKTDown %in% unlist(AllMotifDownList[which(!names(AllMotifDownList) == "NKTDown")])]



CDClustersDown <- intersect(AllMotifDownList$CD4Down, AllMotifDownList$NaiveCD4Down)[which(!intersect(AllMotifDownList$CD4Down, AllMotifDownList$NaiveCD4Down) %in% 
                                                                                             unlist(AllMotifDownList[which(!names(AllMotifDownList) %in% c("CD4Down", "NaiveCD4Down"))]))]

ThreeClustersDown <- Reduce(intersect, list(AllMotifDownList$CD4Down,
                                            AllMotifDownList$CD8Down,
                                            AllMotifDownList$NaiveCD4Down))[which(!Reduce(intersect, list(AllMotifDownList$CD4Down,
                                                                                                          AllMotifDownList$CD8Down,
                                                                                                          AllMotifDownList$NaiveCD4Down)) %in% unlist(AllMotifDownList[which(!names(AllMotifDownList) %in% c("CD4Down", "CD8Down", "NaiveCD4Down"))]))]
for (TF in c(IntersectTFBindUp, CD4SpecificUp, SpecificUpExceptNKT,
             TregSpecificUp, NKTSpecificUp, NKTSpecificDown,
             CDClustersDown, ThreeClustersDown
)) {
  svg(filename = paste("UMAPMotifEmbedding/", TF, ".svg", sep = "")    , width = 6, height = 6)
  print(suppressMessages(InspectEmbedding(TNK.ATAC, QueryGene = TF, 
                                    UseMatrix = getAvailableMatrices(TNK.ATAC)[3],
                                    Embed = "HarmonyUMAP")))
  dev.off()
}

############################ Normalize Signal ##################################
# # As indicated previously, addPeak2GeneLinks require at least 100 cells for each subcluster, therefore, Treg, NKT was filtered
# ValidP2G <- names(which(table(Patient.TNK.ATAC$predictedGroup_Co) > 100))
# 
# ValidP2G <- unlist(lapply(strsplit(ValidP2G, split = " "), paste, collapse = "_"))
# 
# for(i in paste("Patient", ValidP2G, sep = "_")) {
#   Check_normalized_Signal(file = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/ATACorrect/",
#                                        i, "/", i, "_fragments_cl_bamGR_pe_s_corrected.bw", sep = ""),
#                           savefile = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/NormalizedSignal/", i, sep = ""))
# }
# 
# library(motifmatchr)
# # find the motifs in the peaks
# # convert the results of motifmatchr to a GRanges object
# cisbpwm <- chromVARmotifs::human_pwms_v2
# Total_footprint_Motif <- matchMotifs(cisbpwm, getPeakSet(TNK.ATAC), genome = "hg38", out = 'positions', p.cutoff = 5e-05)
# 
# # ArchR function addMotifAnnotations was following the default motifmatchr parameter (size of window for filtration = 7)
# Total_footprint_Motif_GR <- Must_to_GR(Total_footprint_Motif)
# out_all_ext <- data.frame(cbind(names(Total_footprint_Motif), sapply(strsplit(names(cisbpwm), split = "_"), '[[', 3)))
# colnames(out_all_ext) <- c("Motif", "TFs")
# CellTypeFootPrints <- list()
# 
# for(CellType in ValidP2G)  {
#   Fragment_By_CellType <- paste(c("Patient", unlist(strsplit(CellType, split = " "))), collapse = "_")
#   NormalizedSignal <- readRDS(paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/TOBIAS/NormalizedSignal/", Fragment_By_CellType, sep = ""))
#   CellNameNoDash <- paste(unlist(strsplit(CellType, split = "_")), collapse = " ")
#   EnrichedTF <- rownames(DEG[[CellNameNoDash]][rownames(DEG[[CellNameNoDash]]) %in% out_all_ext$TFs,]) ###### This step limited the subsequent analyses to TF that deferentially expressed between patient and control
#   CellTypeFootPrints[[CellNameNoDash]] <- Calculate_footprint_celltypes(Total_footprint_Motif_GR, NormalizedSignal, EnrichedTF, out_all_ext)
#   CellTypeFootPrints[[CellNameNoDash]]$TF <- sapply(strsplit(CellTypeFootPrints[[CellNameNoDash]]$motifs, split = "_"), '[[', 3)
#   CellTypeFootPrints[[CellNameNoDash]] <- Filter_footprints(CellTypeFootPrints[[CellNameNoDash]], delta = 0.1)
# }
# 
# save(CellTypeFootPrints, Total_footprint_Motif, file = "CellTypeFootPrints.Rdata")
load(file = "CellTypeFootPrints.Rdata")

##################### MAGIC inference on gene expression #######################
library(Rmagic)
MagicRes <- list()

for (CellType in names(which(table(Patient.TNK.ATAC$predictedGroup_Co) > 100))) {
  MagicInput <- TNK.RNA[,TNK.RNA$orig.ident == "PatientPBMC" & TNK.RNA$IntermediateGranularity == CellType]
  MagicInput <-  as.matrix(Matrix::t(MagicInput[['RNA']]@data))
  ZeroExpAcrossAllSamples <- names(which(colSums(MagicInput) == 0))
  MagicInput <- MagicInput[,!colnames(MagicInput) %in% ZeroExpAcrossAllSamples]
  MagicOut <- magic(MagicInput, genes = colnames(MagicInput))$result
  MagicRes[[CellType]] <- MagicOut
}

############################# Export to Cytoscape ##############################
Cytoscape <- list()

# Original network
for (CellType in names(which(table(Patient.TNK.ATAC$predictedGroup_Co) > 100))) {
  ########### We only retain CREs regions that overlap with footprints ###########
  k1 <- which(countOverlaps(FilteredPatientP2G[[CellType]], CellTypeFootPrints[[CellType]]) > 0)
  k2 <- which(countOverlaps(CellTypeFootPrints[[CellType]], FilteredPatientP2G[[CellType]]) > 0)
  OverlapP2G <- FilteredPatientP2G[[CellType]][k1]
  OverlapFP <- CellTypeFootPrints[[CellType]][k2]

  ###################### Combine filtered GRange objects #########################
  Res_find <- data.frame(findOverlaps(OverlapFP, OverlapP2G)) # Match query (OverlapFP-footprints) and subject (OverlapP2G-peaks)
  Res_find$Motifs <- OverlapFP$motifs[Res_find$queryHits] # Transfer TF motif information to result through queryHits entry
  Res_find$TF <- OverlapFP$TF[Res_find$queryHits]
  Res_find$FootPrintRange <- as.character(OverlapFP)[Res_find$queryHits]
  Res_find$TargetGene <- OverlapP2G$Correlatedgene[Res_find$subjectHits]  # Transfer target gene information to result through subjectHits entry
  Res_find$ExpStatusTarget <- DEG[[CellType]]$avg_log2FC[match(Res_find$TargetGene, rownames(DEG[[CellType]]))] # ExpStatusTarget indicate whether the RNA expression of a target of TF is up/down regulated within a patient's cell type by comparing control counterpart
  Res_find$ExpStatusTF <- DEG[[CellType]]$avg_log2FC[match(Res_find$TF, rownames(DEG[[CellType]]))] # ExpStatusTF indicate whether the RNA expression of a TF is up/down regulated within a patient's cell type by comparing control counterpart
  Res_find$PeakRange <- as.character(OverlapP2G)[Res_find$subjectHits]
  Res_find$P2GCorrelation <- OverlapP2G$Correlation[Res_find$subjectHits] # Correlation between CRAs and genes
  Res_find$P2GType <- OverlapP2G$P2GType[Res_find$subjectHits]

  for (i in c(1:nrow(Res_find))) {
    TFName <- Res_find$TF[i]
    TargetGeneName <- Res_find$TargetGene[i]

    # TF or TargetGene might be removed due to its low expression during the pre-processing of MAGIC imputation, therefore,
    # these incomplete TF-TargetGene pairs were removed
    if (TFName %in% colnames(MagicRes[[CellType]]) &
        TargetGeneName %in% colnames(MagicRes[[CellType]])) {
      Res_find$T2GCorrelation[i] <- cor(MagicRes[[CellType]][,TFName], MagicRes[[CellType]][,TargetGeneName])
      Res_find$T2GCorrelation[i] <- as.numeric(Res_find$T2GCorrelation[i])
    } else {
      Res_find$T2GCorrelation[i] <- NA
    }
  }
  Res_find <- Res_find[!is.na(Res_find$T2GCorrelation),]
  
  Res_find$T2GCorrelation <- as.numeric(Res_find$T2GCorrelation)
  Res_find <- Res_find[abs(Res_find$T2GCorrelation) > 0.5,]
  Res_find$Regulation[Res_find$T2GCorrelation > 0] <- "Positive"
  Res_find$Regulation[Res_find$T2GCorrelation < 0] <- "Negative"
  Res_find$T2GCorrelation <- abs(Res_find$T2GCorrelation) # Convert all negative correlation to positive value for Cytoscape to draw connection strength

  ############################## Deduplicate #####################################
  # Some TFs have more than 1 motif that overlap with the peak, for these cases we only retain one of them for Cytoscape to visualize
  Combination <- paste(Res_find$TF, Res_find$TargetGene)
  Deduplicate <- !Combination %in% Combination[duplicated(Combination)]
  Res_find <- Res_find[Deduplicate,]

  # Export
  Cytoscape[[CellType]] <- Res_find[,c("TF", "ExpStatusTF", "TargetGene", "ExpStatusTarget", "T2GCorrelation", "Regulation", "P2GType")]
}

# Modify network for visualization
for (CellType in names(which(table(Patient.TNK.ATAC$predictedGroup_Co) > 100))) {
  Origin <- Cytoscape[[CellType]]
  
  if (nrow(Origin) == nrow(Origin[is.na(Origin$ExpStatusTarget),])) {
    print(paste(CellType, "network has 0 row after eliminating NA expression nodes"))
  }  else {
    Origin <- Origin[!is.na(Origin$ExpStatusTarget),]
    
    colnames(Origin) <- c("SourceNodes", "Expression", "TargetNode", "TargetAttribute", "Edge", "Regulation", "P2GType")
    Origin$IntegratedExp <- Origin$Expression
    
    if (length(which(Origin$SourceNodes == Origin$TargetNode)) == length(Origin$SourceNodes)) {
      print(paste(CellType, "all nodes are self loop"))
    } else {
      TargetNodesAsSourceNodes <- Origin
      Origin$Type <- "TF"
      TargetNodesAsSourceNodes$Type <- "Target"
      TargetNodesAsSourceNodes$SourceNodes <- TargetNodesAsSourceNodes$TargetNode
      
      AllTF <- unique(Origin$SourceNodes[Origin$Type == "TF"])
      
      TargetNodesAsSourceNodes <- TargetNodesAsSourceNodes[!TargetNodesAsSourceNodes$SourceNodes %in% AllTF,]
      # The following line adds an empty node necessary for target node to be 
      # considered as source node in Cytoscape (we can not respectively assign two color gradients 
      # for source nodes and target nodes, therefore, we listed target nodes in source node column)
      # the empty node will be manually deleted in Cytoscape after the construction of network.
      TargetNodesAsSourceNodes$TargetNode <- "" 
      
      # Some TF regulate themselves, therefore, these target nodes were set as "TF" in "TargetNodesAsSourceNodes"
      TargetNodesAsSourceNodes$Type[TargetNodesAsSourceNodes$SourceNodes %in% Origin$SourceNodes] <-  "TF"
      TargetNodesAsSourceNodes$IntegratedExp <- TargetNodesAsSourceNodes$TargetAttribute
      RetainColumn <- c("SourceNodes", "TargetNode", "Edge", "Regulation", "P2GType", "IntegratedExp", "Type")
      Cytoscape[[CellType]] <- rbind(Origin[,RetainColumn], TargetNodesAsSourceNodes[,RetainColumn])
      
      ################################ Color bar #####################################
      colmax <- round(max(Cytoscape[[CellType]]$IntegratedExp), digits = 2) + 0.01
      colmin <- round(min(Cytoscape[[CellType]]$IntegratedExp), digits = 2) - 0.01
      
      col_fun <- colorRamp2(c(colmin, 0, colmax), c("#1E64AB", "white", "#B11326"))
      
      lgd = Legend(col_fun = col_fun, title = "Node fill color: Log2 Fold Change", title_position = "topcenter",
                   legend_width = unit(6, "cm"), at = (c(colmin, 0, colmax)), direction = "h")
      
      svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/Cytoscape/", CellType, "legend.svg", sep = ""), height = 3, width = 5)
      plot.new()
      draw(lgd, x = unit(6, "cm"), y = unit(5, "cm"), just = c("left", "bottom"))
      dev.off()
      
      write.csv(Cytoscape[[CellType]], file = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/Cytoscape/", CellType, "Cytoscape.csv", sep = ""), row.names = F, quote = F)
    }
  }
}

# Inspect previously identified critical genes in the network

unique(Cytoscape[["CD4"]]$SourceNodes[Cytoscape[["CD4"]]$SourceNodes %in% AllMotifUpList[["CD4Up"]]])
unique(Cytoscape[["CD4"]]$SourceNodes[Cytoscape[["CD4"]]$SourceNodes %in% AllMotifDownList[["CD4Down"]]])
unique(Cytoscape[["CD4"]]$SourceNodes[Cytoscape[["CD4"]]$SourceNodes %in% GSEACD4Cyc])
unique(Cytoscape[["CD4"]]$SourceNodes[Cytoscape[["CD4"]]$SourceNodes %in% GSEACD4Tox])
unique(Cytoscape[["CD4"]]$SourceNodes[Cytoscape[["CD4"]]$SourceNodes %in% PatientTF])

unique(Cytoscape[["CD4"]]$SourceNodes[Cytoscape[["CD4"]]$SourceNodes %in% as.character(DiffImmuneMarkers)])


# Extract subset from large Cytoscape object

ExtractSubNetWork(CellType = "CD4", CytoscapeNetWork = Cytoscape,
                  InterestedGenes = c(GSEACD4Cyc, GSEACD4Tox, as.character(DiffImmuneMarkers), PatientTF), 
                  OutPutDir = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/Cytoscape/")

############################ Correlation analysis ##############################

CellType <- "CD8"
Gene1 <- "LAG3"
Gene2 <- "DUSP1"
SelectExpandedCells <- rownames(MagicRes[[CellType]]) %in% colnames(Immune)[Immune$orig.ident == "PatientPBMC" & !is.na(Immune$CloneSize)]
ScatterPlotDf <- data.frame(MagicRes[[CellType]][SelectExpandedCells, Gene1], MagicRes[[CellType]][SelectExpandedCells, Gene2])
colnames(ScatterPlotDf) <- c(Gene1, Gene2)

svg(filename = paste("LAG3Cor/", "DUSP1.svg"), width = 4, height = 4)
ggpubr::ggscatter(ScatterPlotDf, 
                  x = Gene1, y = Gene2,
                  add = "reg.line",  
                  add.params = list(color = "blue", fill = "lightgray"),  
                  conf.int = TRUE  
) + ggpubr::stat_cor(method = "pearson", label.x = 0.4, label.y = 0.65)
dev.off()

####

CellType <- "CD8"
Gene1 <- "LAG3"
Gene2 <- "GZMK"
SelectExpandedCells <- rownames(MagicRes[[CellType]]) %in% colnames(Immune)[Immune$orig.ident == "PatientPBMC" & !is.na(Immune$CloneSize)]
ScatterPlotDf <- data.frame(MagicRes[[CellType]][SelectExpandedCells, Gene1], MagicRes[[CellType]][SelectExpandedCells, Gene2])
colnames(ScatterPlotDf) <- c(Gene1, Gene2)

svg(filename = paste("LAG3Cor/", "GZMK.svg"), width = 4, height = 4)
ggpubr::ggscatter(ScatterPlotDf, 
                  x = Gene1, y = Gene2,
                  add = "reg.line",  
                  add.params = list(color = "blue", fill = "lightgray"),  
                  conf.int = TRUE  
) + ggpubr::stat_cor(method = "pearson", label.x = 0.4, label.y = 1.75)
dev.off()

####

CellType <- "CD8"
Gene1 <- "LAG3"
Gene2 <- "PRF1"
SelectExpandedCells <- rownames(MagicRes[[CellType]]) %in% colnames(Immune)[Immune$orig.ident == "PatientPBMC" & !is.na(Immune$CloneSize)]
ScatterPlotDf <- data.frame(MagicRes[[CellType]][SelectExpandedCells, Gene1], MagicRes[[CellType]][SelectExpandedCells, Gene2])
colnames(ScatterPlotDf) <- c(Gene1, Gene2)

svg(filename = paste("LAG3Cor/", "PRF1.svg"), width = 4, height = 4)
ggpubr::ggscatter(ScatterPlotDf, 
                  x = Gene1, y = Gene2,
                  add = "reg.line",  
                  add.params = list(color = "blue", fill = "lightgray"),  
                  conf.int = TRUE  
) + ggpubr::stat_cor(method = "pearson", label.x = 0.4, label.y = 0.75)
dev.off()

####

CellType <- "CD8"
Gene1 <- "LAG3"
Gene2 <- "GZMH"
SelectExpandedCells <- rownames(MagicRes[[CellType]]) %in% colnames(Immune)[Immune$orig.ident == "PatientPBMC" & !is.na(Immune$CloneSize)]
ScatterPlotDf <- data.frame(MagicRes[[CellType]][SelectExpandedCells, Gene1], MagicRes[[CellType]][SelectExpandedCells, Gene2])
colnames(ScatterPlotDf) <- c(Gene1, Gene2)

svg(filename = paste("LAG3Cor/", "GZMH.svg"), width = 4, height = 4)
ggpubr::ggscatter(ScatterPlotDf, 
                  x = Gene1, y = Gene2,
                  add = "reg.line",  
                  add.params = list(color = "blue", fill = "lightgray"),  
                  conf.int = TRUE  
) + ggpubr::stat_cor(method = "pearson", label.x = 0.4, label.y = 0.75)
dev.off()

####

CellType <- "CD8"
Gene1 <- "LAG3"
Gene2 <- "GZMB"
SelectExpandedCells <- rownames(MagicRes[[CellType]]) %in% colnames(Immune)[Immune$orig.ident == "PatientPBMC" & !is.na(Immune$CloneSize)]
ScatterPlotDf <- data.frame(MagicRes[[CellType]][SelectExpandedCells, Gene1], MagicRes[[CellType]][SelectExpandedCells, Gene2])
colnames(ScatterPlotDf) <- c(Gene1, Gene2)

svg(filename = paste("LAG3Cor/", "GZMB.svg"), width = 4, height = 4)
ggpubr::ggscatter(ScatterPlotDf, 
                  x = Gene1, y = Gene2,
                  add = "reg.line",  
                  add.params = list(color = "blue", fill = "lightgray"),  
                  conf.int = TRUE  
) + ggpubr::stat_cor(method = "pearson", label.x = 0.4, label.y = 0.75)
dev.off()

####

CellType <- "CD8"
Gene1 <- "LAG3"
Gene2 <- "GZMA"
SelectExpandedCells <- rownames(MagicRes[[CellType]]) %in% colnames(Immune)[Immune$orig.ident == "PatientPBMC" & !is.na(Immune$CloneSize)]
ScatterPlotDf <- data.frame(MagicRes[[CellType]][SelectExpandedCells, Gene1], MagicRes[[CellType]][SelectExpandedCells, Gene2])
colnames(ScatterPlotDf) <- c(Gene1, Gene2)

svg(filename = paste("LAG3Cor/", "GZMA.svg"), width = 4, height = 4)
ggpubr::ggscatter(ScatterPlotDf, 
                  x = Gene1, y = Gene2,
                  add = "reg.line",  
                  add.params = list(color = "blue", fill = "lightgray"),  
                  conf.int = TRUE  
) + ggpubr::stat_cor(method = "pearson", label.x = 0.4, label.y = 0.75)
dev.off()

####

CellType <- "CD8"
Gene1 <- "LAG3"
Gene2 <- "GZMM"
SelectExpandedCells <- rownames(MagicRes[[CellType]]) %in% colnames(Immune)[Immune$orig.ident == "PatientPBMC" & !is.na(Immune$CloneSize)]
ScatterPlotDf <- data.frame(MagicRes[[CellType]][SelectExpandedCells, Gene1], MagicRes[[CellType]][SelectExpandedCells, Gene2])
colnames(ScatterPlotDf) <- c(Gene1, Gene2)

svg(filename = paste("LAG3Cor/", "GZMM.svg"), width = 4, height = 4)
ggpubr::ggscatter(ScatterPlotDf, 
                  x = Gene1, y = Gene2,
                  add = "reg.line",  
                  add.params = list(color = "blue", fill = "lightgray"),  
                  conf.int = TRUE  
) + ggpubr::stat_cor(method = "pearson", label.x = 0.4, label.y = 0.75)
dev.off()

################################################################################

# Peak2Gene Legend
svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/BroswerTrack/", "P2Glegend.svg", sep = ""), height = 5, width = 5)
plot.new()

LoopMax <- 1
LoopMin <- 0.5
LoopSpan <- LoopMax - LoopMin
lgd <- Legend(col_fun = colorRamp2(c(LoopMin, LoopMin + LoopSpan/3 , LoopMax - LoopSpan/3, LoopMax), c("#E6E7E8", "#3A97FF", 
                                                         "#8816A7", "black")), title = "P2G Cor", title_position = "topcenter",
             legend_width = unit(4, "cm"), at = (c(0.5, 0.6, 0.7, 0.8, 0.9, 1)), direction = "v")
draw(lgd, x = unit(6, "cm"), y = unit(5, "cm"), just = c("left", "bottom"))
dev.off()

# Create a list containing P2G link of each celltype

AllPeak2GeneLink <- list()

for (CellType in names(which(table(Patient.TNK.ATAC$predictedGroup_Co) > 100))) {
  AllPeak2GeneLink[[CellType]] <- addPeak2GeneLinks(Patient.TNK.ATAC[Patient.TNK.ATAC$predictedGroup_Co == CellType,], reducedDims = "PeakIterativeLSI")
}

############################### KLF13-IL6ST ####################################
Gene <- "IL6ST"
TF <- "KLF13"
CellType <- "CD4"
svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/BroswerTrack/", Gene, TF, CellType, ".svg", sep = ""), width = 5.5, height = 6)
PlotImprovedBroswerTrack(AllProject = TNK.ATAC,
                         LoopList = AllPeak2GeneLink, 
                         Gene = Gene,
                         TF = TF,
                         CellType = CellType,
                         upstream = 15000,
                         downstream = 15000, useGroups = c("Patient CD4", "Control CD4", 
                                                           "Patient Naive CD4", "Control Naive CD4",
                                                           "Patient CD8", "Control CD8",
                                                           "Patient Naive CD8", "Control Naive CD8",
                                                           "Patient NK", "Control NK",
                                                           "Patient NKT", "Control NKT",
                                                           "Patient Treg", "Control Treg"),
                         pal = ATACPalette, PlotDAP = "Both", Total_footprint_Motif = Total_footprint_Motif,
                         CellTypeFootPrints = CellTypeFootPrints
)
dev.off()

################################ ETV7-CTLA4 #####################################
Gene <- "CTLA4"
TF <- "ETV7"
CellType <- "CD4"
svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/BroswerTrack/", Gene, TF, CellType, ".svg", sep = ""), width = 5.5, height = 6)
PlotImprovedBroswerTrack(AllProject = TNK.ATAC,
                         LoopList = AllPeak2GeneLink, 
                         Gene = Gene,
                         TF = TF,
                         CellType = CellType,
                         upstream = 15000,
                         downstream =15000, useGroups = c("Patient CD4", "Control CD4", 
                                                           "Patient Naive CD4", "Control Naive CD4",
                                                           "Patient CD8", "Control CD8",
                                                           "Patient Naive CD8", "Control Naive CD8",
                                                           "Patient NK", "Control NK",
                                                          "Patient NKT", "Control NKT",
                                                           "Patient Treg", "Control Treg"),
                         pal = ATACPalette, PlotDAP = "Both", Total_footprint_Motif = Total_footprint_Motif,
                         CellTypeFootPrints = CellTypeFootPrints
)
dev.off()

################################ RORA-CD28 #####################################
Gene <- "CD28"
TF <- "RORA"
CellType <- "CD4"
svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/BroswerTrack/", Gene, TF, CellType, ".svg", sep = ""), width = 5.5, height = 6)
PlotImprovedBroswerTrack(AllProject = TNK.ATAC,
                         LoopList = AllPeak2GeneLink, 
                         Gene = Gene,
                         TF = TF,
                         CellType = CellType,
                         upstream = 15000,
                         downstream =15000, useGroups = c("Patient CD4", "Control CD4", 
                                                          "Patient Naive CD4", "Control Naive CD4",
                                                          "Patient CD8", "Control CD8",
                                                          "Patient Naive CD8", "Control Naive CD8",
                                                          "Patient NK", "Control NK",
                                                          "Patient NKT", "Control NKT",
                                                          "Patient Treg", "Control Treg"),
                         pal = ATACPalette, PlotDAP = "Both", Total_footprint_Motif = Total_footprint_Motif,
                         CellTypeFootPrints = CellTypeFootPrints
)
dev.off()

################################ ETV6-LAG3 #####################################
Gene <- "LAG3"
TF <- "ETV6"
CellType <- "CD4"
svg(filename = paste("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/BroswerTrack/", Gene, TF, CellType, ".svg", sep = ""), width = 5.5, height = 6)
PlotImprovedBroswerTrack(AllProject = TNK.ATAC,
                         LoopList = AllPeak2GeneLink, 
                         Gene = Gene,
                         TF = TF,
                         CellType = CellType,
                         upstream = 15000,
                         downstream =15000, useGroups = c("Patient CD4", "Control CD4", 
                                                          "Patient Naive CD4", "Control Naive CD4",
                                                          "Patient CD8", "Control CD8",
                                                          "Patient Naive CD8", "Control Naive CD8",
                                                          "Patient NK", "Control NK",
                                                          "Patient NKT", "Control NKT",
                                                          "Patient Treg", "Control Treg"),
                         pal = ATACPalette, PlotDAP = "Both", Total_footprint_Motif = Total_footprint_Motif,
                         CellTypeFootPrints = CellTypeFootPrints
)
dev.off()

Footprints <- ImprovedPlotFootPrints(seFoot = getFootprints(
  ArchRProj = TNK.ATAC,
  positions = motifPositions[grep(paste(unique(c(SubNetworkTFs, PatientTF, AllSigTFBindUp, AllSigTFBindDown)), collapse = "|"),
                                  names(motifPositions), value = T)],
  groupBy = "GroupSpecificPseudoBulk"
), pal = ATACPalette, plot = F, force = T)

for (TF in names(Footprints)) {
  OutPutDir <- "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/PBMC/TNK/DownStream/FootPrint/"
  
  TFname <- unlist(strsplit(TF, split = "_"))[1]
  
  if (TFname %in% AllSigTFBindUp) {
    OutPutDir <- paste(OutPutDir, "AllSigTFBindUp/", sep = "")
  } else if (TFname %in% AllSigTFBindDown) {
    OutPutDir <- paste(OutPutDir, "AllSigTFBindDown/", sep = "")
  }  else {
    OutPutDir <- paste(OutPutDir, "Others/", sep = "")
  }

  svg(filename = paste(OutPutDir, TF, ".svg", sep = ""), width = 5, height = 5)
  grid::grid.newpage()
  grid::grid.draw(Footprints[[TF]])
  dev.off()
}

################################### TEST #######################################

DotPlot(TNK.RNA,  features = c("ASCL2"), split.by = "Group") + RotatedAxis() +
  scale_x_discrete("") + scale_y_discrete("")



suppressMessages(InspectEmbedding(TNK.ATAC, QueryGene = "DNMT1", 
                                  UseMatrix = getAvailableMatrices(TNK.ATAC)[3],
                                  Embed = "HarmonyUMAP"))
