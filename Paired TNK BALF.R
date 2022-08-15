######################## Paired Analysis Preprocessing1 ########################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFandPBMC/")
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# 
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/QualityControl/PreprocessedPBMCSamples.Rdata")
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/QualityControl/PreprocessedBALFSamples.Rdata")
# 
# pc.use <- 30
# 
# Paired <- merge(SeuratPatientPBMC, SeuratPatientBALF)
# 
# # Cell cycle scoring
# load("/mnt/TRAINING/Zhongjiacheng/XRY/Results/Seurat/cycle.rda")
# Paired <- NormalizeData(Paired) %>%
#   CellCycleScoring(g2m.features = g2m_genes, s.features = s_genes) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# Paired <- ScaleData(Paired)
# Paired <- RunPCA(Paired)
# 
# # #Considerable differences due to cell cycle phase were observed. Based on this plot, we would regress out the variation due to cell cycle.
# 
# svg(filename = paste("Cell_cycle_scoring.svg", sep = ""), width = 10, height = 4)
# print(DimPlot(Paired,
#               reduction = "pca",
#               group.by = "Phase",
#               split.by = "Phase") +
#         theme(plot.title = element_blank()))
# dev.off()
# 
# Paired <-RunUMAP(Paired, dims = 1:pc.use)
# 
# svg(filename = paste("PairedPBMCBALF_Before_batch_correction.svg", sep = ""), width = 6, height = 5)
# print(DimPlot(Paired, reduction = "umap", group.by = "orig.ident") +
#         theme(plot.title = element_blank()))
# dev.off()
# 
# DefaultAssay(Paired) <- "RNA"
# Paired <- DietSeurat(Paired, assays = "RNA")
# Paired.list <- SplitObject(Paired, split.by = "orig.ident")
# Paired.list <- lapply(Paired.list, FUN = function(x) {
#   x <- NormalizeData(x, verbose = T) # This step was intended to display the original normalized counts in subsequent visualization (SCTransform does not normalize the original RNA assay)
#   x <- SCTransform(x, vars.to.regress = c("percent.mt", "nCount_RNA", "S.Score", "G2M.Score"), verbose = T)
# })
# 
# integ_features <- SelectIntegrationFeatures(object.list = Paired.list,
#                                             nfeatures = 3000)
# 
# Paired.list <- PrepSCTIntegration(object.list = Paired.list,
#                                    anchor.features = integ_features)
# 
# integ_anchors <- FindIntegrationAnchors(object.list = Paired.list,
#                                         normalization.method = "SCT",
#                                         anchor.features = integ_features)
# 
# Paired.integrated <- IntegrateData(anchorset = integ_anchors,
#                                    normalization.method = "SCT")
# 
# Paired.integrated <- RunPCA(object = Paired.integrated)
# 
# save(Paired.integrated, file = "Paired.integrated.Rdata")

######################## Paired Analysis Preprocessing2 ########################
# rm(list = ls())
# setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFandPBMC/")
# library(Seurat)
# library(dplyr)
# library(ggplot2)
# 
# load("Paired.integrated.Rdata")
# 
# svg(filename = paste("Elbow_Plot.svg", sep = ""), width = 7, height = 5)
# ElbowPlot(Paired.integrated, ndims = 30) # to determine number of dimensions for clustering
# dev.off()
# 
# pc.use <- 46
# Paired.integrated <- FindNeighbors(Paired.integrated, dims = 1:pc.use, verbose = T)
# 
# Paired.integrated <- RunUMAP(Paired.integrated, dims = 1:pc.use)
# Paired.integrated <- FindClusters(Paired.integrated, resolution = 5)
# 
# DimPlot(Paired.integrated, label = T, split.by = "orig.ident")
# 
# DefaultAssay(Paired.integrated) <- "RNA"
# Cluster.markers.paired <- FindAllMarkers(Paired.integrated, only.pos = F, test.use = "MAST")
# 
# save(Paired.integrated, Cluster.markers.paired, file = "Paired_reclustering.Rdata")

############################ Cell Communication ################################
rm(list = ls())
setwd("/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFandPBMC/Downstream/")
library(Seurat)
library(dplyr)
library(ggplot2)
library(CellChat)
library(ArchR)
library(emojifont)

load(file = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/CombinedAnalyses/BALF/BALFandPBMC/Paired_reclustering.Rdata")

Cluster.markers.paired$Delta <- Cluster.markers.paired$pct.1 - Cluster.markers.paired$pct.2
top10 <- Cluster.markers.paired %>% dplyr::group_by(cluster) %>% dplyr::top_n(10, avg_log2FC)

DimPlot(Paired.integrated, label = T)

PairedPalette <- c('Naive CD8' = "#A8ABEA",
                   'Naive CD4' = "#B2E6EC",
                   'CD4' = "#77C6DD",
                   'GZMK+ CD4' = "#55C0F5",
                   'CD8' = "#0D97D5",
                   'Treg' = "#4AB4B1",
                   'NK' = "#8672D0",
                   'NKT' = "#6F9DB9",
                   'φ' = "#6BC567",
                   'CD14+ mono' = "#C9EF75",
                   'CD16+ mono' = "#73CF9F",
                   'Dendritic' = "#93E9B6",
                   'B' = "#FF7979",
                   'pDC' = "#BBF785",
                   'EMP1+ cells' = "#B6B6B6",
                   'Airway epithelial' = "#D5D5D5"
)

Paired.integrated$FineGranularity <- ""

Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("8", "10", "20", "30", "41", "46", "53", "56")] <- "Naive CD4"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("28", "40")] <- "Naive CD8"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("34", "35")] <- "Treg"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("19")] <- "GZMK+ CD4"

Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("4", "5", "6", "17", "22", "23", "47", "49", "55", "57")] <- "CD8"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("29", "39", "52")] <- "CD4"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("48")] <- "EMP1+ cells"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("51")] <- "Airway epithelial"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("50")] <- "NK"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("16", "38")] <- "NKT"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("31")] <- "Dendritic"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("0", "1", "2", "3", "7",
                                                                           "9", "11", "14", "18", "26", "27", 
                                                                           "36", "37", "42", "43", "44")] <- "φ"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("12", "25", "45")] <- "CD16+ mono"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("13", "15", "24")] <- "CD14+ mono"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("21", "32", "33")] <- "B"
Paired.integrated$FineGranularity[Paired.integrated$seurat_clusters %in% c("54")] <- "pDC" #LILRA4

Cluster2Unicode <- readxl::read_xlsx(path = "/mnt/TRAINING/Zhongjiacheng/XRY/Results/Dec2UnicodeSymbol.xlsx")
NumCluster <- length(table(Paired.integrated$seurat_clusters)) - 1

Reanno.cluster.ids <- Paired.integrated$FineGranularity[match(as.character(0:NumCluster), table = Paired.integrated$seurat_clusters)]
names(Reanno.cluster.ids) <- levels(Paired.integrated@active.ident)
Paired.integrated <- RenameIdents(Paired.integrated, Reanno.cluster.ids)

svg(filename = "PairedFineGranularity.svg", width = 8.5, height = 6.5)
tmp <- FetchData(Paired.integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident", "FineGranularity"))

tmp$Labels <- match(tmp$FineGranularity, names(PairedPalette))

LabelPositions <- data.frame(UMAP1 = data.frame(UMAP1 = tapply((tmp[,"UMAP_1"]), INDEX = tmp[,"Labels"], FUN = mean)),
                             UMAP2 = data.frame(UMAP2 = tapply((tmp[,"UMAP_2"]), INDEX = tmp[,"Labels"], FUN = mean)))
LabelPositions$color <- "black"

ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = FineGranularity)) +
  geom_point(size = 0.2, alpha = 1) +
  scale_color_manual(values = PairedPalette[names(PairedPalette) %in% tmp$FineGranularity]) +
  annotate(geom = "text", x = LabelPositions$UMAP1, y = LabelPositions$UMAP2, label = rownames(LabelPositions),
           color = LabelPositions$color, size = 7,  fontface = "bold") +
  theme_ArchR() + theme(legend.position = "right", legend.title = element_blank(),
                        legend.key.size = unit(1, 'cm'),
                        legend.text = element_text(size = 17),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 9,
                                                  shape = stringi::stri_unescape_unicode(paste0("\\u", Cluster2Unicode$Unicode))[1:length(table(tmp$FineGranularity))]
  )))
dev.off()


svg(filename = "PairedFineGranularity Split.svg", width = 6, height = 4)
tmp <- FetchData(Paired.integrated, vars = c("UMAP_1", "UMAP_2", "orig.ident", "FineGranularity"))

tmp$Labels <- match(tmp$FineGranularity, names(PairedPalette))

ggplot(tmp, aes(x = UMAP_1, y = UMAP_2, color = FineGranularity)) +
  geom_point(size = 0.2, alpha = 1) + facet_wrap(~orig.ident) +
  scale_color_manual(values = PairedPalette[names(PairedPalette) %in% tmp$FineGranularity]) +
  theme_void() +
  theme(legend.position = "none", legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())
dev.off()

#---------------------> Seurat-Style stacked violin plot <----------------------
svg(filename = paste("Paired_Cell_Violin.svg", sep = ""), width = 6, height = 4.5)
Features <- c("CD3E", "CD4", "CD8A", "GZMK", "LEF1", "FOXP3", "KLRF1", "CD14", "FCGR3A", "MARCO", "CD1C", "MS4A1", "LILRA4",
              "EMP1", "SCGB3A1")

long <- FetchData(Paired.integrated, vars = Features, slot = "data")
long$MainType <- Paired.integrated@active.ident
long$Cell <- rownames(long)
long <- reshape2::melt(long, id.vars = c("Cell","MainType"), measure.vars = Features,
                       variable.name = "Feat", value.name = "Expr")

long$MainType <- factor(x = long$MainType, levels = rev(names(PairedPalette)))
long$Color <- PairedPalette[match(long$MainType, table = names(PairedPalette))]

ggplot(long, aes(Expr, factor(MainType))) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, aes(fill = MainType, color = MainType)) +
  scale_fill_manual("legend", values = PairedPalette) +
  scale_color_manual("legend", values = PairedPalette) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(Feat), scales = "free")  +
  cowplot::theme_cowplot(font_size = 15) +
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

Paired.integrated$GroupCluster <- paste(unlist(lapply(strsplit(Paired.integrated$orig.ident, split = "Patient"), '[[', 2)),
                                        Paired.integrated$FineGranularity)

RetainGroupCluster <- names(which(table(Paired.integrated$GroupCluster) > 100))

################################## Cell Chat ###################################
# cellchat <- createCellChat(object = Paired.integrated[,Paired.integrated$GroupCluster %in% RetainGroupCluster]@assays$RNA@data,
#                            meta = Paired.integrated[,Paired.integrated$GroupCluster %in% RetainGroupCluster]@meta.data,
#                            group.by = "GroupCluster")
# 
# cellchat <- setIdent(cellchat, ident.use = "GroupCluster")
# 
# levels(cellchat@idents)
# 
# groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
# 
# ################# Set the ligand-receptor interaction database #################
# 
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)
# 
# # use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# cellchat@DB <- CellChatDB.use # set the used database in the object
# 
# cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
# 
# future::plan("multiprocess", workers = 16)
# 
# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# 
# # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)
# 
# # Compute the communication probability and infer cellular communication network
# cellchat <- computeCommunProb(cellchat)
# 
# # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# 
# # Infer the cell-cell communication at a signaling pathway level
# cellchat <- computeCommunProbPathway(cellchat)
# 
# # Calculate the aggregated cell-cell communication network
# cellchat <- aggregateNet(cellchat)
# 
# save(cellchat, file = "cellchat.Rdata")
load(file = "cellchat.Rdata")

library(colorspace)

cols1 <- readhex(file = textConnection(paste(PairedPalette, collapse = "\n")),
                 class = "RGB")

#transform to hue/lightness/saturation colorspace

cols1 <- as(cols1, "HLS")
PBMCPalette <- cols1 -> BALFPalette

#multiplicative decrease of lightness

PBMCPalette@coords[, "L"] <- PBMCPalette@coords[, "L"] * 1.1
BALFPalette@coords[, "L"] <- BALFPalette@coords[, "L"] * 0.5

#going via rgb seems to work better

cols1 <- as(cols1, "RGB") %>% hex()
PBMCPalette <- as(PBMCPalette, "RGB") %>% hex()
BALFPalette <- as(BALFPalette, "RGB") %>% hex()

names(BALFPalette) <- paste("BALF", names(PairedPalette))
names(PBMCPalette) <- paste("PBMC", names(PairedPalette))
CellChatPalette <- c(PBMCPalette, BALFPalette)

UseCellNames <- names(table(cellchat@idents))
names(UseCellNames) <- CellChatPalette[match(UseCellNames, names(CellChatPalette))]

# Use the following command to inspect  all available ligand-receptor interaction pathways
cellchat@netP$pathways


svg(filename = paste("Paired_CellCommuniationChordICOS.svg", sep = ""), width = 6, height = 8)
netVisual_aggregate(cellchat, signaling = "ICOS", layout = "chord", color.use = names(UseCellNames))
dev.off()

svg(filename = paste("Paired_CellCommuniationChordMHCI.svg", sep = ""), width = 6, height = 8)
netVisual_aggregate(cellchat, signaling = "MHC-I", layout = "chord", color.use = names(UseCellNames))
dev.off()


svg(filename = paste("Paired_CellCommuniationChordMHCII.svg", sep = ""), width = 6, height = 8)
netVisual_aggregate(cellchat, signaling = "MHC-II", layout = "chord", color.use = names(UseCellNames))
dev.off()




netVisual_bubble(cellchat, sources.use = grep(" B", UseCellNames), 
                 targets.use = grep(paste(c("GZMK+ CD4", "CD4", " CD8", "NKT", "Dendritic", "φ"), collapse = "|"), UseCellNames),
                 remove.isolate = FALSE)







svg(filename = paste("legend.svg", sep = ""), width = 12, height = 8)
netVisual_chord_gene(cellchat, 
                     sources.use = grep(" B", UseCellNames), 
                     targets.use = grep(paste(c("GZMK+ CD4", "CD4", " CD8", "NKT", "Dendritic", "φ"), collapse = "|"), UseCellNames), 
                     lab.cex = 0.3, title.name = "Signaling from B to T", 
                     color.use = names(UseCellNames), legend.pos.x = 220,  legend.pos.y = 120)
dev.off()


svg(filename = paste("BtoT.svg", sep = ""), width = 5, height = 5)
netVisual_chord_gene(cellchat, 
                     sources.use = grep(" B", UseCellNames), 
                     targets.use = grep(paste(c("GZMK+ CD4", "CD4", " CD8", "NKT"), collapse = "|"), UseCellNames), 
                     lab.cex = 0.125, title.name = "",
                     color.use = names(UseCellNames), show.legend = FALSE)
dev.off()

svg(filename = paste("DtoT.svg", sep = ""), width = 5, height = 5)
netVisual_chord_gene(cellchat, 
                     sources.use = grep("Dendritic", UseCellNames), 
                     targets.use = grep(paste(c("GZMK+ CD4", "CD4", " CD8", "NKT"), collapse = "|"), UseCellNames), 
                     lab.cex = 0.125, title.name = "Signaling from BALF to PBMC", 
                     color.use = names(UseCellNames), show.legend = FALSE)
dev.off()

svg(filename = paste("φtoT.svg", sep = ""), width = 5, height = 5)
netVisual_chord_gene(cellchat, 
                     sources.use = grep("φ", UseCellNames), 
                     targets.use = grep(paste(c("GZMK+ CD4", "CD4", " CD8", "NKT"), collapse = "|"), UseCellNames), 
                     lab.cex = 0.125, title.name = "Signaling from BALF to PBMC", 
                     color.use = names(UseCellNames), show.legend = FALSE)
dev.off()





svg(filename = paste("Paired_CellCommuniationChordMHCII.svg", sep = ""), width = 6, height = 8)
netVisual_aggregate(cellchat, signaling = "MHC-II", layout = "chord", color.use = names(UseCellNames))
dev.off()

svg(filename = paste("Paired_CellCommuniationChordTIGIT.svg", sep = ""), width = 6, height = 8)
netVisual_aggregate(cellchat, signaling = "TGFb", layout = "chord", color.use = names(UseCellNames))
dev.off()

svg(filename = paste("S100A8.svg", sep = ""), width = 4, height = 4)
FeaturePlot(Paired.integrated, features = "S100A8", split.by = "orig.ident")
dev.off()



netVisual_aggregate(cellchat, signaling = "IL16", layout = "chord", color.use = names(UseCellNames))



netAnalysis_contribution(cellchat, signaling = "IL16")


FeaturePlot(Paired.integrated, features = "LILRB1")


netVisual_aggregate(cellchat, signaling = "TIGIT", layout = "chord", color.use = names(UseCellNames))





netVisual_aggregate(cellchat, signaling = "MIF", layout = "chord", color.use = names(UseCellNames))



netAnalysis_contribution(cellchat, signaling = "MHC-II")


FeaturePlot(Paired.integrated, features = c("TIGIT", "PVR"))


FeaturePlot(Paired.integrated, features = "IL6ST")






svg(filename = paste("Paired_CellCommuniationBubble.svg", sep = ""), width = 4, height = 8)
netVisual_bubble(cellchat, sources.use = which(UseCellNames == "BALF CD14+ mono"), 
                 targets.use = grep(UseCellNames, pattern = "PBMC"), remove.isolate = FALSE)
dev.off()


save(cellchat, file = "cellchat.Rdata")


cellchat@meta$GroupCluster
