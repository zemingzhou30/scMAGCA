library(ggplot2)
library(Seurat)

# fig1 阿尔兹海默症
#
# --- fig6 a -----
load("/home/zhangpeiru/Rworkspace/scmgca/code/data/Seurat_Anno.RData")
table(obj$subtype)
colors <- c("Astrocyte" = "#009eff",
            "Excitatory neurons-1" = "#8dc691",
            "Excitatory neurons-2" = "#637a50",
            "Excitatory neurons-3" = "#64a83d",
            "Excitatory neurons-4" = "#2CA02C",
            "Grm8+ neurons" = "#9e6975",
            "Inhibitory neurons-1" = "#1d4a9d",
            "Inhibitory neurons-2" = "#205479",
            "Microglia-1" = "#f397c0",
            "Microglia-2" = "#ed5a9e",
            "Oligodendrocyte" = "#9467BD",
            "Oligodendrocyte precursor cells" = "#7200DA",
            "Pericytes endothelial" = "#8db5ce",
            "Spiny projection neurons" = "#cc7816",
            "Type I spiral ganglion neuron" = "#ee4e21")

colors <- c("Astrocyte" = "#AEC7E8",
            "Excitatory neurons-1" = "#C5B0D5",
            "Excitatory neurons-2" = "#9467BD",
            "Excitatory neurons-3" = "#766bd3",
            "Excitatory neurons-4" = "#7200DA",
            "Grm8+ neurons" = "#17BECF",
            "Inhibitory neurons-1" = "#C49C94",
            "Inhibitory neurons-2" = "#8C564B",
            "Microglia-1" = "#FFBB78",
            "Microglia-2" = "#FF7F0E",
            "Oligodendrocyte" = "#2CA02C",
            "Oligodendrocyte precursor cells" = "#E376C2",
            "Pericytes endothelial" = "#009eff",
            "Spiny projection neurons" = "#BCBD22",
            "Type I spiral ganglion neuron" = "#D62728")

pdf("/home/zhangpeiru/Rworkspace/scmgca/code/aer/fig6a.pdf",width = 7.5,height = 5)
DimPlot(
  object = obj,
  group.by = 'subtype',
  label = F,
  cols = colors,
  pt.size = 0.5,
  label.size = 4,
  repel = TRUE) + ggtitle('')
dev.off()


# fig6 f -----
# install.packages("E:/200files/240packages/CNEr_1.42.0.zip",type = "source",repos = NULL)
# devtools::install_github("GreenleafLab/chromVARmotifs")
# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
# install.packages("E:/200files/240packages/hexbin_1.28.5.zip",repos = NULL,type = "source")
# ArchR::installExtraPackages()
proj_sub <- readRDS("/home/zhangpeiru/Rworkspace/scmgca/code/data/proj_sub_2025_03_21_12h.rds")
umap_coords <- getEmbedding(proj_sub, embedding = "UMAP")

# 提取细胞类型信息
cell_types <- getCellColData(proj_sub, select = "celltype_New")

# 创建一个数据框
plot_data <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  cellType = cell_types
)

col <- c("0: Astrocyte" = "#AEC7E8"
         ,"1: Excitatory neurons" = "#766bd3"
         ,"2: Oligodendrocyte" = "#2CA02C"
         ,"3: Spiny projection neurons" = "#BCBD22"
         ,"4: Inhibitory neurons" = "#C49C94"
         ,"5: Grm8+ neurons" = "#17BECF"
         ,"6: Oligodendrocyte precursor cells" = "#E376C2"
         ,"7: Excitatory neurons" = "#766bd3"
         ,"8: Pericytes endothelial" = "#009eff"
         ,"9: Microglia" = "#F7B778"
         ,"10: Type I spiral ganglion neuron" = "#D62728"
         ,"11: Microglia" = "#F7B778"
         ,"12: Excitatory neurons" = "#766bd3"
         ,"13: Excitatory neurons" = "#766bd3"
         ,"14: Grm8+ neurons"  = "#17BECF")


col <- c("0: Astrocyte" = "#AEC7E8"
         ,"1: Excitatory neurons" = "#C5B0D5"
         ,"2: Oligodendrocyte" = "#2CA02C"
         ,"3: Spiny projection neurons" = "#BCBD22"
         ,"4: Inhibitory neurons" = "#C49C94"
         ,"5: Grm8+ neurons" = "#8C564B"
         ,"6: Oligodendrocyte precursor cells" = "#E376C2"
         ,"7: Excitatory neurons" = "#9467BD"
         ,"8: Pericytes endothelial" = "#009eff"
         ,"9: Microglia" = "#FFBB78"
         ,"10: Type I spiral ganglion neuron" = "#D62728"
         ,"11: Microglia" = "#FF7F0E"
         ,"12: Excitatory neurons" = "#766bd3"
         ,"13: Excitatory neurons" = "#7200DA"
         ,"14: Grm8+ neurons"  = "#17BECF")

colors <- c("Astrocyte" = "#AEC7E8",
            "Excitatory neurons-1" = "#C5B0D5",
            "Excitatory neurons-2" = "#9467BD",
            "Excitatory neurons-3" = "#766bd3",
            "Excitatory neurons-4" = "#7200DA",
            "Grm8+ neurons" = "#17BECF",
            "Inhibitory neurons-1" = "#C49C94",
            "Inhibitory neurons-2" = "#8C564B",
            "Microglia-1" = "#FFBB78",
            "Microglia-2" = "#FF7F0E",
            "Oligodendrocyte" = "#2CA02C",
            "Oligodendrocyte precursor cells" = "#E376C2",
            "Pericytes endothelial" = "#009eff",
            "Spiny projection neurons" = "#BCBD22",
            "Type I spiral ganglion neuron" = "#D62728")


pt <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = celltype_New)) +
  geom_point(size = 0.6) +
  theme_minimal() +
  scale_color_manual(values = col) +
  labs(x = "UMAP1", y = "UMAP2")+
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black",linewidth = 1),
    axis.line = element_line(color = "black",linewidth = 1), # 只保留坐标轴
    legend.position = 'right',
    legend.justification = c(0, 1),
    axis.title = element_text(size = 14, face = "bold"), # 美化轴标题
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
pt
ggsave("/home/zhangpeiru/Rworkspace/scmgca/code/aer/fig6f-2.pdf", plot = pt, width = 7.5, height = 4.5)


# fig 6g

markerGenes2 <- c('Gja1','Aqp4', #Astrocyte
                  'Cx3cr1','Csf3r', #Microglia
                  'Mal','Olig1','Neat1','Olig2', #Oligodendrocyte
                  'Pdgfra','Cspg4','Vcan', #oligodendrocyte precursor cells 
                  'Gad1','Gad2',# Inhibitory neurons
                  'Slc17a7','Camk2a','Nrgn','Syn3', # Excitatory neurons
                  'Flt1','Cldn5', #pericytes Endothelial
                  'Penk','Drd2','Adora2a','Tac1', #spiny projection neurons
                  #'Snap25' #spiral ganglion neuron
                  #'Rbfox3', #Neurons
                  'Ntng1','Grm8','Prox1','Pcdh20'# Type I spiral ganglion neuron
)

proj_sub <- addImputeWeights(proj_sub)

custom_pal <- colorRampPalette(c("#008bd0", "#d3d7e7", "#ffa61d"))(100)
for (i in 1:length(markerGenes2)) {
  p <- plotEmbedding(
    ArchRProj = proj_sub, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes2[i], 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj_sub),
    pal = custom_pal
  )+ 
    # guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank()
    ) 
  pdf(paste0("/home/zhangpeiru/Rworkspace/scmgca/code/aer/fig6g_",markerGenes2[i],".pdf"))
  print(p)
  dev.off()
}


# 验证结果
table(proj_sub@cellColData@listData[["celltype_New"]])


### 鉴定Marker基因
markersGS <- getMarkerFeatures(
  ArchRProj = proj_sub, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "celltype_New",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6


markerGenes <- c('Gja1','Aqp4', #Astrocyte
                 'Cx3cr1','Csf3r', #Microglia
                 'Mal','Olig1','Neat1','Olig2', #Oligodendrocyte
                 'Pdgfra','Cspg4','Vcan', #oligodendrocyte precursor cells 
                 'Gad1','Gad2',# Inhibitory neurons
                 'Slc17a7','Camk2a','Nrgn','Syn3', # Excitatory neurons
                 'Flt1','Cldn5', #pericytes Endothelial
                 'Penk','Drd2','Adora2a','Tac1', #spiny projection neurons
                 #'Snap25' #spiral ganglion neuron
                 #'Rbfox3', #Neurons
                 'Ntng1','Grm8','Prox1','Pcdh20'# Type I spiral ganglion neuron
)

custom_pal <- colorRampPalette(c("#008bd0", "#eeeeee", "#ffa61d"))(100)
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  pal = custom_pal,
  limits = c(-3,3),
  transpose = TRUE
)

pdf("/home/zhangpeiru/Rworkspace/scmgca/code/aer/fig6h_heatmap.pdf",width = 7,height = 4.5)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()





# fig h i
## 鉴定marker peaks
markersPeaks_proj <- getMarkerFeatures(
  ArchRProj = proj_sub, 
  useMatrix = "PeakMatrix", 
  groupBy = "celltype_New",
  bias = c("TSSEnrichment", "log10(nFrags)"),  # 消除TSS富集和每个细胞的fragment数对结果的影响
  testMethod = "wilcoxon"
)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks_proj, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)


plotPDF(heatmapPeaks, 
        name = "fig6hi_Peak-Marker-Heatmap.pdf", 
        width = 8, height = 6, 
        ArchRProj = proj_sub, addDOC = FALSE,verbose = TRUE)


# 鉴定marker peaks的motif 
# 选取marker
markerList_sub <- getMarkers(markersPeaks_proj, cutOff = "FDR <= 0.01 & Log2FC >= 0.25")
# 这个差异就比较低
markerList_sub

proj_sub <- addMotifAnnotations(ArchRProj = proj_sub, motifSet = "cisbp", name = "Motif",force = TRUE)


enrichMotifs_sub <- peakAnnoEnrichment(
  seMarker = markersPeaks_proj,
  ArchRProj = proj_sub,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.25"
)

custom_pal <- colorRampPalette(c("#eeeeee", "#FF7F0E"))(100)
heatmapEM <- plotEnrichHeatmap(enrichMotifs_sub, n = 18,pal = custom_pal, transpose = TRUE)
pdf("Motifs-Enriched-Marker-Heatmap.pdf", width = 24, height = 6)
heatmapEM
dev.off()

plotPDF(heatmapEM, name = "fig6_hi_Motifs-Enriched-Marker-Heatmap.pdf", width = 24, height = 6, ArchRProj = proj_sub, addDOC = FALSE)

custom_pal <- colorRampPalette(c("#eeeeee", "#FF7F0E"))(100)
heatmapEM <- plotEnrichHeatmap(enrichMotifs_sub, n = 10,pal = custom_pal, transpose = TRUE)
pdf("Motifs-Enriched-Marker-Heatmap.pdf", width = 20, height = 6)
heatmapEM
dev.off()
plotPDF(heatmapEM, name = "fig6i_Motifs-Enriched-Marker-Heatmap_n10.pdf", width = 12, height = 6, ArchRProj = proj_sub, addDOC = FALSE)

save.image("xxx.RData")
















