
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)

BiocManager::install("EnsDb.Mmusculus.v79")



counts <- Read10X_h5("M_Kidney_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5")
fragpath <- "M_Kidney_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biovizBase")

library("biovizBase")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

dat<-read.csv("new_dat_meta.csv",row.names = 1)

pbmc<-subset(pbmc,cells=rownames(dat))

DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

DefaultAssay(pbmc) <- "RNA"

options(future.globals.maxSize = 1024 * 1024^2)  # 设置为1 GiB


pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

pbmc@meta.data$new_type<-dat$new_type
pbmc@meta.data$tu<-dat$tu
pbmc@meta.data$label<-dat$label

pbmc <- RunUMAP(object = pbmc,, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.atac')

saveRDS(pbmc,"raw_pbmc.rds")

pbmc<-SetIdent(pbmc,value = pbmc@meta.data$new_type)
da_peaks_multiple <- FindAllMarkers(
  object = pbmc,
  test.use = 'LR',  
  latent.vars = 'nCount_ATAC',
  only.pos = TRUE,
  logfc.threshold = 0.4,    
  min.pct = 0.1             
)

head(da_peaks_multiple)

library(dplyr)
top10 <- da_peaks_multiple %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

saveRDS(pbmc,"raw_all.rds")

write.csv(da_peaks_multiple,"Diff_gene_table.csv")

library(Seurat)
library(scMayoMap)

qc_data<-read.csv("mouse_kindy_10.csv") 
dat<-readRDS("mouse_RNA_new.rds")
head(dat)
table(dat@meta.data$label)
dat<-SetIdent(dat,value=dat@meta.data$label)
pbmc.markers <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #寻找差异表达基因

head(pbmc.markers)
table(pbmc.markers$cluster)

scMayoMap.obj <- scMayoMap(data = pbmc.markers, database=scMayoMapDatabase,tissue = "kidney") #进行自动注释
res <- scMayoMap.obj$res
head(res)

scMayoMap.plot(scMayoMap.object = scMayoMap.obj)
table(dat@meta.data$tu,dat@meta.data$label)

pdf("new_annotation_plot.pdf",width=10,height=4.8)
table(dat@meta.data$tu,dat@meta.data$label)
scMayoMap.plot(scMayoMap.object = scMayoMap.obj)
dev.off()

cell_table <-as.data.frame(table(dat@meta.data$tu,dat@meta.data$label))
rna_data$Freq <- rna_data$Freq / ave(rna_data$Freq, rna_data$Var2, FUN = sum)
rna_data$Freq <- dat$Freq / ave(dat$Freq, dat$Var2, FUN = sum)
cell_table$Freq <- cell_table$Freq / ave(cell_table$Freq, cell_table$Var2, FUN = sum)

dat@meta.data$new_type<-ifelse(dat@meta.data$label=="0","Endothelial cell 1",ifelse(dat@meta.data$label=="1","Loop of Henle cell",ifelse(dat@meta.data$label=="2","Proximal tubule cell 1",ifelse(dat@meta.data$label=="3","Endothelial cell 2",ifelse(dat@meta.data$label=="4","Distal convoluted tubule cell",ifelse(dat@meta.data$label=="5","Proximal tubule cell 2",ifelse(dat@meta.data$label=="6","Podocyte",ifelse(dat@meta.data$label=="7","Proximal tubule cell 3",ifelse(dat@meta.data$label=="8","Fibroblast","Intercalated cell"))))))))) #对cluster进行重新注释
head(dat)

saveRDS(dat,"new_mouse_RNA_new.rds")


library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
#options(stringsAsFactors = F)
data_obj <- readRDS("new_mouse_RNA_new.rds")

unique(data_obj$seurat_annotations)

data_obj<-SetIdent(data_obj,value=data_obj$new_type)


high<-subset(data_obj,idents=c("Distal convoluted tubule cell","Fibroblast","Intercalated cell","Loop of Henle cell","Podocyte","Proximal tubule cell 1","Proximal tubule cell 2","Proximal tubule cell 3"))
table(high@meta.data$new_type)

high<-SetIdent(high,value=high$new_type)

CellChatDB <- CellChatDB.mouse
cellchat <- createCellChat(object = high,
                           meta = high@meta.data,
                           group.by = "new_type")
cellchat
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat@data.project[1:4,1:4]
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,"prox_cellchat.rds")



# fig7 a----
library(Seurat)
library(devtools)
#install_github("sajuukLyu/ggunchull", type = "source")  # 从GitHub拉取并安装
#install_github("yikeshu0611/tidydr")  
library(ggunchull)
library(ggplot2)
library(tidydr)
library(dplyr)

umap <- read.csv("mouse_kindy_10.csv")
meta <- read.csv("new_dat_meta.csv")

merged <- merge(umap, meta, by = "X", all = F)
colors <- c("Proximal tubule cell 1" = "#C5B0D5",
            "Loop of Henle cell" = "#AEC7E8",
            "Distal convoluted tubule cell" = "#17BECF",
            "Endothelial cell 1" = "#FFBB78",
            "Endothelial cell 2" = "#FF7F0E",
            "Fibroblast" = "#D62728",
            "Neutrophil" = "#C49C94",
            "Intercalated cell" = "#2CA02C",
            "Principal cell" = "#8C564B",
            "Podocyte" = "#E376C2",
            "Proximal tubule cell 2" = "#9467BD",
            "Macrophage" = "#BCBD22",
            "Proximal tubule cell 3" = "#7200DA")

celltype <- c("Proximal tubule cell 1",
              "Proximal tubule cell 3",
              "Proximal tubule cell 2",
              "Endothelial cell 1",
              "Endothelial cell 2")

p <- ggplot(merged, aes(x = V1, y = V2, 
                 fill = new_type, color = new_type)) +
  lapply(celltype, function(ct) {
    stat_unchull(data = subset(merged, new_type == ct),
                 alpha = 0, 
                 linewidth = 0.75,
                 show.legend = FALSE,
                 color = "black",     
                 linetype = "dashed", 
                 delta = 0.5, 
                 th = 0.2,
                 n = 10)
  }) +
  geom_point(size = 0.5,alpha = 0.5) +
  theme_dr() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1)
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  guides(fill = "none") 
p
dev.off()
ggsave("kidney/fig7_a.pdf",p,width = 8.5,height = 6)


# fig7 b
# 热图----
colors <- c("Proximal tubule cell 1" = "#C5B0D5",
            "Loop of Henle cell" = "#AEC7E8",
            "Distal convoluted tubule cell" = "#17BECF",
            "Endothelial cell 1" = "#FFBB78",
            "Endothelial cell 2" = "#FF7F0E",
            "Fibroblast" = "#D62728",
            "Neutrophil" = "#C49C94",
            "Intercalated cell" = "#2CA02C",
            "Principal cell" = "#8C564B",
            "Podocyte" = "#E376C2",
            "Proximal tubule cell 2" = "#9467BD",
            "Macrophage" = "#BCBD22",
            "Proximal tubule cell 3" = "#7200DA")

new <- readRDS("new_mouse_RNA_new.rds")
new <- SetIdent(new, value = new$new_type)
all.marker <- read.csv("end_all_markers.csv",row.names = 1)
top10 <- all.marker %>% group_by(cluster) %>% top_n(n = 100,wt = avg_log2FC)
#pdf("kidney/fig7b.pdf",width = 6,height = 6)
DoHeatmap(new, features = top10$gene, group.colors = colors,label = F,) +
  scale_fill_gradient2(low = "#008bd0", mid = "#dddddc", high = "#ff5743", midpoint = 0) +
  theme(axis.text.y = element_blank())
#dev.off()


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
atac <- readRDS("raw_all.rds")
all.marker <- read.csv("Diff_gene_table.csv",row.names = 1)
top10 <- all.marker %>% group_by(cluster) %>% top_n(n = 100,wt = avg_log2FC)

atacdata <- atac@assays$ATAC$data 
atacdata_df <- data.frame(atacdata)
atacdata_df$gene <- rownames(atacdata_df)
atacheatmap <- left_join(top10, atacdata_df, by = "gene")

atacdata <- atac@meta.data$new_type
group <- data.frame(group = atac@meta.data$new_type)
rownames(group) <- rownames(atac@meta.data)
rownames(group) <- gsub("-", ".", rownames(group))
length(intersect(colnames(atacheatmap[,8:9656]),rownames(group)))

custom_colors <- colorRampPalette(c("#008bd0", "#eeeeee", "#ffa61d"))(100)
# 修正 annotation_colors
anno_color <- list(
  colors <- c("Proximal tubule cell 1" = "#C5B0D5",
              "Loop of Henle cell" = "#AEC7E8",
              "Distal convoluted tubule cell" = "#17BECF",
              "Endothelial cell 1" = "#FFBB78",
              "Endothelial cell 2" = "#FF7F0E",
              "Fibroblast" = "#D62728",
              "Neutrophil" = "#C49C94",
              "Intercalated cell" = "#2CA02C",
              "Principal cell" = "#8C564B",
              "Podocyte" = "#E376C2",
              "Proximal tubule cell 2" = "#9467BD",
              "Macrophage" = "#BCBD22",
              "Proximal tubule cell 3" = "#7200DA")  
)


# 加载 pheatmap 包
library(pheatmap)

# 绘制热图
pheatmap(
  atacheatmap[,8:9656],
  scale = "none",
  color = custom_colors,
  border_color = NULL,
  cellwidth = 0.03, cellheight = 0.3,
  cluster_rows = F, 
  cluster_cols = F,
  #treeheight_row = 20, treeheight_col = 10,
  #fontsize = 12,
  show_rownames = F, 
  show_colnames = F,
  # annotation_col = group,
  #annotation_colors = anno_color,
  # labels_row = heat$X1,
  # filename = "1110xiu/修图/atac_heatmap.pdf", width = 8, height = 8
)


gc()
DefaultAssay(atac) <- 'ATAC'
atac <- ScaleData(atac,features = rownames(atac))
#pdf("kidney/fig7b_atac_heatmap.pdf",width = 6,height = 6)
DoHeatmap(atac, features = top10$gene, group.colors = colors, label = FALSE) +
  scale_fill_gradientn(
    colors = c("#8d90e3", "#b7b9f7", "#DBDBFF", "#f7f7f7", "#ffb673", "#ffa257","#ff6a13"), # 五种颜色
    values = scales::rescale(c(0, 0.1, 0.3, 0.6, 1, 1.8, 2.5)) # 颜色对应的数值范围
  ) +
  theme(axis.text.y = element_blank())
#dev.off()
#gc()





# fig7 c 左图----

colors <- c("aneuploid"="#ff5743",
            "diploid"="#008bd0",
            "not.defined"="#99b3c6")
p <- ggplot(merged, aes(x = V1, y = V2, 
                        fill = tu, color = tu)) +
  geom_point(size = 1) +
  theme_dr() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1) 
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  guides(fill = "none") 
p
dev.off()
#ggsave("kidney/fig7c_left.pdf",p,width = 6.5,height = 6)

# fig7 c右图
#result <- as.data.frame.matrix(table(merged$new_type, merged$tu))

library(reshape2)
library(ggplot2)
library(ggalluvial)
library(ggh4x)
result <- table(merged$new_type, merged$tu) %>%
  as.data.frame.matrix() %>%
  # 将行名转换为新的一列
  tibble::rownames_to_column(var = "new_type") %>%
  rowwise() %>%
  mutate(across(-new_type, ~ . / sum(c_across(-new_type)))) %>%
  ungroup()
melt.result <- melt(result)
colnames(melt.result) <- c("Part","var","value")
# 绘图
# 用到上述colors
p <- ggplot(melt.result, aes(x = Part, y = value, fill = var,
                             stratum = var, alluvium = var)) +
  geom_col(position = 'stack', width = 0.6) +
  geom_stratum(width = 0.6, color = 'white') +
  geom_alluvium(alpha = 0.4, width = 0.6, color = 'white', linewidth = 1, curve_type = "linear") +
  scale_fill_manual(values = colors) +
  xlab('') + 
  ylab('') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 12) + 
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # 设置 x 轴标签旋转 45 度
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
p
#ggsave("kidney/fig7c_right.pdf", plot = p, height = 4, width = 6)


install.packages("rms")  
library(rms)

ra <- atac
# fig7 g----
Idents(ra) <- "new_type"
object = ra
region = "Hnf4a"
extend.upstream = 3000
extend.downstream = 10000
#object <- SortIdents(object)
#pdf("kidney/fig_7g.pdf",width = 6.5,height = 6)

#install.packages("patchwork")

#.rs.restartR()

#library(Signac)
#library(Seurat)

# 4. 再次运行
CoveragePlot(
  object = object,
  region = region,
  extend.upstream = extend.upstream,
  extend.downstream = extend.downstream,
  annotation = TRUE,
  peaks = TRUE,
  links = TRUE
)
#dev.off()



# fig7 i----
Idents(ra) <- "new_type"
object = ra
cell_types <- levels(Idents(object))
sorted_cell_types <- sort(cell_types)
Idents(object) <- factor(Idents(object), levels = sorted_cell_types)
region = "Pparg"
idents = c("Proximal tubule cell 1", "Proximal tubule cell 2", "Proximal tubule cell 3")
twocol <- c("Proximal tubule cell 1"="#ffa61d", "Proximal tubule cell 2"="#9a99e1", "Proximal tubule cell 3"="#9a99e1")
extend.upstream = 0
extend.downstream = 0


#pdf("kidney/fig7i.pdf",width = 6.5,height = 3.5)
CoveragePlot(
  object = object,
  region = region,
  extend.upstream = extend.upstream,
  extend.downstream = extend.downstream,
  annotation = T,
  idents = idents,
  peaks = T,
  links = T
)
#dev.off()



new <- readRDS("new_mouse_RNA_new.rds")
Idents(new) <- "new_type"
pdf("kidney/fig_7i-2.pdf",width = 6.5,height = 5)
VlnPlot(new,features = region,idents = idents,cols = twocol)
dev.off()


# fig8 b----
Idents(ra) <- "new_type"
object = ra
region = "Pparg"
extend.upstream = 3000

extend.downstream = 10000
#object <- SortIdents(object)
#pdf("kidney/fig_8b.pdf",width = 6.5,height = 6)
CoveragePlot(
  object = object,
  region = region,
  #group.by = "new_type",
  extend.upstream = extend.upstream,
  extend.downstream = extend.downstream,
  annotation = T,
  peaks = T,
  links = T
)
#dev.off()





Idents(ra) <- "new_type"
object = ra
cell_types <- levels(Idents(object))
sorted_cell_types <- sort(cell_types)
Idents(object) <- factor(Idents(object), levels = sorted_cell_types)

region = "chr16-11300000-11800000"
idents = c("Proximal tubule cell 1", "Proximal tubule cell 2", "Proximal tubule cell 3")
twocol <- c("Proximal tubule cell 1"="#ffb69c", "Proximal tubule cell 2"="#e2dda2", "Proximal tubule cell 3"="#6f89ad")
extend.upstream = 0
extend.downstream = 0
  
#pdf(paste0("kidney/",region,"_1.pdf"),width = 6.5,height = 3.8)
CoveragePlot(
  object = object,
  region = region,
  extend.upstream = extend.upstream,
  extend.downstream = extend.downstream,
  annotation = T,
  idents = idents,
  peaks = T,
  links = T
)
#dev.off()


idents = c("Proximal tubule cell 1", "Proximal tubule cell 2", "Proximal tubule cell 3")
twocol <- c("Proximal tubule cell 1"="#ffb69c", "Proximal tubule cell 2"="#e2dda2", "Proximal tubule cell 3"="#6f89ad")
extend.upstream = 0
extend.downstream = 0




# fig8 e f
library(CellChat)
cellchat <- readRDS("prox_cellchat.rds")
custom_colors <- colorRampPalette(c("#019f97", "white", "#9c93e5"))(100)
levels(cellchat@idents)
colors <- c("Distal convoluted tubule cell" = "#17BECF",
            "Fibroblast" = "#D62728",
            "Intercalated cell" = "#2CA02C",
            "Loop of Henle cell" = "#AEC7E8",
            "Podocyte" = "#E376C2",
            "Proximal tubule cell 1" = "#C5B0D5",
            "Proximal tubule cell 2" = "#9467BD",
            "Proximal tubule cell 3" = "#7200DA")  
#pdf("kidney/fig8_cellchat.heatmep.pdf",width = 12,height = 6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = colors)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use = colors)
ht1 + ht2
#dev.off()

#pdf("kidney/fig8cellchat.heatmep1.pdf",width = 12,height = 6)
gg1 <- netVisual_heatmap(cellchat,, color.use = colors)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", color.use = colors)
gg1 + gg2
#dev.off()
 
#pdf("kidney/fig8_F.pdf",width = 8,height = 6)
pathways.show <- "EGF"
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy",color.use = colors)
#dev.off()

#pdf("kidney/fig8e_cellchat.pdf",width = 8,height = 6)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use = colors)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = colors)
#dev.off()


