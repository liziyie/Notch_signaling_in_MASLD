source("./MASLD_utils.R")
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(scCustomize)

# Load Data----
mouse_liver_sc_atlas_macs <- readRDS("./mouse_liver_sc_atlas_macs_release.rds")
RBPJ_new <- mouse_liver_sc_atlas_macs %>% subset(orig.ident == "RBPJ")
RBPJ_new@meta.data$Group <- factor(RBPJ_new@meta.data$Group, levels = c("ND_WT","ND_KO","MCD_WT","MCD_KO"))
RBPJ_new <- SetIdent(RBPJ_new, value = RBPJ_new@meta.data$Cluster)

cluster_color_panel <- c(
  "Ly6Chi Monocytes" = "#C45240", "Patrolling Monocytes" = "#F2C987",
  "Clec4e+ Macs" = "#63B571", "C1qc+ Macs" = "#F8C3B4",
  "Capsule Macs" = "#86C2D7",
  "Spp1+ Macs" = "#C86EA8",
  "moKC" = "#F0CE39", "KCs" = "#60B4F3",
  "Prolif. KCs" = "#F0C1ED"
)

# Figure 4D----
RBPJ_Macro <- RBPJ_new %>% subset(Cluster %in% c("Ly6Chi Monocytes","Patrolling Monocytes","Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC"))

library(monocle3)
library(SeuratWrappers)
library(magrittr)

cells_used_WT <- RBPJ_Macro@meta.data %>% filter(Group == "MCD_WT") %>% row.names()
RBPJ_Macro_WT <- RBPJ_Macro %>% subset(subset = CellName %in% cells_used_WT)
RBPJ_Macro_WT <- RBPJ_Macro_WT %>% NormalizeData() %>% FindVariableFeatures()
hvg <- RBPJ_Macro_WT %>% FindAllMarkers() %>% filter(p_val_adj < 0.05, pct.1 > .3) %>% pull(gene) %>% unique()
RBPJ_Macro_WT@assays$RNA@var.features <- hvg
RBPJ_Macro_WT <- RBPJ_Macro_WT %>% ScaleData() %>% RunPCA() 
RBPJ_Macro_WT <- RBPJ_Macro_WT %>% RunUMAP(dims = 1:15)
DimPlot(RBPJ_Macro_WT,label = T)
RBPJ_Macro_WT.cds <- as.cell_data_set(RBPJ_Macro_WT)
RBPJ_Macro_WT.cds <- cluster_cells(RBPJ_Macro_WT.cds, cluster_method = "louvain")
colData(RBPJ_Macro_WT.cds)$Cluster <- colData(RBPJ_Macro_WT.cds)$Cluster
RBPJ_Macro_WT.cds@clusters$UMAP$clusters <- as.character(RBPJ_Macro_WT.cds@colData$Cluster)
names(RBPJ_Macro_WT.cds@clusters$UMAP$clusters) <- RBPJ_Macro_WT.cds@colData$CellName
RBPJ_Macro_WT.cds <- learn_graph(RBPJ_Macro_WT.cds)
RBPJ_Macro_WT.cds <- order_cells(RBPJ_Macro_WT.cds)

cells_used_KO <- RBPJ_Macro@meta.data %>% filter(Group == "MCD_KO") %>% row.names()
RBPJ_Macro_KO <- RBPJ_Macro %>% subset(subset = CellName %in% cells_used_KO)
RBPJ_Macro_KO <- RBPJ_Macro_KO %>% NormalizeData() %>% FindVariableFeatures()
hvg <- RBPJ_Macro_KO %>% FindAllMarkers() %>% filter(p_val_adj < 0.05, pct.1 > .3) %>% pull(gene) %>% unique()
RBPJ_Macro_KO@assays$RNA@var.features <- hvg
RBPJ_Macro_KO <- RBPJ_Macro_KO %>% ScaleData() %>% RunPCA() 
RBPJ_Macro_KO <- RBPJ_Macro_KO %>% RunUMAP(dims = 1:15, min.dist = .01, spread = 1.2)
DimPlot(RBPJ_Macro_KO,label = T)
RBPJ_Macro_KO.cds <- as.cell_data_set(RBPJ_Macro_KO)
RBPJ_Macro_KO.cds <- cluster_cells(RBPJ_Macro_KO.cds, cluster_method = "louvain")
colData(RBPJ_Macro_KO.cds)$Cluster <- colData(RBPJ_Macro_KO.cds)$Cluster
RBPJ_Macro_KO.cds@clusters$UMAP$clusters <- as.character(RBPJ_Macro_KO.cds@colData$Cluster)
names(RBPJ_Macro_KO.cds@clusters$UMAP$clusters) <- RBPJ_Macro_KO.cds@colData$CellName
RBPJ_Macro_KO.cds <- learn_graph(RBPJ_Macro_KO.cds)
RBPJ_Macro_KO.cds <- order_cells(RBPJ_Macro_KO.cds)

RBPJ_Macro@meta.data[RBPJ_Macro_WT.cds@colData$CellName,"Monocle_pseudotime"] <- RBPJ_Macro_WT.cds@principal_graph_aux@listData$UMAP$pseudotime
RBPJ_Macro@meta.data[RBPJ_Macro_KO.cds@colData$CellName,"Monocle_pseudotime"] <- RBPJ_Macro_KO.cds@principal_graph_aux@listData$UMAP$pseudotime
RBPJ_Macro@meta.data[RBPJ_Macro_WT.cds@colData$CellName,"Monocle_pseudotime_scaled"] <- RBPJ_Macro_WT.cds@principal_graph_aux@listData$UMAP$pseudotime/max(RBPJ_Macro_WT.cds@principal_graph_aux@listData$UMAP$pseudotime[RBPJ_Macro_WT.cds@principal_graph_aux@listData$UMAP$pseudotime != "Inf"])
RBPJ_Macro@meta.data[RBPJ_Macro_KO.cds@colData$CellName,"Monocle_pseudotime_scaled"] <- RBPJ_Macro_KO.cds@principal_graph_aux@listData$UMAP$pseudotime/max(RBPJ_Macro_KO.cds@principal_graph_aux@listData$UMAP$pseudotime[RBPJ_Macro_KO.cds@principal_graph_aux@listData$UMAP$pseudotime != "Inf"])

p1 <- plot_cells(
  RBPJ_Macro_WT.cds, 
  color_cells_by = "Cluster",
  cell_size = 1,
  alpha = .8,
  show_trajectory_graph = T,
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  scale_color_manual(values = cluster_color_panel) +
  theme(legend.position = "none")

p2 <- plot_cells(
  RBPJ_Macro_KO.cds, 
  color_cells_by = "Cluster",
  cell_size = 1,
  alpha = .8,
  show_trajectory_graph = T,
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE) +
  scale_color_manual(values = cluster_color_panel) +
  theme(legend.position = "none")

# Figure 4I and Figure S4H----
RBPJ_subset <- RBPJ_new %>% subset(Cluster %in% c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC"))
RBPJ_markers <- FindAllMarkers(RBPJ_subset)
Clec4e_markers <- RBPJ_markers %>% filter(cluster == "Clec4e+ Macs") %>% arrange(desc(avg_log2FC))
gene.df <- bitr(Clec4e_markers$gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
degenes_entrez <- merge(Clec4e_markers, gene.df, by.x = "gene", by.y = "SYMBOL")
degenes_entrez <- degenes_entrez %>% arrange(desc(avg_log2FC))
geneList <- degenes_entrez$avg_log2FC
names(geneList) <- degenes_entrez$ENTREZID
ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Mm.eg.db,
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE)
# The 18 pathway is "Positive regulation of apoptotic process"
p <- gseaplot(ego, geneSetID = 7, title = ego$Description[18]) +
  geom_text(x = 4000, y = 2.5, 
            label = paste0("NES=",
                           round(ego@result[18,"NES"],2),"\n",
                           "P=",
                           formatC(ego@result[18,"p.adjust"], format = "e", digits = 2)),
            hjust = 0)

enriched_entrezid <- strsplit(ego@result$core_enrichment[18],"/")[[1]]
enriched_symbol <- bitr(enriched_entrezid, fromType = "ENTREZID", toType = c("SYMBOL"), OrgDb = org.Mm.eg.db)
exp_matrix <- aggregate(t(RBPJ_subset@assays$RNA@data[enriched_symbol$SYMBOL,]),list(Cluster = RBPJ_subset@meta.data[,"Cluster"]),mean)
row.names(exp_matrix) <- exp_matrix$Cluster
exp_matrix$Cluster <- c()
plot.matrix <- apply(exp_matrix,2,zscore)
plot.matrix <- t(plot.matrix)
plot.matrix_quantile <- quantile(plot.matrix, c(0.01, 0.99))
plot.matrix <- pmax(plot.matrix, plot.matrix_quantile[1])
plot.matrix <- pmin(plot.matrix, plot.matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(plot.matrix), max(plot.matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
genes_labeled <- c("Nr4a3","Fas","Trem1","S100a8","Cd40","Bcl2a1a","Tnf","Nupr1","Thbs1","Cd24a","Hif1a","Nfkbid","S100a9","Atf3","Ccl3","Cd44","Nlrc4","Dusp1","Cd274","Sik1","Notch2","Mcl1","Traf2","Tspo","Il1b")
p <- Heatmap(plot.matrix,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp") +
  rowAnnotation(link = anno_mark(
    at = which(row.names(plot.matrix) %in% genes_labeled),
    labels = row.names(plot.matrix)[which(row.names(plot.matrix) %in% genes_labeled)],
    labels_gp = gpar(fontsize = 7), padding = unit(1, "mm"))
  )

# Figure S4I----
library(ggbeeswarm)
library(ggsignif)
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_WT"), Cluster %ni% c("Prolif. KCs","Capsule Macs")) %>% pull(CellName)
plot.data <- data.frame(
  Fas = RBPJ_new@assays$RNA@data["Fas",cells_used] + 
    abs(rnorm(length(cells_used))/100000),
  Cluster = RBPJ_new@meta.data[cells_used,"Cluster"],
  Group = RBPJ_new@meta.data[cells_used,"Group"])
p <- 
  ggplot(plot.data, aes(x = Cluster, y = Fas)) +
  geom_quasirandom(aes(col = Cluster), cex = 1.5, width = .25, alpha = 0.5) +
  labs(y = "Fas expression") +
  scale_color_manual(values = cluster_color_panel) +
  geom_violin(aes(color = Cluster, fill = Cluster), alpha = 0.3, scale = "width") +
  geom_signif(comparisons = list(c("C1qc+ Macs","Clec4e+ Macs"),c("Clec4e+ Macs","Spp1+ Macs"),c("Clec4e+ Macs","moKC")), step_increase = .05, tip_length = 0, textsize = 3) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank())