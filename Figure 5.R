source("./MASLD_utils.R")
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(scCustomize)

# Load Data----
mouse_liver_sc_atlas_macs <- readRDS("./mouse_liver_sc_atlas_macs_release.rds")
RBPJ_new <- mouse_liver_sc_atlas_macs %>% subset(orig.ident == "RBPJ")
RBPJ_new@meta.data$Group <- factor(RBPJ_new@meta.data$Group, levels = c("ND_WT","ND_KO","MCD_WT","MCD_KO"))
RBPJ_new <- SetIdent(RBPJ_new, value = RBPJ_new@meta.data$Cluster)

liver_sc_atlas_macs <- readRDS("./liver_sc_atlas_macs_release.rds")
merged_metadata <- merge(
  liver_sc_atlas_macs@meta.data %>% mutate(CellID = paste0(orig.ident,"_",Raw_CellName)),
  mouse_liver_sc_atlas_macs@meta.data %>% mutate(CellID = paste0(orig.ident,"_",Raw_CellName)),
  by = "CellID", all.x = T
)
row.names(merged_metadata) <- merged_metadata$CellName.x
liver_sc_atlas_macs@meta.data$Cluster_mouse <- merged_metadata[row.names(liver_sc_atlas_macs@meta.data),"Cluster.y"]

cluster_color_panel <- c(
  "Ly6Chi Monocytes" = "#C45240", "Patrolling Monocytes" = "#F2C987",
  "Clec4e+ Macs" = "#63B571", "C1qc+ Macs" = "#F8C3B4",
  "Capsule Macs" = "#86C2D7",
  "Spp1+ Macs" = "#C86EA8",
  "moKC" = "#F0CE39", "KCs" = "#60B4F3",
  "Prolif. KCs" = "#F0C1ED"
)

# Figure 5B----
library(clusterProfiler)
library(org.Mm.eg.db)
RBPJ_subset <- RBPJ_new %>% subset(Cluster %in% c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC"))
RBPJ_markers <- FindAllMarkers(RBPJ_subset)
universe_symbol <- row.names(RBPJ_subset@assays$RNA@counts)[rowSums(RBPJ_subset@assays$RNA@counts) > 10]
universe_entrezid <- bitr(universe_symbol, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
universe_entrezid <- unique(universe_entrezid$ENTREZID)
gene_list <- list()
for(cluster_used in as.character(unique(RBPJ_markers$cluster))){
  genes_symbol <- RBPJ_markers %>% filter(cluster == cluster_used, p_val_adj < 0.05, avg_log2FC > 0) %>% top_n(n = 200, wt = avg_log2FC) %>% pull(gene)
  gene.df <- bitr(genes_symbol, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
  gene_list[[cluster_used]] <- unique(gene.df$ENTREZID)
}
cluster_go <- compareCluster(gene_list, fun = "enrichGO",
                             OrgDb = org.Mm.eg.db, ont = "ALL", 
                             universe = universe_entrezid,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = TRUE, pool = TRUE)
description_used <- c(
  "immune effector process",
  "acute inflammatory response",
  "positive regulation of programmed cell death",
  "positive regulation of inflammatory response",
  "positive regulation of apoptotic process",
  "integrin binding",
  "response to lipoprotein particle",
  "fatty acid binding",
  "lipoprotein particle binding",
  "antigen processing and presentation",
  "transmembrane signaling receptor activity",
  "MHC protein complex",
  "extracellular matrix",
  "complement activation",
  "wound healing"
)
description_used <- unique(description_used)
cluster_go_plot <- cluster_go@compareClusterResult %>% 
  filter(Description %in% description_used)
dim_x <- unique(as.character(cluster_go_plot$Cluster))
dim_y <- description_used
map_x <- setNames(seq_along(dim_x), dim_x)
map_y <- setNames(seq_along(dim_y), dim_y)
cluster_go_plot_mat <- sparseMatrix(
  i=map_x[as.character(cluster_go_plot$Cluster)], 
  j=map_y[as.character(cluster_go_plot$Description)], 
  x=-log10(cluster_go_plot$p.adjust), 
  dims=c(length(dim_x), length(dim_y)), 
  dimnames=list(dim_x, dim_y)
)
cluster_go_plot_mat <- as.matrix(cluster_go_plot_mat)
cluster_go_plot_mat <- cluster_go_plot_mat[c("Clec4e+ Macs","Spp1+ Macs","C1qc+ Macs","moKC"),]
cluster_go_plot_mat <- t(cluster_go_plot_mat)
max.exp <- apply(cluster_go_plot_mat,1,which.max)
pathway_order <- c()
for(i in 1:4){
  pathway_order <- c(pathway_order, names(sort(cluster_go_plot_mat[names(max.exp)[max.exp == i],i],decreasing = T)))
}
cluster_go_plot_mat <- cluster_go_plot_mat[pathway_order,]
cluster_go_plot_mat_quantile <- quantile(cluster_go_plot_mat, c(0.05, 0.95))
cluster_go_plot_mat <- pmax(cluster_go_plot_mat, cluster_go_plot_mat_quantile[1])
cluster_go_plot_mat <- pmin(cluster_go_plot_mat, cluster_go_plot_mat_quantile[2])
color_used <- circlize::colorRamp2(seq(min(cluster_go_plot_mat), max(cluster_go_plot_mat), length = 9), rev(RColorBrewer::brewer.pal(11,"RdYlBu")[2:10]))
p <- Heatmap(cluster_go_plot_mat,
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "-log10(adjusted P-value)")

# Figure 5C----
proinflammation_genes <- c("Vcan","Mnda","Lyz2","S100a8","S100a6","S100a4","Il1r2","Tspo","S100a9","Capg","Lgals1","Coro1a","Cd37","S100a10","Il1b","Tnf","Ifng","Ccl2","Ccl3","Ccl4","Ccl5","Csf3r")

cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_KO","MCD_WT"), Cluster %ni% c("Prolif. KCs","Capsule Macs","Ly6Chi Monocytes","Patrolling Monocytes")) %>% pull(CellName)
RBPJ_subset <- RBPJ_new[,cells_used]
RBPJ_subset <- AddModuleScore(RBPJ_subset, 
                              list(proinflammation = proinflammation_genes), 
                              ctrl = length(proinflammation_genes), 
                              name = "proinflammation")
plot.data <- data.frame(
  proinflammation_score = RBPJ_subset@meta.data$proinflammation1,
  Cluster = RBPJ_subset@meta.data$Cluster)
p <- 
  ggplot(plot.data, aes(x = Cluster, y = proinflammation_score)) +
  geom_quasirandom(aes(col = Cluster), cex = 1, width = 0.25, alpha = 0.5) +
  geom_boxplot(aes(color = Cluster), outlier.size = -1, alpha = 0) +
  geom_signif(comparisons = list(c("C1qc+ Macs","Clec4e+ Macs"),c("Clec4e+ Macs","Spp1+ Macs"),c("Clec4e+ Macs","moKC")), step_increase = .05, tip_length = 0, textsize = 3) +
  labs(y = "Pro-inflammation Score") +
  scale_color_manual(values = cluster_color_panel) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text = element_blank(),
        axis.title.x = element_blank()) 

# Figure 5F----
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_KO","MCD_WT"), Cluster %ni% c("Ly6Chi Monocytes","Patrolling Monocytes","Prolif. KCs","Capsule Macs")) %>% pull(CellName)
plot.data <- data.frame(
  Ccl2 = RBPJ_new@assays$RNA@data["Ccl2",cells_used] + abs(rnorm(length(cells_used))/100000),
  Cluster = RBPJ_new@meta.data[cells_used,"Cluster"],
  Group = RBPJ_new@meta.data[cells_used,"Group"])
p <- 
  ggplot(plot.data, aes(x = Group, y = Ccl2)) +
  geom_quasirandom(aes(col = Group), cex = 1, width = .25, alpha = 0.5) +
  geom_boxplot(aes(color = Group), outlier.size = -1, alpha = 0) +
  geom_signif(comparisons = list(c("MCD_WT","MCD_KO")),
              tip_length = 0, textsize = 3) +
  labs(y = "Ccl2 expression") +
  ylim(c(0,6)) +
  facet_wrap(~Cluster, nrow = 1) +
  scale_color_manual(values = c("MCD_WT" = "#e08380", "MCD_KO" = "#82be9c")) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank())  

# Figure 5K----
p1 <- ggplot(liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "human_liver"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = .5, alpha = .8) +
  theme_cowplot(font_size = 7) +
  xlim(c(-8.6,12.7)) + ylim(c(-7.5,6.7)) +
  scale_color_manual(name = "", values = cluster_color_panel) + 
  addLabels(centers = calcCenters(liver_sc_atlas_macs$Cluster,liver_sc_atlas_macs@meta.data[,c("UMAP_1","UMAP_2")])) + 
  ggtitle("Human") +
  theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")
p2 <- ggplot(liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "RBPJ"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster_mouse), size = .5, alpha = .8) +
  theme_cowplot(font_size = 7) +
  xlim(c(-8.6,12.7)) + ylim(c(-7.5,6.7)) +
  scale_color_manual(name = "", values = cluster_color_panel)  + 
  addLabels(centers = calcCenters(liver_sc_atlas_macs$Cluster_mouse,liver_sc_atlas_macs@meta.data[,c("UMAP_1","UMAP_2")])) + 
  ggtitle("RBPJ") +
  theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none")

# Figure 5L----
library(scibetR)
human_liver_cells <- liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "human_liver", Cluster %in% c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC")) %>% pull(CellName)
RBPJ_cells <- liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "RBPJ", Cluster %in% c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC")) %>% pull(CellName)
train_set <- as.matrix(t(liver_sc_atlas_macs@assays$RNA@data[,RBPJ_cells]))
train_set <- data.frame(train_set, label = liver_sc_atlas_macs@meta.data[RBPJ_cells,"Cluster"] %>% droplevels())
test_set <- as.matrix(t(liver_sc_atlas_macs@assays$RNA@data[,human_liver_cells]))
test_set <- data.frame(test_set, label = liver_sc_atlas_macs@meta.data[human_liver_cells,"Cluster"] %>% droplevels())
prd <- SciBet_R(train_set, test_set[,-ncol(test_set)])
confusion_matrix <- table(as.character(test_set$label),prd)
confusion_matrix <- t(t(confusion_matrix) / colSums(confusion_matrix))
confusion_matrix <- as.data.frame(confusion_matrix)
confusion_matrix <- dcast(confusion_matrix, Var1~prd, value.var = "Freq")
row.names(confusion_matrix) <- confusion_matrix$Var1
confusion_matrix$Var1 <- c()
confusion_matrix <- as.matrix(confusion_matrix)
test_order <- c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC") 
train_order <- c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC")
color_used <- 
  circlize::colorRamp2(c(seq(0,0.2,length.out = 50),0.7), 
                       c(colorRampPalette(rev(brewer.pal(6,'Blues')))(100)[51:100],"red2"), space = "RGB")
p <- Heatmap(
  confusion_matrix[test_order, train_order],
  name = 'Cell Proportion',
  column_title = "RBPJ (train data)", row_title = "Human Liver (test data)",
  col = color_used,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_side = "left", row_dend_side = "right",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7), 
                              labels_gp = gpar(fontsize = 7),
                              legend_height = unit(3, "cm"))
)

# Figure 5M----
proinflammation_genes <- c("Vcan","Mnda","Lyz2","S100a8","S100a6","S100a4","Il1r2","Tspo","S100a9","Capg","Lgals1","Coro1a","Cd37","S100a10","Il1b","Tnf","Ifng","Ccl2","Ccl3","Ccl4","Ccl5","Csf3r")
cells_used <- liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "human_liver", Cluster %ni% c("Prolif. KCs","Ly6Chi Monocytes","Patrolling Monocytes")) %>% pull(CellName)
human_liver_subset <- liver_sc_atlas_macs[,cells_used]
human_liver_subset <- AddModuleScore(human_liver_subset, 
                                     list(inflammation = proinflammation_genes), 
                                     ctrl = length(inflammation_genes), 
                                     name = "Inflammation")
plot.data <- data.frame(
  Inflammation_score = human_liver_subset@meta.data$Inflammation1,
  Cluster = human_liver_subset@meta.data$Cluster)
p <- 
  ggplot(plot.data, aes(x = Cluster, y = Inflammation_score)) +
  ggbeeswarm::geom_quasirandom(aes(col = Cluster), cex = 1, width = 0.25, alpha = 0.5) +
  geom_boxplot(aes(color = Cluster), outlier.size = -1, alpha = 0) +
  ggsignif::geom_signif(comparisons = list(c("C1qc+ Macs","Clec4e+ Macs"),c("Clec4e+ Macs","Spp1+ Macs"),c("Clec4e+ Macs","moKC")), step_increase = .05, tip_length = 0, textsize = 3) +
  labs(y = "Pro-inflammation Score") +
  scale_color_manual(values = cluster_color_panel) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text = element_blank(),
        axis.title.x = element_blank()) 

# Figure S5A----
inflammation_genes <- c("Il1b","Tnf","Ifng","Ccl2","Ccl3","Ccl4","Ccl5","Csf3r")
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_KO","MCD_WT"), Cluster %ni% c("Prolif. KCs","Capsule Macs","Ly6Chi Monocytes","Patrolling Monocytes")) %>% pull(CellName)
plot.matrix <- aggregate(t(RBPJ_new@assays$RNA@data[inflammation_genes,cells_used]),list(RBPJ_new@meta.data[cells_used,"Cluster"]),mean)
row.names(plot.matrix) <- plot.matrix$Group.1
plot.matrix$Group.1 <- c()
plot.matrix <- apply(plot.matrix, 2, zscore)
plot.matrix_quantile <- quantile(plot.matrix, c(0.05, 0.95))
plot.matrix <- pmax(plot.matrix, plot.matrix_quantile[1])
plot.matrix <- pmin(plot.matrix, plot.matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(plot.matrix), max(plot.matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
p <- Heatmap(plot.matrix,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 7),
             row_names_gp = gpar(fontsize = 7),
             col = color_used,
             name = "Exp")

# Figure S5B----
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_KO","MCD_WT"), Cluster %ni% c("Prolif. KCs","Capsule Macs","Ly6Chi Monocytes","Patrolling Monocytes")) %>% pull(CellName)
RBPJ_subset <- RBPJ_new[,cells_used]
RBPJ_subset <- AddModuleScore(RBPJ_subset, 
                              list(inflammation = inflammation_genes), 
                              ctrl = length(inflammation_genes), 
                              name = "Inflammation")
plot.data <- data.frame(
  Inflammation_score = RBPJ_subset@meta.data$Inflammation1,
  Cluster = RBPJ_subset@meta.data$Cluster,
  Group = RBPJ_subset@meta.data$Group)
p <- 
  ggplot(plot.data %>% filter(Cluster == "Clec4e+ Macs"), aes(x = Group, y = Inflammation_score)) +
  geom_quasirandom(aes(col = Group), cex = 1, width = 0.25, alpha = 0.5) +
  geom_boxplot(aes(color = Group), outlier.size = -1, alpha = 0) +
  geom_signif(comparisons = list(c("MCD_WT","MCD_KO")),
              tip_length = 0, textsize = 3) +
  labs(y = "Pro-inflammation cytokines in Clec4e+ Macs") +
  scale_color_manual(values = c("MCD_WT" = "#e08380", "MCD_KO" = "#82be9c")) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank())

# Figure S5C----
lipid_genes <- c("Trem2", "Cd9", "Spp1", "Gpnmb", "Fabp5", "Cd36")
lipid_genes <- lipid_genes[lipid_genes %in% row.names(RBPJ_new)]
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_KO","MCD_WT"), Cluster %ni% c("Prolif. KCs","Capsule Macs","Ly6Chi Monocytes","Patrolling Monocytes")) %>% pull(CellName)
RBPJ_subset <- RBPJ_new[,cells_used]
RBPJ_subset <- AddModuleScore(RBPJ_subset, 
                              list(lipid = lipid_genes), 
                              ctrl = length(lipid_genes), 
                              name = "lipid")
plot.data <- data.frame(
  lipid_score = RBPJ_subset@meta.data$lipid1,
  Cluster = RBPJ_subset@meta.data$Cluster)
p <- 
  ggplot(plot.data, aes(x = Cluster, y = lipid_score)) +
  geom_quasirandom(aes(col = Cluster), cex = 1, width = 0.25, alpha = 0.5) +
  geom_boxplot(aes(color = Cluster), outlier.size = -1, alpha = 0) +
  geom_signif(comparisons = list(c("C1qc+ Macs","Spp1+ Macs"),c("Clec4e+ Macs","Spp1+ Macs"),c("Spp1+ Macs","moKC")), step_increase = .05, tip_length = 0, textsize = 3) +
  labs(y = "Pro-lipid Score") +
  scale_color_manual(values = cluster_color_panel) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text = element_blank(),
        axis.title.x = element_blank())

# Figure S5D----
lipid_genes <- c("Trem2", "Cd9", "Spp1", "Gpnmb", "Fabp5", "Cd36", "Vegfa", "Pdgfb")
inflammation_genes <- c("Vcan","Mnda","Lyz2","S100a8","S100a6","S100a4","Il1r2","Tspo","S100a9","Capg","Lgals1","Coro1a","Cd37","S100a10")
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_KO","MCD_WT"), Cluster %ni% c("Prolif. KCs","Capsule Macs","Ly6Chi Monocytes","Patrolling Monocytes")) %>% pull(CellName)
p <- 
  RBPJ_new@assays$RNA@data[c(inflammation_genes,lipid_genes),cells_used] %>% 
  as.data.frame() %>% 
  mutate(Gene = c(inflammation_genes,lipid_genes)) %>% 
  melt(id.vars = "Gene") %>% 
  mutate(Gene = factor(Gene, levels = c(inflammation_genes,lipid_genes))) %>% 
  merge(., RBPJ_new@meta.data[cells_used,c("CellName","Cluster")], by.x = "variable", by.y = "CellName") %>% 
  mutate(Cluster = factor(Cluster, levels = rev(levels(Cluster)))) %>%
  ggplot(aes(x = Cluster, y = value, fill = Cluster)) +
  geom_violin(scale = "width", color = NA) +
  coord_flip() +
  facet_grid(~Gene, switch = "x", scales = "free_y") +
  scale_fill_manual(values = cluster_color_panel) +
  theme_cowplot(font_size = 7) +
  scale_y_continuous(position = "right") +
  theme(legend.position = "null",
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size = .25),
        strip.background = element_rect(fill=NA, color=NA),
        strip.text.x = element_text(size = 7, angle = 270, hjust = 0))

# Figure S5K----
cells_used <- RBPJ_new@meta.data %>% filter(Cluster == "KCs", Group %in% c("ND_KO","ND_WT")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Group"]
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "KO"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "WT"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol),deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol),"Cd36")
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val + 1e-303))) +
  geom_point(aes(color = Sig), size = 2) +
  scale_color_manual(values = c("WT" = "#e08380", "KO" = "#82be9c", "FALSE" = "lightgrey")) +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc'))) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none")

# Figure S5L----
cells_used <- RBPJ_new@meta.data %>% filter(Cluster %in% c("moKC","KCs","Prolif. KCs")) %>% pull(CellName)
plot.data <- data.frame(
  Cd36 = RBPJ_new@assays$RNA@data["Cd36",cells_used] + abs(rnorm(length(cells_used))/100000),
  Cluster = RBPJ_new@meta.data[cells_used,"Cluster"],
  Group = RBPJ_new@meta.data[cells_used,"Group"])
p <- 
  ggplot(plot.data, aes(x = Group, y = Cd36)) +
  geom_quasirandom(aes(col = Group), cex = 1, width = .25, alpha = 0.5) +
  geom_boxplot(aes(color = Group), outlier.size = -1, alpha = 0) +
  geom_signif(comparisons = list(c("MCD_WT","MCD_KO"),c("ND_WT","ND_KO")),
              tip_length = 0, textsize = 3) +
  labs(y = "Ccl2 expression") +
  ylim(c(0,4.5)) +
  facet_wrap(~Cluster, nrow = 1) +
  scale_color_manual(values = c("ND_WT" = "#e08380", "ND_KO" = "#82be9c", "MCD_WT" = "#e08380", "MCD_KO" = "#82be9c")) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank())

# Figure S5O----
alluvial_cluster_color_panel <- cluster_color_panel
names(alluvial_cluster_color_panel) <- paste0("Atlas_", names(alluvial_cluster_color_panel))
library(ggalluvial)
p <- liver_sc_atlas_macs@meta.data %>% 
  filter(orig.ident %in% c("human_liver")) %>% 
  select(Raw_Cluster,  Cluster) %>% table() %>% as.data.frame() %>% filter(Freq > 0) %>%
  mutate(Raw_Cluster = factor(paste0("Human_",Raw_Cluster), 
                              levels = paste0("Human_",
                                              intersect(names(raw_cluster_color_panel),
                                                        as.character(unique(Raw_Cluster)))))) %>%
  mutate(Cluster = factor(paste0("Atlas_",Cluster), 
                          levels = paste0("Atlas_",
                                          intersect(names(cluster_color_panel),
                                                    as.character(unique(Cluster)))))) %>%
  ggplot(aes(axis1 = Cluster, axis2 = Raw_Cluster, y = Freq)) +
  geom_alluvium(aes(fill = Cluster), curve_type = "sigmoid", width = 1/8) +
  geom_stratum(width = 1/8) +
  ggtitle("Human Liver - Human Liver Raw") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = alluvial_cluster_color_panel) +
  theme_void(base_size = 7) +
  theme(legend.position = "none")

# Figure S5P----
library(scMetabolism)
library(VISION)
library(rsvd)
gmtFile <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
scMetabolism_mouse <- list()
for(dataset in unique(mouse_liver_sc_atlas_macs@meta.data$orig.ident)){
  cells_used <- mouse_liver_sc_atlas_macs@meta.data %>% filter(orig.ident == dataset) %>% pull(CellName)
  countexp <- mouse_liver_sc_atlas_macs@assays$RNA@counts[,cells_used]
  countexp <- countexp[rowSums(countexp) > 50,]
  gene_symbol <- convertGeneList(row.names(countexp), species = "mouse")
  gene_symbol <- gene_symbol[!duplicated(gene_symbol$MGI.symbol),]
  row.names(gene_symbol) <- gene_symbol$MGI.symbol
  countexp <- countexp[gene_symbol$MGI.symbol,]
  row.names(countexp) <- gene_symbol$HGNC.symbol
  result.completed <- alra(as.matrix(countexp))
  countexp2 <- result.completed[[3]]
  row.names(countexp2) <- row.names(countexp)
  n.umi <- colSums(countexp2)
  scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
  vis <- Vision(scaled_counts, signatures = gmtFile)
  options(mc.cores = 12)
  vis <- analyze(vis)
  scMetabolism_mouse[[dataset]] <- data.frame(t(vis@SigScores))
}

cells_used <- liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "human_liver") %>% pull(CellName)
countexp <- liver_sc_atlas_macs@assays$RNA@counts[,cells_used]
countexp <- countexp[rowSums(countexp) > 50,]
gene_symbol <- convertGeneList(row.names(countexp), species = "mouse")
gene_symbol <- gene_symbol[!duplicated(gene_symbol$MGI.symbol),]
row.names(gene_symbol) <- gene_symbol$MGI.symbol
countexp <- countexp[gene_symbol$MGI.symbol,]
row.names(countexp) <- gene_symbol$HGNC.symbol
result.completed <- alra(as.matrix(countexp))
countexp2 <- result.completed[[3]]
row.names(countexp2) <- row.names(countexp)
n.umi <- colSums(countexp2)
scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
vis <- Vision(scaled_counts, signatures = gmtFile)
options(mc.cores = 12)
vis <- analyze(vis)
scMetabolism_human <- data.frame(t(vis@SigScores))

scMetabolism_mean.list <- list(
  mouse_NAFLD = aggregate(
    t(scMetabolism_mouse$mouse_liver_NAFLD),
    list(Cluster = mouse_liver_sc_atlas_macs@meta.data[gsub("[.]","-",colnames(scMetabolism_mouse$mouse_liver_NAFLD)),"Cluster"]),
    mean),
  mouse_WD = aggregate(
    t(scMetabolism_mouse$mouse_liver_WD),
    list(Cluster = mouse_liver_sc_atlas_macs@meta.data[gsub("[.]","-",colnames(scMetabolism_mouse$mouse_liver_WD)),"Cluster"]),
    mean),
  mouse_RBPJ = aggregate(
    t(scMetabolism_mouse$RBPJ),
    list(Cluster = mouse_liver_sc_atlas_macs@meta.data[gsub("[.]","-",colnames(scMetabolism_mouse$RBPJ)),"Cluster"]),
    mean),
  human_liver = aggregate(
    t(scMetabolism_human),
    list(Cluster = liver_sc_atlas_macs@meta.data[gsub("[.]","-",colnames(scMetabolism_human)),"Cluster"]),
    mean)
)
for(i in names(scMetabolism_mean.list)){
  temp_rowname <- scMetabolism_mean.list[[i]]$Cluster
  row.names(scMetabolism_mean.list[[i]]) <- temp_rowname
  scMetabolism_mean.list[[i]]$Cluster <- c()
  scMetabolism_mean.list[[i]] <- t(scMetabolism_mean.list[[i]])
  scMetabolism_mean.list[[i]] <- apply(scMetabolism_mean.list[[i]],1,scale)
  row.names(scMetabolism_mean.list[[i]]) <- temp_rowname
  scMetabolism_mean.list[[i]] <- t(scMetabolism_mean.list[[i]])
  colnames(scMetabolism_mean.list[[i]]) <- paste0(i,"_",colnames(scMetabolism_mean.list[[i]]))
}
pathway_used <- Reduce(intersect, list(row.names(scMetabolism_mean.list$mouse_NAFLD),
                                       row.names(scMetabolism_mean.list$mouse_WD),
                                       row.names(scMetabolism_mean.list$mouse_RBPJ),
                                       row.names(scMetabolism_mean.list$human_liver)))
scMetabolism_mean_all <- cbind(
  scMetabolism_mean.list$mouse_NAFLD[pathway_used,],
  scMetabolism_mean.list$mouse_WD[pathway_used,],
  scMetabolism_mean.list$mouse_RBPJ[pathway_used,],
  scMetabolism_mean.list$human_liver[pathway_used,])
scMetabolism_mean_macs <- cbind(
  scMetabolism_mean.list$mouse_NAFLD[pathway_used,paste0("mouse_NAFLD_",c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC","KCs"))],
  scMetabolism_mean.list$mouse_WD[pathway_used,paste0("mouse_WD_",c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC","KCs"))],
  scMetabolism_mean.list$mouse_RBPJ[pathway_used,paste0("mouse_RBPJ_",c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC","KCs"))],
  scMetabolism_mean.list$human_liver[pathway_used,paste0("human_liver_",c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","moKC","KCs"))])

plot.matrix <- scMetabolism_mean_all
plot.matrix <- t(apply(plot.matrix, 1, scale))
colnames(plot.matrix) <- colnames(scMetabolism_mean_all)
plot.matrix_quantile <- quantile(plot.matrix, c(0.01, 0.99))
plot.matrix <- pmax(plot.matrix, plot.matrix_quantile[1])
plot.matrix <- pmin(plot.matrix, plot.matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(plot.matrix), max(plot.matrix), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
ha_column <-  HeatmapAnnotation(
  df = data.frame(
    Dataset = apply(stringr::str_split_fixed(colnames(plot.matrix),"_",3)[,c(1,2)],1,function(x) paste0(x[1],"_",x[2])),
    Cluster = stringr::str_split_fixed(colnames(plot.matrix),"_",3)[,3]
  ),
  col = list(
    Dataset = c("mouse_NAFLD" = "#fb942a", "mouse_WD" = "#d945a0", "mouse_RBPJ" = "#415f96", "human_liver" = "#c8ebe8"),
    Cluster = cluster_color_panel),
  annotation_name_gp = gpar(fontsize = 7)
)
p <- Heatmap(
  plot.matrix,
  col = color_used,
  cluster_columns = F,
  column_split = factor(apply(stringr::str_split_fixed(colnames(plot.matrix),"_",3)[,c(1,2)],1,function(x) paste0(x[1],"_",x[2]))),
  top_annotation = ha_column,
  row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = 
    list(title_gp = gpar(fontsize = 7),
         labels_gp = gpar(fontsize = 7),
         legend_height = unit(3, "cm"))
)

