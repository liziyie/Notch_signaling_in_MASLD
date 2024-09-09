setwd("/Volumes/Macintosh HD/Users/liziyi/Desktop/")
source("./MASLD_utils.R")
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(grid)
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(scCustomize)

# Load data----
mouse_liver_sc_atlas_macs <- readRDS("./mouse_liver_sc_atlas_macs_release.rds")
NAFLD_cells_used <- mouse_liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "mouse_liver_NAFLD") %>% pull(CellName)
WD_cells_used <- mouse_liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "mouse_liver_WD") %>% pull(CellName)
RBPJ_cells_used <- mouse_liver_sc_atlas_macs@meta.data %>% filter(orig.ident == "RBPJ") %>% pull(CellName)

RBPJ_new <- mouse_liver_sc_atlas_macs %>% subset(orig.ident == "RBPJ")
RBPJ_new@meta.data$Group <- factor(RBPJ_new@meta.data$Group, levels = c("ND_WT","ND_KO","MCD_WT","MCD_KO"))
RBPJ_new <- SetIdent(RBPJ_new, value = RBPJ_new@meta.data$Cluster)

# Figure 3A----
ggplot(mouse_liver_sc_atlas_macs@meta.data[RBPJ_cells_used,], aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = .5, alpha = .8) +
  scale_color_manual(values = cluster_color_panel) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  theme_cowplot(font_size = 7) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1)) +
  ggtitle("RBPJ") +
  addLabels(centers = calcCenters(mouse_liver_sc_atlas_macs@meta.data[RBPJ_cells_used,"Cluster"],
                                  mouse_liver_sc_atlas_macs@meta.data[RBPJ_cells_used,c("UMAP_1","UMAP_2")]))

# Figure 3B----
genes_used <- c("Mki67","Top2a","Ube2c","Stmn1","Mcm6","Mcm2","Mcm5","Cd209f","Cd209g","Cxcl13","Marco","Cd163","Lyve1","Cd36","Lyz2","Folr2","Timd4","Il1b","Mertk","Mrc1","Clec1b","Cd5l","Slc40a1","C1qa","C1qb","Pltp","Adgre1","Cdh5","H2-Eb1","H2-Aa","Id3","Maf","Tmem176a","Ly86","Marcksl1","Clec4e","Il1rn","Cxcl2","Cd14","Ccl4","Ehd1","Ace","Klf4","Ear2","Cd300ld","S100a8","Ly6c2","S100a9","Sell","Vcan","Thbs1","Mmp12","Gpnmb","Cd63","Fabp5","Spp1","Cd9")
gene.mean.matrix <- 
  aggregate(t(as.matrix(RBPJ_new@assays$RNA@data[genes_used,])),
            list(Cluster = RBPJ_new@meta.data[,"Cluster"]), 
            mean)
gene.mean.zscore <- apply(gene.mean.matrix[,2:ncol(gene.mean.matrix)], 2, function(x){(x - mean(x)) / sd(x)})
row.names(gene.mean.zscore) <- gene.mean.matrix[,1]
gene.mean.zscore.df <- 
  data.frame(Group = gene.mean.matrix[,1], gene.mean.zscore, check.names = F) %>% 
  melt(id.vars = "Group", variable.name = "Gene", value.name = "Exp")
gene.mean.zscore <- t(gene.mean.zscore)
max.avg <- apply(gene.mean.zscore, 1, which.max)
gene_order <- c()
for(i in 1:ncol(gene.mean.zscore)){
  if(sum(max.avg == i) == 1){
    temp <- data.frame(Gene = names(max.avg)[max.avg == i], Gene.Group = colnames(gene.mean.zscore)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
  if(sum(max.avg == i) > 1){
    temp <- data.frame(Gene = names(sort(gene.mean.zscore[names(max.avg)[max.avg == i],i], decreasing = T)), Gene.Group = colnames(gene.mean.zscore)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
}
gene.per <- 
  aggregate(t(as.matrix(RBPJ_new@assays$RNA@data[genes_used,])),
            list(Group = RBPJ_new@meta.data[,"Cluster"]),
            function(x){sum(x > 2) / length(x)}) %>% 
  melt(id.vars = "Group", variable.name = "Gene", value.name = "Per")
plot.data <- merge(merge(gene.mean.zscore.df, gene.per), gene_order)
plot.data$Gene <- as.character(plot.data$Gene)
plot.data[plot.data$Gene %in% c("C1qb"), "Gene.Group"] <- "C1qc+ Macs"
plot.data$Group <- factor(plot.data$Group, levels = rev(c("Ly6Chi Monocytes","Patrolling Monocytes","Clec4e+ Macs","C1qc+ Macs","Capsule Macs","Spp1+ Macs","moKC","KCs","Prolif. KCs")))
plot.data$Gene.Group <- factor(plot.data$Gene.Group, levels = c("Ly6Chi Monocytes","Patrolling Monocytes","Clec4e+ Macs","C1qc+ Macs","Capsule Macs","Spp1+ Macs","moKC","KCs","Prolif. KCs"))
plot.data.order <- plot.data[order(plot.data$Gene.Group, plot.data$Exp),]
plot.data$Gene <- factor(plot.data$Gene, levels = unique(plot.data.order$Gene))
plot.data$Gene.Group <- as.character(plot.data$Gene.Group)
plot.data$Gene.Group <- plyr::revalue(
  plot.data$Gene.Group,
  replace = c("Patrolling Monocytes" = "Mono", "Ly6Chi Monocytes" = "Mono",
              "Clec4e+ Macs" = "Macro", "Spp1+ Macs" = "Macro", "C1qc+ Macs" = "Macro",
              "Capsule Macs" = "Macro",
              "moKC" = "KC", "KCs" = "KC", "Prolif. KCs" = "KC(pro)")
)
plot.data$Gene.Group <- factor(plot.data$Gene.Group, levels = c("KC(pro)","KC","Macro","Mono"))
myColorPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
ggplot(plot.data, aes(x = Group, y = Gene)) +
  geom_point(aes(size = Per, fill = Exp, color = Exp), shape = 21) +
  facet_grid(Gene.Group~., scales = "free", space = "free") +
  theme_minimal() + 
  scale_x_discrete(limits = rev(levels(plot.data$Group))) +
  theme(
    axis.ticks = element_line(size = .5),
    axis.ticks.length = unit(0.05,"in"),
    strip.text = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", angle = 90, 
                               hjust = 1, vjust = 0.5, size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "lightgrey",
      linetype = "dashed",
      size = 0.2
    )) + 
  labs(x = "", y = "") +
  scale_fill_gradientn("Exp",
                       colours = myColorPalette(100), 
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black", barwidth = 0.8)) +
  scale_color_gradientn("Exp",
                        colours = myColorPalette(100), 
                        guide = guide_colorbar(frame.colour = "black",
                                               ticks.colour = "black", barwidth = 0.8)) +
  scale_size_continuous("Percentage",
                        breaks = seq(0, 0.8, 0.2), range = c(1,3))

# Figure 3C----
p1 <- ggplot(RBPJ_new@meta.data %>% filter(Group == "ND_WT"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 1, alpha = .8) +
  theme_cowplot(font_size = 7) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  scale_color_manual(name = "", values = cluster_color_panel) +
  ggforce::theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none") + 
  ggtitle("ND WT")
p2 <- ggplot(RBPJ_new@meta.data %>% filter(Group == "ND_KO"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 1, alpha = .8) +
  theme_cowplot(font_size = 7) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  scale_color_manual(name = "", values = cluster_color_panel) +
  ggforce::theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none") + 
  ggtitle("ND KO")
p3 <- ggplot(RBPJ_new@meta.data %>% filter(Group == "MCD_WT"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 1, alpha = .8) +
  theme_cowplot(font_size = 7) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  scale_color_manual(name = "", values = cluster_color_panel) +
  ggforce::theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none") + 
  ggtitle("MCD WT")
p4 <- ggplot(RBPJ_new@meta.data %>% filter(Group == "MCD_KO"), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = 1, alpha = .8) +
  theme_cowplot(font_size = 7) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  scale_color_manual(name = "", values = cluster_color_panel) +
  ggforce::theme_no_axes() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none") + 
  ggtitle("MCD KO")
plot_grid(p1, p2, p3, p4, nrow = 2)

# Figure 3D----
myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:6])
RBPJ_new@meta.data[,c("Cluster","Group")] %>%
  mutate(Cluster = forcats::fct_rev(Cluster)) %>%
  table() %>% ROIE() %>% 
  ROIE_plot(font.size = 7, text.size = 2.5) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 12))

# Figure 3E----
gene.A <- "Timd4"
gene.B <- "Clec4f"
plot.data <- data.frame(
  Cluster = RBPJ_new@meta.data$Cluster,
  Timd4 = RBPJ_new@assays$RNA@data[gene.A,],
  Clec4f = RBPJ_new@assays$RNA@data[gene.B,]
)
plot.data <- plot.data %>% filter(Cluster %in% c("Ly6Chi Monocytes","Patrolling Monocytes","C1qc+ Macs","Clec4e+ Macs","Spp1+ Macs","moKC","KCs"))
plot.data$Cluster <- as.character(plot.data$Cluster)
plot.data$Cluster[plot.data$Cluster %in% c("Ly6Chi Monocytes","Patrolling Monocytes","C1qc+ Macs","Clec4e+ Macs","Spp1+ Macs")] <- "Mono/Mac"
plot.data$Cluster <- factor(plot.data$Cluster, levels = c("Mono/Mac","moKC","KCs"))
x.cutoff = 1.5
y.cutoff = 1.5
plot.data[,gene.A] <- plot.data[,gene.A] + abs(rnorm(nrow(plot.data))/1000)
plot.data[,gene.B] <- plot.data[,gene.B] + abs(rnorm(nrow(plot.data))/1000)
color.by = "Cluster"
color_panel <- c("Mono/Mac" = "#73458a", "moKC" = "#F0CE39","KCs" = "#60B4F3")
myColorPalette <- colorRampPalette(rev(brewer.pal(8, "RdBu")))
pmain <- 
  ggplot(plot.data, aes_string(x = gene.A, y = gene.B)) +
  geom_density_2d(data = plot.data[plot.data[,gene.A] > 0.1 & plot.data[,gene.B] > 0.1,], aes(color = Cluster)) +
  geom_point(size = -1) +
  scale_color_manual("", values = color_panel) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.A]) + 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.B]) + 0.5)) +
  facet_wrap(~Cluster) +
  geom_hline(yintercept = y.cutoff, linetype = "dashed") +  
  geom_vline(xintercept = x.cutoff, linetype = "dashed") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  ) +
  guides(alpha = FALSE, fill = FALSE)

gene.A <- "Cd74"
gene.B <- "Ly6c2"
plot.data <- data.frame(
  Cluster = RBPJ_new@meta.data$Cluster,
  Cd74 = RBPJ_new@assays$RNA@data[gene.A,],
  Ly6c2 = RBPJ_new@assays$RNA@data[gene.B,]
)
plot.data <- plot.data %>% filter(Cluster %in% c("Ly6Chi Monocytes","C1qc+ Macs","Clec4e+ Macs","Spp1+ Macs"))
plot.data$Cluster <- as.character(plot.data$Cluster)
plot.data$Cluster[plot.data$Cluster %in% c("C1qc+ Macs","Spp1+ Macs")] <- "Other mononuclear"
plot.data$Cluster <- factor(plot.data$Cluster, levels = c("Ly6Chi Monocytes","Clec4e+ Macs","Other mononuclear"))
x.cutoff = 1.5
y.cutoff = 1.5
plot.data[,gene.A] <- plot.data[,gene.A] + abs(rnorm(nrow(plot.data))/1000)
plot.data[,gene.B] <- plot.data[,gene.B] + abs(rnorm(nrow(plot.data))/1000)
color.by = "Cluster"
color_panel <- c("Ly6Chi Monocytes" = "#C45240", "Clec4e+ Macs" = "#63B571", "Other mononuclear" = "#748198")
myColorPalette <- colorRampPalette(rev(brewer.pal(8, "RdBu")))
pmain <- 
  ggplot(plot.data, aes_string(x = gene.A, y = gene.B)) +
  geom_density_2d(data = plot.data[plot.data[,gene.A] > 0.1 & plot.data[,gene.B] > 0.1,], aes(color = Cluster)) +
  geom_point(size = -1) +
  scale_color_manual("", values = color_panel) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.A]) + 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.B]) + 0.5)) +
  facet_wrap(~Cluster) +
  geom_hline(yintercept = y.cutoff, linetype = "dashed") +  
  geom_vline(xintercept = x.cutoff, linetype = "dashed") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  ) +
  guides(alpha = FALSE, fill = FALSE)

gene.A <- "Cd36"
gene.B <- "Cd63"
plot.data <- data.frame(
  Cluster = RBPJ_new@meta.data$Cluster,
  Cd36 = RBPJ_new@assays$RNA@data[gene.A,],
  Cd63 = RBPJ_new@assays$RNA@data[gene.B,]
)
plot.data <- plot.data %>% filter(Cluster %in% c("C1qc+ Macs","Spp1+ Macs"))
plot.data$Cluster <- as.character(plot.data$Cluster)
x.cutoff = 1.5
y.cutoff = 1.5
plot.data[,gene.A] <- plot.data[,gene.A] + abs(rnorm(nrow(plot.data))/1000)
plot.data[,gene.B] <- plot.data[,gene.B] + abs(rnorm(nrow(plot.data))/1000)
color.by = "Cluster"
color_panel <- c("C1qc+ Macs" = "#F8C3B4", "Spp1+ Macs" = "#C86EA8")
myColorPalette <- colorRampPalette(rev(brewer.pal(8, "RdBu")))
pmain <- 
  ggplot(plot.data, aes_string(x = gene.A, y = gene.B)) +
  geom_density_2d(data = plot.data[plot.data[,gene.A] > 0.1 & plot.data[,gene.B] > 0.1,], aes(color = Cluster)) +
  scale_color_manual("", values = color_panel) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.A]) + 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot.data[,gene.B]) + 0.5)) +
  geom_hline(yintercept = y.cutoff, linetype = "dashed") +  
  geom_vline(xintercept = x.cutoff, linetype = "dashed") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  ) +
  guides(alpha = FALSE, fill = FALSE)
xdens <- axis_canvas(pmain, axis = "x") +
  geom_density(data = plot.data[plot.data[,gene.A] > 0.1 & plot.data[,gene.B] > 0.1,], 
               aes_string(x = gene.A, fill = color.by),
               color = NA, alpha = 0.7, size = 0.2) +
  scale_fill_manual(values = color_panel)
p0 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE) +
  geom_density(data = plot.data[plot.data[,gene.A] > 0.1 & plot.data[,gene.B] > 0.1,], 
               aes_string(x = gene.B, fill = color.by),
               color = NA, alpha = 0.7, size = 0.2) +
  coord_flip() +
  scale_fill_manual(values = color_panel)
p <- insert_yaxis_grob(p0, ydens, grid::unit(.2, "null"), position = "right")
grid.draw(p)

# Figure 3I and Figure S3H----
momac_bulk <- read.csv("./momac_sorting_tpm.csv", row.names = 1)
colnames(momac_bulk) <- c("Clec4e+ Macs","Spp1+ Macs","C1qc+ Macs","Tim4+ moKC","Tim4- moKC","ResKCs")
momac_bulk <- momac_bulk[,c("Clec4e+ Macs","C1qc+ Macs","Spp1+ Macs","Tim4- moKC","Tim4+ moKC","ResKCs")]
sc_cells_used <- RBPJ_new@meta.data %>% filter(Cluster %ni% c("Patrolling Monocytes","Prolif. KCs", "Capsule Macs", "Ly6Chi Monocytes")) %>% pull(CellName)
RBPJ_new_subset_markers <- FindAllMarkers(RBPJ_new %>% subset(CellName %in% sc_cells_used))

genes_used <- intersect(row.names(RBPJ_new),row.names(momac_bulk))
genes_used <- genes_used[!grepl("^Ig|^Gm|^Tcr|Rik|-|[.]",genes_used)]
genes_used1 <- RBPJ_new_subset_markers %>% filter(gene %in% genes_used, avg_log2FC > 0, pct.1 > .3) %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 100) %>% pull(gene) %>% unique()
genes_used2 <- names(sort(apply(momac_bulk[genes_used,], 1, sd), decreasing = T))[1:500]
genes_used <- unique(c(genes_used1, genes_used2))
sc_bulk_comparison <- Compare2Clusters(
  exp1 = RBPJ_new@assays$RNA@data[genes_used,sc_cells_used],
  cluster1 = RBPJ_new@meta.data[sc_cells_used,"Cluster"] %>% 
    plyr::revalue(replace = c("Spp1+ Macs" = "MDMs", "Clec4e+ Macs" = "MDMs",
                              "C1qc+ Macs" = "MDMs")),
  exp2 = momac_bulk[genes_used,],
  cluster2 = colnames(momac_bulk) %>% 
    plyr::revalue(replace = c("Spp1+ Macs" = "MDMs", "Clec4e+ Macs" = "MDMs",
                              "C1qc+ Macs" = "MDMs")),
  name1 = "sc",
  name2 = "bulk"
)

genes_labeled1 <- RBPJ_new_subset_markers %>% filter(gene %in% genes_used, avg_log2FC > 0, pct.1 > .3) %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 3) %>% pull(gene) %>% unique()
genes_labeled2 <- names(sort(apply(momac_bulk[genes_used,], 1, sd), decreasing = T))[1:15]
genes_labeled <- unique(c(genes_labeled1, genes_labeled2,"Cd209f","Cxcl13","Cd163","Lyve1","Cd36","Lyz2","Folr2","Timd4","Il1b","Mertk","Clec4f","Slc40a1","Pltp","Clec4e","Ly6c2","Vcan","Gpnmb","Cd63","Fabp5","Spp1","Cd9"))
plot.matrix <- apply(sc_bulk_comparison$merged_exp,1,scale)
plot.matrix <- t(plot.matrix)
colnames(plot.matrix) <- colnames(sc_bulk_comparison$merged_exp)
plot.matrix_quantile <- quantile(plot.matrix, c(0.01, 0.99))
plot.matrix <- pmax(plot.matrix, plot.matrix_quantile[1])
plot.matrix <- pmin(plot.matrix, plot.matrix_quantile[2])
color_used <- circlize::colorRamp2(seq(min(plot.matrix), max(plot.matrix), length = 9), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:10]))
Heatmap(plot.matrix,
             cluster_rows = T, 
             cluster_columns = as.dendrogram(sc_bulk_comparison$hc_fit),
             show_row_names = T,
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 5),
             col = color_used,
             name = "Exp") +
  rowAnnotation(link = anno_mark(
    at = which(row.names(plot.matrix) %in% genes_labeled),
    labels = row.names(plot.matrix)[which(row.names(plot.matrix) %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

dend <- sc_bulk_comparison$hc_fit %>% as.dendrogram %>%
  set("labels_cex", .6) %>%
  set("branches_lwd", 2.5) %>%
  set("labels_colors", 
      c("#de767b", "#4c8ed3","#de767b", "#4c8ed3","#4c8ed3","#de767b", "#4c8ed3"))
plot(dend)

# Figure S3A----
RBPJ_raw <- readRDS("./RBPJ.rds")
RBPJ_raw@meta.data$Cluster2 <- as.character(RBPJ_raw@meta.data$MajorCluster)
RBPJ_raw@meta.data[RBPJ_raw@meta.data$MajorCluster %in% c("DC","Neutrophil","Lymphocyte","CD45-"),"Cluster2"] <- as.character(RBPJ_raw@meta.data[RBPJ_raw@meta.data$MajorCluster %in% c("DC","Neutrophil","Lymphocyte","CD45-"),"Cluster"])
RBPJ_raw_color_panel <- c("CD4 T" = "#86C2D7", "CD8 T" = "#2D75B7", "NK" = "#007FB7", "B cell" = "#EA95CE", "Plasma B" = "#50C2FB", "Plasma B(pro)" = "#568CC1", "Neutrophil" = "#C6D9DF", "Neutrophil(pro)" = "#918FBD", "cDC1" = "#F0BED7", "cDC2" = "#EEE1E8", "mregDC" = "#FC7542", "cDC(pro)" = "#E69E8E", "Monocyte" = "#F2C987", "Macrophage" = "#63B571", "KC" = "#60B4F3", "Epithelium" = "#EDB889", "Endothelium" = "#D4B696")
RBPJ_raw@meta.data$Cluster2 <- factor(RBPJ_raw@meta.data$Cluster2, levels = names(RBPJ_raw_color_panel))

p <- ggplot(RBPJ_raw@meta.data %>% arrange(Cluster2), aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster2), size = .5, alpha = .8) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = RBPJ_raw_color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

p <- ggplot(RBPJ_raw@meta.data %>% filter(Cluster %in% c("KC(pro)_G2M","KC(pro)_S","Mertk+ KC","Mertk- KC","Lyve1+ KC")), aes(x = KC_UMAP_1, y = KC_UMAP_2)) +
  geom_point(aes(color = Cluster), size = .5, alpha = .8) +
  theme_cowplot(font_size = 7) +
  scale_color_manual(name = "", values = color_panel) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Figure S3B----
RBPJ_raw <- SetIdent(RBPJ_raw, value = RBPJ_raw@meta.data$Cluster2)
RBPJ_all_markers <- FindAllMarkers(RBPJ_raw)
genes_used <- RBPJ_all_markers %>% group_by(cluster) %>% filter(pct.2 < 0.8, pct.1 > 0.5) %>% slice_max(n = 20, order_by = avg_log2FC) %>% pull(gene) %>% as.character() %>% unique()
genes_used <- unique(c(genes_used, "Cd8a", "Ly6c2", "Lyve1", "Epcam"))
markers_mean_exp <- 
  aggregate(
    as.matrix(t(RBPJ_raw@assays$RNA@data[genes_used,])),
    list(Cluster = RBPJ_raw$Cluster2),
    mean)
row.names(markers_mean_exp) <- markers_mean_exp$Cluster
markers_mean_exp$Cluster <- c()
markers_plot_matrix <- apply(markers_mean_exp, 2, scale) %>% t()
colnames(markers_plot_matrix) <- row.names(markers_mean_exp)
color_used <- circlize::colorRamp2(seq(min(markers_plot_matrix), max(markers_plot_matrix), length = 9), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:10]))
max.avg <- apply(markers_plot_matrix, 1, which.max)
gene_order <- c()
for(i in 1:ncol(markers_plot_matrix)){
  if(sum(max.avg == i) == 1){
    temp <- data.frame(gene = names(max.avg)[max.avg == i], cluster = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
  if(sum(max.avg == i) > 1){
    temp <- data.frame(gene = names(sort(markers_plot_matrix[names(max.avg)[max.avg == i],i], decreasing = T)), cluster = colnames(markers_plot_matrix)[i], stringsAsFactors = F)
    gene_order <- rbind(gene_order, temp)
  }
}
gene_order <- gene_order[rev(1:nrow(gene_order)),]
markers_plot_matrix_quantile <- quantile(markers_plot_matrix, c(0.01, 0.99))
markers_plot_matrix <- pmax(markers_plot_matrix, markers_plot_matrix_quantile[1])
markers_plot_matrix <- pmin(markers_plot_matrix, markers_plot_matrix_quantile[2])
genes_labeled <- c("Cd3d","Cd8a","Ncr1","Klrb1c","Nkg7","Cd19","Cd79a","Jchain","Ighg2b","Ube2c","Mki67","Retnlg","S100a9","S100a8","Elane","Ms4a3","Clec9a","Flt3","Irf8","Cd209a","Itgax","Ccr7","Fscn1","Stmn1","Top2a","Nr4a1","Cx3cr1","Cebpb","Vcan","S100a6","Ly6c2","S100a6","Clec4e","Cd14","Gpnmb","Spp1","Marco","Cd163","Lyve1","Folr2","Clec4f","Timd4","Mcm6","Krt18","Epcam","Ptprb","Ehd3","Gpr182")
Heatmap(markers_plot_matrix[gene_order$gene,],
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 5),
             col = color_used,
             name = "Exp") +
  rowAnnotation(link = anno_mark(
    at = which(gene_order$gene %in% genes_labeled),
    labels = gene_order$gene[which(gene_order$gene %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

# Figure S3C----
p <- ggplot(mouse_liver_sc_atlas@meta.data[NAFLD_cells_used,],
            aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = .5, alpha = .8) +
  scale_color_manual(values = cluster_color_panel) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  ggtitle("GSE192742") +
  theme_cowplot(font_size = 7) +
  addLabels(centers = calcCenters(mouse_liver_sc_atlas@meta.data[NAFLD_cells_used,"Cluster"],
                                  mouse_liver_sc_atlas@meta.data[NAFLD_cells_used,c("UMAP_1","UMAP_2")])) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

p <- ggplot(mouse_liver_sc_atlas@meta.data[WD_cells_used,], 
            aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = .5, alpha = .8) +
  scale_color_manual(values = cluster_color_panel) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  ggtitle("GSE156059") +
  theme_cowplot(font_size = 7) +
  addLabels(centers = calcCenters(mouse_liver_sc_atlas@meta.data[WD_cells_used,"Cluster"],
                                  mouse_liver_sc_atlas@meta.data[WD_cells_used,c("UMAP_1","UMAP_2")])) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

p <- ggplot(mouse_liver_sc_atlas@meta.data[RBPJ_cells_used,], 
            aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Cluster), size = .5, alpha = .8) +
  scale_color_manual(values = cluster_color_panel) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  ggtitle("RBPJ") +
  theme_cowplot(font_size = 7) +
  addLabels(centers = calcCenters(mouse_liver_sc_atlas@meta.data[RBPJ_cells_used,"Cluster"],
                                  mouse_liver_sc_atlas@meta.data[RBPJ_cells_used,c("UMAP_1","UMAP_2")])) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Figure S3D----
p <- ggplot(mouse_liver_sc_atlas@meta.data[NAFLD_cells_used,], 
            aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Raw_Cluster), size = .5, alpha = .8) +
  scale_color_manual(values = raw_cluster_color_panel) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  ggtitle("GSE192742") +
  theme_cowplot(font_size = 7) +
  addLabels(centers = calcCenters(mouse_liver_sc_atlas@meta.data[NAFLD_cells_used,"Raw_Cluster"],
                                  mouse_liver_sc_atlas@meta.data[NAFLD_cells_used,c("UMAP_1","UMAP_2")])) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

p <- ggplot(mouse_liver_sc_atlas@meta.data[WD_cells_used,], 
            aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = Raw_Cluster), size = .5, alpha = .8) +
  scale_color_manual(values = raw_cluster_color_panel) +
  xlim(c(-8.2,6.3)) + ylim(c(-4,3.5)) +
  ggtitle("GSE156059") +
  theme_cowplot(font_size = 7) +
  addLabels(centers = calcCenters(mouse_liver_sc_atlas@meta.data[WD_cells_used,"Raw_Cluster"],
                                  mouse_liver_sc_atlas@meta.data[WD_cells_used,c("UMAP_1","UMAP_2")])) +
  guides(color = guide_legend(override.aes = list(size = 3, height = 1), ncol = 1))

# Figure S3E----
alluvial_cluster_color_panel <- cluster_color_panel
names(alluvial_cluster_color_panel) <- paste0("Atlas_", names(alluvial_cluster_color_panel))
library(ggalluvial)
p <- mouse_liver_sc_atlas_macs@meta.data %>% 
  filter(orig.ident %in% c("mouse_liver_NAFLD"), Raw_Cluster %ni% c("cDC1s","cDC2s","Mig cDCs")) %>% 
  select(Raw_Cluster,  Cluster) %>% table() %>% as.data.frame() %>% filter(Freq > 0) %>%
  mutate(Raw_Cluster = factor(paste0("NAFLD_",Raw_Cluster), 
                              levels = paste0("NAFLD_",
                                              intersect(names(raw_cluster_color_panel),
                                                        as.character(unique(Raw_Cluster)))))) %>%
  mutate(Cluster = factor(paste0("Atlas_",Cluster), 
                          levels = paste0("Atlas_",
                                          intersect(names(cluster_color_panel),
                                                    as.character(unique(Cluster)))))) %>%
  ggplot(aes(axis1 = Cluster, axis2 = Raw_Cluster, y = Freq)) +
  geom_alluvium(aes(fill = Cluster), curve_type = "sigmoid", width = 1/8) +
  geom_stratum(width = 1/8) +
  ggtitle("RBPJ - GSE192742") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = alluvial_cluster_color_panel) +
  theme_void(base_size = 7) +
  theme(legend.position = "none")

p <- mouse_liver_sc_atlas_macs@meta.data %>% 
  filter(orig.ident %in% c("mouse_liver_WD"), Raw_Cluster %ni% c("cDC1s","cDC2s","Mig cDCs")) %>% 
  select(Raw_Cluster,  Cluster) %>% table() %>% as.data.frame() %>% filter(Freq > 0) %>%
  mutate(Raw_Cluster = factor(paste0("WD_",Raw_Cluster), 
                              levels = paste0("WD_",
                                              intersect(names(raw_cluster_color_panel),
                                                        as.character(unique(Raw_Cluster)))))) %>%
  mutate(Cluster = factor(paste0("Atlas_",Cluster), 
                          levels = paste0("Atlas_",
                                          intersect(names(cluster_color_panel),
                                                    as.character(unique(Cluster)))))) %>%
  ggplot(aes(axis1 = Cluster, axis2 = Raw_Cluster, y = Freq)) +
  geom_alluvium(aes(fill = Cluster), curve_type = "sigmoid", width = 1/8) +
  geom_stratum(width = 1/8) +
  ggtitle("RBPJ - GSE156059") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Survey", "Response"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = alluvial_cluster_color_panel) +
  theme_void(base_size = 7) +
  theme(legend.position = "none")

# Figure S3F----
p1 <- mouse_liver_sc_atlas_macs@meta.data[mouse_liver_sc_atlas_macs@meta.data$orig.ident == "mouse_liver_NAFLD",c("Cluster","Group")] %>%
  mutate(Cluster = forcats::fct_rev(Cluster)) %>%
  table() %>% ROIE() %>%
  ROIE_plot(font.size = 7, text.size = 2.5) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 12))
p2 <- mouse_liver_sc_atlas_macs@meta.data[mouse_liver_sc_atlas_macs@meta.data$orig.ident == "mouse_liver_WD",c("Cluster","Group")] %>%
  mutate(Cluster = forcats::fct_rev(Cluster)) %>%
  table() %>% ROIE() %>%
  ROIE_plot(font.size = 7, text.size = 2.5) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 12))

# Figure S3G----
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_WT"), Cluster %in% c("Clec4e+ Macs","Spp1+ Macs")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Cluster"]
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "Clec4e+ Macs"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "Spp1+ Macs"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol),deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol))
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val+1e-294))) +
  geom_point(aes(color = Sig), size = 1) +
  scale_color_manual(values = c(cluster_color_panel[c("Clec4e+ Macs","Spp1+ Macs")], "FALSE" = "lightgrey")) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none") +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc')))

cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_WT"), Cluster %in% c("Clec4e+ Macs","C1qc+ Macs")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Cluster"]
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "C1qc+ Macs"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "Clec4e+ Macs"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol),deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol),"Ly6c2")
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val+1e-294))) +
  geom_point(aes(color = Sig), size = 1) +
  scale_color_manual(values = c(cluster_color_panel[c("Clec4e+ Macs","C1qc+ Macs")], "FALSE" = "lightgrey")) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none") +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc')))

cells_used <- RBPJ_new@meta.data %>% filter(Group == "MCD_WT", Cluster %in% c("Spp1+ Macs","C1qc+ Macs")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Cluster"]
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "C1qc+ Macs"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "Spp1+ Macs"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol),deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol),"Ly6c2")
genes_sig <- genes_sig[genes_sig %ni% c("C1qc")]
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val+1e-318))) +
  geom_point(aes(color = Sig), size = 1) +
  scale_color_manual(values = c(cluster_color_panel[c("Spp1+ Macs","C1qc+ Macs")], "FALSE" = "lightgrey")) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none") +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc')))

cells_used <- RBPJ_new@meta.data %>% filter(Group == "MCD_WT", Cluster %in% c("Ly6Chi Monocytes","Clec4e+ Macs")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Cluster"]
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "Ly6Chi Monocytes"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "Clec4e+ Macs"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol),deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol),"Cd74")
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val+1e-318))) +
  geom_point(aes(color = Sig), size = 1cells_used <- RBPJ_new@meta.data %>% filter(Group == "MCD_WT", Cluster %in% c("moKC","KCs")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Cluster"]
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "moKC"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "KCs"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt|^Rp|^Gm",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol),deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol), "Timd4")
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val+1e-318))) +
  geom_point(aes(color = Sig), size = 2) +
  scale_color_manual(values = c(cluster_color_panel[c("moKC","KCs")], "FALSE" = "lightgrey")) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none")
p.5 <- ggAIplot(volcano_plot) +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc')))) +
  scale_color_manual(values = c(cluster_color_panel[c("Clec4e+ Macs","Ly6Chi Monocytes")], "FALSE" = "lightgrey")) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none") +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc')))

cells_used <- RBPJ_new@meta.data %>% filter(Group == "MCD_WT", Cluster %in% c("moKC","KCs")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Cluster"]
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "moKC"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "KCs"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt|^Rp|^Gm",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol),deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol), "Timd4")
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val+1e-318))) +
  geom_point(aes(color = Sig), size = 1) +
  scale_color_manual(values = c(cluster_color_panel[c("moKC","KCs")], "FALSE" = "lightgrey")) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none") +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc')))

cells_used <- RBPJ_new@meta.data %>% filter(Group == "MCD_WT", Cluster %in% c("moKC","Spp1+ Macs","Clec4e+ Macs","C1qc+ Macs")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Cluster"] == "moKC"
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "MDM"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "MoKC"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt|^Rp|^Gm",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol),deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol))
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val+1e-318))) +
  geom_point(aes(color = Sig), size = 1) +
  scale_color_manual(values = c("MoKC" = "#F0CE39", "MDM" = "black", "FALSE" = "lightgrey")) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none") +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc')))