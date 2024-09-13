source("./MASLD_utils.R")
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(scCustomize)

# Load Data----
mouse_liver_sc_atlas_macs <- readRDS("./mouse_liver_sc_atlas_macs_release.rds")
RBPJ_new <- mouse_liver_sc_atlas_macs %>% subset(orig.ident == "RBPJ")
RBPJ_new@meta.data$Group <- factor(RBPJ_new@meta.data$Group, levels = c("ND_WT","ND_KO","MCD_WT","MCD_KO"))
RBPJ_new <- SetIdent(RBPJ_new, value = RBPJ_new@meta.data$Cluster)

group_color_panel <- c("ND_WT" = "#e08380", "ND_KO" = "#82be9c", "MCD_WT" = "#e08380", "MCD_KO" = "#82be9c")
cluster_color_panel <- c(
  "Ly6Chi Monocytes" = "#C45240", "Patrolling Monocytes" = "#F2C987",
  "Clec4e+ Macs" = "#63B571", "C1qc+ Macs" = "#F8C3B4",
  "Capsule Macs" = "#86C2D7",
  "Spp1+ Macs" = "#C86EA8",
  "moKC" = "#F0CE39", "KCs" = "#60B4F3",
  "Prolif. KCs" = "#F0C1ED"
)

RBPJ_raw <- readRDS("./RBPJ.rds")
RBPJ_raw@meta.data$Cluster2 <- as.character(RBPJ_raw@meta.data$MajorCluster)
RBPJ_raw@meta.data[RBPJ_raw@meta.data$MajorCluster %in% c("DC","Neutrophil","Lymphocyte","CD45-"),"Cluster2"] <- as.character(RBPJ_raw@meta.data[RBPJ_raw@meta.data$MajorCluster %in% c("DC","Neutrophil","Lymphocyte","CD45-"),"Cluster"])
RBPJ_raw_color_panel <- c("CD4 T" = "#86C2D7", "CD8 T" = "#2D75B7", "NK" = "#007FB7", "B cell" = "#EA95CE", "Plasma B" = "#50C2FB", "Plasma B(pro)" = "#568CC1", "Neutrophil" = "#C6D9DF", "Neutrophil(pro)" = "#918FBD", "cDC1" = "#F0BED7", "cDC2" = "#EEE1E8", "mregDC" = "#FC7542", "cDC(pro)" = "#E69E8E", "Monocyte" = "#F2C987", "Macrophage" = "#63B571", "KC" = "#60B4F3", "Epithelium" = "#EDB889", "Endothelium" = "#D4B696")
RBPJ_raw@meta.data$Cluster2 <- factor(RBPJ_raw@meta.data$Cluster2, levels = names(RBPJ_raw_color_panel))

# Figure 6C----
patro_markers <- c("Nr4a1","Ace","Ear2","Spn")
cells_used <- RBPJ_new@meta.data %>% filter(Group == "MCD_KO") %>% pull(CellName)
p <- 
  RBPJ_new@assays$RNA@data[patro_markers,cells_used] %>% 
  as.data.frame() %>% 
  mutate(Gene = patro_markers) %>% 
  melt(id.vars = "Gene") %>% 
  mutate(Gene = factor(Gene, levels = rev(patro_markers))) %>% 
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

# Figure 6G----
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(clusterProfiler)
library(msigdbr)
mouse_db <- msigdbr(species = "Mus musculus")
mouse_db_used = mouse_db %>% dplyr::select(gs_name,gene_symbol) %>% split(x = .$gene_symbol, f = .$gs_name)
exp <- as.matrix(RBPJ_new@assays$RNA@data[,RBPJ_new@meta.data$Cluster == "Patrolling Monocytes"])
gsva.matrix <- gsva(exp, mouse_db_used, parallel.sz=2L)
saveRDS(gsva.matrix, file = "./02.processed_data/RBPJ_new_gsva.rds")
gsva.matrix <- readRDS("./02.processed_data/RBPJ_new_gsva.rds")
gsva.matrix.mean <- aggregate(t(gsva.matrix), list(as.character(RBPJ_new@meta.data[colnames(gsva.matrix),"Group"])), mean)
row.names(gsva.matrix.mean) <- gsva.matrix.mean$Group.1
gsva.matrix.mean$Group.1 <- c()
gsva.matrix.mean <- t(gsva.matrix.mean)
pathways_used <- c(
  "GOBP_EXTRACELLULAR_MATRIX_CELL_SIGNALING",
  "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
  "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT",
  "GOBP_NEGATIVE_REGULATION_OF_MONOCYTE_CHEMOTAXIS",
  "GOBP_METANEPHRIC_GLOMERULUS_VASCULATURE_DEVELOPMENT",
  "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX",
  "GOMF_MHC_CLASS_II_RECEPTOR_ACTIVITY",
  "BIOCARTA_INTEGRIN_PATHWAY",
  "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX",
  "GOCC_FIBRINOGEN_COMPLEX"
)
mat_for_plot <- gsva.matrix.mean[pathways_used,]
mat_for_plot <- apply(mat_for_plot,1,scale)
row.names(mat_for_plot) <- colnames(gsva.matrix.mean)
mat_for_plot <- t(mat_for_plot)
row.names(mat_for_plot) <- tolower(row.names(mat_for_plot))
color_used <- circlize::colorRamp2(seq(min(mat_for_plot, na.rm = T), max(mat_for_plot, na.rm = T), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
mat_for_plot <- mat_for_plot[,c("ND_WT","ND_KO","MCD_WT","MCD_KO")]
row.names(mat_for_plot) <- stringr::str_split_fixed(row.names(mat_for_plot),"_",2)[,2]
row.names(mat_for_plot) <- gsub("_"," ",row.names(mat_for_plot))
ha_column <-  HeatmapAnnotation(
  df = data.frame(Group = factor(colnames(mat_for_plot),levels = c("ND_WT","ND_KO","MCD_WT","MCD_KO"))),
  col = list(Group = c("ND_WT" = "#e08380", "ND_KO" = "#82be9c", "MCD_WT" = "#e08380", "MCD_KO" = "#82be9c")),
  annotation_name_gp= gpar(fontsize = 7)
)
p <- Heatmap(
  mat_for_plot,
  col = color_used,
  cluster_columns = F,
  column_split = factor(c(1,1,2,2)),
  top_annotation = ha_column,
  row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = 
    list(title_gp = gpar(fontsize = 7),
         labels_gp = gpar(fontsize = 7),
         legend_height = unit(3, "cm"))
)

# Figure 6J and Figure S6K----
cells_used <- RBPJ_raw@meta.data %>% filter(Cluster %in% c("Endothelium")) %>% pull(CellName)
inflammation_genes <- c("Il1b","Ccl6","Il6","Il27","Tnf","Il18","Cxcl2","Socs3","Marco","Ccl2","Ccl3","Ccl5","Spp1")
plot.matrix <- aggregate(t(RBPJ_raw@assays$RNA@data[inflammation_genes,cells_used]),list(RBPJ_raw@meta.data[cells_used,"Group"]),mean)
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

RBPJ_subset <- RBPJ_raw[,cells_used]
RBPJ_subset <- AddModuleScore(RBPJ_subset, 
                              list(inflammation = inflammation_genes), 
                              ctrl = length(inflammation_genes), 
                              name = "Inflammation")
p <- 
  ggplot(RBPJ_subset@meta.data[cells_used,], aes(x = Group, y = Inflammation1)) +
  geom_quasirandom(aes(col = Group), cex = 1, width = 0.25, alpha = 0.5) +
  geom_boxplot(aes(color = Group), outlier.size = -1, alpha = 0) +
  geom_signif(comparisons = list(c("MCD_WT","MCD_KO"),c("MCD_WT","ND_WT"),c("MCD_KO","ND_KO")), step_increase = .05, tip_length = 0, textsize = 2) +
  labs(y = "Pro-inflammatory cytokines and chemokines") +
  scale_color_manual(values = group_color_panel) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "null",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text = element_blank(),
        axis.title.x = element_blank())

# Figure S6A----
gene.A <- "Fcgr4"
gene.B <- "Spn"
plot.data <- data.frame(
  Cluster = RBPJ_new@meta.data$Cluster,
  Group = RBPJ_new@meta.data$Group,
  Fcgr4 = RBPJ_new@assays$RNA@data["Fcgr4",],
  Spn = RBPJ_new@assays$RNA@data["Spn",]
)
plot.data <- plot.data %>% filter(Group == "MCD_KO")
plot.data$Cluster <- as.character(plot.data$Cluster)
plot.data$Cluster[plot.data$Cluster != "Patrolling Monocytes"] <- "Other mononuclear"
x.cutoff = 1.5
y.cutoff = 1.25
plot.data[,gene.A] <- plot.data[,gene.A] + abs(rnorm(nrow(plot.data))/1000)
plot.data[,gene.B] <- plot.data[,gene.B] + abs(rnorm(nrow(plot.data))/1000)
color.by = "Cluster"
color_panel <- c("Patrolling Monocytes" = "#e08380", "Other mononuclear" = "#82be9c")
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
plot_grid(p)

# Figure S6B----
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_KO","MCD_WT"), Cluster %in% c("Patrolling Monocytes")) %>% pull(CellName)
deg <- LIMMA(
  expression_matrix = RBPJ_new@assays$RNA@data[,cells_used],
  groupid = RBPJ_new@meta.data[cells_used,"Group"]
)
deg$Sig <- FALSE
deg[deg$logFC > 0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "MCD_KO"
deg[deg$logFC < -0.25 & deg$adj.P.Val < 0.05 ,"Sig"] <- "MCD_WT"
deg$label <- c()
genes_sig <- c(deg %>% filter(Sig != FALSE, !grepl("^mt",Symbol)) %>% top_n(20, logFC) %>% pull(Symbol), deg %>% filter(Sig != FALSE) %>% top_n(20, -logFC) %>% pull(Symbol))
deg[deg$Symbol %in% genes_sig,"label"] <- deg[deg$Symbol %in% genes_sig,"Symbol"]
volcano_plot <- 
  ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val+1e-318))) +
  geom_point(aes(color = Sig), size = 1) +
  scale_color_manual(values = c(group_color_panel[c("MCD_KO","MCD_WT")], "FALSE" = "lightgrey")) +
  labs(x = "log2 fold change", y = "-log10 Adjusted P-value") +
  geom_vline(xintercept = c(-0.25,0.25), color = "grey", linetype = "dashed", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed", lwd = 1) +
  annotate("text", x = min(deg$logFC) + .1, y = 0, label = paste0(deg$Grp[nrow(deg)]," <--"), vjust = 1.5, size = 3) +
  annotate("text", x = max(deg$logFC) - .1, y = 0, label = paste0("--> ",deg$Grp[1]), vjust = 1.5, size = 3) +
  geom_text_repel(data = deg, aes(label = label), size = 4, max.overlaps = 100, segment.alpha = 0.8, force = 2, segment.size = 0.5, arrow = arrow(length = unit(0.005, 'npc'))) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none")

# Figure S6C----
cells_used <- RBPJ_new@meta.data %>% filter(Group %in% c("MCD_KO","MCD_WT"), Cluster %in% c("Patrolling Monocytes")) %>% pull(CellName)
plot.data <- data.frame(
  t(RBPJ_new@assays$RNA@data[c("Cx3cr1","Lyz1","S1pr5","Ccr2","Cd74","H2-Aa","H2-Ab1","H2-Eb1","H2-K1","Fabp4"),cells_used]),
  Cluster = RBPJ_new@meta.data[cells_used,"Cluster"],
  Group = RBPJ_new@meta.data[cells_used,"Group"]) %>% melt()
plot.list <- list()
for(i in unique(plot.data$variable)){
  plot.list[[i]] <- 
    ggplot(plot.data %>% filter(variable == i), aes(x = Group, y = value)) +
    geom_quasirandom(aes(col = Group), cex = 1.5, width = .25, alpha = 0.5) +
    geom_boxplot(aes(color = Group), outlier.size = -1, alpha = 0) +
    geom_signif(comparisons = list(c("MCD_WT","MCD_KO")),
                tip_length = 0, textsize = 3) +
    ggtitle(i) +
    labs(y = paste0("Expression of ",i)) +
    scale_color_manual(values = group_color_panel) +
    theme_cowplot(font_size = 7) +
    theme(legend.position = "null",
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()) 
}
p <- plot_grid(plotlist = plot.list, nrow = 2)

# Figure S6D----
cells_used <- RBPJ_raw@meta.data %>% filter(Diet == "MCD", Cluster %in% c("Endothelium")) %>% pull(CellName)
plot.data <- data.frame(
  t(RBPJ_raw@assays$RNA@data[c("Icam1","Icam2","Stab2","Itgb1"),cells_used]),
  Cluster = RBPJ_raw@meta.data[cells_used,"Cluster"],
  Group = RBPJ_raw@meta.data[cells_used,"Group"]) %>% melt()
plot.list <- list()
for(i in unique(plot.data$variable)){
  plot.list[[i]] <- 
    ggplot(plot.data %>% filter(variable == i), aes(x = Group, y = value)) +
    geom_quasirandom(aes(col = Group), cex = 3, width = .25, alpha = 0.5) +
    geom_boxplot(aes(color = Group), outlier.size = -1, alpha = 0) +
    geom_signif(comparisons = list(c("MCD_WT","MCD_KO")),
                tip_length = 0, textsize = 3) +
    ggtitle(i) +
    labs(y = paste0("Expression of ",i)) +
    scale_color_manual(values = group_color_panel) +
    theme_cowplot(font_size = 7) +
    theme(legend.position = "null",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.title.x = element_blank())
}
p <- plot_grid(plotlist = plot.list, nrow = 1)

# Figure S6E----
cells_used <- RBPJ_raw@meta.data %>% filter(Diet == "MCD", Cluster %in% c("Endothelium")) %>% pull(CellName)
Kit_status <- c("Neg","Pos")[as.numeric(RBPJ_raw@assays$RNA@data["Kit",cells_used] > 0) + 1]
plot.data <- data.frame(
  t(RBPJ_raw@assays$RNA@data[c("Icam1","Icam2","Esam","Pecam1","Itgb1","Stab2"),cells_used]),
  Kit = Kit_status,
  KO = RBPJ_raw@meta.data[cells_used,"KO"]) %>% melt()
plot.data <- plot.data %>% 
  mutate(Group = paste0(KO," ",Kit)) %>% 
  mutate(Group = plyr::mapvalues(
    Group, 
    from = c("KO Neg","KO Pos","WT Neg","WT Pos"),
    to = c("KO Kit-","KO Kit+","WT Kit-","WT Kit+")))%>%
  mutate(Group = factor(Group, levels = c("WT Kit+","KO Kit+","WT Kit-","KO Kit-")))
plot.list <- list()
for(i in unique(plot.data$variable)){
  plot.list[[i]] <- 
    ggplot(plot.data %>% filter(variable == i), 
           aes(x = Group, y = value)) +
    geom_quasirandom(aes(col = Group), cex = 3, width = .25, alpha = 0.5) +
    geom_boxplot(aes(color = Group), outlier.size = -1, alpha = 0) +
    geom_signif(comparisons = list(c("WT Kit+","KO Kit+"),
                                   c("WT Kit+","WT Kit-"),
                                   c("KO Kit+","KO Kit-"),
                                   c("WT Kit-","KO Kit-")),
                tip_length = 0, step_increase = .1, textsize = 3) +
    ggtitle(i) +
    labs(y = paste0("Expression of ",i)) +
    scale_color_manual(values = c("WT Kit+" = "#de827f","KO Kit+" = "#81be9c","WT Kit-" = "#e69e55","KO Kit-" = "#6895ca")) +
    theme_cowplot(font_size = 7) +
    theme(legend.position = "null",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.title.x = element_blank()) 
}
p <- plot_grid(plotlist = plot.list, nrow = 2)

# Figure S6L----
library(CellChat)
library(patchwork)
RBPJ_raw@meta.data$Cluster_new <- as.character(RBPJ_raw@meta.data$Cluster2)
metadata <- RBPJ_new@meta.data
row.names(metadata) <- metadata$Raw_CellName
RBPJ_raw@meta.data[row.names(metadata),"Cluster_new"] <- as.character(metadata$Cluster)
RBPJ_raw <- RBPJ_raw %>% subset(Cluster_new %ni% c("KC","Macrophage","Monocyte"))

cluster_used <- c("CD4 T","CD8 T","NK","B cell","Plasma B","Neutrophil","cDC1","cDC2","mregDC","Patrolling Monocytes","Ly6Chi Monocytes","C1qc+ Macs","Clec4e+ Macs","Spp1+ Macs","moKC","KCs","Epithelium","Endothelium")
for(group in c("MCD_WT","MCD_KO")){
  cat(group,"\n")
  cells_used <- RBPJ_raw@meta.data %>% filter(Cluster_new %in% cluster_used, Group == group) %>% pull(CellName)
  cellchat <- createCellChat(object = RBPJ_raw %>% subset(subset = CellName %in% cells_used), group.by = "Cluster_new")
  cellchat@DB <- CellChatDB.mouse
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 20)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat <- computeNetSimilarity(cellchat, type = "structural")
  cellchat <- netEmbedding(cellchat, type = "structural")
  cellchat <- netClustering(cellchat, type = "structural", do.parallel = F)
  saveRDS(cellchat, file = paste0("/cellchat_",group,".rds"))
}

cellchat_MCD_KO <- readRDS("./cellchat_MCD_KO.rds")
cellchat_MCD_WT <- readRDS("./cellchat_MCD_WT.rds")
MCD_object.list <- list(WT = cellchat_MCD_WT, KO = cellchat_MCD_KO)
cellchat_MCD <- mergeCellChat(MCD_object.list, add.names = names(MCD_object.list))
pairLR.use <- data.frame(interaction_name = c("THBS1_SDC4","APP_CD74","ANXA1_FPR1","CCL6_CCR2","CCL6_CCR1","CSF1_CSF1R","CXCL2_CXCR2","DLL4_NOTCH2","IL1B_IL1R2","LGALS9_HAVCR2","DLL4_NOTCH1","ITGB2_ICAM1"), stringsAsFactors = F)
netVisual_bubble(cellchat_MCD, pairLR.use = pairLR.use, sources.use = 8, targets.use = c(11,16,7,2,18,12,10), comparison = c(1, 2), angle.x = 90)