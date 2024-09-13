source("./MASLD_utils.R")

# Load Data----
Apoe_bulk <- read.csv("./WD_ApoeKO_LysMCre_Rbpjflox_bulk_tpm.csv")
Apoe_sample <- data.frame(
  row.names = colnames(Apoe_bulk)[2:16],
  SampleID = colnames(Apoe_bulk)[2:16],
  Group = as.vector(t(Apoe_bulk[1,2:16]))
)
Apoe_sample$Group <- factor(Apoe_sample$Group, levels = c("WT","Rbpj","RIN1"))
Apoe_bulk <- Apoe_bulk[-1,]
Apoe_bulk <- Apoe_bulk[!duplicated(Apoe_bulk$X),]
row.names(Apoe_bulk) <- Apoe_bulk$X
Apoe_bulk$X <- c()
Apoe_bulk[] <- lapply(Apoe_bulk, function(x) if(is.character(x)) as.numeric(x) else x)
Apoe_sample <- Apoe_sample %>% filter(SampleID %ni% c("L_8","L_15","L_85"))
Apoe_bulk <- Apoe_bulk[,Apoe_sample$SampleID]

# Figure 7I----
genes_used <- readRDS("./bulk_genes_used.rds")
inflammasome_genes <- unique(c(intersect(row.names(Apoe_bulk), genes_used %>% filter(Pathway == "GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY") %>% pull(GeneName) %>% stringr::str_to_title()),"Il1b","Tnf","Ifng","Ccl2","Ccl3","Ccl4","Ccl5","Csf3r"))
neg_inflammasome_genes <- unique(c(genes_used %>% filter(Pathway == "negativeregulationofNLRP3inflammasomecomplexassembly") %>% pull(GeneName),"Hspa8","Irgm1","Igtp","Lamp2","Trim30a","Trim31","Nlrc3","Irgm2","Mefv","Sirt2","Zdhhc12","Fbxl2","Cptp","Trem2","Csnk1a1","Trim11")) %>% stringr::str_to_title()
infla_gene <- inflammasome_genes[inflammasome_genes %ni% neg_inflammasome_genes]
Apoe_sample$Infla_Score <- colMeans(Apoe_bulk[infla_gene,])
p <- ggplot(Apoe_sample, aes(x = Group, y = Infla_Score)) +
  geom_boxplot(aes(fill = Group), alpha = .6, lwd = .3) +
  labs(y = "Inflammasome Score", x = "") +
  scale_fill_manual(values = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")) +
  ggsignif::geom_signif(comparisons = list(c("Rbpj","RIN1"),
                                           c("Rbpj","WT"),
                                           c("RIN1","WT")),
                        step_increase = .1, tip_length = 0,
                        test = "t.test", textsize = 2) +
  theme_cowplot(font_size = 7)

# Figure 7J----
fibrotic_gene <- c("Acta2","Col1a1","Col1a2","Ccn2","Tgfb1","Fgf2","Lox","Loxl2","Timp1","Timp2","Timp3","Stat1")
Apoe_sample$Fibrotic_Score <- colMeans(Apoe_bulk[fibrotic_gene,])
p <- ggplot(Apoe_sample, aes(x = Group, y = Fibrotic_Score)) +
  geom_boxplot(aes(fill = Group), alpha = .6, lwd = .3) +
  labs(y = "Fibrotic Score", x = "") +
  scale_fill_manual(values = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")) +
  ggsignif::geom_signif(comparisons = list(c("Rbpj","RIN1"),
                                           c("Rbpj","WT"),
                                           c("RIN1","WT")),
                        step_increase = .1, tip_length = 0,
                        test = "t.test", textsize = 2) +
  theme_cowplot(font_size = 7)

# Figure 7K----
chol_biosyn_gene <- intersect(row.names(Apoe_bulk), gsub("\t","",genes_used %>% filter(Pathway == "Biosynthesis") %>% pull(GeneName)) %>% convertGeneList(species = "human") %>% pull(MGI.symbol))
Apoe_sample$Chol_biosyn_Score <- colMeans(Apoe_bulk[chol_biosyn_gene,])
chol_trans_gene <- intersect(row.names(Apoe_bulk), gsub("\t","",genes_used %>% filter(Pathway == "Transport") %>% pull(GeneName)) %>% convertGeneList(species = "human") %>% pull(MGI.symbol))
Apoe_sample$Chol_trans_Score <- colMeans(Apoe_bulk[chol_trans_gene,])
chol_efflux_gene <- intersect(row.names(Apoe_bulk), gsub("\t","",genes_used %>% filter(Pathway == "Efflux") %>% pull(GeneName)) %>% convertGeneList(species = "human") %>% pull(MGI.symbol))
Apoe_sample$Chol_efflux_Score <- colMeans(Apoe_bulk[chol_efflux_gene,])
p1 <- ggplot(Apoe_sample, aes(x = Group, y = Chol_biosyn_Score)) +
  geom_boxplot(aes(fill = Group), alpha = .6, lwd = .3) +
  labs(y = "Cholesterol Biosysnthesis Score", x = "") +
  scale_fill_manual(values = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")) +
  ggsignif::geom_signif(comparisons = list(c("Rbpj","RIN1"),
                                           c("Rbpj","WT"),
                                           c("RIN1","WT")),
                        step_increase = .1, tip_length = 0,
                        test = "t.test", textsize = 2) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none")
p2 <- ggplot(Apoe_sample, aes(x = Group, y = Chol_trans_Score)) +
  geom_boxplot(aes(fill = Group), alpha = .6, lwd = .3) +
  labs(y = "Cholesterol Transport Score", x = "") +
  scale_fill_manual(values = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")) +
  ggsignif::geom_signif(comparisons = list(c("Rbpj","RIN1"),
                                           c("Rbpj","WT"),
                                           c("RIN1","WT")),
                        step_increase = .1, tip_length = 0,
                        test = "t.test", textsize = 2) +
  theme_cowplot(font_size = 7) +
  theme(legend.position = "none")
p3 <- ggplot(Apoe_sample, aes(x = Group, y = Chol_efflux_Score)) +
  geom_boxplot(aes(fill = Group), alpha = .6, lwd = .3) +
  labs(y = "Cholesterol Efflux Score", x = "") +
  scale_fill_manual(values = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")) +
  ggsignif::geom_signif(comparisons = list(c("Rbpj","RIN1"),
                                           c("Rbpj","WT"),
                                           c("RIN1","WT")),
                        step_increase = .1, tip_length = 0,
                        test = "t.test", textsize = 2) +
  theme_cowplot(font_size = 7)

hdl2ldl_gene <- intersect(row.names(Apoe_bulk), genes_used %>% filter(Pathway == "GenerationofHDLtoLDL") %>% pull(GeneName))
Apoe_sample$HDL2LDL_Score <- colMeans(Apoe_bulk[hdl2ldl_gene,])
p4 <- ggplot(Apoe_sample, aes(x = Group, y = HDL2LDL_Score)) +
  geom_boxplot(aes(fill = Group), alpha = .6, lwd = .3) +
  labs(y = "HDL to LDL Score", x = "") +
  scale_fill_manual(values = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")) +
  ggsignif::geom_signif(comparisons = list(c("Rbpj","RIN1"),
                                           c("Rbpj","WT"),
                                           c("RIN1","WT")),
                        step_increase = .1, tip_length = 0,
                        test = "t.test", textsize = 2) +
  theme_cowplot(font_size = 7)

p <- p1 + p2 + p3 + p4 + plot_layout(nrow = 1)

# Figure 7L----
genes_used <- unique(c(infla_gene, fibrotic_gene, chol_biosyn_gene, chol_trans_gene, chol_efflux_gene, hdl2ldl_gene))
mean_exp <- aggregate(t(Apoe_bulk[genes_used,]),list(Group = Apoe_sample$Group),mean)
row.names(mean_exp) <- mean_exp$Group
mean_exp$Group <- c()
matrix_for_heatmap <- t(mean_exp)
matrix_for_heatmap <- t(apply(matrix_for_heatmap,1,scale))
colnames(matrix_for_heatmap) <- row.names(mean_exp)
matrix_for_heatmap_quantile <- quantile(matrix_for_heatmap, c(0.01, 0.99))
matrix_for_heatmap <- pmax(matrix_for_heatmap, matrix_for_heatmap_quantile[1])
matrix_for_heatmap <- pmin(matrix_for_heatmap, matrix_for_heatmap_quantile[2])
color_used <- circlize::colorRamp2(seq(min(matrix_for_heatmap), max(matrix_for_heatmap), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
ha_column <-  HeatmapAnnotation(
  df = data.frame(Group = colnames(matrix_for_heatmap)),
  col = list(Group = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")))
p <- Heatmap(matrix_for_heatmap,
             cluster_rows = F, 
             cluster_columns = F,
             top_annotation = ha_column,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 1),
             col = color_used,
             name = "Exp")

# Figure S7P and Figure S7Q----
deg_list <- list()
for(clustergroup in ClustComb(unique(Apoe_sample$Group))){
  sample_used <- Apoe_sample %>% filter(Group %in% clustergroup) %>% pull(SampleID)
  deg_list[[paste0(clustergroup[1],"_",clustergroup[2])]] <- LIMMA(
    Apoe_bulk[,sample_used],
    Apoe_sample[sample_used,"Group"]
  )
}

genes_used <- c(
  deg_list[["WT_Rbpj"]] %>% group_by(Grp) %>% slice_max(order_by = abs(logFC), n = 200) %>% pull(Symbol),
  deg_list[["WT_RIN1"]] %>% group_by(Grp) %>% slice_max(order_by = abs(logFC), n = 200) %>% pull(Symbol),
  deg_list[["RIN1_Rbpj"]] %>% group_by(Grp) %>% slice_max(order_by = abs(logFC), n = 200) %>% pull(Symbol)
) %>% unique()
pca_plot <- data.frame(
  prcomp(t(Apoe_bulk[genes_used,]))$x,
  Apoe_sample
)
p <- ggplot(pca_plot, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group)) +
  geom_text_repel(aes(label = SampleID)) +
  scale_color_manual(values = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")) +
  theme_cowplot()

matrix_for_heatmap <- t(apply(Apoe_bulk[genes_used,],1,scale))
colnames(matrix_for_heatmap) <- colnames(Apoe_bulk)
matrix_for_heatmap_quantile <- quantile(matrix_for_heatmap, c(0.01, 0.99))
matrix_for_heatmap <- pmax(matrix_for_heatmap, matrix_for_heatmap_quantile[1])
matrix_for_heatmap <- pmin(matrix_for_heatmap, matrix_for_heatmap_quantile[2])
genes_labeled <- sample(genes_used,50)
color_used <- circlize::colorRamp2(seq(min(matrix_for_heatmap), max(matrix_for_heatmap), length = 8), rev(RColorBrewer::brewer.pal(11,"RdBu")[2:9]))
ha_column <-  HeatmapAnnotation(
  df = data.frame(Group = Apoe_sample$Group),
  col = list(Group = c("WT" = "#415f96", "Rbpj" = "#fb942a", "RIN1" = "#d945a0")))
p <- Heatmap(matrix_for_heatmap,
             cluster_rows = T, 
             cluster_columns = T,
             top_annotation = ha_column,
             show_row_names = T,
             column_names_gp = gpar(fontsize = 8),
             row_names_gp = gpar(fontsize = 5),
             col = color_used,
             name = "Exp") +
  rowAnnotation(link = anno_mark(
    at = which(genes_used %in% genes_labeled),
    labels = genes_used[which(genes_used %in% genes_labeled)],
    labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

# Figure S7R----
library(clusterProfiler)
library(org.Mm.eg.db)
genes_used <- list(
  WT = deg_list[["WT_RIN1"]] %>% filter(Grp == "WT", P.Value < 0.05) %>% slice_max(order_by = abs(logFC), n = 200) %>% pull(Symbol),
  RIN1 = deg_list[["WT_RIN1"]] %>% filter(Grp == "RIN1", P.Value < 0.05) %>% slice_max(order_by = abs(logFC), n = 200) %>% pull(Symbol)
)
Apoe_RIN1_pathway <- compareCluster(geneCluster = genes_used, fun = enrichGO, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL")
Apoe_RIN1_description_used <- c(
  "sterol biosynthetic process",
  "cholesterol biosynthetic process",
  "secondary alcohol biosynthetic process",
  "sister chromatid segregation",
  "nuclear division",
  "long-chain fatty acid metabolic process",
  "fatty acid metabolic process",
  "olefinic compound metabolic process",
  "xenobiotic metabolic process",
  "steroid metabolic process"
)
Apoe_RIN1_pathway_used <- Apoe_RIN1_pathway@compareClusterResult %>% filter(Description %in% Apoe_RIN1_description_used) %>% dplyr::select(Cluster, Description, GeneRatio, p.adjust)
Apoe_RIN1_pathway_used$GeneRatio <- sapply(Apoe_RIN1_pathway_used$GeneRatio, FUN = function(x){eval(parse(text = x))})
Apoe_RIN1_pathway_used$GeneRatio[Apoe_RIN1_pathway_used$Cluster == "WT"] <- -Apoe_RIN1_pathway_used$GeneRatio[Apoe_RIN1_pathway_used$Cluster == "WT"]
Apoe_RIN1_pathway_used <- Apoe_RIN1_pathway_used %>% arrange(desc(Cluster), desc(GeneRatio))
Apoe_RIN1_pathway_used$Description <- factor(Apoe_RIN1_pathway_used$Description, levels = unique(Apoe_RIN1_pathway_used$Description))
Apoe_RIN1_pathway_used$logP <- -log10(Apoe_RIN1_pathway_used$p.adjust)
Apoe_RIN1_pathway_used$logP <- pmin(Apoe_RIN1_pathway_used$logP,15)
p <- ggplot(Apoe_RIN1_pathway_used, aes(x = Description, y = GeneRatio)) +
  geom_bar(aes(fill = logP), stat = "identity") +
  scale_fill_gradientn(name = "-log10 (P-value)", limits = c(0,15), colors = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)) +
  theme_bw(base_size = 12) +
  labs(y = "Gene Ratio", x = "") +
  coord_flip() +
  theme(strip.text.x = element_blank(),strip.background = element_rect(fill=NA, color=NA))

genes_used <- list(
  WT = deg_list[["WT_Rbpj"]] %>% filter(Grp == "WT", P.Value < 0.05) %>% slice_max(order_by = abs(logFC), n = 200) %>% pull(Symbol),
  Rbpj1 = deg_list[["WT_Rbpj"]] %>% filter(Grp == "Rbpj", P.Value < 0.05) %>% slice_max(order_by = abs(logFC), n = 200) %>% pull(Symbol)
)
Apoe_Rbpj_pathway <- compareCluster(geneCluster = genes_used, fun = enrichGO, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "ALL")
Apoe_Rbpj_description_used <- c(
  "chromosome segregation",
  "mitotic nuclear division",
  "spindle checkpoint signaling",
  "extracellular matrix structural constituent",
  "collagen binding",
  "oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen",
  "monooxygenase activity",
  "heme binding",
  "steroid metabolic process",
  "fatty acid metabolic process",
  "cholesterol metabolic process"
)
Apoe_Rbpj_pathway_used <- Apoe_Rbpj_pathway@compareClusterResult %>% filter(Description %in% Apoe_Rbpj_description_used) %>% dplyr::select(Cluster, Description, GeneRatio, p.adjust)
Apoe_Rbpj_pathway_used$Description[9] <- "oxidoreductase activity"
Apoe_Rbpj_pathway_used$GeneRatio <- sapply(Apoe_Rbpj_pathway_used$GeneRatio, FUN = function(x){eval(parse(text = x))})
Apoe_Rbpj_pathway_used$GeneRatio[Apoe_Rbpj_pathway_used$Cluster == "WT"] <- -Apoe_Rbpj_pathway_used$GeneRatio[Apoe_Rbpj_pathway_used$Cluster == "WT"]
Apoe_Rbpj_pathway_used <- Apoe_Rbpj_pathway_used %>% arrange(desc(Cluster), desc(GeneRatio))
Apoe_Rbpj_pathway_used$Description <- factor(Apoe_Rbpj_pathway_used$Description, levels = unique(Apoe_Rbpj_pathway_used$Description))
Apoe_Rbpj_pathway_used$logP <- -log10(Apoe_Rbpj_pathway_used$p.adjust)
Apoe_Rbpj_pathway_used$logP <- pmin(Apoe_Rbpj_pathway_used$logP,15)
p <- ggplot(Apoe_Rbpj_pathway_used, aes(x = Description, y = GeneRatio)) +
  geom_bar(aes(fill = logP), stat = "identity") +
  scale_fill_gradientn(name = "-log10 (P-value)", limits = c(0,15), colors = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)) +
  theme_bw(base_size = 12) +
  labs(y = "Gene Ratio", x = "") +
  coord_flip() +
  theme(strip.text.x = element_blank(),strip.background = element_rect(fill=NA, color=NA))
