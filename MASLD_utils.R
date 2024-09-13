library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(grid)
library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
library(patchwork)

# Color panel----
Binomial_color_panel <- c("TRUE" = "#E64B35", "FALSE" = "lightgrey")
group_color_panel <- c("ND_WT" = "#F0CE39", "ND_KO" = "#9E6BAB", "MCD_WT" = "#EC7000", "MCD_KO" = "#019FAA", "SD" = "#9E6BAB", "WD" = "#EC7000")
cluster_color_panel <- c(
  "Ly6Chi Monocytes" = "#C45240", "Patrolling Monocytes" = "#F2C987",
  "Clec4e+ Macs" = "#63B571", "C1qc+ Macs" = "#F8C3B4",
  "Capsule Macs" = "#86C2D7",
  "Spp1+ Macs" = "#C86EA8",
  "moKC" = "#F0CE39", "KCs" = "#60B4F3",
  "Prolif. KCs" = "#F0C1ED"
)
raw_cluster_color_panel <- c(
  "cDC1" = "#F0BED7", "cDC1s" = "#F0BED7",
  "cDC2" = "#EEE1E8", "cDC2s" = "#EEE1E8",
  "mregDC" = "#FC7542", "Mig cDCs" = "#FC7542", "cDC(pro)" = "#E69E8E", 
  "Ly6c+ Mono" = "#C45240", "Ly6Chi monocytes" = "#C45240", "Monocytes" = "#C45240", "Mono" = "#C45240",
  "Ly6c- Mono" = "#F2C987", "Patrolling Monocytes" = "#F2C987", "Patrolling\nmonocytes" = "#F2C987", "Pat.Mono" = "#F2C987",
  "Trans Monocytes 1" = "#2D75B7", "Trans Monocytes 2" = "#007FB7", "Transitioning monocytes" = "#2D75B7", "MoMac1" = "#007FB7",
  "Clec4e+ Macro" = "#63B571", "immLAMs" = "#63B571", "C1qc+ Macro" = "#F8C3B4", 
  "Capsule macs" = "#86C2D7", "CV and Capsule" = "#86C2D7",
  "Spp1+ Macro" = "#C86EA8", "LAMs" = "#C86EA8", "mac1" = "#C86EA8", "matLAMs" = "#C86EA8",
  "moKCs" = "#F0CE39", "Pre-MoKC and MoKC" = "#F0CE39", "Pre-moKCs and moKCs" = "#F0CE39",
  "KCs" = "#6E8AD4", "Resident KCs" = "#6E8AD4", "resKCs" = "#6E8AD4",
  "Mertk- KC" = "#60B4F3", "Mertk+ KC" = "#6E8AD4", "Lyve1+ KC" = "#F17591",
  "KC(pro)_S" = "#F0C1ED", "KC(pro)_G2M" = "#7371B0", "Prolif. macs" = "#F0C1ED",
  "Peritoneal macs" = "#D4B696"
)

# >>Not include in----
`%ni%` <- Negate(`%in%`)

# >>Add labels----
calcCenters <- function(cluster, reduction, filter = T, filter.num = 10) {
  df <- data.frame(Cluster = as.factor(cluster), Dim1 = reduction[,1], Dim2 = reduction[,2])
  centers <- df %>%
    group_by(Cluster) %>%
    summarise(mean_x = median(Dim1),
              mean_y = median(Dim2))
  if(filter){
    cluster_used <- data.frame(table(cluster)) %>% filter(Freq > filter.num) %>% pull(cluster) %>% as.character()
    centers <- centers %>% filter(Cluster %in% cluster_used)
  }
  return(centers)
  
}

addLabels <- function(centers, label_size = 3, label_short = FALSE) {
  if (label_short) centers <- suppressWarnings(
    tidyr::separate(centers, Cluster, into = c("Cluster", "Cluster_long"), extra = "drop"))
  ggrepel::geom_text_repel(data = centers,
                           aes(x = mean_x, y = mean_y),
                           label = centers$Cluster,
                           size = label_size,
                           alpha = 0.8,
                           segment.alpha = 0.8,
                           force = 2,
                           segment.size = 0.5,
                           arrow = arrow(length = unit(0.01, 'npc')))
}

# >>LIMMA----
LIMMA <- function(expression_matrix, groupid, doAUC = FALSE) {
  ## Differential expressed genes in two groups by Limma.
  ##
  ## Args:
  #' @expression_matrix: Gene*cell matrix.
  #' @groupid: The groupid of each cell, there should be only two groups.
  ##
  ## Returns:
  ## A dataframe with the output of limma and expression percentage.
  library(limma)
  library(ROCR)
  expression_matrix <- as.matrix(expression_matrix)
  groupid <- as.character(groupid)
  groupid_raw <- unique(groupid)
  groupid[groupid == groupid_raw[1]] <- "GroupA"
  groupid[groupid == groupid_raw[2]] <- "GroupB"
  contrast <<- paste0(levels(factor(groupid)), collapse = "-")
  design <- model.matrix( ~ 0 + factor(groupid))
  colnames(design) <- levels(factor(groupid))
  rownames(design) <-
    colnames(expression_matrix)  # design data used in limma
  contrast.matrix <- makeContrasts(contrast, levels = design)
  fit <- lmFit(expression_matrix, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput <-  topTable(fit2, coef = 1, n = Inf)
  nrDEG <-  na.omit(tempOutput)
  nrDEG$Symbol <- row.names(nrDEG)
  positive_group <-
    row.names(fit2$contrasts)[fit2$contrasts == 1]  # high expression when logFC > 0
  negative_group <-
    row.names(fit2$contrasts)[fit2$contrasts == -1]  # low expression when logFC < 0
  nrDEG$Grp <-
    c(negative_group, positive_group)[as.numeric(nrDEG$logFC > 0) + 1]
  nrDEG$Grp[nrDEG$Grp == "GroupA"] <- groupid_raw[1]
  nrDEG$Grp[nrDEG$Grp == "GroupB"] <- groupid_raw[2]
  cell.Grp1 <- which(groupid == levels(as.factor(groupid))[1])
  cell.Grp2 <- which(groupid == levels(as.factor(groupid))[2])
  Exp.Mean.Grp1 <-
    rowMeans(expression_matrix[nrDEG$Symbol, cell.Grp1])  # mean expression in the group
  Exp.Mean.Grp2 <-
    rowMeans(expression_matrix[nrDEG$Symbol, cell.Grp2])  # mean expression out the group
  Exp.Per.Grp1 <-
    apply(expression_matrix[nrDEG$Symbol, cell.Grp1], 1, Expression.Per)  # expression percentage in the group
  Exp.Per.Grp2 <-
    apply(expression_matrix[nrDEG$Symbol, cell.Grp2], 1, Expression.Per)  # expression percentage out the group
  de.genes.all <- cbind(nrDEG, Exp.Mean.Grp1, Exp.Per.Grp1, Exp.Mean.Grp2, Exp.Per.Grp2)
  colnames(de.genes.all)[9:12] <-
    c(paste0(c("Exp.Mean.", "Exp.Per."), groupid_raw[1]),
      paste0(c("Exp.Mean.", "Exp.Per."), groupid_raw[2]))
  if (!is.na(de.genes.all[1, 1])) {
    if(doAUC){
      for (k in 1:nrow(de.genes.all)) {
        category <- as.numeric(groupid_new == de.genes.all$Grp[k])
        pred <-
          prediction(expression_matrix[de.genes.all$Symbol[k], ], category)
        pauc <- performance(pred, measure = "auc")
        de.genes.all[k, "AUC"] <- pauc@y.values[[1]]
      }  # calculate the AUC for each gene
    }
    de.genes.all <- de.genes.all %>% arrange(desc(logFC))
    de.genes.all
    return(de.genes.all)
  } else{
    print("No significant genes!")
  }
  
}

Expression.Per <- function(x, cutoff = 0.1) {
  # percent of gene-expressed cell, cutoff = 0.1
  return(sum(x > cutoff) / length(x))
}

# >>BatchEntropy----
BatchEntropy <- function(input.data, group.id, k.used = 30, dimension.used = "tSNE") {
  ## Calculate the cell entropy according to the cell group.
  ##
  ## Args:
  #' @input.data: A matrix with each cell in a row. The column could be
  #' tSNE coordinate, expression matrix or even a dist format.
  #' @group.id: A vector with the same length as the row of input.data, or
  #' a list of vectors.
  #' @k.used: The k used to build a kNN-graph.
  #' @dimension.used: The method to reduce the dimension, tSNE by default,
  #' could also be PCA or raw.
  ##
  ## Returns:
  ## A vector with each cell's entropy.
  library(dbscan)
  library(entropy)
  if (dimension.used == "raw") {
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = input.data, k = k.used, sort = FALSE, search = "dist")
  }
  if (dimension.used == "tSNE") {
    cat(paste("Reduce the dimension using", dimension.used, "at", Sys.time(), "\n"))
    tSNE.coor <- Rtsne::Rtsne(input.data)
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = tSNE.coor$Y, k = k.used, sort = FALSE, search = "dist")
  }
  if (dimension.used == "PCA") {
    cat(paste("Reduce the dimension using", dimension.used, "at", Sys.time(), "\n"))
    PCA.coor <- prcomp(input.data, center = FALSE)
    PCA.cumsd <- cumsum(PCA.coor$sdev) / sum(PCA.coor$sdev)
    nPCs.used <- which(PCA.cumsd > 0.9)[1]
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = PCA.coor$x[, 1:nPCs.used], k = k.used, sort = FALSE, search = "dist")
  }
  if (!is.list(group.id)) {
    group.id <- list(default = group.id)
  }
  cell_entropy <- list()
  for (i in names(group.id)) {
    knn_group <- matrix(group.id[[i]][knn_graph$id],
                        nrow = nrow(input.data),
                        byrow = FALSE)
    row.names(knn_group) <- row.names(input.data)
    colnames(knn_group) <- 1:k.used
    cat(paste("Calculate the cell entropy of", i, "at", Sys.time(), "\n"))
    cell_entropy[[i]] <- apply(knn_group, 1, function(x) {
      entropy(table(x))
    })
  }
  return(cell_entropy)
}

# >>Confusion heatmap----
Confusion_heatmap_new <- function (ori, prd, color = NULL) 
{
  cross.validation.filt <- tibble(ori = ori, prd = prd) %>% 
    dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
  cross.validation.filt[is.na(cross.validation.filt)] = 0
  cross.validation.filt[, -1] <- round(cross.validation.filt[, 
                                                             -1]/rowSums(cross.validation.filt[, -1]), 2)
  cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", 
                                                                   value = "Prob", -ori)
  p <- cross.validation.filt %>% ggplot(aes(ori, prd, fill = Prob)) + 
    geom_tile() + theme(axis.title = element_text(size = 0)) + 
    theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
    theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                         panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                         axis.ticks = element_blank(), axis.title = element_blank()) + 
    theme(axis.text.y = element_text(color = "black"), 
          axis.text.x = element_text(color = "black", 
                                     angle = 45, hjust = 1))
  if(is.null(color)){
    p <- p + scale_fill_viridis()
  }else{
    p <- p + scale_fill_gradientn(colors = color)
  }
  return(p)
}

# Chi-square test plot---
CrossTabPlot <- function(crosstab){
  pvalue <- fisher.test(crosstab)
  ct <- crosstab %>% data.frame()
  ct$Prop <- round(ct$Freq/sum(ct$Freq) * 100,2)
  item1 <- unique(ct[,1])[1]
  item2 <- unique(ct[,2])[1]
  ct[(ct[,1] == item1) == (ct[,2] == item2),"Group"] <- "Group1"
  ct[(ct[,1] == item1) != (ct[,2] == item2),"Group"] <- "Group2"
  p <- ggplot(ct, aes_string(x = colnames(ct)[1], y = colnames(ct)[2])) +
    geom_tile(aes(fill = Group)) +
    ggtitle(paste0("p-value = ",format(pvalue$p.value,digits=3))) +
    geom_text(aes(label = paste0(Freq,"\n(",Prop,"%)"))) +
    scale_fill_manual(values = c("Group1" = "#9bc0cc", "Group2" = "#e1e2e3")) +
    theme_minimal() +
    theme(legend.position = "none")
  return(p)
}

# >>RO/E----
ROIE <- function(crosstab, filter = NULL){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  if(is.null(filter)){filter = 10}
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  rowsum.matrix[rowsum.matrix <= filter] <- 0
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  colsum.matrix[colsum.matrix <= filter] <- 0
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  roie <- roie[rowSums(roie)!=0,]
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        if(m2[i,j] == 0){
          div.result[i,j] <- 0
        }else{
          div.result[i,j] <- m1[i,j] / m2[i,j]
        }
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

# >>ROIE_plot----
ROIE_plot <- function(ROIE_table, max = 2.5, font.size = 12, text.size = 4){
  ROIE_table <- melt(ROIE_table)
  ROIE_table$text <- round(ROIE_table$value,2)
  ROIE_table$value[ROIE_table$value > max] <- max
  ROIE_table$value[ROIE_table$value == 0] <- NA
  ROIE_table$text[ROIE_table$text == 0] <- "-"
  p <- 
    ggplot(data = ROIE_table, aes(Var2, Var1, fill = value)) +
    geom_tile() + 
    scale_fill_gradientn(name = "Ro/e", colours = colorRampPalette(brewer.pal(4, "YlOrRd"))(100), na.value = "lightgrey") +
    geom_text(aes(label = text), size = text.size) +
    labs(x = "", y = "") +
    theme_cowplot(font_size = font.size) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          axis.line = element_blank(), 
          axis.ticks = element_blank()) +
    guides(fill = guide_colourbar(barwidth = .5, barheight = 6))
  return(p)
}

## Compare clusters in different datasets----
Compare2Clusters <- function(exp1, cluster1, exp2, cluster2, ngenes = 1000, upper_gene = F, homo_gene = F, name1 = "D1", name2 = "D2"){
  if(homo_gene == T){
    gene_homo <- read.csv("/work/lzy/project/utils/homomuris.csv", row.names = 3, stringsAsFactors = F)
    exp1_gene <- row.names(exp1)
    exp2_gene <- row.names(exp2)
    exp1_gene[exp1_gene %in% row.names(gene_homo)] <- gene_homo[exp1_gene[exp1_gene %in% row.names(gene_homo)] ,"symbol"]
    exp2_gene[exp2_gene %in% row.names(gene_homo)] <- gene_homo[exp2_gene[exp2_gene %in% row.names(gene_homo)] ,"symbol"]
    exp1 <- SetRowNames(exp1,exp1_gene)
    exp2 <- SetRowNames(exp2,exp2_gene)
  }
  if(homo_gene == T | upper_gene == T){
    exp1 <- SetRowNames(exp1,toupper(row.names(exp1)))
    exp2 <- SetRowNames(exp2,toupper(row.names(exp2)))
  }
  genes_used <- intersect(row.names(exp1), row.names(exp2))
  exp1_avg <- aggregate(t(exp1[genes_used,]), list(Cluster = cluster1), mean)
  row.names(exp1_avg) <- exp1_avg$Cluster
  exp1_avg$Cluster <- c()
  exp1_avg <- t(exp1_avg)
  colnames(exp1_avg) <- paste0(name1,"_",colnames(exp1_avg))
  exp2_avg <- aggregate(t(exp2[genes_used,]), list(Cluster = cluster2), mean)
  row.names(exp2_avg) <- exp2_avg$Cluster
  exp2_avg$Cluster <- c()
  exp2_avg <- t(exp2_avg)
  colnames(exp2_avg) <- paste0(name2,"_",colnames(exp2_avg))
  merged_exp <- cbind(exp1_avg, exp2_avg)
  merged_exp <- merged_exp[rowSums(merged_exp) > 0,]
  merged_exp_norm <- preprocessCore::normalize.quantiles(as.matrix(merged_exp))
  row.names(merged_exp_norm) <- row.names(merged_exp)
  colnames(merged_exp_norm) <- c(colnames(exp1_avg),colnames(exp2_avg))
  merged_exp_norm <- sva::ComBat(as.matrix(merged_exp_norm), c(rep(1,ncol(exp1_avg)), rep(2,ncol(exp2_avg))), mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  if(ngenes > nrow(merged_exp_norm)){
    genes_used <- row.names(merged_exp_norm)
  }else{
    genes_used <- names(sort(apply(merged_exp_norm,1,sd),decreasing = T))[1:ngenes] 
  }
  pca_analysis <- prcomp(t(merged_exp_norm[genes_used,]))
  pca_analysis <- summary(pca_analysis)
  fit <- pvclust::pvclust(t(pca_analysis$x), method.hclust = "ward.D", method.dist = "euclidean")
  return(list(
    merged_exp = merged_exp_norm, 
    genes_used = genes_used, 
    pca_analysis = pca_analysis, 
    hc_fit = fit)
  )
}
