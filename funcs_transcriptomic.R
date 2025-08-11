# functions for transcriptomic analysis
library(pheatmap)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(Rgb)
library(dplyr)
library(stringr)


# sample size before(always known) and after each QC step
## for each T-type (?) & for PV/SST cells -- later when we have different set of parameters
func_sample_size <- function(raw_count, basic_seurat, basic_CS_seurat, meta_pv, meta_sst, group_name = "planA") {
  # raw count: the raw count matrix (data frame), which should also be the input of func_basic_QC
  # basic_count: the output of func_basic_QC, which should be a seurat object with filtered raw count
  # basic_CS: the output of func_SC_QC, which should be a list of colnames (cell id) indicates the cells that passed the basic + CS QC
  
  ss_raw <- ncol(raw_count)
  ss_basic <- ncol(basic_seurat)
  ss_basic_cs <- ncol(basic_CS_seurat)
  
  ss_df <- data.frame(sample_size = c(ss_raw, ss_basic, ss_basic_cs), 
                      QC_step = c("raw count", "basic QC", "basic + CS"), 
                      para_group = rep(group_name, 3))
  
  return(ss_df)
}







# Differential expression analysis
# a function to apply DE analysis to the input 2 raw count matrix (for T-type, make up the input for the two T-type)
func_DE <- function(basic_seurat, CSlist_a, CSlist_b, cell_type_a, cell_type_b) {
  # cell_type_a/b: string describe the cell type (PV/SST or T-types)
  raw_count_a <- subset(basic_seurat, cells = CSlist_a)
  raw_count_a <- as.data.frame(raw_count_a@assays$RNA@counts)
  raw_count_b <- subset(basic_seurat, cells = CSlist_b)
  raw_count_b <- as.data.frame(raw_count_b@assays$RNA@counts)
  
  # merge two raw count data frame
  raw_count_ab <- merge(raw_count_a, raw_count_b, by = "row.names", all = T)
  rownames(raw_count_ab) <- raw_count_ab[, 1]
  raw_count_ab <- raw_count_ab[, -1]
  
  # generate colData to describe design
  ab_colData <- as.data.frame(matrix(ncol = 2, nrow = ncol(raw_count_ab)))
  colnames(ab_colData) <- c("sample_name", "cell_type")
  ab_colData[, "sample_name"] <- c(colnames(raw_count_a), colnames(raw_count_b))
  ab_colData[, "cell_type"] <- factor(rep(c(cell_type_a, cell_type_b), c(ncol(raw_count_a), ncol(raw_count_b))), 
                                      levels = c(cell_type_a, cell_type_b))
  rownames(ab_colData) <- ab_colData$sample_name
  
  # generate DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = raw_count_ab,
                                colData = ab_colData, 
                                design = ~ cell_type)
  dds <- DESeq(dds)
  res <- results(dds)
  
  res <- na.omit(res)
  # res <- res[which(res$padj <= 0.05), ] # fold change cut off not applied here for later volcano plot
  
  return(res)
  
}


# for the given normalized count matrix and gene list, generate heatmap (unclustered) with given color labels (? not now)...
# func_heatmap <- function(gene_name_list, gene_type_list, norm_count_df) {
#   
# }


# volcano plot & label some cell type marker genes to indicate the performance of DE analysis
func_volcano <- function(DE_output, padj_cutoff = 0.05, abs_logFC_cutoff = 0.5) {
  # input is the output of func_DE, a large DESeqResults
  DE_df <- as.data.frame(DE_output)
  # add 'enriched in PV', 'enriched in SST', 'not DE'
  DE_df$enriched_in <- as.factor(ifelse(DE_df$padj<padj_cutoff & abs(DE_df$log2FoldChange)>abs_logFC_cutoff, 
                                        ifelse(DE_df$log2FoldChange > abs_logFC_cutoff, "enriched in SST", "enriched in PV"), "not DE"))
  
  vol_plot <- ggplot(data = DE_df, aes(x = log2FoldChange, y = -log10(padj), color = enriched_in)) + 
    geom_point(alpha = 0.8, size = 1) + theme_bw(base_size = 15) +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + 
    scale_color_manual(name = "", values = c("red", "green", "black"), limits = c("enriched in SST", "enriched in PV", "not DE"))
    
  
  return(vol_plot)
  
}






# merge normalized data with ephy feature, then set threshold to remove low exprssion genes
# a function to transpose cpm matrix for latter table merging ##################
t_cpm_matrix <- function(cpm_matrix, cell_col_name) {
  t_cpm <- t(cpm_matrix)
  t_cpm <- cbind(rownames(t_cpm), t_cpm)
  colnames(t_cpm)[1] <- cell_col_name
  t_cpm <- as.data.frame(t_cpm, stringsAsFactors = F)
  ncol_t <- ncol(t_cpm)
  t_cpm[, 2:ncol_t] <- as.data.frame(lapply(t_cpm[, 2:ncol_t], as.numeric))
  return(t_cpm)
}


func_ephy_Ttype_logtNorm <- function(ephy_Ttype_df, Norm_df, merge_name = "transcriptomics_sample_id") {
  # substract cells to make the transpose easier
  # directly get the intersection of two dfs' cells
  cells_intersection <- intersect(colnames(Norm_df), ephy_Ttype_df[, merge_name])
  # Norm_df_sub <- Norm_df[, colnames(Norm_df) %in% ephy_Ttype_df[, merge_name]]
  Norm_df_sub <- Norm_df[, cells_intersection]
  
  # t_cpm_matrix required
  tNorm_df <- t_cpm_matrix(Norm_df_sub, cell_col_name = merge_name)
  logtNorm_df <- tNorm_df
  logtNorm_df[, 2:ncol(logtNorm_df)] <- log2(logtNorm_df[, 2:ncol(logtNorm_df)] + 1)
  
  
  merged_df <- merge(ephy_Ttype_df, logtNorm_df, by = merge_name, all = F)
  
  return(merged_df)
  
}

# a function to select genes based on the gene's count of non-zero-cell. quantile_cutoff is the lower bound
func_non_zero_quantile <- function(ephy_log_cpm_df, ncol_ephy = 6, quantile_cutoff = 0.1) {
  # ncol_ephy describeds the first certain number of columns that are not genes
  # for merged single cell ephy and logcpm, ncol_ephy = 6
  log_cpm_df <- ephy_log_cpm_df[, -c(1:ncol_ephy)]
  num_cutoff <- as.numeric(quantile(colSums(log_cpm_df > 0), quantile_cutoff))
  log_cpm_df <- log_cpm_df[, colSums(log_cpm_df > 0) >= num_cutoff]
  
  filtered_ephy_log_cpm_df <- cbind(ephy_log_cpm_df[, 1:ncol_ephy], log_cpm_df)
  return(filtered_ephy_log_cpm_df)
  
}

# an extra function to select genes that have good expression level in all cell types
func_each_type_non_zero_quantile <- function(ephy_log_cpm_df, ncol_ephy = 6, quantile_cutoff = 0.1, type_col_name = "cell_type") {
  # ncol_ephy describeds the first certain number of columns that are not genes
  # for merged single cell ephy and logcpm, ncol_ephy = 6
  type_list <- unique(ephy_log_cpm_df[, type_col_name])
  num_type <- length(type_list)
  
  # put the first cell type out of the for loop so that we could have the initial cpm df for later rbind
  log_cpm_df <- ephy_log_cpm_df[which(ephy_log_cpm_df[, type_col_name] == type_list[1]),  -c(1:ncol_ephy)]
  num_cutoff <- as.numeric(quantile(colSums(log_cpm_df > 0), quantile_cutoff))
  log_cpm_df <- log_cpm_df[, colSums(log_cpm_df > 0) >= num_cutoff]
  gene_list <- colnames(log_cpm_df)
  
  # for loop to append the following cell types' filtered cpm df
  for (i in 2:num_type) {
    i_log_cpm_df <- ephy_log_cpm_df[which(ephy_log_cpm_df[, type_col_name] == type_list[i]),  -c(1:ncol_ephy)]
    i_num_cutoff <- as.numeric(quantile(colSums(i_log_cpm_df > 0), quantile_cutoff))
    i_log_cpm_df <- i_log_cpm_df[, colSums(i_log_cpm_df > 0) >= i_num_cutoff]
    # get intersect of colnames of cpm df so that all remained genes will not have low expression level in either cell type
    gene_list <- intersect(gene_list, colnames(i_log_cpm_df))
  }
  
  filtered_log_cpm_df <- ephy_log_cpm_df[, gene_list]
  
  # log_cpm_df <- ephy_log_cpm_df[, -c(1:ncol_ephy)]
  # num_cutoff <- as.numeric(quantile(colSums(log_cpm_df > 0), quantile_cutoff))
  # log_cpm_df <- log_cpm_df[, colSums(log_cpm_df > 0) >= num_cutoff]
  
  filtered_ephy_log_cpm_df <- cbind(ephy_log_cpm_df[, 1:ncol_ephy], filtered_log_cpm_df)
  return(filtered_ephy_log_cpm_df)
  
}








# T-type aggregate for trimmed mean, then merge with sample size for each T-type for later weighted lm
func_trim_mean <- function(x) {
  return(mean(x, trim = 0.25))
}

# recover single-cell log(CPM+1) to CPM, then trim mean, then log norm
func_recover_trim_mean <- function(x) {
  # first recover x from log2(CPM+1) to CPM
  x_cpm <- 2^x-1
  # trim mean of CPM
  x_trim_cpm <- mean(x_cpm, trim = 0.25)
  # log norm
  x_log_of_trim_cpm <- log2(x_trim_cpm+1)
  return(x_log_of_trim_cpm)
  # return(mean(x, trim = 0.25))
}

# for ephy columns only
func_Ttype_trim_mean <- function(AIHW_Ttype_log, remove_col_index = c(1, 4:6), group_col_name = "T-type_Label") {
  Ttype_trim_ephy_log_df <- aggregate(AIHW_Ttype_log[, -remove_col_index], list(AIHW_Ttype_log[, group_col_name]), func_trim_mean)
  colnames(Ttype_trim_ephy_log_df)[1] <- "T-type_Label"
  # T-type sample size
  Ttype_sample_size <- aggregate(AIHW_Ttype_log[, group_col_name], list(AIHW_Ttype_log[, group_col_name]), length)
  colnames(Ttype_sample_size) <- c("T-type_Label", "sample_size")
  # merge sample size, set a column of cell_type (not T-type)
  Ttype_num_trim_ephy_log_df <- merge(Ttype_sample_size, Ttype_trim_ephy_log_df, by = "T-type_Label")
  rownames(Ttype_num_trim_ephy_log_df) <- Ttype_num_trim_ephy_log_df$`T-type_Label`
  Ttype_num_trim_ephy_log_df[, 1] <- str_split_fixed(Ttype_num_trim_ephy_log_df$`T-type_Label`, " ", n = 2)[, 1]
  colnames(Ttype_num_trim_ephy_log_df)[1] <- "cell_type"
  
  return(Ttype_num_trim_ephy_log_df)
}

# for single-cell logCPM part of the data, remove_col_index should include ephy columns
func_Ttype_recover_trim_mean <- function(AIHW_Ttype_log, remove_col_index = c(1, 4:6), group_col_name = "T-type_Label") {
  Ttype_trim_ephy_log_df <- aggregate(AIHW_Ttype_log[, -remove_col_index], list(AIHW_Ttype_log[, group_col_name]), func_recover_trim_mean)
  colnames(Ttype_trim_ephy_log_df)[1] <- "T-type_Label"
  # T-type sample size
  Ttype_sample_size <- aggregate(AIHW_Ttype_log[, group_col_name], list(AIHW_Ttype_log[, group_col_name]), length)
  colnames(Ttype_sample_size) <- c("T-type_Label", "sample_size")
  # merge sample size, set a column of cell_type (not T-type)
  Ttype_num_trim_ephy_log_df <- merge(Ttype_sample_size, Ttype_trim_ephy_log_df, by = "T-type_Label")
  rownames(Ttype_num_trim_ephy_log_df) <- Ttype_num_trim_ephy_log_df$`T-type_Label`
  Ttype_num_trim_ephy_log_df[, 1] <- str_split_fixed(Ttype_num_trim_ephy_log_df$`T-type_Label`, " ", n = 2)[, 1]
  colnames(Ttype_num_trim_ephy_log_df)[1] <- "cell_type"
  
  return(Ttype_num_trim_ephy_log_df)
}



func_trimed_df_select_genes <- function(Ttype_num_trim_ephy_log_df, gene_list, kept_col_names = c("cell_type", "sample_size", "AI", "HW")) {
  # the select genes may be filtered out during the previous quality control steps
  # if so, print the warning and return the subset data frame with the remained genes
  # also, delete genes whose trimmed mean of all T-types are 0, and print their names
  remained_gene_list <- gene_list[gene_list %in% colnames(Ttype_num_trim_ephy_log_df)]
  select_Ttype_num_trim_ephy_log_df <- Ttype_num_trim_ephy_log_df[, remained_gene_list]
  # remove genes whose all trimmed mean are 0
  remained_Ttype_num_trim_ephy_log_df <- select_Ttype_num_trim_ephy_log_df[, colSums(select_Ttype_num_trim_ephy_log_df)>0]
  
  if(length(remained_gene_list) != length(gene_list)) {
    print("Warning! Following genes are not available: ")
    print(gene_list[!(gene_list %in% colnames(Ttype_num_trim_ephy_log_df))])
  }
  
  if(ncol(remained_Ttype_num_trim_ephy_log_df) != ncol(select_Ttype_num_trim_ephy_log_df)) {
    print("Warning! Folloing genes are removed due to their 0 trimmed mean of all T-types: ")
    print(colnames(select_Ttype_num_trim_ephy_log_df)[colSums(select_Ttype_num_trim_ephy_log_df)<=0])
  }
  
  remained_Ttype_num_trim_ephy_log_df <- cbind(Ttype_num_trim_ephy_log_df[, kept_col_names], 
                                               remained_Ttype_num_trim_ephy_log_df)
  
  return(remained_Ttype_num_trim_ephy_log_df)
}


func_filtered_scdf_select_genes <- function(filtered_ephy_Ttype_scdf, gene_list) {
  # select genes from the single-cell (sc) data frame (low expression genes already filtered)
  # some genes in the given list may be filtered out during the previous quality control steps
  # if so, print the warning and return the subset data frame with the remained genes
  remained_gene_list <- gene_list[gene_list %in% colnames(filtered_ephy_Ttype_scdf)]
  
  if(length(remained_gene_list) != length(gene_list)) {
    print("Warning! Following genes are nor available: ")
    print(gene_list[!(gene_list %in% colnames(filtered_ephy_Ttype_scdf))])
  }
  
  remained_filtered_scdf <- filtered_ephy_Ttype_scdf[, c("transcriptomics_sample_id", "AI", "HW", "cell_type", "T-type_Label", "structure", remained_gene_list)]
  return(remained_filtered_scdf)
  
  
}












func_Ttype_weighted_lm_plot <- function(Ttype_num_med_df, gene_name, ephy = "HW/AI", save_plot = F, plot_title_annotation = "...") {
  # plot the scatter plot and fitted line, return the p-value, R-squred of lm result
  w_df <- Ttype_num_med_df[, c("cell_type", "sample_size", ephy, gene_name)]
  
  # lm for the given ephy ~ gene
  col_ephy <- w_df[,ephy]
  col_gene <- w_df[,gene_name]
  col_ss <- w_df$sample_size
  
  # lm_formula <- as.formula(paste(ephy, "~", gene_name, sep = ""))
  Ttype_w_lm <- lm(col_ephy ~ col_gene, weights = col_ss)
  sum_Ttype_w_lm <- summary(Ttype_w_lm)
  # print(sum_weighted_lm)
  
  
  lm_df <- as.data.frame(matrix(nrow = 1, ncol = 3))
  colnames(lm_df) <- c("p_value", "R_squared", "slope")
  rownames(lm_df) <- gene_name
  lm_df[1, ] <- c(sum_Ttype_w_lm$coefficients[2,4], sum_Ttype_w_lm$r.squared, sum_Ttype_w_lm$coefficients[2,1])
  
  if (save_plot) {
    w_df$i_pred <- predict(Ttype_w_lm)
    lm_ggplot <- ggplot(data = w_df, aes(x = .data[[gene_name]], y = .data[[ephy]])) + 
    geom_point(aes(color = cell_type)) + 
    geom_line(aes(y = i_pred), colour = "#000000") + 
    ggtitle(paste(ephy, "~", gene_name, ", ", plot_title_annotation, sep = ""), 
            subtitle = paste("p-value=", signif(lm_df[1,1],4), "  R.squared=", signif(lm_df[1,2],4), " slope=", signif(lm_df[1,3],4), sep = ""))
    ggsave(filename = paste("Output_Data/Ttype_weighted/", ephy, "_", gene_name, ".png", sep = ""), plot = lm_ggplot, width = 6, height = 4)
  }
  
  return(lm_df)
  
}




# a function to wlm all input genes in a df/sub-df and save the result of wlm for latter bubble plot
func_genes_wlm_df <- function(sub_Ttype_num_trim_ephy_logcpm_df, ephy_name = "AI") {
  # the first 4 columns shoud be: cell_type, sample_size, AI, HW
  # the 5th to the last columns are genes
  if (ephy_name == "AI") { # include HW in the 4th column
    gene_list <- colnames(sub_Ttype_num_trim_ephy_logcpm_df)[4:ncol(sub_Ttype_num_trim_ephy_logcpm_df)]
  } else {
    gene_list <- colnames(sub_Ttype_num_trim_ephy_logcpm_df)[5:ncol(sub_Ttype_num_trim_ephy_logcpm_df)]
  }
  
  # generate an empty data frame, each row for a gene's wlm results
  wlm_df <- as.data.frame(matrix(nrow = length(gene_list), ncol = 5))
  colnames(wlm_df) <- c("gene_name", "p_value", "R_squared", "slope", "n_log_p")
  rownames(wlm_df) <- gene_list
  wlm_df$gene_name <- gene_list
  
  # for loop to get p_value, R_squared and slope
  for (i in 1:length(gene_list)) {
    i_gene_name <- gene_list[i]
    i_ephy <- sub_Ttype_num_trim_ephy_logcpm_df[, ephy_name]
    i_trim_cpm <- sub_Ttype_num_trim_ephy_logcpm_df[, i_gene_name]
    
    # if sum(i_trim_cpm) > 0, run the wlm; else, skip to next one
    if(sum(i_trim_cpm) > 0) {
      i_wlm <- lm(i_ephy ~ i_trim_cpm, weights = sub_Ttype_num_trim_ephy_logcpm_df$sample_size)
      i_sum_wlm <- summary(i_wlm)
      wlm_df[i, 2:4] <- c(i_sum_wlm$coefficients[2,4], i_sum_wlm$r.squared, i_sum_wlm$coefficients[2,1])
      wlm_df[i, 5] <- -log10(wlm_df[i,2])
    } else {
      print(paste(i_gene_name, "ZERO trimmed cpm in total, skip"))
    }
    
  }
  
  return(wlm_df)
  
}




















