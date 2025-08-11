# functions for transcriptomic quality control

# library packages
library(Seurat)
library(ggplot2)
library(ggpubr)
library(Rgb)
library(dplyr)
library(DESeq2)
library(stringr)
library(limma)
library(edgeR)





# functions to filter low count cells and low count genes by Seurat

# use R package Seurat
func_basic_QC <- function(raw_count, nCount = 5e5, nFeature = 7500, output_type = "seurat_object", plot_violin = F) {
  # output_type should be "seurat_object" or "raw_count_df"
  
  # create Seurat object
  seurat_patchseq <- CreateSeuratObject(counts = raw_count, project = "basic_QC", min.cells = 100, min.features = 200)
  
  if(plot_violin) {
    png(filename = "basic_QC.png")
    VlnPlot(seurat_patchseq, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    dev.off()
  }
  
  seurat_patchseq <- subset(seurat_patchseq, subset = nFeature_RNA > nFeature & nFeature_RNA < 15000 & nCount_RNA > nCount)
  
  count_pass_basic_QC <- as.data.frame(seurat_patchseq@assays$RNA@counts)
  
  output <- NA
  if(output_type == "seurat_object") {
    output <- seurat_patchseq
  } else if (output_type == "raw_count_df") {
    output <- count_pass_basic_QC
  } else {
    print("wrong output type reauired. Should be seurat_object or raw_count_df")
  }
  
  return(output)
}

# normalization: CPM, TPM or TMM
## but TMM is kind of tricky for single-cell data because of its dependency of the reference column/sample


# function for CPM, called by the later summary function (func_normalization)
func_CPM <- function(input_seurat) {
  # called by the summary function
  input_seurat <- NormalizeData(input_seurat, normalization.method = "RC", scale.factor = 1e6)
  cpm_df <- as.data.frame(input_seurat@assays$RNA@data)
  
  return(cpm_df)
}

# function for TPM, called by the later summary normalization function
func_TPM <- function(input_seurat) {
  # transcript length data frame and col cuntion required
  ## merge transcript length with the raw count data frame
  raw_count_df <- as.data.frame(input_seurat@assays$RNA@counts)
  raw_count_len_df <- merge(raw_count_df, mm10_gene_transcript_length_df, by = "row.names", all = F, sort = F)
  rownames(raw_count_len_df) <- raw_count_len_df$Row.names
  
  raw_count_df["sample_name", ] <- colnames(raw_count_df)
  ## a function to calculate TPM for each column, would be called by apply later
  func_col_TPM <- function(length_read_col) {
    # use the raw_count_len_df generated above in this function
    i_sample_name <- length_read_col["sample_name"]
    
    # sample name is added as the last row, so delete it before calculating TPM
    i_length_read_col <- length_read_col[1:(length(length_read_col)-1)]
    i_read_over_len <- as.numeric(i_length_read_col)/raw_count_len_df$gene_transcript_length
    sum_read_over_len <- sum(i_read_over_len)
    i_TPM <- (i_read_over_len/sum_read_over_len)*10^6
    i_TPM_col <- data.frame(sample_TPM = i_TPM,
                            row.names = raw_count_len_df$gene_name)
    colnames(i_TPM_col) <- i_sample_name
    
    return(i_TPM_col)
  }
  
  tpm_df <- bind_cols(apply(raw_count_df[rownames(raw_count_df) %in% rownames(raw_count_len_df), ], 2, func_col_TPM))
  
  return(tpm_df)
}

# TMM
func_TMM <- function(input_seurat) {
  # get the matrix/data frame out of the seurat object
  count_pass_basic_QC <- as.data.frame(input_seurat@assays$RNA@counts)
  dge <- calcNormFactors(count_pass_basic_QC, method = "TMM")
  tmm <- cpm(dge)
  
  return(tmm)
}


# the summary function for normalization
func_normalization <- function(input_data, input_type = "seurat_object", method) {
  # usually, the input data should be a seurate object that passed the basic QC about nFeature and nCount
  # input_type default to be "seurat_object. There might be some time that the input is a data frame of the subset of raw counts
  # current avaliable method are CPM and TPM. TMM may be added later.
  
  # prepare the seurat_object for later normalization
  input_seurat <- NA
  if(input_type == "raw_count_df") {
    input_seurat <- CreateSeuratObject(counts = input_data, project = method, min.cells = 100, min.features = 200)
  } else if(input_type == "seurat_object") {
    input_seurat <- input_data
  } else {
    print("Error. Wrong input data type.")
    break
  }
  
  output_df <- NA
  # normalization (CPM or TPM)
  if(method == "CPM") {
    output_df <- func_CPM(input_seurat)
    
  } else if(method == "TPM") {
    output_df <- func_TPM(input_seurat)
    
  } else if(method == "TMM") {
    output_df <- func_TMM(input_seurat)
    
  } else{
    print("Error. Normalization method not abvaliable.")
    break
  }
  
  return(output_df)
  
}


# function to filter high contaminated cells based on a given threshold
func_CS_QC <- function(CS_df, cs_threshold = 0.5) {
  csqc_colname <- colnames(CS_df)[CS_df <= cs_threshold]
  return(csqc_colname)
}


# function to subset output seurat object of func_basic_QC absed on CS filtering results
func_CS_seurat <- function(seurat_object, CS_lists) {
  CS_seurat_object <- subset(seurat_object, cells = CS_lists)
  
  return(CS_seurat_object)
}






################################################################################

# a function for sample size bar plot of a single set of parameters
func_single_set_sample_size <- function(raw_count, basic_QC, basic_QC_CS, ephy_count, set_name = "default") {
  ## calculate sample size for PV and SST separately at each QC step
  ## 2 * 4 = 8 rows, each row a number of sample size
  ## 4 columns: sample_size, QC_step, cell_type and parameter_set
  
  # select pv/sst cells based on meta_patchseq_pv/sst$transcriptomics_sample_id
  ss_df <- as.data.frame(matrix(nrow = 8, ncol = 3))
  colnames(ss_df) <- c("sample_size", "QC_step", "cell_type_set_name")
  ss_df$QC_step <- rep(c("raw", "basic_QC", "basic_QC_CS", "ephy_count"),2)
  ss_df$QC_step <- factor(ss_df$QC_step, levels = c("raw", "basic_QC", "basic_QC_CS", "ephy_count") )
  ss_df$cell_type_set_name <- rep(c(paste("Pvalb", set_name, sep = "_"), paste("Sst", set_name, sep = "_")), each = 4)
  ## raw
  ss_df[1,1] <- sum(colnames(raw_count) %in% meta_patchseq_pv$transcriptomics_sample_id)
  ss_df[5,1] <- sum(colnames(raw_count) %in% meta_patchseq_sst$transcriptomics_sample_id)
  ## basic_QC
  ss_df[2,1] <- sum(colnames(basic_QC) %in% meta_patchseq_pv$transcriptomics_sample_id)
  ss_df[6,1] <- sum(colnames(basic_QC) %in% meta_patchseq_sst$transcriptomics_sample_id)
  ## basic_QC_CS
  ss_df[3,1] <- sum(colnames(basic_QC_CS) %in% meta_patchseq_pv$transcriptomics_sample_id)
  ss_df[7,1] <- sum(colnames(basic_QC_CS) %in% meta_patchseq_sst$transcriptomics_sample_id)
  ## ephy_count
  ss_df[4,1] <- sum(ephy_count$transcriptomics_sample_id %in% meta_patchseq_pv$transcriptomics_sample_id)
  ss_df[8,1] <- sum(ephy_count$transcriptomics_sample_id %in% meta_patchseq_sst$transcriptomics_sample_id)
  
  return(ss_df)
}


func_sample_size_bar <- function(ss_df) {
  bar_plot <- ggplot(ss_df, aes(fill = QC_step, y = sample_size, x = cell_type_set_name)) + 
    geom_bar(position = "dodge", stat = "identity")
  
  return(bar_plot)
}


















