# Load data
mydata <- read.csv("processed_data.csv", stringsAsFactors = FALSE)
head(mydata)
ref = read.csv("reference.csv")



#############################
#Add ref_values, not exatcly same with reference in RG processing
#############################

# Extract replicate numbers from column names
rep_ids <- as.numeric(sub(".*_(\\d+)$", "\\1", colnames(mydata)[-1]))

# Compute the reference values
ref_values <- sapply(sort(unique(rep_ids)), function(rep) {
  cols <- which(rep_ids == rep)
  tmp <- rowMeans(mydata[,-1][, cols, drop = FALSE], na.rm = TRUE)
  round(tmp, 0)
})

colnames(ref) <- c("Gene", paste0("ref_value_", sort(unique(rep_ids))))

# Combine with the reference values
mydata <- cbind(mydata, ref_values)

all.equal(ref[,2:5], ref_values)


# DESeq2
# install.packages("BiocManager")
# BiocManager::install('DESeq2')

# Load necessary libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)

# Set rownames to gene names
rownames(mydata) <- mydata$Gene
mydata$Gene <- NULL

# Identify reference columns
all_columns <- colnames(mydata)
ref_cols <- grep("^ref_value", all_columns, value = TRUE)

# Identify tissue columns
tissue_cols <- setdiff(all_columns, ref_cols)

# Extract tissue names
getTissueName <- function(colname) {
  sub("_[0-9]+$", "", colname)
}
tissue_names <- sapply(tissue_cols, getTissueName)
unique_tissues <- unique(tissue_names)


# Loop over each tissue and perform DESeq2 analysis against the reference
for (tissue in unique_tissues) {
  
  # Subset data
  tissue_subset <- mydata[, which(tissue_names == tissue), drop = FALSE]
  ref_subset <- mydata[, ref_cols, drop = FALSE]
  counts_subset <- cbind(tissue_subset, ref_subset)
  
  #counts_subset[is.na(counts_subset)] <- 0
  
  # Condition vector, tissue name vs "ref"
  condition <- factor(c(rep(tissue, ncol(tissue_subset)),
                        rep("ref", ncol(ref_subset))))
  coldata <- data.frame(condition = condition)
  rownames(coldata) <- colnames(counts_subset)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                                colData = coldata,
                                design = ~ condition)
  
  # Remove rows with very low counts
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # Run the DESeq2 analysis
  dds <- DESeq(dds)
  
  # Obtain significant genes at 0.05
  res <- results(dds, contrast = c("condition", tissue, "ref"))
  resOrdered <- res[order(res$padj), ]
  sig_genes <- subset(as.data.frame(resOrdered), padj < 0.05)
  
  # Save the significant genes
  outfile <- paste0("DESeq2_Results/DESeq2_sig_", tissue, ".csv")
  write.csv(sig_genes, file = outfile, row.names = TRUE)
  
  # Save all genes
  outfile <- paste0("DESeq2_Results/DESeq2_all_", tissue, ".csv")
  write.csv(as.data.frame(resOrdered), file = outfile, row.names = TRUE)
  
  # Plot volcano plot
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

  top_genes <- rownames(res_df)[order(res_df$padj)][1:10]  # Select the top 20 by p-value


  p <- EnhancedVolcano(res_df,
                  lab = rownames(res_df),
                  selectLab = top_genes,
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = paste0(tissue, " vs. Reference"),
                  pCutoff = .05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 5.0)


  outfile <- paste0("DESeq2_Results/DESeq2_volcano_", tissue, ".png")
  ggsave(outfile, plot = p, width = 10, height = 10, dpi = 300)
  
  message("Finished analysis for tissue: ", tissue, " -> Results saved to: ", outfile)
}
