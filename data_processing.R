library(tidyverse)
library(data.table)

#load in data
dat = read.csv(file = "25_03_16_Parsed_Revised_full.csv")

load("raw_data.rda")

#make a dataframe with just the tissues with genes as row names:
dat <- dat[,c(1,8,10:101)]

#rename first NFx4 tissues with same naming conventions as the rest of the tissues
colnames(dat)[3:6] <- c("NFx4_1", "NFx4_2", "NFx4_3", "NFx4_4")

# Convert to data.table
setDT(dat)

# Extract unique tissue prefixes
tissue_cols <- names(dat)[-(1:2)]
tissue_names <- unique(gsub("_[1-4]$", "", tissue_cols))

# Function to count valid tissues 
count_valid_tissues <- function(dt, tissue_prefixes) {
  valid_counts <- rowSums(sapply(tissue_prefixes, function(tissue) {
    rowSums(dt[, grep(paste0("^", tissue, "_[1-4]$"), names(dt)), with = FALSE] > 0.125) >= 3
  }))
  return(valid_counts)
}

# Function to count gene categories
count_gene_categories <- function(dt) {
  gene_status <- dt[, .(spliced = any(Splicing == "Spliced"), 
                        unspliced = any(Splicing == "unspliced")), 
                    by = Gene]
  
  only_spliced <- sum(gene_status$spliced & !gene_status$unspliced)
  only_unspliced <- sum(!gene_status$spliced & gene_status$unspliced)
  both <- sum(gene_status$spliced & gene_status$unspliced)
  
  return(list(only_spliced = only_spliced, only_unspliced = only_unspliced, both = both))
}

# Count gene categories in raw data
raw_gene_counts <- count_gene_categories(dat)

# Filter out rows with less than two valid tissues
dat$valid_tissues = count_valid_tissues(dat, tissue_names)
dat <- dat[valid_tissues >= 2]

filtered_gene_counts <- count_gene_categories(dat)

dat <- as.data.frame(dat)

#Calculate Splicing Efficiency
dat <- dat %>%
  group_by(Gene, Splicing) %>%
  summarise(across(NFx4_1:lfBr_4, sum)) %>%
  pivot_wider(names_from = Splicing, values_from = c(NFx4_1:lfBr_4)) %>%
  mutate(across(contains("_Spliced"), 
                .fns = ~ {
                  spliced_TPM <- .
                  unspliced_TPM <- get(sub("_Spliced", "_unspliced", cur_column()))
                  total_TPM <- sum(spliced_TPM, unspliced_TPM, na.rm = TRUE)
                  
                  if (is.na(spliced_TPM)) return(0) 
                  if (is.na(unspliced_TPM)) return(1) 
                  if (total_TPM == 0) return(NA)       
                  
                  return(spliced_TPM / total_TPM)      # Default calculation
                },
                .names = "{.col}_efficiency")) %>%
  ungroup() %>%
  select(Gene, contains("_efficiency"))
    
#Multiply by 100 and round to nearest integer
dat <- dat %>%
  mutate(across(contains("_efficiency"), ~ round(. * 100))) # Multiply by 100 and round
  
#Calculate reference values for each gene in each replicate
colnames(dat)[2:93] <- tissue_cols

#Create empty reference dataframe
reference <- data.frame(matrix(NA, nrow = nrow(dat), ncol = 4))
colnames(reference) <- c("rep_1", "rep_2", "rep_3", "rep_4")

#loop through the 4 replicates and find reference values for each gene
for (i in 1:4){
  #Form temporary dataframe with all replicate i's
  temp = dat %>%
    select(contains(paste0("_",i)))
  
  #Calculate mean across that row and add to reference dataframe
  temp = rowMeans(temp, na.rm = T)
  reference[,i] = temp
  rm(temp)
}

#add in a gene name column to the reference dataframe and round reference values
reference <- round(reference)
reference <- add_column(reference, Gene = dat$Gene, .before = "rep_1")

write.csv(reference, file = "reference.csv", quote = F, row.names = F)
write.csv(dat, file = "processed_data_NAs_included.csv", quote = F, row.names = F)

#Replace NA values in the data with corresponding reference values
gene_col <- "Gene"  
replicate_cols <- grep("_\\d+$", colnames(dat), value = TRUE)

for (col in replicate_cols) {
  # Extract the replicate number (e.g., 1, 2, 3, 4)
  rep_num <- as.integer(gsub("^.*_(\\d+)$", "\\1", col))
  
  # Get the corresponding reference column name
  ref_col <- paste0("rep_", rep_num)
  
  # Replace NA values with the reference value for that replicate and gene
  dat[[col]] <- ifelse(
    is.na(dat[[col]]),
    reference[match(dat[[gene_col]], reference[[gene_col]]), ref_col],
    dat[[col]]
  )
}

length(which(is.na(dat)))

write.csv(dat, file = "processed_data.csv", quote = F, row.names = F)

