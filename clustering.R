### This script assesses within-colony clustering structure in expression space.
### It loads RNA-seq metadata and counts, performs DESeq2 normalization and VST,
### selects the top 500 most variable genes, and then quantifies sample clustering
### by colony using two complementary approaches:
### (1) silhouette widths based on sample–sample distances, and
### (2) within-colony distances in PCA space (PC1–PC3).
### Differences between colonies are tested using Wilcoxon rank-sum tests.

library(DESeq2)
library(dplyr)
library(cluster)

# Retrieve metadata
metadata <- read.csv("./metadata.csv", row.names=1)

# Read header row to get sample IDs for column names
rawcounts <- read.table("./counts.tsv", nrows=1)
colnames <- as.character(rawcounts)

# Retrieve raw counts for all genes for all samples
# row.names: use the first column as the row names  
# col.names: use the column names we retrieved above as column names  
# check.names: in this instance prevents R from prepending an "X" to the column names  
# skip: skip the first row because it's just the sample ID  
rawcounts <- read.table("./counts.tsv", row.names=1, col.names=colnames,
                        check.names=FALSE,skip=1)

# Arrange columns to match metadata row order (required by DESeq2)
rawcounts <- dplyr::select(rawcounts, rownames(metadata))


# Create the DESeq2 object specifying colony location as the condition of interest
dds <- DESeqDataSetFromMatrix(countData = rawcounts, 
                              colData = metadata,
                              design = ~ colony_location)

# Estimate size factors using median ratio method
dds <- estimateSizeFactors(dds)

# Extract the normalised counts  
normalised_counts <- counts(dds, normalized=TRUE)

# Variance stabilising transformation (blind to design)
vsd <- vst(dds, blind=TRUE)

# Extract VST-transformed expression matrix
vsd_matrix <- assay(vsd)

# Compute pairwise correlation values
vsd_correlation <- cor(vsd_matrix)


# Calculate variance across samples for each gene
gene_vars <- apply(vsd_matrix, 1, var)

# Get top 500 most variable genes across samples
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:500]

# Subset the VSD matrix
vsd_matrix <- vsd_matrix[top_genes, ]



# Recode Colony: cruickshank → B, newburgh → A

colony_recode <- ifelse(metadata$colony_location == "cruickshank", "B",
                        ifelse(metadata$colony_location == "newburgh", "A", NA))
metadata$Colony <- factor(colony_recode, levels = c("A", "B"))




## Silhouette Width Calculation + Test

# Compute distance matrix across samples (using all genes or filtered)
dist_matrix <- dist(t(vsd_matrix))  # Transpose so samples = rows

# Convert colony labels to numeric (1 = A, 2 = B)
group_numeric <- as.numeric(metadata$Colony)

# Calculate silhouette widths for each sample
sil <- silhouette(group_numeric, dist_matrix)

# Extract silhouette widths and colony labels
sil_widths <- sil[, "sil_width"]
colony_labels <- metadata$Colony

# Average silhouette width per colony
mean_sil_A <- mean(sil_widths[colony_labels == "A"])
mean_sil_B <- mean(sil_widths[colony_labels == "B"])

cat("Average Silhouette Width:\n")
cat("Colony A:", round(mean_sil_A, 3), "\n")
cat("Colony B:", round(mean_sil_B, 3), "\n\n")

# Wilcoxon test to compare silhouette widths between colonies
wilcox_test_sil <- wilcox.test(sil_widths[colony_labels == "A"],
                               sil_widths[colony_labels == "B"])

cat("Wilcoxon test for Silhouette Widths:\n")
print(wilcox_test_sil)

## PCA Distance Calculation + Test

# Perform PCA on VSD matrix
pca_res <- prcomp(t(vsd_matrix))
pca_scores <- as.data.frame(pca_res$x[, 1:3])  # Use PC1–PC3
pca_scores$Colony <- metadata$Colony

# Pairwise distances within each colony (Euclidean)
dist_A_vals <- as.numeric(dist(pca_scores[pca_scores$Colony == "A", 1:3]))
dist_B_vals <- as.numeric(dist(pca_scores[pca_scores$Colony == "B", 1:3]))

# Average distance within each colony
mean_dist_A <- mean(dist_A_vals)
mean_dist_B <- mean(dist_B_vals)

cat("Average Within-Colony Distance in PCA Space:\n")
cat("Colony A:", round(mean_dist_A, 3), "\n")
cat("Colony B:", round(mean_dist_B, 3), "\n\n")

# Wilcoxon test to compare PCA distances
wilcox_test_dist <- wilcox.test(dist_A_vals, dist_B_vals)

cat("Wilcoxon test for Within-Colony PCA Distances:\n")
print(wilcox_test_dist)


## Amount of variance explained by PC1-PC3
pca_var <- pca_res$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)

# Sum for PC1–PC3
sum(pca_var_perc[1:3])

