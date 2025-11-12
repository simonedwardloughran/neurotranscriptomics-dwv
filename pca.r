### Loads metadata and count data, applies VST, and produces a PCA plot (PC1 vs PC2).
### Expects 'metadata.csv' and 'counts.tsv' to be present in the working directory.

# for reproducibility
set.seed(123)

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

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
vsd_matrix_top <- vsd_matrix[top_genes, ]

# Perform PCA on Filtered Data
pca_filtered <- prcomp(t(vsd_matrix_top))

# Calculate % variance explained for each PC
pca_var <- pca_filtered$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)

## Get PCA scores (1 and 2)
pca_scores <- as.data.frame(pca_filtered$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$sample <- colnames(vsd_matrix_top)

# Map behavioural_test ranks to labels ("poor","good")
# Recode colony_location to A/B for plotting
pca_scores$Performance <- factor(metadata$behavioural_test,
                                 levels = c("rank 0", "rank 21"),
                                 labels = c("poor", "good"))

colony_recode <- ifelse(metadata$colony_location == "cruickshank", "B",
                        ifelse(metadata$colony_location == "newburgh", "A", NA))
pca_scores$Colony <- factor(colony_recode, levels = c("A", "B"))

# Prepare ellipses
ellipse_data_1_2 <- pca_scores[, c("PC1", "PC2", "Colony")]
ellipse_data_1_2 <- ellipse_data_1_2[complete.cases(ellipse_data_1_2), ]


# Plot PC1 vs PC2
x_lab_1_2 <- paste0("PC1 (", pca_var_perc[1], "%)")
y_lab_1_2 <- paste0("PC2 (", pca_var_perc[2], "%)")

plot_1_2 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Colony, shape = Performance)) +
    stat_ellipse(data = ellipse_data_1_2,
                 aes(x = PC1, y = PC2, fill = Colony, group = Colony),
                 type = "norm",
                 geom = "polygon",
                 alpha = 0.1,
                 color = NA,
                 inherit.aes = FALSE) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = sample), size = 6, max.overlaps = 10) +
    labs(x = x_lab_1_2, y = y_lab_1_2,
         color = "Colony",
         fill = "Colony",
         shape = "Performance") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")

plot_1_2 <- plot_1_2 +
    theme(
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22)
    )

# Save hi-res plot
ggsave("PCA_PC1_PC2_Top500.png", plot = plot_1_2,
       width = 16, height = 12, dpi = 600, units = "in", bg = "white")

