### Builds a DESeq2 object, normalises counts (median-ratio) and applies a 
### variance-stabilising transform, computes the sample–sample correlation matrix, 
### then plots a heatmap of those correlations with samples ordered by a fixed 
### vector and columns annotated by colony_location


set.seed(123) # for reproducibility

library("DESeq2")
library("dplyr")
library("pheatmap")

## All samples

### File retrieval and preparation

# Retrieve metadata
metadata <- read.csv("./metadata.csv", row.names=1)

# Retrieve first row of the counts file, which will be the sample ID
rawcounts <- read.table("./counts.tsv", nrows=1)
colnames <- as.character(rawcounts)

# Retrieve raw counts for all genes for all samples

# row.names: use the first column as the row names  
# col.names: use the column names we retrieved above as column names  
# check.names: in this instance prevents R from prepending an "X" to the column names  
# skip: skip the first row because it's just the sample ID  
rawcounts <- read.table("./counts.tsv", row.names=1, col.names=colnames,
                        check.names=FALSE,skip=1)

# Arrange rawcounts columns in the same 
# order as the metadata rows. This is required by DESeq2
rawcounts <- dplyr::select(rawcounts, rownames(metadata))


# Create the DESeq2 object specifying colony location as the condition of interest
dds <- DESeqDataSetFromMatrix(countData = rawcounts, 
                              colData = metadata,
                              design = ~ colony_location)

# Estimate size factors using median ratio method
dds <- estimateSizeFactors(dds)

# Extract the normalised counts  
normalised_counts <- counts(dds, normalized=TRUE)

# Variance stabilising transformation
vsd <- vst(dds, blind=TRUE)

# Extract matrix of normalised counts
vsd_matrix <- assay(vsd)

# Compute pairwise correlation values
vsd_correlation <- cor(vsd_matrix)


# Ensure reproducibility of output for tie-break cases
ord <- c("7", "5", "21", "30", "13", "14", "4", "31", "3", "15", "18",
         "23", "29", "11", "27", "12", "16", "9", "6", "19", "22", "32",
         "1", "10", "33", "2", "24", "28", "20", "34", "17", "25", "26")

# Group by colony location
pheatmap(
  vsd_correlation[ord, ord, drop = FALSE],
  annotation_col = dplyr::select(metadata, colony_location),
  annotation_names_col = FALSE
)



