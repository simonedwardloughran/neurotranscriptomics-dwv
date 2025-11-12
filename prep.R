### Loads metadata and counts, merges gene ID/symbol/description lookups,
### constructs SummarizedExperiment objects (full and subset by colony/behaviour),
### defines viral-load contrast, and performs DESeq2 differential expression analyses.
### Includes helper function to extract annotated significant DE genes.
### Assumes metadata.csv and counts.tsv are present in the current working directory.
### Assumes a database table mapping Apis mellifera IDs, symbols, and descriptions is available.
### See am_details.R for details on how the lookup table is populated.

set.seed(123) # for reproducibility


# Load libraries
library(SummarizedExperiment)
library(EnrichmentBrowser)
library(DESeq2)
library(dplyr)
library(stringr)
library(DBI)
library(RMariaDB)


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


# Retrieve Apis mellifera ID, symbol, description  from database
# See am_details.R
con <- dbConnect(RMariaDB::MariaDB(), username="", password="", dbname="")
query  <- "SELECT am_id, am_symbol, am_description FROM am_ids_symbols_descriptions"
res <- dbSendQuery(con, query)
am_ids_symbols_desc <- dbFetch(res)

# Add the current rownames (symbols) as a column, then merge in ID/description
rawcounts$apis_symbol <- rownames(rawcounts)
rawcounts <- merge(rawcounts, am_ids_symbols_desc, by.x="apis_symbol", by.y="am_symbol")
rownames(rawcounts) <- rawcounts$am_id

# Prepare rowData (symbols) for the SummarizedExperiment
rowdata <- rawcounts[, c("apis_symbol")]

# Keep only the count columns (drop helper/annotation fields)
rawcounts <- as.matrix(rawcounts[, c(2:34)])

# Get list of all the ids
am_ids <- am_ids_symbols_desc$am_symbol

# Create SummarizedExperiment object
se <- SummarizedExperiment(assay = list(count=rawcounts), colData = metadata, rowData = rowdata)
metadata(se)$annotation <- "ame"

# Colony subsets
se.newburgh.viral <- subset(se, , colony_location == "newburgh")
se.cruickshank.viral <- subset(se, , colony_location == "cruickshank")

# Behavioural test subsets (poor = rank 0; good = rank 21)
se.poor <- subset(se, , behavioural_test == "rank 0")
se.good <- subset(se, , behavioural_test == "rank 21")

# Blocking factor by colony
se.poor$BLOCK <- se.poor$colony_location
se.good$BLOCK <- se.good$colony_location

# Define contrasts: high vs low viral load (1 = high, 0 = low)
se$GROUP <- ifelse(se$viral_load == 'high', 1, 0)
se.newburgh.viral$GROUP <- ifelse(se.newburgh.viral$viral_load == 'high', 1, 0)
se.cruickshank.viral$GROUP <- ifelse(se.cruickshank.viral$viral_load == 'high', 1, 0)
se.poor$GROUP <- ifelse(se.poor$viral_load == 'high', 1, 0)
se.good$GROUP <- ifelse(se.good$viral_load == 'high', 1, 0)

# DESEq2 analyses
se <- deAna(se, de.method = 'DESeq2', filter.by.expr = TRUE)
se.newburgh.viral <- deAna(se.newburgh.viral, de.method = 'DESeq2', filter.by.expr = TRUE)
se.cruickshank.viral <- deAna(se.cruickshank.viral, de.method = 'DESeq2', filter.by.expr = TRUE)
se.poor <- deAna(se.poor, de.method = 'DESeq2', filter.by.expr = TRUE)
se.good <- deAna(se.good, de.method = 'DESeq2', filter.by.expr = TRUE)


# Helper:
# Extract significantly DE genes from a SummarizedExperiment and merge symbol/description
getSigs <- function(group, direction, sigp, sigl2fc){
    sumEx <- get(group)
    rd <- na.omit(rowData(sumEx))
    if(direction == 'up'){
        sig <- rd[rd$ADJ.PVAL < sigp & rd$FC >= sigl2fc,]
    }else{
        sig <- rd[rd$ADJ.PVAL < sigp & rd$FC <= sigl2fc*-1,]
    }    
    if(nrow(sig)>0){
        # Get significant IDs and join to lookup
        
        ids <- rownames(sig)
        sub <- am_ids_symbols_desc[am_ids_symbols_desc$am_id %in% ids,]
        
        # Order by absolute fold change (descending)
        merged <- merge(sig, sub, by.x="X", by.y="am_symbol")        
        
        # order by decreasing fold change
        merged <- as.data.frame(merged[order(abs(merged$FC), decreasing=T),])    
        return(merged)    
    }else{
        return('np')
    }
}

# Global significance thresholds
sigp <- 0.05
sigl2fc <- 1
