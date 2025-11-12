### Build Apis mellifera gene lookup (Entrez ID, symbol, description) from symbols in counts.tsv
### Source: HymenopteraMine template "gene_symbol_gene_id" via InterMineR
### Cleans descriptions (custom replacements, trims trailing commas) and writes to SQL table
### Output table: am_ids_symbols_descriptions
### Expects 'metadata.csv' and 'counts.tsv' to be present in the working directory.

library(tidyverse)
library(InterMineR)
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

# Arrange columns to match metadata row order
rawcounts <- dplyr::select(rawcounts, rownames(metadata))

# Connect to db (add personal credentials)
con <- dbConnect(RMariaDB::MariaDB(), username="", password="", dbname="")

# Get character vector of all Apis mellifera gene symbols
geneSymbols <- rownames(rawcounts)

# Load Intermine, specifically HymenopteraMine and get all the templates
im <- initInterMine(mine=listMines()["HymenopteraMine"])
templates = getTemplates(im)

# Create query from template to get apis mellifera gene id from gene symbol
queryG = getTemplateQuery(
  im = im, 
  name = "gene_symbol_gene_id"
)

# Gene.description manually added to template
queryG$select[4] = 'Gene.description'

# Query a single symbol; returns (am_id, am_symbol, am_description)
symbolToId <- function(x){
  queryG$where[[1]][["value"]] <- x
  res <- runQuery(im, queryG)
  return(cbind(
    res$Gene.primaryIdentifier,
    res$Gene.symbol,
    res$Gene.description
  ))      
}

# Retrieve (takes a while)
res <- lapply(geneSymbols, FUN=symbolToId)

# Create dataframe from results
am_details <- data.frame(matrix(unlist(res), ncol=3, byrow=TRUE))

# Add column names
colnames(am_details) <- c("am_id", "am_symbol", "am_description")

# Drop rows with duplicate IDs
am_details <- am_details |>
  distinct(am_id, .keep_all = TRUE)

# Manual replacement of unorthodox/out-of-date description values
replacements = c(
  "uncharacterized LOC113218549" = "LOC113218549 (lncRNA)",
  "uncharacterized LOC100578156" = "LOC100578156 (lncRNA)",
  "uncharacterized LOC107965291" = "LOC107965291 (lncRNA)",
  "uncharacterized LOC113218947" = "LOC113218947 (lncRNA)",
  "uncharacterized LOC102654076" = "LOC102654076 (lncRNA)",
  "uncharacterized LOC100577198" = "LOC100577198 (lncRNA)",
  "uncharacterized LOC102655375" = "LOC102655375 (lncRNA)",
  "uncharacterized LOC726505" = "LOC726505",
  "uncharacterized LOC411159" = "LOC411159",
  "uncharacterized LOC100578156" = "LOC100578156",
  "uncharacterized LOC102654257" = "sosie oogenesis-related protein sosie"
)

am_details <- am_details %>%
  mutate(am_description = str_replace_all(am_description, replacements))

# Remove trailing commas from description field
am_details <- am_details %>%
  mutate(am_description = str_replace(am_description, ",\\s*$", ""))


# Write id -> symbols look up table to database
dbWriteTable(con, "am_ids_symbols_descriptions", am_details)





























