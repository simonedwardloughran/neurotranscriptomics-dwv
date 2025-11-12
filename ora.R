### Directional GO and KEGG ORA on VST-normalized gene sets,
### plus annotation of enriched terms with GO names/definitions,
### and export of both gene-level memberships and term-level summaries to the database.

# Loads se objects and am_ids_symbols_desc lookup into the environment
source("prep.R")

library(biomaRt)
library(stringr)
library(qusage)

# Load Apis mellifera gene set file
am.go.gs <- read.gmt('./am.gmt')


# Normalise SummarizedExperiments
se.newburgh.viral <- normalize(se.newburgh.viral, norm.method = "vst")
se.poor <- normalize(se.poor, norm.method = "vst")

# Split each into up and down regulated genes
se.poor.up <- subset(se.poor, FC > 0, )
se.poor.down <- subset(se.poor, FC < 0, )
se.newburgh.viral.up <- subset(se.newburgh.viral, FC > 0, )
se.newburgh.viral.down <- subset(se.newburgh.viral, FC < 0, )

# Carry out ORA
ora.se.poor.up <- sbea(method = "ora", se = se.poor.up, gs = am.go.gs, perm=0, sig.stat="fc", padj.method="bonferroni")
ora.se.poor.down <- sbea(method = "ora", se = se.poor.down, gs = am.go.gs, perm=0, padj.method="bonferroni")
ora.se.newburgh.viral.up <- sbea(method = "ora", se = se.newburgh.viral.up, gs = am.go.gs, perm=0, padj.method="bonferroni")
ora.se.newburgh.viral.down <- sbea(method = "ora", se = se.newburgh.viral.down, gs = am.go.gs, perm=0, padj.method="bonferroni")


# Creates a database table based on the sbea_result name and creates a row
# for each GO term:Gene combination
writeGenesToDB <- function(sbea_result){
  
  table_name <- paste(gsub('\\.', '', sbea_result), '_genes', sep='')
  
  cols <- c('varchar(100)', 'varchar(100)', 'DOUBLE', 'DOUBLE', 'DOUBLE', 'VARCHAR(100)')
  names(cols) <- c('go_term', 'symbol', 'l2fc', 'p', 'padj', 'id')
  dbCreateTable(con, table_name, cols)    
  
  sbea_result = get(sbea_result)
  
  goTerms <- sbea_result$res.tbl$GENE.SET    
  
  for(i in 1:length(goTerms)){
    goTerm <- goTerms[i]
    ids <- sbea_result$gs[[goTerm]]
    rd <- rowData(sbea_result$se)
    rd$id <- rownames(rd)
    rows <- rd[rd$id %in% ids,]
    rows <- rows[order(rows$PVAL, decreasing = F),]
    
    for(j in 1:nrow(rows)){
      df <- data.frame(goTerm, rows[j,]$X, rows[j,]$FC, rows[j,]$PVAL, rows[j,]$ADJ.PVAL, rows[j,]$id)
      colnames(df) <- c('go_term', 'symbol', 'l2fc', 'p', 'padj', 'id')
      dbAppendTable(con, table_name, df)
    }
  }
}



writeGenesToDB("ora.se.poor.up")
writeGenesToDB("ora.se.poor.down")
writeGenesToDB("ora.se.newburgh.viral.up")
writeGenesToDB("ora.se.newburgh.viral.down")



# Connect to the Ensembl Metazoa BioMart
ensembl_metazoa = useEnsemblGenomes(biomart = "metazoa_mart")

# Select the Apis mellifera gene dataset within Ensembl Metazoa
ensembl_apis_mellifera <- useEnsemblGenomes(
  biomart = "metazoa_mart", 
  dataset = "amellifera_eg_gene")


# GO annotation fields to retrieve for each gene
attributes <- c("ensembl_gene_id", "go_id", "name_1006", "definition_1006")

# Retrieve GO annotations for Apis mellifera genes
annot <- getBM(
  attributes = attributes,
  mart = ensembl_apis_mellifera,
  uniqueRows = T
)

# Sort by go_id
annot <- annot[order(annot$go_id, decreasing = F),]

# Eemove rows with empty go_id column
annot <- annot[!(is.na(annot$go_id) | annot$go_id == ''),]

# Remove 'LOC' and 'GeneID_' prefixes
for(i in 1:nrow(annot)){
  annot[i,1] <- str_replace(annot[i,1], "LOC", "") %>% 
    str_replace("GeneID_", "")
}




### For each enriched GO term, look up the GO name and definition retrieved from BioMart,
### append these fields to the ORA results table, rename columns for clarity,
### and write the annotated results to the database.

gr <- gsRanking(ora.se.poor.up)
for(x in 1:nrow(gr)){
  gs <- gr[x,]$GENE.SET
  name <- annot[annot$go_id == gs, ][1,]$name_1006
  def <- annot[annot$go_id == gs, ][1,]$definition_1006
  gr[x,5] <- name
  gr[x,6] <- def
}

colnames(gr) <- c("go_id", "num_genes", "num_sig_genes", "p", "go_name", "go_def")
dbWriteTable(con, "ora_deseq2_poor_up", data.frame(gr))

gr


gr <- gsRanking(ora.se.poor.down)
for(x in 1:nrow(gr)){
  gs <- gr[x,]$GENE.SET
  name <- annot[annot$go_id == gs, ][1,]$name_1006
  def <- annot[annot$go_id == gs, ][1,]$definition_1006
  gr[x,5] <- name
  gr[x,6] <- def
}
colnames(gr) <- c("go_id", "num_genes", "num_sig_genes", "p", "go_name", "go_def")
dbWriteTable(con, "ora_deseq2_poor_down", data.frame(gr))
gr


gr <- gsRanking(ora.se.newburgh.viral.up)
for(x in 1:nrow(gr)){
  gs <- gr[x,]$GENE.SET
  name <- annot[annot$go_id == gs, ][1,]$name_1006
  def <- annot[annot$go_id == gs, ][1,]$definition_1006
  gr[x,5] <- name
  gr[x,6] <- def
}
colnames(gr) <- c("go_id", "num_genes", "num_sig_genes", "p", "go_name", "go_def")
dbWriteTable(con, "_ora_deseq2_newburgh_viral_up", data.frame(gr))
gr


gr <- gsRanking(ora.se.newburgh.viral.down)
for(x in 1:nrow(gr)){
  gs <- gr[x,]$GENE.SET
  name <- annot[annot$go_id == gs, ][1,]$name_1006
  def <- annot[annot$go_id == gs, ][1,]$definition_1006
  gr[x,5] <- name
  gr[x,6] <- def
}
colnames(gr) <- c("go_id", "num_genes", "num_sig_genes", "p", "go_name", "go_def")
dbWriteTable(con, "_ora_deseq2_newburgh_viral_down", data.frame(gr))
gr




# Fetch KEGG gene sets for A. mellifera and run directional ORA per subset.
am.kegg.gs <- getGenesets(org = "ame", db = "kegg", onto="BP", gene.id.type="ENTREZID")
ora.se.poor.kegg.up <- sbea(method = "ora", se = se.poor.up, gs = am.kegg.gs, perm = 0, beta="1", sig.stat="fc",  padj.method="none")
ora.se.poor.kegg.down <- sbea(method = "ora", se = se.poor.down, gs = am.kegg.gs, perm = 0, beta="1", sig.stat="fc", padj.method="none")
ora.se.newburgh.viral.kegg.up <- sbea(method = "ora", se = se.newburgh.viral.up, gs = am.kegg.gs, perm = 0, beta="1", sig.stat="fc", padj.method="none")
ora.se.newburgh.viral.kegg.down <- sbea(method = "ora", se = se.newburgh.viral.down, gs = am.kegg.gs, perm = 0, beta="1", sig.stat="fc",  padj.method="none")


# Parse KEGG gene set IDs/names into separate columns for readability

gr <- gsRanking(ora.se.poor.kegg.up)
gr <- data.frame(gr)
for(x in 1:nrow(gr)){
    gs <- gr[x,]$GENE.SET
    parts <- str_split(gs, "_", n=2)
    id <- parts[[1]][1]
    id <- str_replace(id, "ame", "")
    name <- parts[[1]][2]
    gr[x,5] <- id
    gr[x,6] <- name
}
colnames(gr) <- c("gene_set", "num_genes", "num_sig_genes", "p", "id", "name")
gr


gr <- gsRanking(ora.se.poor.kegg.down)
gr <- data.frame(gr)
for(x in 1:nrow(gr)){
    gs <- gr[x,]$GENE.SET
    parts <- str_split(gs, "_", n=2)
    id <- parts[[1]][1]
    id <- str_replace(id, "ame", "")
    name <- parts[[1]][2]
    gr[x,5] <- id
    gr[x,6] <- name
}
colnames(gr) <- c("gene_set", "num_genes", "num_sig_genes", "p", "id", "name")
gr

gr <- gsRanking(ora.se.newburgh.viral.kegg.up)
gr <- data.frame(gr)
for(x in 1:nrow(gr)){
    gs <- gr[x,]$GENE.SET
    parts <- str_split(gs, "_", n=2)
    id <- parts[[1]][1]
    id <- str_replace(id, "ame", "")
    name <- parts[[1]][2]
    gr[x,5] <- id
    gr[x,6] <- name
}
colnames(gr) <- c("gene_set", "num_genes", "num_sig_genes", "p", "id", "name")
gr

gr <- gsRanking(ora.se.newburgh.viral.kegg.down)
gr <- data.frame(gr)
for(x in 1:nrow(gr)){
    gs <- gr[x,]$GENE.SET
    parts <- str_split(gs, "_", n=2)
    id <- parts[[1]][1]
    id <- str_replace(id, "ame", "")
    name <- parts[[1]][2]
    gr[x,5] <- id
    gr[x,6] <- name
}

colnames(gr) <- c("gene_set", "num_genes", "num_sig_genes", "p", "id", "name")
gr
