### Generates annotated volcano plots for selected DESeq2 result sets.
### It loads the pre-processed SummarizedExperiment objects (via prep.R),
### merges significant genes with their Apis mellifera symbol/description annotations,
### and visualises differential expression (log2 fold change vs adjusted p-value),
### highlighting up- and down-regulated genes and labelling significant points.

# Loads se objects and am_ids_symbols_desc lookup into the environment
source("prep.R")

library(ggplot2)
library(ggrepel)
library(grid)

# Merge rowData from a SummarizedExperiment with gene annotation
mergeRowDataAmDetails <- function(group){
  
  sig <- get(group)
  
  # IDs of significant genes in this SE object
  ids <- rownames(sig)
  
  # Look up annotation entries for these IDs
  sub <- am_ids_symbols_desc[am_ids_symbols_desc$am_id %in% ids,]
  
  # Extract rowData and drop NAs
  rd <- na.omit(rowData(sig))
  
  # Merge gene-level statistics with annotation
  merged <- merge(rd, sub, by.x="X", by.y="am_symbol")        
  
  # Order by absolute fold change (largest first)
  merged <- as.data.frame(merged[order(abs(merged$FC), decreasing=T),])    
  return(merged)
  
}


## Quick inspection of significant genes at different thresholds
getSigs('se.poor', 'up', 0.05, 0)
getSigs('se.poor', 'down', 0.05, 0)
getSigs('se.poor', 'up', 0.05, 1)
getSigs('se.poor', 'down', 0.05, 1)

getSigs('se.newburgh.viral', 'up', 0.05, 0)
getSigs('se.newburgh.viral', 'down', 0.05, 0)
getSigs('se.newburgh.viral', 'up', 0.05, 1)
getSigs('se.newburgh.viral', 'down', 0.05, 1)

getSigs('se.cruickshank.viral', 'up', 0.05, 0)
getSigs('se.cruickshank.viral', 'down', 0.05, 0)

getSigs('se.good', 'up', 0.05, 0)
getSigs('se.good', 'down', 0.05, 0)



# Volcano plot function
volcanoPlot <- function(df, label, xfrom, xto){
  
  # Find the greatest of the minimum and maximum log2FoldChange column
  # so that we can make the plot symmetrical
  xlimit <-
    max(
      abs(floor(min(na.omit(df$FC)))),
      ceiling(max(na.omit(df$FC)))
    )
  
  
  # Find greatest y-limit
  ylimit <- ceiling(max(na.omit(-log10(df$ADJ.PVAL))))
  
  df <- na.omit(df)
  
  # Add significance type column
  df <- df %>%
    dplyr::mutate(gene_type = case_when(FC >= 1 & ADJ.PVAL <= sigp ~ "up",
                                        FC <= -1*sigl2fc & ADJ.PVAL <= sigp ~ "down",
                                        TRUE ~ "ns")) 
  
  # Add gene_id as its own column to the dataframe
  df$gene_id <- rownames(df)
  
  
  # Add -log10 calculation to the dataframe as a column to make it clearer what we're plotting
  df$minusLog10padj <- -log10(df$ADJ.PVAL)
  
  # Specify colour, alphas and sizes for the plot
  cols <- c("up" = "#ffcd93", "down" = "#86d3ff", "ns" = "grey") 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  
  # Subset significant genes
  sig_genes <- subset(df, gene_type == "down" | gene_type == "up")
  
  # draw the plot
  volcano_plot_viral_load <- 
    ggplot(df,
           aes(x = FC, 
               y = minusLog10padj,
               fill = gene_type,
               size = gene_type,
               alpha = gene_type)
    ) + 
    geom_point(shape=21) +
    geom_hline(yintercept = -log10(sigp), col = "red", linetype = "dashed") +
    geom_vline(xintercept = c(-1*sigl2fc, sigl2fc), col = "red", linetype = "dashed") +
    geom_label_repel(
      data = sig_genes,
      aes(label = am_description),
      size = 3.3,
      force = 2,
      nudge_y = 5,
      max.overlaps = nrow(sig_genes),
      show.legend = FALSE,
      label.padding = unit(0.5, "lines")
    ) +
    geom_point(
      aes(fill = gene_type),
      shape = 21,
      size = 2,
      alpha = 0.8
    ) +
    scale_fill_manual(
      values = cols,
      name = "Expression",
      guide = guide_legend(override.aes = list(shape = 22, size = 8)) # 22 = filled square
    ) +
    guides(
      size = "none",
      alpha = "none"
    ) +
    scale_size_manual(values = sizes) +
    scale_alpha_manual(values = alphas) + 
    scale_x_continuous(breaks = c(seq(-xfrom, xto, 1)),
                       limits = c(-xfrom, xto)) +
    scale_y_continuous(breaks = c(seq(0, ylimit, 2))) +
    labs(title = label,
         x = expression(paste("Log"[2], " fold change")),
         y = expression(paste("-Log"[10], "P"))) +    
    
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid   = element_line(color = "#eeeef5", size = 0.3),
      
      # Axis titles
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      
      # Axis tick labels
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      
      # Legend text
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 15)
    )
  
  return(volcano_plot_viral_load)
}

# Get volcano plots
volcanoPlot(mergeRowDataAmDetails('se.poor'), "Gene expression changes in poorly performing bees: High vs low viral load", 7, 7)
volcanoPlot(mergeRowDataAmDetails('se.newburgh.viral'), "Gene expression changes in Newburgh colony bees. High vs low viral load", 5, 7)
