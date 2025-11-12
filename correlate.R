# Computes gene-wise Spearman correlations between log10-normalised expression 
# and DWV load for each differential-expression contrast, and generates 
# correlation plots for selected genes

# Loads se objects and am_ids_symbols_desc lookup into the environment
source("prep.R")

library(ggplot2)
library(cowplot)
library(grid)

# Runs a Spearman correlation between counts and viral_loads, 
# then returns just the p-value and the correlation coefficient in a named vector
calcSpearman <- function(counts, viral_loads){
  res <- cor.test(viral_loads, counts, method='spearman', exact=FALSE)
  bits <- c(res$p.value, res$estimate)
  names(bits) <- c("corr_p", "rho")
  bits
}


# Takes a SummarizedExperiment (se) and a direction (up/down).
# Pulls out the significant genes for that direction using getSigs(). 
# Then it pulls out the underlying data object s and rebuilds a DESeq2 dataset 
# from it so it can recompute normalised counts
# Calls calcSpearman to get p-value and correlation coefficient and merges these with 
# significance results and the gene annoation table
mergeCors <- function(se, direction){
  sigs <- getSigs(se, direction, sigp, sigl2fc)  
  s <- get(se)
  dds <- DESeqDataSetFromMatrix(countData = assay(s), 
                                colData = colData(s),
                                design = ~ viral_load)
  dds <- estimateSizeFactors(dds)
  normalised_counts <- counts(dds, normalized=TRUE)
  
  # Extract numeric DWV load vector
  viral_load_values <- as.numeric(colData(s)["viral_load_value"][,1])
  
  # Gene-wise Spearman correlations; transpose result and order by rho
  corrs <- apply(as.data.frame(normalised_counts), 1, calcSpearman, viral_load_values)
  corrs <- as.data.frame(t(corrs[, order(corrs[which(rownames(corrs) == 'rho'),]) ]))
  corrs$am_id <- rownames(corrs)
  
  # Merge: significant genes x correlations x AM lookup
  merged <- merge(sigs, corrs, by='am_id')
  moremerged <- merge(merged, am_ids_symbols_desc, by='am_id')  
  moremerged
}

# Build merged tables for each contrast/direction (used later for plotting)
mergeduppoor <- mergeCors('se.poor', 'up');
mergeddownpoor <- mergeCors('se.poor', 'down');
mergedupnewburgh <- mergeCors('se.newburgh.viral', 'up');
mergeddownnewburgh <- mergeCors('se.newburgh.viral', 'down');

# Prepare DESeq2 objects and normalized counts for specific cohorts
dds_newburgh <- DESeqDataSetFromMatrix(countData = assay(se.newburgh.viral), 
                                       colData = colData(se.newburgh.viral),
                                       design = ~ viral_load)
dds_newburgh <- estimateSizeFactors(dds_newburgh)


dds_poor <- DESeqDataSetFromMatrix(countData = assay(se.poor), 
                                   colData = colData(se.poor),
                                   design = ~ viral_load)
dds_poor <- estimateSizeFactors(dds_poor)

normalised_counts_newburgh <- counts(dds_newburgh, normalized=TRUE)
normalised_counts_poor <- counts(dds_poor, normalized=TRUE)

# Numeric DWV loads (order must match columns in the count matrices)
viral_load_newburgh <- colData(se.newburgh.viral)["viral_load_value"]
viral_loads_newburgh <- as.numeric(viral_load_newburgh[,1])

viral_load_poor <- colData(se.poor)["viral_load_value"]
viral_loads_poor <- as.numeric(viral_load_poor[,1])


# Given (df = data.frame with columns x = log10 DWV load, y = log10 norm counts)
# and a 'row' containing am_id, rho, and corr_p, produce a plot
get_correlation_plot <- function(df, row, p_digits = 3){
  
  # Values for the annotation
  am_id <- row$am_id
  rho   <- row$rho
  pval  <- row$corr_p
  
  # Build a plotmath-ready p-value
  .p_plotmath <- function(p, digits = 3){
    if (length(p) == 0 || is.na(p)) return('"NA"')
    if (p >= 0.001) {
      return(paste0('"', formatC(p, format = "f", digits = digits), '"'))
    } else {
      sc <- formatC(p, format = "e", digits = digits - 1)
      parts <- strsplit(sc, "e", fixed = TRUE)[[1]]
      coeff <- parts[1]
      expo  <- as.integer(parts[2])
      coeff <- sub("0+$", "", coeff)
      coeff <- sub("[.]$", "", coeff)
      return(paste0(coeff, " %*% 10^", expo))
    }
  }
  
  # Compose the right-hand-side label "rho = X, p = Y" in plotmath syntax
  rhs <- paste0(
    'paste("rho = ", "', sprintf("%.2f", rho),
    '", ", ", italic(p), " = ", ', .p_plotmath(pval, p_digits), ')'
  )
  
  # Base scatter with clean theming and margin space for the annotation
  base <- ggplot(df, aes(x, y)) +
    geom_point(color="black", size=4, shape=21, fill="lightblue") +
    labs(
      x = expression(paste("Log"[10], " DWV load")),
      y = expression(paste("Log"[10], " normalised count"))
    ) +
    theme_bw() +
    theme(
      plot.margin = margin(36, 24, 10, 24, "pt"),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 1),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 22, margin = unit(c(3,0,0,0), "mm")),
      axis.text.y = element_text(size = 22, margin = unit(c(0,3,0,0), "mm")),
      axis.title.x = element_text(size = 24, margin = unit(c(4,0,0,0), "mm")),
      axis.title.y = element_text(size = 24, margin = unit(c(0,4,0,0), "mm")),
      axis.ticks.length = unit(-0.25, "cm"),
      axis.ticks = element_line(size = 1)
    ) +
    annotate("text", x = Inf, y = Inf, label = rhs, parse = TRUE,
             hjust = 1.02, vjust = -0.2, size = 8) +
    coord_cartesian(clip = "off")
  
  # Add a left-top corner label "LOC<am_id>"
  cowplot::ggdraw(base) +
    cowplot::draw_label(
      paste0("LOC", am_id),
      x = 0.028, y = 0.985, hjust = 0, vjust = 1, size = 24
    )
}


# PLOTS

merged = mergeduppoor

amId = 113218549
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 406140
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 406142
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)




merged = mergeddownpoor

amId = 100577198
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 102654257 
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 102655185 
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 107965291 
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 406131
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 410326    
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 411159    
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 412829    
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 724252    
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 724642    
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 726505
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 726935    
counts <- normalised_counts_poor[as.character(amId), ]
df <- data.frame(x = viral_loads_poor, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)


merged = mergedupnewburgh

amId = 100578156  
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 102654076 
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 113218549 
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 113218947 
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 406140        
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 406142    
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 406144    
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 551263    
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)


merged = mergeddownnewburgh

amId = 102655319 
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 102655375 
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 406155          
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 408544    
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)

amId = 724899    
counts <- normalised_counts_newburgh[as.character(amId), ]
df <- data.frame(x = viral_loads_newburgh, y = log10(counts))
row <- merged[merged$am_id == as.character(amId), ]
get_correlation_plot(df, row)
