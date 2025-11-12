### Infers gene regulatory networks (GRNs) for the poor learners and Newburgh 
### colony subsets using BioNERO, based on predefined lncRNA candidate regulators

set.seed(123) # for reproducibility

# Loads se objects and am_ids_symbols_desc lookup into the environment
source("prep.R")

# Load required packages
library(BioNERO)

# Potential transcription factors
# These are the lncRNAs that were identified by DESeq2
am_tf_poor <- as.matrix(read.table("./tf_poor.csv"))
am_tf_newburgh <- as.matrix(read.table("./tf_newburgh.csv"))

colnames(am_tf_poor) <- c("x")
colnames(am_tf_newburgh) <- c("x")

# Preprocessing
se.poor.pp <- exp_preprocess(
  se.poor, min_exp = 1, 
  variance_filter = TRUE, 
  percentile=.6,
  vstransform = TRUE
)
se.newburgh.pp <- exp_preprocess(
    se.newburgh.viral, min_exp = 1, 
    variance_filter = TRUE, 
    percentile=.6, 
    vstransform = TRUE
)

# Infer GRNs
grn_poor <- exp2grn(
  exp = se.poor.pp,
  regulators = am_tf_poor
)

grn_newburgh <- exp2grn(
  exp = se.newburgh.pp,
  regulators = am_tf_newburgh
)

# View GRNs
grn_poor
grn_newburgh





