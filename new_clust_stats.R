library(dplyr)
library(scales)
#install.packages("GGally")
library(GGally)
library(ggplot2)
library(FactoMineR)
# install.packages("factoextra")
library(factoextra)

# ----------
# Functions
# ----------
heaD <- function(x, ...){
  head(x, ...)
}

my_dens <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_density(..., alpha = 0.6, color = NA) 
}

numeric_cols <- function(data){
  data[, sapply(data, is.numeric)]
}


round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

round_place <- 3

# ----------
# Load data
# ----------
stats_file <- "~/Documents/SNP_analysis/all_lins/csv/New_clusters_stats_table.csv"
stats_table <- read.csv(stats_file, stringsAsFactors = F)

# ---------------
# Clean & select
# ---------------
# Set ok to factor
stats_table$ok <- as.factor(stats_table$ok)
# Round
stats_table <- round_df(stats_table, round_place)

# data <- select(stats_table, branch_ratio, abs_log10_branch_test_stat, med_ratio, abs_log10_t_test_stat, n_SNP_0.99, n_SNP_1, ok)
data <- select(stats_table, branch_ratio, branch_wilcox_stat, med_ratio, dist_mwu_stat, n_SNP_1, ok)
# data <- select(stats_table, log10_branch_p_val, dist_log10_p, n_SNP_1, ok)
# data <- select(stats_table, branch_cohen_d, dist_cohen_d, n_SNP_1, ok)

# data <- data.frame(scale(numeric_cols(data)))
# data$ok <- stats_table$ok

# -----
# Plot
# -----

# Setup
my_cols <- c("#E7B800", "#00AFBB")
Alpha <- 0.6

# Base R
# pairs(data, 
#       col = alpha(my_cols[stats_table$ok], Alpha),
#       pch = 19, cex = 2,
#       upper.panel=NULL)

# GGpairs
pairs <- ggpairs(data, mapping = aes(colour = ok) 
                 # , diag = list(continuous = my_dens)
                 )+
  theme_bw()

for(i in 1:pairs$nrow){
  for(j in 1:pairs$ncol){
    pairs[i, j] <- pairs[i, j] +
      scale_fill_manual(values = alpha(my_cols, Alpha))+
      scale_color_manual(values = alpha(my_cols, Alpha))
  }
}
pairs

# ----
# PCA
# ----

data_pca <- PCA(subset(data, select = -ok), graph = F)

fviz_pca_biplot(data_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = data$ok, # colour by groups
             palette = my_cols,
             label = "var",
             repel = TRUE,
             addEllipses = TRUE, 
             legend.title = "Groups"
)
  

























