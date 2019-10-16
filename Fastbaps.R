
# Fastbaps 

# Read in packages:
Packages <- "https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Packages.R"
source(Packages)

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
#                                       Functions
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

heaD <- function(x,...){
  head(x, ...)
}

h_col <- function(x, nrow = 6){
  head(x[, 1:10], nrow) 
}

# # Set and reset plots
reset_plot <- function(){
  par(mfrow = c(1, 1))
}

Plot_dims <- function(Rows, Cols){
  par(mfrow=c(Rows, Cols))
}

N_clusters <- function(Tree_in, Height, Min_clust_sz){
  length(unique(cutreeDynamic(as.hclust(Tree_in), 
                              cutHeight = Height,
                              minClusterSize=Min_clust_sz
  )))
}

Numeric <- function(x){
  x[, sapply(x, is.numeric)]
}

# Function to take out outliers from each group (lineage in this case) -
# adapt for your group of interest
No_outliers <- function(Data, SDs){
  # SDs <- 3.5 # How many standard devs from mean to take out
  Out <- function(Data, SDs){
    Vals <-c(mean(Data)+(sd(Data)*SDs),mean(Data)-(sd(Data)*SDs))
    # abs(Data - abs(mean(Data) - Val)) > Val
    Data > max(Vals) | Data < min(Vals)
  }
  
  # Need to order or indexing won't work
  Data <- Data[order(Data$main_lineage, Data$sub_lineage), ] 
  Outliers <- unique(c(
    as.vector(which(unlist(sapply(split(Data$V1, as.character(Data$sub_lineage)),
                                  function(x){Out(x, SDs)})))), 
    as.vector(which(unlist(sapply(split(Data$V2, as.character(Data$sub_lineage)),
                                  function(x){Out(x, SDs)}))))
  )) # unique(c(
  
  Data <- Data[-(Outliers), ]
}

# Function to subset rows of data by group (sub_lineage in this case),
# first proportional to amount of data in the group vs other groups, then
# retaining x% of rows (Retain arg - set to 0-1). Group arg - srtring of group col names
Sub_group <- function(x, Group, Retain){
  N_each_group <- plyr::count(x, eval(substitute(Group))) # Nb need library(dplyr)
  N_each_group$retain <- round(N_each_group$freq - 
                                 (N_each_group$freq*
                                    (N_each_group$freq/sum(N_each_group$freq)) ))
  
  N_each_group$retain <- round(N_each_group$retain*Retain)
  
  x_split <- split(x, as.character(x[,Group]))
  do.call("rbind", lapply(seq(x_split), function(i){
    x_split[[i]][sample(1:nrow(x_split[[i]]), 
                        N_each_group[N_each_group[, Group] == N_each_group[i,Group],
                                     "retain"]) , ]
  }))
}

Read_files <- function(File_names, Header = F, Sep = " ", Dec = "."){
  data.frame(do.call("rbind", lapply(File_names, 
                                     fread, header = Header, 
                                     sep = Sep, dec = Dec
  )))
}

# Function to get hierachy of lineage names. 
# e.g. get all lineages belonging to "1" will pull "1.1", "1.1.2", "1.2" etc. 
# i.e. pull all lineage names that contain the first element of 'subj'. 
# Returns INDICES of heirarchy in a list of lineages, so use like this:
# Lineage_vector[Tree("1.1", Lineage_vector)]
Tree <- function(Subj, Lins){
  N_subj <- nchar(Subj)
  Subj_split <- strsplit(Lins, split = "")
  Lins_test <- unlist(lapply(Subj_split, function(x){ paste(x[1:N_subj], collapse = "") }) )
  if (N_subj == 1){
    return(which(Lins == Subj))
  }else{
    return(grep(Subj, Lins_test))
  }
}

Optimise_prior <- function(sparse.data, grid.interval = c(5e-04, 10),
                            type = "optimise.baps", 
                            hc.method = "ward", n.cores = 1){
  
  if (!is.list(sparse.data)) 
    stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if (!(class(sparse.data$snp.matrix) == "dgCMatrix")) 
    stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if (!is.numeric(sparse.data$consensus)) 
    stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if (!is.matrix(sparse.data$prior)) 
    stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if (!(type %in% c("optimise.symmetric", "symmetric", "hc", 
                    "optimise.baps", "baps"))) 
    stop("Invalid value for type. Must be one of 'optimise.symmetric','symmetric', 'hc', 'optimise.baps' or 'baps'")
  if (!all(grid.interval > 0)) 
    stop("grid values must greater than 0")
  if (!(hc.method %in% c("ward", "genie"))) 
    stop("Invalid hc.method!")
  h <- fastbaps:::get_hclust(sparse.data, TRUE, hc.method, n.cores)
  if (type == "hc") {
    initial.prior <- matrix(c(rep(ncol(sparse.data$snp.matrix), 
                                  nrow(sparse.data$snp.matrix)), rowSums(sparse.data$snp.matrix == 
                                                                           1), rowSums(sparse.data$snp.matrix == 2), rowSums(sparse.data$snp.matrix == 
                                                                                                                               3), rowSums(sparse.data$snp.matrix == 4)), nrow = 5, 
                            byrow = TRUE)
    initial.prior[1, ] <- initial.prior[1, ] - colSums(initial.prior[2:5, 
                                                                     ])
    initial.prior <- initial.prior + 1
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
    initial.prior <- ceiling(initial.prior * 1000)/1000
  }
  else if (type == "optimise.baps") {
    initial.prior <- matrix(c(rep(ncol(sparse.data$snp.matrix), 
                                  nrow(sparse.data$snp.matrix)), rowSums(sparse.data$snp.matrix == 
                                                                           1), rowSums(sparse.data$snp.matrix == 2), rowSums(sparse.data$snp.matrix == 
                                                                                                                               3), rowSums(sparse.data$snp.matrix == 4)), nrow = 5, 
                            byrow = TRUE)
    initial.prior[1, ] <- initial.prior[1, ] - colSums(initial.prior[2:5, 
                                                                     ])
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
    initial.prior <- initial.prior > 0
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
  }
  else if (type == "baps") {
    initial.prior <- matrix(c(rep(ncol(sparse.data$snp.matrix), 
                                  nrow(sparse.data$snp.matrix)), rowSums(sparse.data$snp.matrix == 
                                                                           1), rowSums(sparse.data$snp.matrix == 2), rowSums(sparse.data$snp.matrix == 
                                                                                                                               3), rowSums(sparse.data$snp.matrix == 4)), nrow = 5, 
                            byrow = TRUE)
    initial.prior[1, ] <- initial.prior[1, ] - colSums(initial.prior[2:5, 
                                                                     ])
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
    initial.prior <- initial.prior > 0
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
  }
  else {
    initial.prior <- matrix(1, nrow = nrow(sparse.data$prior), 
                            ncol = ncol(sparse.data$prior))
  }
  sparse.data$prior <- initial.prior
  if ((type == "baps") || (type == "symmetric")) 
    return(sparse.data)
  opt <- stats::optimise(fastbaps:::calc_prior_prob, grid.interval, sparse.data, 
                         initial.prior, h, maximum = TRUE, tol = 0.001)
  cc <- round(opt$maximum, digits = 3)
  sparse.data$prior <- initial.prior * cc
  sparse.data$prior[sparse.data$prior < 0.001] <- 0.001
  if (any(abs(grid.interval - cc) < 0.005)) {
    warning("Inferred hyperparameter is very close to interval boundries! Consider changing the interval.")
  }
  # print(paste("Optimised hyperparameter:", cc))
  sparse.data$prior.type = "optimised"
  sparse.data$hclust <- h
  sparse.data$opt_hyperparam <- cc
  return(sparse.data)
}

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
#                                         Setup
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

WD <- "/home/gnapier/Documents/SNP_analysis/"
setwd(WD)

Lin_nums <- 1:4

# Plot params
Line_width <- 0.2
Collapse <- 0.0007
Tick_sz <- 0.7
Round <- 5
Alpha <- 0.4
Ncol <- 2
Width <- 100
Height <- Width*0.4
Scale <- 3
Units <- "mm"
Labs <- labs(x = "PC1", y = "PC2")
par(mar = c(0.5, 1.5, 0.5, 1.5))
# Get lineage names of interest
Lin_names <- c("1", "2", "3", "4")
Lin_names_first_dec <- c("1.1", "1.2", "2.1", "2.2", "3.0", "3.1", "4.0", "4.1", "4.2",                          "4.3", "4.4", "4.5", "4.6", "4.7", "4.8", "4.9")
Lin_names_all <- c(Lin_names, Lin_names_first_dec)
# Lin_names <- subset(Lin_names, subset = !(Lin_names %in% c("5", "6", "7", "BOV")))

# Read in sample names and (sub)lineage to merge with PCA data
# Lins <- "lin_1/txt/lin_1_samps_lins.txt"
# Lins <- data.frame(fread(Lins, header = T, sep = "\t", dec = "."))
Lins <- Read_files("all_lins/txt/lineages_clean.txt", Header = T, Sep = "\t")

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
#                                       Fastbaps
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

#Paper: https://www.biorxiv.org/content/biorxiv/early/2018/10/26/454355.full.pdf
#Github: https://github.com/gtonkinhill/fastbaps

# https://github.com/gtonkinhill/fastbaps

# Read in fasta files & do fastbaps
Fasta_files <- sprintf("lin_%s/results/lin_%s_results.filt.snps.fa", 
                       Lin_names_first_dec, Lin_names_first_dec)
sparse.data <- lapply(Fasta_files, import_fasta_sparse_nt, prior = "baps")
# sparse.data <- lapply(sparse.data, Optimise_prior, type = "optimise.symmetric")

Tree_files <- sprintf("lin_%s/results/lin_%s_results.filt.snps.fa.raxml.support",
                     Lin_names_first_dec, Lin_names_first_dec)

Trees <- lapply(Tree_files, function(x){
  tryCatch({phytools::read.newick(x)}, error = function(e){})
  })

Trees_rooted <- lapply(Trees, function(x){
  tryCatch({phytools::midpoint.root(x)}, error = function(e){})
})

clusters_orig <- lapply(seq(sparse.data), function(i){
  tryCatch({best_baps_partition(sparse.data[[i]], Trees_rooted[[i]])}, 
           error = function(e){})
  })

# baps.hc <- lapply(sparse.data, fast_baps, quiet = T)

# clusters <- lapply(seq(sparse.data), function(i){
#   best_baps_partition(sparse.data[[i]], as.phylo(baps.hc[[i]]), quiet = T)
# })

# Put results into data frames
clusters <- lapply(clusters_orig, function(x){
  data.frame(id = row.names(data.frame(x)), clusters = x)
})

# Display results in markdown
FB_df <- data.frame(Lin = Lin_names_first_dec,
                    N_clusters = do.call("rbind",
                                         lapply(clusters, 
                                                function(x){ 
                                                  tryCatch({max(unique(Numeric(x)))}, 
                                                           error = function(e){"-"})
                                                  })) )

FB_df
save(FB_df, file = "all_lins/R_data/Fastbap_phylo_n_clusts.rda")

# Put all samples and cluster prediction together
All_clusters <- do.call("rbind", clusters)

# # Read in lineages_master.txt and merge with clusters
Lins_master <- read.delim("all_lins/txt/lineages_master.txt",
                          header = TRUE, sep = "\t", dec = ".")

All_clusters <- merge(Lins_master, All_clusters, by.x = "sample", by.y = "id")
# Sort to keep neat
All_clusters <- All_clusters[with(All_clusters, order(main_lineage,
                                                      first_dec,
                                                      sub_lineage,
                                                      clusters)), ]

Alpha <- 0.5

All_clusters_split <- split(All_clusters, All_clusters$first_dec)
  
Cluster_col_lookup <- do.call("rbind", lapply(All_clusters_split, function(x) {
  data.frame(clusters = unique(x$clusters),
             col = rainbow(length(unique(x$clusters)), alpha = Alpha))
}))
 
Cluster_col_lookup$new_sublin <- row.names(Cluster_col_lookup)

# Create csv of colour keys (second csv file to pass to csv2itol.py). 
n <- length(unique(All_clusters$clusters))
write.csv(data.frame(Type = rep("clusters", n),
                     Value = sort(unique(All_clusters$clusters)),
                     Colour = rainbow(n, alpha = Alpha) ),
          file = "all_lins/itol/fastbaps_key.csv", row.names = F, quote = F)

All_clusters <- rename(All_clusters, id = sample) # Need to rename for Python script
write.csv(All_clusters[, c("id", "clusters")],
          file = "all_lins/itol/fastbap_clusts.csv", 
          row.names = F, quote = F)

# Run python script to convert above two csvs to itol file - should 
# then be in all_lins/itol/ folder
Python_path <- "/home/gnapier/anaconda3/bin/python"
Python_script_path <- sprintf("%s/python_scripts/csv2itol.py", WD)
Fastbap_results <- sprintf("%s/all_lins/itol/fastbap_clusts.csv", WD)
Key_path <- sprintf("%s/all_lins/itol/fastbaps_key.csv", WD)
Save_as <- sprintf("%s/all_lins/itol/fastbap_phylo_", WD)
system(paste(Python_path, Python_script_path, Fastbap_results, 
             Key_path, Save_as))

# system("which python", intern = T)
# See all_lins/itol/cluster.meta.itol.txt file


N <- 1

plot.df <- data.frame(id = Trees_rooted[[N]]$tip.label, 
                      fastbaps = clusters_orig[[N]], 
                      stringsAsFactors = FALSE)

gg <- ggtree(Trees_rooted[[N]])
f2 <- facet_plot(gg, panel = "fastbaps", 
                 data = plot.df, geom = geom_tile, aes(x = fastbaps), 
                 color = "blue")
f2


















