suppressMessages(library(stringr))
suppressMessages(library(reshape2))
suppressMessages(library(tidytext))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggthemes))
suppressMessages(library(gghalves))
suppressMessages(library(ggseqlogo))
suppressMessages(library(pheatmap))
suppressMessages(library(patchwork))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(umap))
options(stringsAsFactors = FALSE)
suppressMessages(source("scover_helper_functions.R"))
# Amount of motifs
d <- 600
outdir <- "output/tm/"
if(!dir.exists(outdir)) { dir.create(outdir) }
curr_sce <- readRDS("data/tm/tm_pooled_sce.RDS")
curr_colData <- as.data.frame(colData(curr_sce))
# Load data
category_means_aggregates <- read.csv("data/tm/tm_subfams_category_influence.csv", row.names = 1)
subfam_reproducibility <- read.csv("data/tm/tm_subfams_reproducibility.csv", row.names = 1)
category_means_aggregates <- category_means_aggregates[rownames(subfam_reproducibility)[subfam_reproducibility$Reproducibility >= 0.5],,drop=FALSE]
category_means_aggregates <- category_means_aggregates[rownames(category_means_aggregates) != "N/A",]
subfams_loo_df <- read.csv("data/tm/tm_subfams_influence.csv", row.names = 1)
subfams_loo_df <- subfams_loo_df[rownames(subfam_reproducibility)[subfam_reproducibility$Reproducibility >= 0.5],,drop=FALSE]
subfams_loo_df <- subfams_loo_df[rownames(subfams_loo_df) != "N/A",]
motif_family_annotation <- read.csv("data/tm/tm_motif_annotation.csv", row.names=1)
curr_category_annot <- data.frame(row.names=rownames(category_means_aggregates), 
                                  amount_motifs = rep("", nrow(category_means_aggregates)))
curr_category_annot$amount_motifs <- sapply(rownames(category_means_aggregates), 
                                            FUN=function(x){
                                              sum(motif_family_annotation$subfamily == x)
                                            })
colnames(curr_category_annot) <- "Motifs in cluster"
curr_category_annot$`Summed mean influence` <- rowMeans(subfams_loo_df)[rownames(curr_category_annot)]
# Fig 2a =====