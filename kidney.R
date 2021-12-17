suppressMessages(library(stringr))
suppressMessages(library(reshape2))
suppressMessages(library(tidytext))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggthemes))
suppressMessages(library(gghalves))
suppressMessages(library(pheatmap))
suppressMessages(library(patchwork))
suppressMessages(library(SingleCellExperiment))
options(stringsAsFactors = FALSE)
suppressMessages(source("scover_helper_functions.R"))
# Amount of motifs
d <- 600
outdir <- "output/kidney/"
if(!dir.exists(outdir)) { dir.create(outdir) }
curr_sce <- readRDS("data/kidney/kidney_pooled_sce.RDS")
curr_colData <- as.data.frame(colData(curr_sce))
# Load data
category_means_aggregates <- read.csv("data/kidney/kidney_subfams_category_influence.csv", row.names = 1)
subfam_reproducibility <- read.csv("data/kidney/kidney_subfams_reproducibility.csv", row.names = 1)
category_means_aggregates <- category_means_aggregates[rownames(subfam_reproducibility)[subfam_reproducibility$Reproducibility >= 0.5],,drop=FALSE]
category_means_aggregates <- category_means_aggregates[rownames(category_means_aggregates) != "N/A",]
subfams_loo_df <- read.csv("data/kidney/kidney_subfams_influence.csv", row.names = 1)
subfams_loo_df <- subfams_loo_df[rownames(subfam_reproducibility)[subfam_reproducibility$Reproducibility >= 0.5],,drop=FALSE]
subfams_loo_df <- subfams_loo_df[rownames(subfams_loo_df) != "N/A",]
motif_family_annotation <- read.csv("data/kidney/kidney_motif_annotation.csv", row.names=1)
curr_category_annot <- data.frame(row.names=rownames(category_means_aggregates), 
                                  amount_motifs = rep("", nrow(category_means_aggregates)))
curr_category_annot$amount_motifs <- sapply(rownames(category_means_aggregates), 
                                            FUN=function(x){
                                              sum(motif_family_annotation$subfamily == x)
                                            })
colnames(curr_category_annot) <- "Motifs in cluster"
curr_category_annot$`Summed mean influence` <- rowMeans(subfams_loo_df)[rownames(curr_category_annot)]
# Fig 2a =====
heatmap_colors <- colorRampPalette(c("magenta", "black", "yellow"))(100)
motif_frequency_scores <- read.csv("data/kidney/kidney_motif_freq.csv",
                                    row.names = 1)
category_freq_means <- c()
for(mot_cat in rownames(category_means_aggregates)){
  curr_mots <- motif_family_annotation[motif_family_annotation$subfamily == mot_cat,]$motif
  category_freq_means <- c(category_freq_means, mean(motif_frequency_scores$freq[motif_frequency_scores$motif %in% curr_mots]))
}
curr_category_annot$`Motif frequency` <- category_freq_means
pheatmap(to_z(category_means_aggregates), color=heatmap_colors, 
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, 
         filename = paste0(outdir, "2a.pdf"),
         useDingbats=FALSE,
         annotation_row = curr_category_annot)
