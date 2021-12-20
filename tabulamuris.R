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
# Fig 4a =====
heatmap_colors <- colorRampPalette(c("magenta", "black", "yellow"))(100)
motif_frequency_scores <- read.csv("data/tm/tm_motif_freq.csv",
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
         filename = paste0(outdir, "4a.pdf"),
         useDingbats=FALSE,
         annotation_row = curr_category_annot)
# See object `curr_category_annot` for numbers
all_LOO_mat <- as.matrix(read.csv("data/tm/tm_influence.csv.gz", row.names = 1))
all_LOO_mat_selection <- all_LOO_mat[rownames(motif_family_annotation[motif_family_annotation$subfamily %in% rownames(subfams_loo_df),]),]
LOO_mat_melted <- melt(all_LOO_mat_selection) %>% magrittr::set_colnames(c("Motif", "Pool", "Weight"))
LOO_mat_melted$Celltype <- curr_colData[LOO_mat_melted$Pool,1]
LOO_mat_melted$`Motif cluster annotation` <- motif_family_annotation[as.character(LOO_mat_melted$Motif),]$subfamily
all_LOO_mat_selection_aggregates <- matrix(nrow=length(unique(LOO_mat_melted$`Motif cluster annotation`)), 
                                           ncol=ncol(all_LOO_mat_selection))
for(i in 1:nrow(all_LOO_mat_selection_aggregates)){
  curr_motif_family <- unique(LOO_mat_melted$`Motif cluster annotation`)[i]
  all_LOO_mat_selection_aggregates[i,] <- colSums(all_LOO_mat_selection[rownames(all_LOO_mat_selection) 
                                                                        %in% rownames(motif_family_annotation[motif_family_annotation$subfamily == curr_motif_family,]),])
}
rownames(all_LOO_mat_selection_aggregates) <- unique(LOO_mat_melted$`Motif cluster annotation`)
colnames(all_LOO_mat_selection_aggregates) <- colnames(all_LOO_mat_selection)
all_LOO_mat_selection_aggregates_melt <- read.csv("data/tm/tm_influence_aggregates.csv")
all_LOO_mat_selection_aggregates_melt$`Aggregate of motif weights` <- all_LOO_mat_selection_aggregates_melt$sum
all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` <- all_LOO_mat_selection_aggregates_melt$subfamily
all_LOO_mat_selection_aggregates_melt$Category <- all_LOO_mat_selection_aggregates_melt$cat
all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` <- factor(x=as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`),
                                                                           levels=unique(as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))[order(unique(as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)))])
all_LOO_mat_selection_aggregates_melt_dupl <- all_LOO_mat_selection_aggregates_melt
all_LOO_mat_selection_aggregates_melt_dupl$subfamily[all_LOO_mat_selection_aggregates_melt_dupl$subfamily %in% c("ETS/2:ETS", "ETS/1:ETS")] <- "ETS/1/2:ETS"
all_LOO_mat_selection_aggregates_melt_dupl$`Motif cluster annotation` <- all_LOO_mat_selection_aggregates_melt_dupl$subfamily
ggplot(all_LOO_mat_selection_aggregates_melt_dupl, 
       aes(x=reorder_within(Category,`Aggregate of motif weights`,`Motif cluster annotation`), y=`Aggregate of motif weights`, 
           fill=`Category`)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_x_reordered() +
  scale_fill_stata() + theme_bw(base_size=14) + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x = "Category", y="Aggregate of motif influence scores") +
  theme_Nice() +
  facet_wrap(~`Motif cluster annotation`, scales="free", ncol = 3)
ggsave(paste0(outdir, "/tm_supp.pdf"),
       width=7,height=9, useDingbats=FALSE)
graphics.off()
all_LOO_mat_selection_aggregates_melt$`Aggregate of motif influence scores` <- all_LOO_mat_selection_aggregates_melt$`Aggregate of motif weights`
pheatmap(cor(t(all_LOO_mat_selection_aggregates), method="spearman"), 
         color=heatmap_colors,
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, useDingbats=FALSE,
         filename = paste0(outdir, "tm_supp2.pdf"))
graphics.off()
all_motifs <- read_meme("data/tm/tm_motifs.meme")
# NR/1
ggplot() + geom_logo(all_motifs$`9_122`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_NR1.eps"), width=3.5,height=2)
# IRF/1
ggplot() + geom_logo(all_motifs$`3_433`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_IRF1.eps"), width=3.5,height=2)
# SPI
ggplot() + geom_logo(all_motifs$`9_437`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_SPI.eps"), width=3.5,height=2)
# ETS/2
ggplot() + geom_logo(all_motifs$`3_524`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_ETS2.eps"), width=3.5,height=2)
# ETS/1
ggplot() + geom_logo(all_motifs$`5_581`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_ETS1.eps"), width=3.5,height=2)
# KLF/SP/2
ggplot() + geom_logo(all_motifs$`0_222`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_KLF.eps"), width=3.5,height=2)
# ZFX
ggplot() + geom_logo(all_motifs$`2_477`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_ZFX.eps"), width=3.5,height=2)
# CREB/ATF/1
ggplot() + geom_logo(all_motifs$`3_52`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_CREB1.eps"), width=3.5,height=2)
# CREB/ATF/2
ggplot() + geom_logo(all_motifs$`9_104`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_CREB2.eps"), width=3.5,height=2)
# YY1
ggplot() + geom_logo(all_motifs$`7_8`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_YY1.eps"), width=3.5,height=2)
# eBOX
ggplot() + geom_logo(all_motifs$`3_266`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_EBOX.eps"), width=3.5,height=2)
# e2f
ggplot() + geom_logo(all_motifs$`0_487`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_E2F.eps"), width=3.5,height=2)
# tbx
ggplot() + geom_logo(all_motifs$`5_355`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir,"tm_TBX.eps"), width=3.5,height=2)
# GO term analysis: correlations of motif influence scores with expression of GO-term related genes:
GO_terms <- readLines("data/tm/go_scfind.tsv")
names(GO_terms) <- sapply(GO_terms, FUN=function(x){return(str_split(x, "\t")[[1]][1])})
GO_terms <- sapply(GO_terms, FUN=function(x){
  str_split(str_split(x, "\t")[[1]][2], ",")[[1]]
})
curr_GO_lengths <- c()
pb <- txtProgressBar(0, length(GO_terms), style=3)
for(i in 1:length(GO_terms)){
  curr_GO_name <- names(GO_terms)[i]
  curr_GO_genes <- GO_terms[[i]]
  curr_GO_genes <- unique(curr_GO_genes[curr_GO_genes %in% rownames(curr_sce)])
  curr_GO_length <- length(curr_GO_genes)
  curr_GO_lengths <- c(curr_GO_lengths, curr_GO_length)
  setTxtProgressBar(pb, i)
}
# Select GO terms with at least 1 gene and fewer than 50 genes in set (that are found
# in the current experiment)
GO_selection <- GO_terms[which(curr_GO_lengths < 50 & curr_GO_lengths > 0)]
curr_GO_lengths_selection <- curr_GO_lengths[curr_GO_lengths < 50 & curr_GO_lengths > 0]
GO_selection_corrs <- list()
GO_selection_pvals <- list()
pb <- txtProgressBar(0, length(GO_selection), style=3)
for(j in 1:length(GO_selection)){ # Correlate expression of GO term genes to motif influence scores. This can take a while
  curr_GO_name <- names(GO_selection)[j]
  curr_GO_genes <- GO_selection[[j]]
  curr_GO_genes <- unique(curr_GO_genes[curr_GO_genes %in% rownames(curr_sce)])
  curr_GO_expression <- colMeans(logcounts(curr_sce[curr_GO_genes,]))
  
  curr_GO_corrs <- c()
  curr_GO_pvals <- c()
  for(i in 1:length(unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))){
    curr_cluster_annot <- as.character(unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))[i]
    curr_melt_aggregates <- all_LOO_mat_selection_aggregates_melt[all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` == curr_cluster_annot,]
    curr_GO_corrs <- c(curr_GO_corrs, cor(curr_melt_aggregates$`Aggregate of motif influence scores`, curr_GO_expression, 
                                          method="spearman"))
    curr_GO_pvals <- c(curr_GO_pvals, p.adjust(cor.test(curr_melt_aggregates$`Aggregate of motif influence scores`, 
                                                        curr_GO_expression, 
                                                        method="spearman", 
                                                        exact = FALSE)$p.value, method="fdr"))
  }
  names(curr_GO_corrs) <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)
  names(curr_GO_pvals) <- names(curr_GO_corrs)
  GO_selection_corrs[[j]] <- curr_GO_corrs
  GO_selection_pvals[[j]] <- curr_GO_pvals
  setTxtProgressBar(pb, j)
}
names(GO_selection_corrs) <- names(GO_selection)
names(GO_selection_pvals) <- names(GO_selection)
GO_melt <- melt(GO_selection_corrs)
GO_melt$motif_family <- names(GO_selection_corrs[[1]])
GO_melt_p <- melt(GO_selection_pvals)
GO_melt$p_corr <- GO_melt_p$value
colnames(GO_melt) <- c("Spearman R", "GO term", "Motif cluster annotation", "Corrected p-value")
GO_melt_cast <- dcast(GO_melt, formula=`GO term`~`Motif cluster annotation`, value.var="Spearman R")
rownames(GO_melt_cast) <- GO_melt_cast[,1]
GO_melt_cast <- GO_melt_cast[,-c(1)]
GO_melt_cast_pval <- dcast(GO_melt, formula=`GO term`~`Motif cluster annotation`, value.var="Corrected p-value")
GO_melt_cast_pval <- as.data.frame(GO_melt_cast_pval)
rownames(GO_melt_cast_pval) <- GO_melt_cast_pval[,1]
GO_melt_cast_pval <- GO_melt_cast_pval[,-c(1)]
bottom_top_quantiles <- quantile(as.numeric(unlist(GO_melt_cast)), c(0.01, 0.99))
has_no_significant <- apply(GO_melt_cast, 1, FUN=function(x){
  return((sum(x < bottom_top_quantiles[1]) + 
            sum(x > bottom_top_quantiles[2])) 
         == 0)
})
sum(!has_no_significant)
# Cluster annotation
k <- 3
clust_annot <- cutree(hclust(dist(GO_melt_cast[!has_no_significant,])), k=k)
GO_row_annot <- data.frame(row.names=names(clust_annot),
                           "GO term cluster"=as.character(clust_annot))
colnames(GO_row_annot) <- "GO term cluster"
GO_row_annot_col <- list(
  "GO term cluster" = stata_pal()(k)
)
names(GO_row_annot_col$`GO term cluster`) <- as.character(1:k)
# Fig 5b =====
pheatmap::pheatmap(GO_melt_cast[!has_no_significant,],
                   show_rownames = FALSE, border_color = NA,
                   color=heatmap_colors, cellwidth = 12,
                   angle_col = 45, cellheight = .14, 
                   annotation_row = GO_row_annot, 
                   annotation_colors = GO_row_annot_col,
                   filename = paste0(outdir, "/Fig5b.pdf"),
                   useDingbats=FALSE)
graphics.off()
# # The annotation for GO terms was added in Illustrator, and it was found with these commands:
# significant_cluster <- cutree(hclust(dist(GO_melt_cast[!has_no_significant,])), k=k)
# cat(names(significant_cluster[significant_cluster == 1])[sample(length(significant_cluster[significant_cluster == 1]), size=10)], sep="\n")
# cat(names(significant_cluster[significant_cluster == 2])[sample(length(significant_cluster[significant_cluster == 2]), size=10)], sep="\n")
# cat(names(significant_cluster[significant_cluster == 3])[sample(length(significant_cluster[significant_cluster == 3]), size=10)], sep="\n")

