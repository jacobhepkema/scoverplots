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
outdir <- "output/snare/"
if(!dir.exists(outdir)) { dir.create(outdir) }
curr_sce_ad <- readRDS("data/ad/ad_pooled_sce.RDS")
curr_colData_ad <- as.data.frame(colData(curr_sce_ad))
category_means_aggregates_z_ad <- read.csv("data/ad/ad_subfams_category_influence.csv", row.names = 1)
subfam_reproducibility_ad <- read.csv("data/ad/ad_subfams_reproducibility.csv", row.names = 1)
category_means_aggregates_z_ad <- category_means_aggregates_z_ad[subfam_reproducibility_ad$Reproducibility >= 0.5,]
category_means_aggregates_z_ad <- category_means_aggregates_z_ad[rownames(category_means_aggregates_z_ad) != "N/A",]
subfams_loo_df_ad <- read.csv("data/ad/ad_subfams_influence.csv", row.names = 1)
subfams_loo_df_ad <- subfams_loo_df_ad[subfam_reproducibility_ad$Reproducibility >= 0.5,]
subfams_loo_df_ad <- subfams_loo_df_ad[rownames(subfams_loo_df_ad) != "N/A",]
motif_family_annotation_ad <- read.csv("data/ad/ad_motif_annotation.csv", row.names=1)
curr_category_annot_ad <- data.frame(row.names=rownames(category_means_aggregates_z_ad), 
                                     amount_motifs = rep("", nrow(category_means_aggregates_z_ad)))
curr_category_annot_ad$amount_motifs <- sapply(rownames(category_means_aggregates_z_ad), 
                                               FUN=function(x){
                                                 sum(motif_family_annotation_ad$subfamily == x)
                                               })
colnames(curr_category_annot_ad) <- "Motifs in cluster"
curr_category_annot_ad$`Summed mean influence` <- rowMeans(subfams_loo_df_ad)[rownames(curr_category_annot_ad)]
motif_frequency_scores_ad <- read.csv("data/ad/ad_motif_freq.csv",
                                   row.names = 1)
category_freq_means_ad <- c()
for(mot_cat in rownames(curr_category_annot_ad)){
  curr_mots <- motif_family_annotation_ad[motif_family_annotation_ad$subfamily == mot_cat,]$motif
  category_freq_means_ad <- c(category_freq_means_ad, mean(motif_frequency_scores_ad$freq[motif_frequency_scores_ad$motif %in% curr_mots]))
}
curr_category_annot_ad$`Motif frequency` <- category_freq_means_ad
# Fig 2a =====
category_means_aggregates_z_ad <- category_means_aggregates_z_ad[,order(colnames(category_means_aggregates_z_ad))]
colnames(category_means_aggregates_z_ad) <- sapply(colnames(category_means_aggregates_z_ad), FUN=function(x){str_replace_all(x, "\\.", "/")}) 
category_means_aggregates_z_ad <- category_means_aggregates_z_ad[order(sapply(rownames(category_means_aggregates_z_ad), 
             FUN=function(x){
               str_split(x, ":")[[1]][2]
             })),]
heatmap_colors <- colorRampPalette(c("magenta", "black", "yellow"))(100)
pheatmap(to_z(category_means_aggregates_z_ad), color=heatmap_colors, 
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, cluster_cols = FALSE,
         cluster_rows = FALSE,
         filename = paste0(outdir, "6a_ad.pdf"),
         useDingbats=FALSE,
         annotation_row = curr_category_annot_ad)
category_means_aggregates_z_ad_ex <- to_z(category_means_aggregates_z_ad)[,2:8]
curr_annot_col <- data.frame(row.names=colnames(category_means_aggregates_z_ad_ex),
                             Category=c("L2/3",
                                        "L3/4/5",
                                        "L4",
                                        "L4/5",
                                        "L5",
                                        "L5/6",
                                        "L6"))
curr_annot_color <- list(
  Category = stata_pal()(length(unique(curr_annot_col$Category)))
)
names(curr_annot_color$Category) <- unique(curr_annot_col$Category)
pheatmap(category_means_aggregates_z_ad_ex, color=heatmap_colors, 
         angle_col = 45, cellwidth = 12, cellheight = 9,
         border_color = NA, cluster_cols = FALSE, cluster_rows = TRUE,
         filename = paste0(outdir, "6a_ad_ExOnly.pdf"),
         useDingbats=FALSE, annotation_col = curr_annot_col,
         annotation_colors = curr_annot_color,
         annotation_row = curr_category_annot_ad)
all_LOO_mat_ad <- as.matrix(read.csv("data/ad/ad_influence.csv.gz", row.names = 1))
all_LOO_mat_selection_ad <- all_LOO_mat_ad[rownames(motif_family_annotation_ad[motif_family_annotation_ad$subfamily %in% rownames(subfams_loo_df_ad),]),]
LOO_mat_melted_ad <- melt(all_LOO_mat_selection_ad) %>% magrittr::set_colnames(c("Motif", "Pool", "Weight"))
LOO_mat_melted_ad$Celltype <- curr_colData_ad[LOO_mat_melted_ad$Pool,1]
LOO_mat_melted_ad$`Motif cluster annotation` <- motif_family_annotation_ad[as.character(LOO_mat_melted_ad$Motif),]$subfamily
LOO_mat_melted_ad$Category <- curr_colData_ad[LOO_mat_melted_ad$Pool,]$cell_type_category
all_LOO_mat_selection_aggregates_ad <- matrix(nrow=length(unique(LOO_mat_melted_ad$`Motif cluster annotation`)), 
                                           ncol=ncol(all_LOO_mat_selection_ad))
for(i in 1:nrow(all_LOO_mat_selection_aggregates_ad)){
  curr_motif_family <- unique(LOO_mat_melted_ad$`Motif cluster annotation`)[i]
  all_LOO_mat_selection_aggregates_ad[i,] <- colSums(all_LOO_mat_selection_ad[rownames(all_LOO_mat_selection_ad) 
                                                                        %in% rownames(motif_family_annotation_ad[motif_family_annotation_ad$subfamily ==
                                                                                                                curr_motif_family,,drop=FALSE]),])
}
rownames(all_LOO_mat_selection_aggregates_ad) <- unique(LOO_mat_melted_ad$`Motif cluster annotation`)
colnames(all_LOO_mat_selection_aggregates_ad) <- colnames(all_LOO_mat_selection_ad)
all_LOO_mat_selection_aggregates_melt_ad <- melt(all_LOO_mat_selection_aggregates_ad) %>% 
  magrittr::set_colnames(c("Motif cluster annotation", "Pool", "Aggregate of motif weights"))
all_LOO_mat_selection_aggregates_melt_ad$Celltype <- curr_colData_ad[all_LOO_mat_selection_aggregates_melt_ad$Pool,1]
all_LOO_mat_selection_aggregates_melt_ad$Category <- curr_colData_ad[all_LOO_mat_selection_aggregates_melt_ad$Pool,]$cell_type_category
# Order:
all_LOO_mat_selection_aggregates_melt_ad <- all_LOO_mat_selection_aggregates_melt_ad[all_LOO_mat_selection_aggregates_melt_ad$`Motif cluster annotation` != "NA",]
all_LOO_mat_selection_aggregates_melt_ad$`Motif cluster annotation` <- factor(x=as.character(all_LOO_mat_selection_aggregates_melt_ad$`Motif cluster annotation`),
                                                                           levels=unique(as.character(all_LOO_mat_selection_aggregates_melt_ad$`Motif cluster annotation`))[order(unique(as.character(all_LOO_mat_selection_aggregates_melt_ad$`Motif cluster annotation`)))])
graphics.off()
ggplot(all_LOO_mat_selection_aggregates_melt_ad, 
       aes(x=Category, y=`Aggregate of motif weights`, 
           fill=Category)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_x_reordered() +
  scale_fill_stata() + theme_bw(base_size=14) + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x = "Category", y="Aggregate of motif influence scores") +
  theme_Nice() +
  facet_wrap(~`Motif cluster annotation`, scales="free", ncol = 6)
ggsave(paste0(outdir, "/ad_supp.pdf"),
       width=14,height=8, useDingbats=FALSE)
expression_hits_df_ad <- read.csv("data/ad/ad_expression_corrs.csv", header=TRUE)
expression_hits_df_ad$expression_2 <- expression_hits_df_ad$mean_expression
expression_hits_df_ad$expression_2[expression_hits_df_ad$pval > 0.05] <- NA
expression_hits_df_ad$tf_labels <- NA
top_n <- 3
for(i in 1:length(unique(expression_hits_df_ad$family))){
  curr_cluster <- unique(expression_hits_df_ad$family)[i]
  curr_corr_tf_df <- expression_hits_df_ad[expression_hits_df_ad$family == curr_cluster,,drop=FALSE]
  curr_corr_tf_df <- curr_corr_tf_df[order(abs(curr_corr_tf_df$correlation), decreasing = TRUE),,drop=FALSE]
  annot_tfs <- curr_corr_tf_df$gene[1:(ifelse(nrow(curr_corr_tf_df) < top_n, nrow(curr_corr_tf_df), top_n))]
  expression_hits_df_ad$tf_labels[expression_hits_df_ad$family == curr_cluster &
                                 expression_hits_df_ad$gene %in% annot_tfs] <- 
    expression_hits_df_ad$gene[expression_hits_df_ad$family == curr_cluster &
                              expression_hits_df_ad$gene %in% annot_tfs]
}
expression_hits_df_ad$tf_labels[expression_hits_df_ad$pval > 0.05] <- ""
graphics.off()
ggplot(expression_hits_df_ad, aes(x=family, y=correlation, color=expression_2)) +
  geom_jitter(width = 0) + 
  geom_label_repel(label=expression_hits_df_ad$tf_labels, size=4) +
  theme_bw(base_size=14) + 
  theme_Nice() + theme(legend.position = "right") + 
  labs(x="Motif cluster name", 
       y="Spearman R", color="Mean TF expression across pools") 
ggsave(filename=paste0(outdir, "/6c.pdf"), 
       width = 9, height=6, useDingbats=FALSE)