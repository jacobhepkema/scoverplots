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
curr_sce_p0 <- readRDS("data/p0/p0_pooled_sce.RDS")
curr_colData_p0 <- as.data.frame(colData(curr_sce_p0))
category_means_aggregates_z_p0 <- read.csv("data/p0/p0_subfams_category_influence.csv", row.names = 1)
subfam_reproducibility_p0 <- read.csv("data/p0/p0_subfams_reproducibility.csv", row.names = 1)
category_means_aggregates_z_p0 <- category_means_aggregates_z_p0[subfam_reproducibility_p0$Reproducibility >= 0.5,]
category_means_aggregates_z_p0 <- category_means_aggregates_z_p0[rownames(category_means_aggregates_z_p0) != "N/A",]
subfams_loo_df_p0 <- read.csv("data/p0/p0_subfams_influence.csv", row.names = 1)
subfams_loo_df_p0 <- subfams_loo_df_p0[subfam_reproducibility_p0$Reproducibility >= 0.5,]
subfams_loo_df_p0 <- subfams_loo_df_p0[rownames(subfams_loo_df_p0) != "N/A",]
motif_family_annotation_p0 <- read.csv("data/p0/p0_motif_annotation.csv", row.names=1)
curr_category_annot_p0 <- data.frame(row.names=rownames(category_means_aggregates_z_p0), 
                                     amount_motifs = rep("", nrow(category_means_aggregates_z_p0)))
curr_category_annot_p0$amount_motifs <- sapply(rownames(category_means_aggregates_z_p0), 
                                               FUN=function(x){
                                                 sum(motif_family_annotation_p0$subfamily == x)
                                               })
colnames(curr_category_annot_p0) <- "Motifs in cluster"
curr_category_annot_p0$`Summed mean influence` <- rowMeans(subfams_loo_df_p0)[rownames(curr_category_annot_p0)]
motif_frequency_scores_p0 <- read.csv("data/p0/p0_motif_freq.csv",
                                   row.names = 1)
category_freq_means_p0 <- c()
for(mot_cat in rownames(curr_category_annot_p0)){
  curr_mots <- motif_family_annotation_p0[motif_family_annotation_p0$subfamily == mot_cat,]$motif
  category_freq_means_p0 <- c(category_freq_means_p0, mean(motif_frequency_scores_p0$freq[motif_frequency_scores_p0$motif %in% curr_mots]))
}
curr_category_annot_p0$`Motif frequency` <- category_freq_means_p0
# Fig 2a =====
category_means_aggregates_z_p0 <- category_means_aggregates_z_p0[,order(colnames(category_means_aggregates_z_p0))]
colnames(category_means_aggregates_z_p0) <- sapply(colnames(category_means_aggregates_z_p0), FUN=function(x){str_replace_all(x, "\\.", "/")}) 
category_means_aggregates_z_p0 <- category_means_aggregates_z_p0[order(sapply(rownames(category_means_aggregates_z_p0), 
             FUN=function(x){
               str_split(x, ":")[[1]][2]
             })),]
heatmap_colors <- colorRampPalette(c("magenta", "black", "yellow"))(100)
pheatmap(to_z(category_means_aggregates_z_p0), color=heatmap_colors, 
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, cluster_cols = FALSE,
         cluster_rows = FALSE,
         filename = paste0(outdir, "6a_P0.pdf"),
         useDingbats=FALSE,
         annotation_row = curr_category_annot_p0)
category_means_aggregates_z_p0_ex <- to_z(category_means_aggregates_z_p0)[,2:8]
curr_annot_col <- data.frame(row.names=colnames(category_means_aggregates_z_p0_ex),
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
pheatmap(category_means_aggregates_z_p0_ex, color=heatmap_colors, 
         angle_col = 45, cellwidth = 12, cellheight = 9,
         border_color = NA, cluster_cols = FALSE, cluster_rows = TRUE,
         filename = paste0(outdir, "6a_P0_ExOnly.pdf"),
         useDingbats=FALSE, annotation_col = curr_annot_col,
         annotation_colors = curr_annot_color,
         annotation_row = curr_category_annot_p0)
all_LOO_mat_p0 <- as.matrix(read.csv("data/p0/p0_influence.csv.gz", row.names = 1))
all_LOO_mat_selection_p0 <- all_LOO_mat_p0[rownames(motif_family_annotation_p0[motif_family_annotation_p0$subfamily %in% rownames(subfams_loo_df_p0),]),]
LOO_mat_melted_p0 <- melt(all_LOO_mat_selection_p0) %>% magrittr::set_colnames(c("Motif", "Pool", "Weight"))
LOO_mat_melted_p0$Celltype <- curr_colData_p0[LOO_mat_melted_p0$Pool,1]
LOO_mat_melted_p0$`Motif cluster annotation` <- motif_family_annotation_p0[as.character(LOO_mat_melted_p0$Motif),]$subfamily
LOO_mat_melted_p0$Category <- curr_colData_p0[LOO_mat_melted_p0$Pool,]$cell_type_category
all_LOO_mat_selection_aggregates_p0 <- matrix(nrow=length(unique(LOO_mat_melted_p0$`Motif cluster annotation`)), 
                                           ncol=ncol(all_LOO_mat_selection_p0))
for(i in 1:nrow(all_LOO_mat_selection_aggregates_p0)){
  curr_motif_family <- unique(LOO_mat_melted_p0$`Motif cluster annotation`)[i]
  all_LOO_mat_selection_aggregates_p0[i,] <- colSums(all_LOO_mat_selection_p0[rownames(all_LOO_mat_selection_p0) 
                                                                        %in% rownames(motif_family_annotation_p0[motif_family_annotation_p0$subfamily ==
                                                                                                                curr_motif_family,,drop=FALSE]),])
}
rownames(all_LOO_mat_selection_aggregates_p0) <- unique(LOO_mat_melted_p0$`Motif cluster annotation`)
colnames(all_LOO_mat_selection_aggregates_p0) <- colnames(all_LOO_mat_selection_p0)
all_LOO_mat_selection_aggregates_melt_p0 <- melt(all_LOO_mat_selection_aggregates_p0) %>% 
  magrittr::set_colnames(c("Motif cluster annotation", "Pool", "Aggregate of motif weights"))
all_LOO_mat_selection_aggregates_melt_p0$Celltype <- curr_colData_p0[all_LOO_mat_selection_aggregates_melt_p0$Pool,1]
all_LOO_mat_selection_aggregates_melt_p0$Category <- curr_colData_p0[all_LOO_mat_selection_aggregates_melt_p0$Pool,]$cell_type_category
# Order:
all_LOO_mat_selection_aggregates_melt_p0 <- all_LOO_mat_selection_aggregates_melt_p0[all_LOO_mat_selection_aggregates_melt_p0$`Motif cluster annotation` != "NA",]
all_LOO_mat_selection_aggregates_melt_p0$`Motif cluster annotation` <- factor(x=as.character(all_LOO_mat_selection_aggregates_melt_p0$`Motif cluster annotation`),
                                                                           levels=unique(as.character(all_LOO_mat_selection_aggregates_melt_p0$`Motif cluster annotation`))[order(unique(as.character(all_LOO_mat_selection_aggregates_melt_p0$`Motif cluster annotation`)))])
graphics.off()
ggplot(all_LOO_mat_selection_aggregates_melt_p0, 
       aes(x=Category, y=`Aggregate of motif weights`, 
           fill=Category)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_x_reordered() +
  scale_fill_stata() + theme_bw(base_size=14) + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x = "Category", y="Aggregate of motif influence scores") +
  theme_Nice() +
  facet_wrap(~`Motif cluster annotation`, scales="free", ncol = 6)
ggsave(paste0(outdir, "/p0_supp.pdf"),
       width=14,height=8, useDingbats=FALSE)
expression_hits_df_p0 <- read.csv("data/p0/p0_expression_corrs.csv", header=TRUE)
expression_hits_df_p0$expression_2 <- expression_hits_df_p0$mean_expression
expression_hits_df_p0$expression_2[expression_hits_df_p0$pval > 0.05] <- NA
expression_hits_df_p0$tf_labels <- NA
top_n <- 3
for(i in 1:length(unique(expression_hits_df_p0$family))){
  curr_cluster <- unique(expression_hits_df_p0$family)[i]
  curr_corr_tf_df <- expression_hits_df_p0[expression_hits_df_p0$family == curr_cluster,,drop=FALSE]
  curr_corr_tf_df <- curr_corr_tf_df[order(abs(curr_corr_tf_df$correlation), decreasing = TRUE),,drop=FALSE]
  annot_tfs <- curr_corr_tf_df$gene[1:(ifelse(nrow(curr_corr_tf_df) < top_n, nrow(curr_corr_tf_df), top_n))]
  expression_hits_df_p0$tf_labels[expression_hits_df_p0$family == curr_cluster &
                                 expression_hits_df_p0$gene %in% annot_tfs] <- 
    expression_hits_df_p0$gene[expression_hits_df_p0$family == curr_cluster &
                              expression_hits_df_p0$gene %in% annot_tfs]
}
expression_hits_df_p0$tf_labels[expression_hits_df_p0$pval > 0.05] <- ""
graphics.off()
ggplot(expression_hits_df_p0, aes(x=family, y=correlation, color=expression_2)) +
  geom_jitter(width = 0) + 
  geom_label_repel(label=expression_hits_df_p0$tf_labels, size=4) +
  theme_bw(base_size=14) + 
  theme_Nice() + theme(legend.position = "right") + 
  labs(x="Motif cluster name", 
       y="Spearman R", color="Mean TF expression across pools") 
ggsave(filename=paste0(outdir, "/6c.pdf"), 
       width = 9, height=6, useDingbats=FALSE)