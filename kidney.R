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
all_LOO_mat <- as.matrix(read.csv("data/kidney/kidney_influence.csv", row.names = 1))
all_LOO_mat_selection <- all_LOO_mat[rownames(motif_family_annotation[motif_family_annotation$subfamily %in% rownames(subfams_loo_df),]),]
LOO_mat_melted <- melt(all_LOO_mat_selection) %>% magrittr::set_colnames(c("Motif", "Pool", "Weight"))
LOO_mat_melted$Celltype <- curr_colData[LOO_mat_melted$Pool,1]
LOO_mat_melted$`Motif cluster annotation` <- motif_family_annotation[as.character(LOO_mat_melted$Motif),]$subfamily
all_LOO_mat_selection_aggregates <- matrix(nrow=length(unique(LOO_mat_melted$`Motif cluster annotation`)), 
                                           ncol=ncol(all_LOO_mat_selection))
for(i in 1:nrow(all_LOO_mat_selection_aggregates)){
  curr_motif_family <- unique(LOO_mat_melted$`Motif cluster annotation`)[i]
  all_LOO_mat_selection_aggregates[i,] <- colSums(all_LOO_mat_selection[rownames(all_LOO_mat_selection) 
                                                                        %in% rownames(motif_family_annotation[motif_family_annotation$subfamily ==
                                                                                                                curr_motif_family,,drop=FALSE]),])
}
rownames(all_LOO_mat_selection_aggregates) <- unique(LOO_mat_melted$`Motif cluster annotation`)
colnames(all_LOO_mat_selection_aggregates) <- colnames(all_LOO_mat_selection)
all_LOO_mat_selection_aggregates_melt <- melt(all_LOO_mat_selection_aggregates) %>% 
  magrittr::set_colnames(c("Motif cluster annotation", "Pool", "Aggregate of motif weights"))
all_LOO_mat_selection_aggregates_melt$Celltype <- curr_colData[all_LOO_mat_selection_aggregates_melt$Pool,1]
all_LOO_mat_selection_aggregates_melt$Category <- curr_colData[all_LOO_mat_selection_aggregates_melt$Pool,]$cell_type_category
# Order:
all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` <- factor(x=as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`),
                                                                           levels=unique(as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))[order(unique(as.character(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)))])
graphics.off()
# Fig 3a =====
ggplot(all_LOO_mat_selection_aggregates_melt, 
       aes(x=reorder_within(Category,`Aggregate of motif weights`,`Motif cluster annotation`), y=`Aggregate of motif weights`, 
           fill=`Category`)) +
  geom_half_boxplot() + geom_half_violin(side="r") +
  scale_x_reordered() +
  scale_fill_stata() + theme_bw(base_size=14) + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  labs(x = "Category", y="Aggregate of motif influence scores") +
  theme_Nice() +
  facet_wrap(~`Motif cluster annotation`, scales="free", ncol = 4)
ggsave(paste0(outdir, "/3a.pdf"),
       width=16,height=9, useDingbats=FALSE)
dev.off()
