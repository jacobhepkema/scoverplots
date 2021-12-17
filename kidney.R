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
all_LOO_mat <- as.matrix(read.csv("data/kidney/kidney_influence.csv.gz", row.names = 1))
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
prcomp_mat <- prcomp(t(all_LOO_mat_selection))
eigs <- prcomp_mat$sdev^2
variances_expl <- eigs/sum(sum(eigs))
prcomp_mat <- prcomp_mat$x
prcomp_mat <- as.data.frame(prcomp_mat)
curr_colData$Origin <- sapply(curr_colData$most_abundant_dataset, FUN=function(x){ 
  if(x == "1") { return("Mature") } else { return("Fetal") }  
  })
curr_colData$Category <- curr_colData$cell_type_category
curr_colData$Celltype <- curr_colData$most_abundant_celltype
prcomp_mat$Celltype <- curr_colData$Celltype
prcomp_mat$Category <- curr_colData$Category
prcomp_mat$Origin <- as.character(curr_colData[rownames(prcomp_mat),]$Origin)
graphics.off()
# Fig 2c =====
ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=Category, shape=Origin)) + 
  geom_point(size=2) + theme_bw(base_size=14) +
  labs(x=paste0("PC1 (", round(variances_expl[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl[2]*100, 2), "%)")) +
  coord_fixed() +
  theme_Nice(angled=FALSE) + theme(legend.position="right") +
  scale_color_stata()
ggsave(paste0(outdir, "2c.pdf"), width=7, height=7,
       useDingbats=FALSE)
all_motifs <- read_meme("data/kidney/kidney_motifs.meme")
# CREB
ggplot() + geom_logo(all_motifs$`4_285`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_CREB.eps"), width=3.5,height=2)
# ETS/1
ggplot() + geom_logo(all_motifs$`4_394`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_ETS1.eps"), width=3.5,height=2)
# ETS/2
ggplot() + geom_logo(all_motifs$`4_475`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_ETS2.eps"), width=3.5,height=2)
# KLF/SP/2
ggplot() + geom_logo(all_motifs$`2_297`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_KLF.eps"), width=3.5,height=2)
# E2F
ggplot() + geom_logo(all_motifs$`0_55`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_E2F.eps"), width=3.5,height=2)
# YY1
ggplot() + geom_logo(all_motifs$`2_426`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_YY1.eps"), width=3.5,height=2)
# NRF1
ggplot() + geom_logo(all_motifs$`4_89`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_NRF1.eps"), width=3.5,height=2)
# ZNF143
ggplot() + geom_logo(all_motifs$`2_474`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_ZNF143.eps"), width=3.5,height=2)
# KAISO
ggplot() + geom_logo(all_motifs$`3_195`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_KAISO.eps"), width=3.5,height=2)
# GC-tract
ggplot() + geom_logo(all_motifs$`2_190`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_GC-tract.eps"), width=3.5,height=2)
# ZFX
ggplot() + geom_logo(all_motifs$`4_189`, method = "bits") + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks = c(0,1,2)) + theme_Nice(angled = FALSE)
ggsave(filename=paste0(outdir, "kidney_ZFX.eps"), width=3.5,height=2)
pseudotime_proliferation_corr_df <- read.csv("data/kidney/kidney_correlations.csv", 
                                             row.names=1)
pseudotime_proliferation_corr_df <- t(pseudotime_proliferation_corr_df)
rownames(pseudotime_proliferation_corr_df) <- c("Correlation with proximal tubule pseudotime", 
                                                "Correlation with proliferation",
                                                "Correlation with GO:proximal tubule development")
graphics.off()
# Fig 2e =====
pheatmap::pheatmap(pseudotime_proliferation_corr_df, cluster_rows = FALSE, 
                   color=heatmap_colors, border_color = NA, cellwidth = 12, cellheight = 12,
                   angle_col=45, filename = paste0(outdir, "/2e.pdf"), width=7, height=3)
cat("Done\n")
CREM_df <- read.csv("data/kidney/kidney_CREB_CREM.csv", row.names=1)
ELF1_df <- read.csv("data/kidney/kidney_ETS1_ELF1.csv", row.names=1)
ETV3_df <- read.csv("data/kidney/kidney_ETS1_ETV3.csv", row.names=1)
FLI1_df <- read.csv("data/kidney/kidney_ETS1_FLI1.csv", row.names=1)
YY1_df <- read.csv("data/kidney/kidney_YY1_YY1.csv", row.names=1)
YY2_df <- read.csv("data/kidney/kidney_YY1_YY2.csv", row.names=1)
NRF1_df <- read.csv("data/kidney/kidney_NRF1_NRF1.csv", row.names=1)
BANP_df <- read.csv("data/kidney/kidney_NRF1_BANP.csv", row.names=1)
alphas <- 0.6
CREB_plot <- ggplot(CREM_df, aes(x=agg_scores, y=exp_scores, color=cell_type_category)) +
  geom_point(alpha=alphas) + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="CREB/ATF/1:bZIP", 
       y="CREM expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0, 0.02, 0.04)) + theme(aspect.ratio=1)
YY1_plot <- ggplot(YY1_df, aes(x=agg_scores, y=exp_scores, color=cell_type_category)) +
  geom_point(alpha=alphas) + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE)+
  labs(color="Category", 
       x="YY1:C2H2", 
       y="YY1 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0.25,0.45)) + theme(aspect.ratio=1)
YY2_plot <- ggplot(YY2_df, aes(x=agg_scores, y=exp_scores, color=cell_type_category)) +
  geom_point(alpha=alphas) + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE)+
  labs(color="Category", 
       x="YY1:C2H2", 
       y="YY2 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0.25,0.45)) + theme(aspect.ratio=1)
NRF1_plot <- ggplot(NRF1_df, aes(x=agg_scores, y=exp_scores, color=cell_type_category)) +
  geom_point(alpha=alphas) + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE)+
  labs(color="Category", 
       x="NRF1:CNC-bZIP", 
       y="NRF1 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0, 0.25,0.35,0.45)) + theme(aspect.ratio=1)
BANP_plot <- ggplot(BANP_df, aes(x=agg_scores, y=exp_scores, color=cell_type_category)) +
  geom_point(alpha=alphas) + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE)+
  labs(color="Category", 
       x="NRF1:CNC-bZIP", 
       y="BANP expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0, 0.25,0.35,0.45)) + theme(aspect.ratio=1)
ELF1_plot <- ggplot(ELF1_df, aes(x=agg_scores, y=exp_scores, color=cell_type_category)) +
  geom_point(alpha=alphas) + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE)+
  labs(color="Category", 
       x="ETS/1:ETS", 
       y="ELF1 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0.25,0.45)) + theme(aspect.ratio=1)
FLI1_plot <- ggplot(FLI1_df, aes(x=agg_scores, y=exp_scores, color=cell_type_category)) +
  geom_point(alpha=alphas) + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE)+
  labs(color="Category", 
       x="ETS/1:ETS", 
       y="FLI1 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0.25,0.45)) + theme(aspect.ratio=1)
ETV3_plot <- ggplot(ETV3_df, aes(x=agg_scores, y=exp_scores, color=cell_type_category)) +
  geom_point(alpha=alphas) + theme_bw(base_size=14) +
  theme_Nice(angled=FALSE)+
  labs(color="Category", 
       x="ETS/1:ETS", 
       y="ETV3 expression") + 
  scale_color_stata() + scale_x_continuous(breaks=c(0.25,0.45)) + theme(aspect.ratio=1) +
  theme(legend.position="right")
graphics.off()
# Fig 3b =====
(YY1_plot | NRF1_plot | ELF1_plot) / (YY2_plot | BANP_plot  | FLI1_plot)
ggsave(paste0(outdir, "/3b.pdf"), width=8, height=4, useDingbats=FALSE)
# Correlations were obtained using e.g. cor(ETV5_corr_df$agg_LOO, ETV5_corr_df$exp, method="spearman")