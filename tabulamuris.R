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
category_means_aggregates <- read.csv("data/tm/20221111_tm_subfamilies_cat_loo_df.csv", row.names = 1)
subfam_reproducibility <- read.csv("data/tm/20221111_tm_subfamilies_reproducbility.csv", row.names = 1)
category_means_aggregates <- category_means_aggregates[rownames(subfam_reproducibility)[subfam_reproducibility$Reproducibility >= 0.5],,drop=FALSE]
category_means_aggregates <- category_means_aggregates[rownames(category_means_aggregates) != "N/A",]
subfams_loo_df <- read.csv("data/tm/20221111_tm_subfamilies_loo_df.csv", row.names = 1)
subfams_loo_df <- subfams_loo_df[rownames(subfam_reproducibility)[subfam_reproducibility$Reproducibility >= 0.5],,drop=FALSE]
subfams_loo_df <- subfams_loo_df[rownames(subfams_loo_df) != "N/A",]
motif_family_annotation <- read.csv("data/tm/20221111_tm_motif_families.csv", row.names=1)
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
pheatmap(to_z(category_means_aggregates), color=heatmap_colors, 
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, 
         filename = paste0(outdir, "4a.pdf"),
         useDingbats=FALSE,
         annotation_row = curr_category_annot)


# See object `curr_category_annot` for numbers
all_LOO_mat <- as.matrix(read.csv("data/tm/20221111_tm_all_loo_scores_df.csv", row.names = 1))
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
all_LOO_mat_selection_aggregates_melt <- read.csv("data/tm/20221111_tm_all_loo_scores_df_aggregate_melt.csv")
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

all_motifs <- read_meme("data/tm/20221111_tm_all_MEME_motifs.txt")
plot_motif <- function(mat){
  return(
    ggplot() + geom_logo(mat, method='bits') + theme_logo() +
      scale_y_continuous(limits=c(0,2), breaks=c(0,1,2),
                         expand = c(0, 0)) + theme_Nice(angled=FALSE) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line.x=element_blank(),
            plot.title=element_text(hjust=0.5))
  )
}
klf_plot <- plot_motif(all_motifs$`0_437`) + ggtitle('KLF/SP/2'); print(klf_plot)
ets_plot <- plot_motif(all_motifs$`9_467`) + ggtitle('ETS/1'); print(ets_plot)
yy1_plot <- plot_motif(all_motifs$`7_402`) + ggtitle('YY1'); print(yy1_plot)
zfx_plot <- plot_motif(all_motifs$`2_202`) + ggtitle('ZFX'); print(zfx_plot)
ebox_caccgt_plot <- plot_motif(all_motifs$`6_21`) + ggtitle('Ebox/CACCGT'); print(ebox_caccgt_plot)
gc_tract_plot <- plot_motif(all_motifs$`8_317`) + ggtitle('GC-tract'); print(gc_tract_plot)
creb_plot <- plot_motif(all_motifs$`3_478`) + ggtitle('CREB/ATF/1'); print(creb_plot)
ebox_cacgtg_1_plot <- plot_motif(all_motifs$`8_322`) + ggtitle('Ebox/CACGTG/1'); print(ebox_cacgtg_1_plot)
e2f_plot <- plot_motif(all_motifs$`0_103`) + ggtitle('E2F/2'); print(e2f_plot)
nr_plot <- plot_motif(all_motifs$`1_353`) + ggtitle('NR/1'); print(nr_plot)

(klf_plot | ets_plot | yy1_plot) / (zfx_plot | ebox_caccgt_plot | ebox_cacgtg_1_plot) / (e2f_plot | creb_plot | nr_plot)
ggsave(filename=paste0(outdir, '/4b.pdf'), width=7, height=4, useDingbats=FALSE)



prcomp_mat <- prcomp(t(all_LOO_mat_selection))
eigs <- prcomp_mat$sdev^2
variances_expl <- eigs/sum(sum(eigs))
prcomp_mat <- prcomp_mat$x
prcomp_mat <- as.data.frame(prcomp_mat)
# Add labels randomly:
prcomp_mat$Category <- curr_colData$Category
prcomp_mat$Celltype <- curr_colData$Celltype
prcomp_mat <- prcomp_mat[sample(1:nrow(prcomp_mat), replace = FALSE),]
prcomp_mat$ggrepel_labels <- prcomp_mat$Category
prcomp_mat$ggrepel_labels[duplicated(prcomp_mat$ggrepel_labels)] <- ""
# Fig 4c =====
ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=Category)) + geom_point() +
  scale_color_stata() + theme_Nice(angled = FALSE) +
  labs(x=paste0("PC1 (", round(variances_expl[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl[2]*100, 2), "%)")) +
  geom_label_repel(label=prcomp_mat$ggrepel_labels, size=4, show.legend = FALSE) + 
  coord_fixed()
ggsave(paste0(outdir, "4c.pdf"), width=7, height=7,
       useDingbats=FALSE)
plot(variances_expl[1:10]*100, main='Variance explained for PCA on influence scores',
     xlab='PC', ylab='Variance explained (%)', ylim=c(0,70))
ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=Category)) + geom_point() +
  scale_color_stata() + theme_Nice(angled = FALSE) + theme(legend.position = "right") + 
  labs(x=paste0("PC1 (", round(variances_expl[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl[2]*100, 2), "%)")) +
  geom_label_repel(label=prcomp_mat$ggrepel_labels, size=4, show.legend = FALSE) + 
  coord_fixed()
ggsave(paste0(outdir, "4c_with_legend.pdf"), width=7, height=7,
       useDingbats=FALSE)

# umaps =====
aligned_motif_patterns <- read.csv("data/tm/20230126_tm_repr_aligned_motif_patterns_in_genes.csv", row.names = 1,
                                   header = TRUE)
colnames(aligned_motif_patterns) <- str_split_fixed(colnames(aligned_motif_patterns), "X", n=2)[,2]
motif_families <- read.csv("data/tm/20221111_tm_motif_families.csv", row.names = 1)
motif_families_repro <- read.csv("data/tm/20221111_tm_subfamilies_reproducbility.csv", row.names=1)
prom_mot_fam_scores <- read.csv("data/tm/20230126_tm_prom_mot_fam_scores.csv", row.names=1)
# run umap
custom.config <- umap.defaults
prom_umap <- umap::umap(aligned_motif_patterns, random_state = 42)
prom_umapp <- prom_umap$layout
colnames(prom_umapp) <- c("UMAP1", "UMAP2")
prom_umapp <- cbind(prom_umapp, prom_mot_fam_scores)

ps <- .5
umap_melt <- melt(prom_umapp[,c('UMAP1', 'UMAP2', 
                                'YY1.C2H2', 'E2F.2.E2F',
                                'KLF.SP.2.C2H2', 'Ebox.CACGTG.2.bHLH',
                                'ZFX.C2H2', 'ETS.1.ETS',
                                'NR.1.nuclearreceptor',
                                'CREB.ATF.1.bZIP')],
                  id.vars = c('UMAP1', 'UMAP2'))
umap_melt %>% 
  group_split(variable) %>% 
  map(
    ~ggplot(., aes(UMAP1, UMAP2, color = value)) + 
      geom_point(size=.3) +
      scale_colour_gradient2(
        low = "magenta", 
        mid = "black", 
        high = "yellow", 
        midpoint = 0
      ) +
      theme_Nice() +
      theme(legend.position='right') + 
      coord_equal() +
      facet_grid(~ variable, labeller = function(x) label_value(x, multi_line = FALSE))
  ) %>% 
  plot_grid(plotlist = ., align = 'hv', ncol = 2)
ggsave(paste0(outdir, "5a.pdf"), width = 6, height=10, useDingbats=FALSE)


# # Data for the GO term plot was calculated as follows: 
# GO_terms <- readLines("data/tm/go_scfind.tsv")
# names(GO_terms) <- sapply(GO_terms, FUN=function(x){return(str_split(x, "\t")[[1]][1])})
# GO_terms <- sapply(GO_terms, FUN=function(x){
#   str_split(str_split(x, "\t")[[1]][2], ",")[[1]]
# })
# curr_GO_lengths <- c()
# pb <- txtProgressBar(0, length(GO_terms), style=3)
# for(i in 1:length(GO_terms)){
#   curr_GO_name <- names(GO_terms)[i]
#   curr_GO_genes <- GO_terms[[i]]
#   curr_GO_genes <- unique(curr_GO_genes[curr_GO_genes %in% rownames(curr_sce)])
#   curr_GO_length <- length(curr_GO_genes)
#   curr_GO_lengths <- c(curr_GO_lengths, curr_GO_length)
#   setTxtProgressBar(pb, i)
# }
# # Select GO terms with at least 1 gene and fewer than 50 genes in set (that are found
# # in the current experiment)
# GO_selection <- GO_terms[which(curr_GO_lengths < 50 & curr_GO_lengths > 0)]
# curr_GO_lengths_selection <- curr_GO_lengths[curr_GO_lengths < 50 & curr_GO_lengths > 0]
# GO_selection_corrs <- list()
# GO_selection_pvals <- list()
# pb <- txtProgressBar(0, length(GO_selection), style=3)
# for(j in 1:length(GO_selection)){ # Correlate expression of GO term genes to motif influence scores. This can take a while
#   curr_GO_name <- names(GO_selection)[j]
#   curr_GO_genes <- GO_selection[[j]]
#   curr_GO_genes <- unique(curr_GO_genes[curr_GO_genes %in% rownames(curr_sce)])
#   curr_GO_expression <- colMeans(logcounts(curr_sce[curr_GO_genes,]))
#   curr_GO_corrs <- c()
#   curr_GO_pvals <- c()
#   for(i in 1:length(unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))){
#     curr_cluster_annot <- as.character(unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`))[i]
#     curr_melt_aggregates <- all_LOO_mat_selection_aggregates_melt[all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation` == curr_cluster_annot,]
#     curr_GO_corrs <- c(curr_GO_corrs, cor(curr_melt_aggregates$`Aggregate of motif influence scores`, curr_GO_expression,
#                                           method="spearman"))
#     curr_GO_pvals <- c(curr_GO_pvals, p.adjust(cor.test(curr_melt_aggregates$`Aggregate of motif influence scores`,
#                                                         curr_GO_expression,
#                                                         method="spearman",
#                                                         exact = FALSE)$p.value, method="fdr"))
#   }
#   names(curr_GO_corrs) <- unique(all_LOO_mat_selection_aggregates_melt$`Motif cluster annotation`)
#   names(curr_GO_pvals) <- names(curr_GO_corrs)
#   GO_selection_corrs[[j]] <- curr_GO_corrs
#   GO_selection_pvals[[j]] <- curr_GO_pvals
#   setTxtProgressBar(pb, j)
# }
# names(GO_selection_corrs) <- names(GO_selection)
# names(GO_selection_pvals) <- names(GO_selection)
# GO_melt <- melt(GO_selection_corrs)
# GO_melt$motif_family <- names(GO_selection_corrs[[1]])
# GO_melt_p <- melt(GO_selection_pvals)
# GO_melt$p_corr <- GO_melt_p$value
# colnames(GO_melt) <- c("Spearman R", "GO term", "Motif cluster annotation", "Corrected p-value")
# GO_melt_cast <- dcast(GO_melt, formula=`GO term`~`Motif cluster annotation`, value.var="Spearman R")
# rownames(GO_melt_cast) <- GO_melt_cast[,1]
# GO_melt_cast <- GO_melt_cast[,-c(1)]
# GO_melt_cast_pval <- dcast(GO_melt, formula=`GO term`~`Motif cluster annotation`, value.var="Corrected p-value")
# GO_melt_cast_pval <- as.data.frame(GO_melt_cast_pval)
# rownames(GO_melt_cast_pval) <- GO_melt_cast_pval[,1]
# GO_melt_cast_pval <- GO_melt_cast_pval[,-c(1)]
# bottom_top_quantiles <- quantile(as.numeric(unlist(GO_melt_cast)), c(0.01, 0.99))
# has_no_significant <- apply(GO_melt_cast, 1, FUN=function(x){
#   return((sum(x < bottom_top_quantiles[1]) +
#             sum(x > bottom_top_quantiles[2]))
#          == 0)
# })
# sum(!has_no_significant)
# # Cluster annotation
# k <- 3
# GO_melt_cast_no_sig <- GO_melt_cast[!has_no_significant,]
# clust_annot <- cutree(hclust(dist(GO_melt_cast_no_sig)), k=k)
# GO_row_annot <- data.frame(row.names=names(clust_annot),
#                            "GO term cluster"=as.character(clust_annot))
# colnames(GO_row_annot) <- "GO term cluster"
# GO_row_annot_col <- list(
#   "GO term cluster" = stata_pal()(k)
# )
# names(GO_row_annot_col$`GO term cluster`) <- as.character(1:k)
# saveRDS(GO_row_annot, "data/tm/tm_go_term_annotation.RDS")
# saveRDS(GO_row_annot_col, "data/tm/tm_go_term_colours.RDS")
# saveRDS(GO_melt_cast[!has_no_significant,], "data/tm/tm_go_terms.RDS")

GO_melt_cast_no_sig <- readRDS("data/tm/tm_go_terms.RDS")
GO_row_annot_col <- readRDS("data/tm/tm_go_term_colours.RDS")
GO_row_annot <- readRDS("data/tm/tm_go_term_annotation.RDS")
pheatmap::pheatmap(GO_melt_cast_no_sig,
                   show_rownames = FALSE, border_color = NA,
                   color=heatmap_colors, cellwidth = 12,
                   angle_col = 45, cellheight = .14,
                   annotation_row = GO_row_annot,
                   annotation_colors = GO_row_annot_col,
                   filename = paste0(outdir, "/5b.pdf"),
                   useDingbats=FALSE)


# fig 5c ========

expression_hits_df <- read.csv("data/tm/20221111_tm_expression_corrs.csv", header=TRUE)
expression_hits_df$expression_2 <- expression_hits_df$mean_expression
expression_hits_df$expression_2[expression_hits_df$pval > 0.05] <- NA
expression_hits_df$tf_labels <- NA
top_n <- 3
for(i in 1:length(unique(expression_hits_df$family))){
  curr_cluster <- unique(expression_hits_df$family)[i]
  curr_corr_tf_df <- expression_hits_df[expression_hits_df$family == curr_cluster,,drop=FALSE]
  curr_corr_tf_df <- curr_corr_tf_df[order(abs(curr_corr_tf_df$correlation), decreasing = TRUE),,drop=FALSE]
  annot_tfs <- curr_corr_tf_df$gene[1:(ifelse(nrow(curr_corr_tf_df) < top_n, nrow(curr_corr_tf_df), top_n))]
  expression_hits_df$tf_labels[expression_hits_df$family == curr_cluster &
                                 expression_hits_df$gene %in% annot_tfs] <- 
    expression_hits_df$gene[expression_hits_df$family == curr_cluster &
                              expression_hits_df$gene %in% annot_tfs]
}
expression_hits_df$tf_labels[expression_hits_df$pval > 0.05] <- ""
# Fig 3c =====
ggplot(expression_hits_df, aes(x=family, y=correlation, color=expression_2)) +
  geom_jitter(width = 0) + 
  geom_label_repel(label=expression_hits_df$tf_labels, size=4) +
  theme_bw(base_size=14) + 
  theme_Nice() + theme(legend.position = "right") + 
  labs(x="Motif cluster name", 
       y="Spearman R", color="Mean TF expression across pools") 
ggsave(filename=paste0(outdir, "/5c.pdf"), 
       width = 12, height=6, useDingbats=FALSE)



# fig 5d ========
alphas=.5
nr1h2_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Nr1h2',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['NR/1:nuclearreceptor',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="NR/1:nuclearreceptor", 
       y="Nr1h2 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.02,0.07)) +
  theme(aspect.ratio=1); plot(nr1h2_plot)
round(cor(as.numeric(logcounts(curr_sce['Nr1h2',])), 
          as.numeric(all_LOO_mat_selection_aggregates['NR/1:nuclearreceptor',]), method = "spearman"), 2)
nr2f6_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Nr2f6',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['NR/1:nuclearreceptor',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="NR/1:nuclearreceptor", 
       y="Nr2f6 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.02,0.07)) +
  theme(aspect.ratio=1); plot(nr2f6_plot)
round(cor(as.numeric(logcounts(curr_sce['Nr2f6',])), 
          as.numeric(all_LOO_mat_selection_aggregates['NR/1:nuclearreceptor',]), method = "spearman"), 2)

atf1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Atf1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['CREB/ATF/1:bZIP',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="CREB/ATF/1:bZIP", 
       y="Atf1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.04,0.16)) +
  theme(aspect.ratio=1); plot(atf1_plot)
round(cor(as.numeric(logcounts(curr_sce['Atf1',])), 
          as.numeric(all_LOO_mat_selection_aggregates['CREB/ATF/1:bZIP',]), method = "spearman"), 2)
creb1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Creb1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['CREB/ATF/1:bZIP',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="CREB/ATF/1:bZIP", 
       y="Creb1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.04,0.16)) +
  theme(aspect.ratio=1); plot(creb1_plot)

clock_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Clock',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['Ebox/CACGTG/2:bHLH',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="Ebox/CACGTG/2:bHLH", 
       y="Clock expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.05,0.05)) +
  theme(aspect.ratio=1); plot(clock_plot)
round(cor(as.numeric(logcounts(curr_sce['Clock',])), 
          as.numeric(all_LOO_mat_selection_aggregates['Ebox/CACGTG/2:bHLH',]), method = "spearman"), 2)
hes1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Hes1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['Ebox/CACGTG/2:bHLH',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="Ebox/CACGTG/2:bHLH", 
       y="Hes1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.05,0.05)) +
  theme(aspect.ratio=1); plot(hes1_plot)
round(cor(as.numeric(logcounts(curr_sce['Hes1',])), 
          as.numeric(all_LOO_mat_selection_aggregates['Ebox/CACGTG/2:bHLH',]), method = "spearman"), 2)
elf1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Elf1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['ETS/1:ETS',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="ETS/1:ETS", 
       y="Elf1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(1,3)) +
  theme(aspect.ratio=1); plot(elf1_plot)
round(cor(as.numeric(logcounts(curr_sce['Elf1',])), 
          as.numeric(all_LOO_mat_selection_aggregates['ETS/1:ETS',]), method = "spearman"), 2)
elf4_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Elf4',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['ETS/1:ETS',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="ETS/1:ETS", 
       y="Elf4 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(1,3)) +
  theme(aspect.ratio=1); plot(elf4_plot)
round(cor(as.numeric(logcounts(curr_sce['Elf4',])), 
          as.numeric(all_LOO_mat_selection_aggregates['ETS/1:ETS',]), method = "spearman"), 2)
zfx_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Zfx',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['ZFX:C2H2',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="ZFX:C2H2", 
       y="Zfx expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.075,0.150)) +
  theme(aspect.ratio=1); plot(zfx_plot)
round(cor(as.numeric(logcounts(curr_sce['Zfx',])), 
          as.numeric(all_LOO_mat_selection_aggregates['ZFX:C2H2',]), method = "spearman"), 2)
sp1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Sp1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['KLF/SP/2:C2H2',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="KLF/SP/2:C2H2", 
       y="Sp1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.8,1.6)) +
  theme(aspect.ratio=1); plot(sp1_plot)
round(cor(as.numeric(logcounts(curr_sce['Sp1',])), 
          as.numeric(all_LOO_mat_selection_aggregates['KLF/SP/2:C2H2',]), method = "spearman"), 2)
vezf1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Vezf1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['KLF/SP/2:C2H2',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="KLF/SP/2:C2H2", 
       y="Vezf1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.8,1.6)) +
  theme(aspect.ratio=1); plot(vezf1_plot)
round(cor(as.numeric(logcounts(curr_sce['Vezf1',])), 
          as.numeric(all_LOO_mat_selection_aggregates['KLF/SP/2:C2H2',]), method = "spearman"), 2)
tcf4_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Tcf4',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['Ebox/CACCTG:bHLH',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="Ebox/CACCTG:bHLH", 
       y="Tcf4 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.7,-0.2)) +
  theme(aspect.ratio=1); plot(tcf4_plot)
round(cor(as.numeric(logcounts(curr_sce['Tcf4',])), 
          as.numeric(all_LOO_mat_selection_aggregates['Ebox/CACCTG:bHLH',]), method = "spearman"), 2)
snai2_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['Snai2',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['Ebox/CACCTG:bHLH',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="Ebox/CACCTG:bHLH", 
       y="Snai2 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.7,-0.2)) +
  theme(aspect.ratio=1); plot(snai2_plot)
(elf1_plot | elf4_plot | clock_plot) / (sp1_plot | vezf1_plot | nr1h2_plot) / (zfx_plot | atf1_plot | tcf4_plot)
ggsave(paste0(outdir, "/5d.pdf"), width=6.5, height=6.5, useDingbats=FALSE)

cat("Done\n")