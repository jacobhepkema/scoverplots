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
suppressMessages(library(cowplot))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
options(stringsAsFactors = FALSE)
suppressMessages(source("scover_helper_functions.R"))
# Amount of motifs
d <- 600
outdir <- "output/humankidney/"
if(!dir.exists(outdir)) { dir.create(outdir) }
curr_sce <- readRDS("data/humankidney/kidney_pooled_sce.RDS")
curr_colData <- as.data.frame(colData(curr_sce))
# Load data
category_means_aggregates <- read.csv("data/humankidney/20230123_kidney_subfamilies_cat_loo_df.csv", row.names = 1)
subfam_reproducibility <- read.csv("data/humankidney//20230123_kidney_subfamilies_reproducbility.csv", row.names = 1)
category_means_aggregates <- category_means_aggregates[rownames(subfam_reproducibility)[subfam_reproducibility$Reproducibility >= 0.5],,drop=FALSE]
category_means_aggregates <- category_means_aggregates[rownames(category_means_aggregates) != "N/A",]
subfams_loo_df <- read.csv("data/humankidney//20230123_kidney_subfamilies_loo_df.csv", row.names = 1)
subfams_loo_df <- subfams_loo_df[rownames(subfam_reproducibility)[subfam_reproducibility$Reproducibility >= 0.5],,drop=FALSE]
subfams_loo_df <- subfams_loo_df[rownames(subfams_loo_df) != "N/A",]
motif_family_annotation <- read.csv("data/humankidney//20230123_kidney_motif_family_annotation.csv", row.names=1)
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
pheatmap(to_z(category_means_aggregates), color=heatmap_colors, 
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, 
         filename = paste0(outdir, "2a.pdf"),
         useDingbats=FALSE,
         annotation_row = curr_category_annot)
all_LOO_mat <- as.matrix(read.csv("data/humankidney/20230123_kidney_all_loo_scores_df.csv", row.names = 1))
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
prcomp_mat$`E2F/2` <- as.numeric(all_LOO_mat_selection_aggregates['E2F/2:E2F',])
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

prcomp_mat$`E2F/2` <- as.numeric(all_LOO_mat_selection_aggregates['E2F/2:E2F',])
prcomp_mat$`KLF/SP/2` <- as.numeric(all_LOO_mat_selection_aggregates['KLF/SP/2:C2H2',])
prcomp_mat$`ETS/1` <- as.numeric(all_LOO_mat_selection_aggregates['ETS/1:ETS',])
e2f_plot <- ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=`E2F/2`, shape=Origin)) + 
  geom_point(size=2) + theme_bw(base_size=14) +
  labs(x=paste0("PC1 (", round(variances_expl[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl[2]*100, 2), "%)")) +
  coord_fixed() +
  theme_Nice(angled=FALSE) + theme(legend.position="right")
ets_plot <- ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=`ETS/1`, shape=Origin)) + 
  geom_point(size=2) + theme_bw(base_size=14) +
  labs(x=paste0("PC1 (", round(variances_expl[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl[2]*100, 2), "%)")) +
  coord_fixed() +
  theme_Nice(angled=FALSE) + theme(legend.position="right")
klf_plot <- ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=`KLF/SP/2`, shape=Origin)) + 
  geom_point(size=2) + theme_bw(base_size=14) +
  labs(x=paste0("PC1 (", round(variances_expl[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl[2]*100, 2), "%)")) +
  coord_fixed() +
  theme_Nice(angled=FALSE) + theme(legend.position="right")



all_motifs <- read_meme("data/humankidney/20230122_all_MEME_motifs_kidney.txt")
# NRF1
p1 <- ggplot() + geom_logo(all_motifs$`1_282`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('NRF1') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_NRF1.eps'), width=3.5, height=2)
# ZFX
p2 <- ggplot() + geom_logo(all_motifs$`1_586`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('ZFX') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_motif_ZFX.eps'), width=3.5, height=2)
# YY1
p3 <- ggplot() + geom_logo(all_motifs$`7_575`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('YY1') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_motif_YY1.eps'), width=3.5, height=2)
# E2F/2
p4 <- ggplot() + geom_logo(all_motifs$`8_336`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('E2F/2') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_motif_E2F_2.eps'), width=3.5, height=2)
# KLF/SP/2
p5 <- ggplot() + geom_logo(all_motifs$`3_467`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('KLF/SP/2') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_motif_KLF_SP_2.eps'), width=3.5, height=2)
# Ebox/CACCTG
p6 <- ggplot() + geom_logo(all_motifs$`9_69`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('Ebox/CACCTG') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_motif_Ebox_CACCTG.eps'), width=3.5, height=2)
# CREB/ATF/1
p7 <- ggplot() + geom_logo(all_motifs$`1_340`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('CREB/ATF/1') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_motif_CREB_ATF_1.eps'), width=3.5, height=2)
# GLI
p8 <- ggplot() + geom_logo(all_motifs$`2_143`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('GLI') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_motif_GLI.eps'), width=3.5, height=2)
# ETS/1
p9 <- ggplot() + geom_logo(all_motifs$`6_491`, method='bits') + theme_logo() +
  scale_y_continuous(limits=c(0,2), breaks=c(0,1,2)) + theme_Nice(angled=FALSE) +
  ggtitle('ETS/1') + theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0(outdir, 'kidney_motif_ETS_1.eps'), width=3.5, height=2)


(p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9)
ggsave(filename=paste0(outdir, 'kidney_motifs.eps'), width=7, height=5)
ggsave(filename=paste0(outdir, 'kidney_motifs.pdf'), width=7, height=5)


pseudotime_proliferation_corr_df <- read.csv("data/humankidney/20230124_kidney_pseudotime_corrs.csv", 
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


alphas=.6
yy1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['YY1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['YY1:C2H2',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="YY1:C2H2", 
       y="YY1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.6,1.2)) +
  theme(aspect.ratio=1); plot(yy1_plot)
round(cor(as.numeric(logcounts(curr_sce['YY1',])), 
          as.numeric(all_LOO_mat_selection_aggregates['YY1:C2H2',]), method = "spearman"), 2)
yy2_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['YY2',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['YY1:C2H2',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="YY1:C2H2", 
       y="YY2 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.6,1.2)) +
  theme(aspect.ratio=1); plot(yy2_plot)
round(cor(as.numeric(logcounts(curr_sce['YY2',])), 
          as.numeric(all_LOO_mat_selection_aggregates['YY1:C2H2',]), method = "spearman"), 2)
nrf1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['NRF1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['NRF1:CNC-bZIP',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="NRF1:CNC-bZIP", 
       y="NRF1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.1,0.4)) +
  theme(aspect.ratio=1); plot(nrf1_plot)
round(cor(as.numeric(logcounts(curr_sce['NRF1',])), 
          as.numeric(all_LOO_mat_selection_aggregates['NRF1:CNC-bZIP',]), method = "spearman"), 2)
banp_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['BANP',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['NRF1:CNC-bZIP',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="NRF1:CNC-bZIP", 
       y="BANP expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(0.1,0.4)) +
  theme(aspect.ratio=1); plot(banp_plot)
round(cor(as.numeric(logcounts(curr_sce['BANP',])), 
          as.numeric(all_LOO_mat_selection_aggregates['NRF1:CNC-bZIP',]), method = "spearman"), 2)
e2f6_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['E2F6',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['E2F/2:E2F',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="E2F/2:E2F", 
       y="E2F6 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.05,0.15)) +
  theme(aspect.ratio=1); plot(e2f6_plot)
round(cor(as.numeric(logcounts(curr_sce['E2F6',])), 
          as.numeric(all_LOO_mat_selection_aggregates['E2F/2:E2F',]), method = "spearman"), 2)
e2f4_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['E2F4',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['E2F/2:E2F',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="E2F/2:E2F", 
       y="E2F4 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.05,0.15)) +
  theme(aspect.ratio=1); plot(e2f4_plot)
round(cor(as.numeric(logcounts(curr_sce['E2F4',])), 
          as.numeric(all_LOO_mat_selection_aggregates['E2F/2:E2F',]), method = "spearman"), 2)
graphics.off()

snai2_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['SNAI2',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['SNAI2:C2H2',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="SNAI2:C2H2", 
       y="SNAI2 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.1,-0.025)) +
  theme(aspect.ratio=1); plot(snai2_plot)
round(cor(as.numeric(logcounts(curr_sce['SNAI2',])), 
          as.numeric(all_LOO_mat_selection_aggregates['SNAI2:C2H2',]), method = "spearman"), 2)
zeb1_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['ZEB1',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['SNAI2:C2H2',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="SNAI2:C2H2", 
       y="ZEB1 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.1,-0.025)) +
  theme(aspect.ratio=1); plot(zeb1_plot)
round(cor(as.numeric(logcounts(curr_sce['ZEB1',])), 
          as.numeric(all_LOO_mat_selection_aggregates['SNAI2:C2H2',]), method = "spearman"),
      2)
tcf12_plot <- ggplot(data.frame(
  exp=as.numeric(logcounts(curr_sce['TCF12',])),
  inf=as.numeric(all_LOO_mat_selection_aggregates['Ebox/CAGCTG:bHLH',]),
  cat=curr_colData$Category
), aes(x=inf, y=exp, color=cat)) + geom_point(alpha=alphas) + 
  theme_bw(base_size=14) +
  theme_Nice(angled=FALSE) +
  labs(color="Category", 
       x="Ebox/CAGCTG:bHLH", 
       y="TCF12 expression") +
  scale_color_stata() + scale_x_continuous(breaks=c(-0.25,-0.05)) +
  theme(aspect.ratio=1); plot(tcf12_plot)
round(cor(as.numeric(logcounts(curr_sce['TCF12',])), as.numeric(all_LOO_mat_selection_aggregates['Ebox/CAGCTG:bHLH',]), method = "spearman"),
      2)
(yy1_plot | nrf1_plot | e2f6_plot) / (yy2_plot | banp_plot  | e2f4_plot) / (snai2_plot | zeb1_plot | tcf12_plot)
ggsave(paste0(outdir, "/3c.pdf"), width=8, height=6, useDingbats=FALSE)



aligned_motif_patterns <- read.csv("data/humankidney/20230125_kidney_repr_aligned_motif_patterns_in_genes.csv", row.names = 1,
                                   header = TRUE)
colnames(aligned_motif_patterns) <- str_split_fixed(colnames(aligned_motif_patterns), "X", n=2)[,2]

motif_families <- read.csv("data/humankidney/20230123_kidney_motif_family_annotation.csv", row.names = 1)
motif_families_repro <- read.csv("data/humankidney/20230123_kidney_subfamilies_reproducbility.csv", row.names=1)
prom_mot_fam_scores <- read.csv("data/humankidney/20230125_kidney_prom_mot_fam_scores.csv", row.names=1)
# run umap
custom.config <- umap.defaults
prom_umap <- umap::umap(aligned_motif_patterns, random_state = 42)
prom_umapp <- prom_umap$layout
colnames(prom_umapp) <- c("UMAP1", "UMAP2")
prom_umapp <- cbind(prom_umapp, prom_mot_fam_scores)

ps <- .5
m_klf <- ggplot(prom_umapp, aes(x=UMAP1, y=UMAP2, color=KLF.SP.2.C2H2)) + 
  geom_point(size=ps) + theme_bw(base_size=14) +
  coord_fixed() +
  theme_Nice(angled=FALSE) + theme(legend.position="right")
nrow(prom_umapp[prom_umapp$UMAP2 > 40,])
umap_outlier_filter <- prom_umapp$UMAP2 < 40
almot2 <- aligned_motif_patterns[umap_outlier_filter,]
prom_umap <- umap::umap(almot2, random_state = 1337)
prom_umapp <- prom_umap$layout
colnames(prom_umapp) <- c("UMAP1", "UMAP2")
prom_umapp <- cbind(prom_umapp, prom_mot_fam_scores[umap_outlier_filter,])

umap_melt <- melt(prom_umapp[,c('UMAP1', 'UMAP2', 
                                'YY1.C2H2', 'E2F.2.E2F',
                                'Ebox.CAGCTG.bHLH', 'NRF1.CNC.bZIP',
                                'ZFX.C2H2', 'ETS.1.ETS')],
                  id.vars = c('UMAP1', 'UMAP2'))
umap_melt$variable <- str_replace_all(str_replace_all(str_replace_all(str_replace_all(str_replace_all(str_replace_all(umap_melt$variable, 
                                                                                                                      '\\.C2H2', ':C2H2'),
                                                                                                      '\\.ETS', ':ETS'),
                                                                                      '\\.bHLH', ':bHLH'),
                                                                      '\\.CNC\\.bZIP', ':CNC-bZIP'),
                                                      '\\.E2F', ':E2F'),
                                      '\\.', '/')
umap_melt %>% 
  group_split(variable) %>% 
  map(
    ~ggplot(., aes(UMAP1, UMAP2, color = value)) + 
      geom_point(size=.3) +
      labs(color='Promoter activation') +
      scale_colour_gradient2(
        low = "magenta", 
        mid = "black", 
        high = "yellow", 
        midpoint = 0
      ) +
      theme_Nice(angled = FALSE) +
      theme(legend.position='right') + 
      # coord_equal() +
      facet_grid(~ variable, labeller = function(x) label_value(x, multi_line = FALSE))
  ) %>% 
  plot_grid(plotlist = ., align = 'hv', ncol = 3)
ggsave(paste0(outdir, "3a.pdf"), width = 15, height=7)



# Correlations with expression values
expression_hits_df <- read.csv("data/humankidney/20230123_kidney_expression_hits_df.csv", header=TRUE)
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
graphics.off()
ggplot(expression_hits_df, aes(x=family, y=correlation, color=expression_2)) +
  geom_jitter(width = 0) + 
  geom_label_repel(label=expression_hits_df$tf_labels, size=4) +
  theme_bw(base_size=14) + 
  theme_Nice() + theme(legend.position = "right") + 
  labs(x="Motif cluster name", 
       y="Spearman R", color="Mean TF expression across pools") 
ggsave(filename=paste0(outdir, "/3b.pdf"), 
       width = 12, height=6, useDingbats=FALSE)

cat("Done\n")