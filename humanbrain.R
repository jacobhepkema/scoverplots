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
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
options(stringsAsFactors = FALSE)
suppressMessages(source("scover_helper_functions.R"))
# Amount of motifs
d <- 600
outdir <- "output/humanbrain/"
if(!dir.exists(outdir)) { dir.create(outdir) }

curr_sce_trevino <- readRDS("data/humanbrain/humanbrain_pooled_sce.RDS")
curr_colData_trevino <- as.data.frame(colData(curr_sce_trevino))
category_means_aggregates_z_trevino <- read.csv("data/humanbrain/20230110_trevino_100neighbours_subfamilies_cat_loo_df.csv", row.names = 1)
subfam_reproducibility_trevino <- read.csv("data/humanbrain/20230110_trevino_100neighbours_subfamilies_reproducbility.csv", row.names = 1)
category_means_aggregates_z_trevino <- category_means_aggregates_z_trevino[subfam_reproducibility_trevino$Reproducibility >= 0.5,]
category_means_aggregates_z_trevino <- category_means_aggregates_z_trevino[rownames(category_means_aggregates_z_trevino) != "N/A",]
subfams_loo_df_trevino <- read.csv("data/humanbrain/20230110_trevino_100neighbours_subfamilies_loo_df.csv", row.names = 1)
subfams_loo_df_trevino <- subfams_loo_df_trevino[subfam_reproducibility_trevino$Reproducibility >= 0.5,]
subfams_loo_df_trevino <- subfams_loo_df_trevino[rownames(subfams_loo_df_trevino) != "N/A",]
motif_family_annotation_trevino <- read.csv("data/humanbrain/20230110_trevino_100neighbours_motif_family_annotation.csv", row.names=1)
curr_category_annot_trevino <- data.frame(row.names=rownames(category_means_aggregates_z_trevino), 
                                          amount_motifs = rep("", nrow(category_means_aggregates_z_trevino)))
curr_category_annot_trevino$amount_motifs <- sapply(rownames(category_means_aggregates_z_trevino), 
                                                    FUN=function(x){
                                                      sum(motif_family_annotation_trevino$subfamily == x)
                                                    })
colnames(curr_category_annot_trevino) <- "Motifs in cluster"
curr_category_annot_trevino$`Summed mean influence` <- rowMeans(subfams_loo_df_trevino)[rownames(curr_category_annot_trevino)]

category_means_aggregates_z_trevino <- category_means_aggregates_z_trevino[,order(colnames(category_means_aggregates_z_trevino))]
colnames(category_means_aggregates_z_trevino) <- sapply(colnames(category_means_aggregates_z_trevino), FUN=function(x){str_replace_all(x, "\\.", "/")}) 
category_means_aggregates_z_trevino <- category_means_aggregates_z_trevino[order(sapply(rownames(category_means_aggregates_z_trevino), 
                                                                                        FUN=function(x){
                                                                                          str_split(x, ":")[[1]][2]
                                                                                        })),]
heatmap_colors <- colorRampPalette(c("magenta", "black", "yellow"))(100)
colnames(category_means_aggregates_z_trevino) <- c(
  'Cycling progenitor','Inhibitory','L1/Progenitor','L2','L3','L4','L5','Microglia','Progenitor','Radial glia','Subplate'
)
pheatmap(to_z(category_means_aggregates_z_trevino), color=heatmap_colors, 
         angle_col = 45, cellwidth = 14, cellheight = 14,
         border_color = NA, cluster_cols = FALSE,
         cluster_rows = TRUE,
         filename = paste0(outdir, "heatmap_humanbrain.pdf"),
         useDingbats=FALSE,
         annotation_row = curr_category_annot_trevino)


celltype_to_category = c('IN1'= 'Inhibitory',
                         'IN2' = 'Inhibitory',
                         'GluN5' = 'L5',
                         'GluN3' = 'L3',
                         'IN3' = 'Inhibitory',
                         'nIPC/GluN1' = 'L1 / Progenitor',
                         'RG' = 'Radial glia',
                         'Cyc. Prog.' = 'Cycling progenitor',
                         'GluN2' = 'L2',
                         'SP' = 'Subplate',
                         'MG' = 'Microglia',
                         'oIPC/OPC' = 'Progenitor',
                         'GluN4' = 'L4'
)
curr_colData_trevino$cell_type_category <- celltype_to_category[curr_colData_trevino$most_abundant_celltype]


all_LOO_mat_trevino <- as.matrix(read.csv("data/humanbrain/20230110_trevino_100neighbours_all_loo_scores_df.csv", row.names = 1))
all_LOO_mat_selection_trevino <- all_LOO_mat_trevino[rownames(motif_family_annotation_trevino[motif_family_annotation_trevino$subfamily %in% rownames(subfams_loo_df_trevino),]),]
LOO_mat_melted_trevino <- melt(all_LOO_mat_selection_trevino) %>% magrittr::set_colnames(c("Motif", "Pool", "Weight"))
LOO_mat_melted_trevino$Celltype <- curr_colData_trevino[LOO_mat_melted_trevino$Pool,1]
LOO_mat_melted_trevino$`Motif cluster annotation` <- motif_family_annotation_trevino[as.character(LOO_mat_melted_trevino$Motif),]$subfamily
LOO_mat_melted_trevino$Category <- curr_colData_trevino[LOO_mat_melted_trevino$Pool,]$cell_type_category
all_LOO_mat_selection_aggregates_trevino <- matrix(nrow=length(unique(LOO_mat_melted_trevino$`Motif cluster annotation`)), 
                                                   ncol=ncol(all_LOO_mat_selection_trevino))
for(i in 1:nrow(all_LOO_mat_selection_aggregates_trevino)){
  curr_motif_family <- unique(LOO_mat_melted_trevino$`Motif cluster annotation`)[i]
  all_LOO_mat_selection_aggregates_trevino[i,] <- colSums(all_LOO_mat_selection_trevino[rownames(all_LOO_mat_selection_trevino) 
                                                                                        %in% rownames(motif_family_annotation_trevino[motif_family_annotation_trevino$subfamily ==
                                                                                                                                        curr_motif_family,,drop=FALSE]),])
}
rownames(all_LOO_mat_selection_aggregates_trevino) <- unique(LOO_mat_melted_trevino$`Motif cluster annotation`)
colnames(all_LOO_mat_selection_aggregates_trevino) <- colnames(all_LOO_mat_selection_trevino)


all_motifs <- read_meme("data/humanbrain/20230110_all_MEME_motifs_trevino_100neighbours.txt")
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
nfi_plot <- plot_motif(all_motifs$`2_226`) + ggtitle('NFI/3'); print(nfi_plot)
mef_plot <- plot_motif(all_motifs$`3_597`) + ggtitle('MEF2'); print(mef_plot)
hd_plot <- plot_motif(all_motifs$`7_567`) + ggtitle('HD/2'); print(hd_plot)
tbx_plot <- plot_motif(all_motifs$`5_291`) + ggtitle('TBX/4'); print(tbx_plot)
ebox_cagctg_plot <- plot_motif(all_motifs$`0_537`) + ggtitle('Ebox/CAGCTG'); print(ebox_cagctg_plot)
ebox_cacctg_plot <- plot_motif(all_motifs$`1_465`) + ggtitle('Ebox/CACCTG'); print(ebox_cacctg_plot)
rfx_plot <- plot_motif(all_motifs$`4_329`) + ggtitle('RFX/3'); print(rfx_plot)
creb_plot <- plot_motif(all_motifs$`8_4`) + ggtitle('CREB/ATF/1'); print(creb_plot)
pou_plot <- plot_motif(all_motifs$`7_255`) + ggtitle('POU/2'); print(pou_plot)
ctcf_plot <- plot_motif(all_motifs$`1_248`) + ggtitle('CTCF'); print(ctcf_plot)
sox_plot <- plot_motif(all_motifs$`3_144`) + ggtitle('SOX/1'); print(sox_plot)
mies_plot <- plot_motif(all_motifs$`0_148`) + ggtitle('MIES'); print(mies_plot)
(nfi_plot | hd_plot | ebox_cagctg_plot) / (mef_plot | tbx_plot | ebox_cacctg_plot) /
  (rfx_plot | creb_plot | pou_plot) / (ctcf_plot | sox_plot | mies_plot)
ggsave(filename=paste0(outdir, '/6b.pdf'), width=7, height=6, useDingbats=FALSE)



prcomp_mat <- prcomp(t(all_LOO_mat_selection_trevino))
eigs2 <- prcomp_mat$sdev^2
variances_expl2 <- eigs2/sum(sum(eigs2))
prcomp_mat <- prcomp_mat$x
prcomp_mat <- as.data.frame(prcomp_mat)
# Add labels randomly:
prcomp_mat$Category <- curr_colData_trevino$cell_type_category
prcomp_mat$Celltype <- curr_colData_trevino$most_abundant_celltype
prcomp_mat <- cbind(prcomp_mat, t(all_LOO_mat_selection_aggregates_trevino))
prcomp_mat <- prcomp_mat[sample(1:nrow(prcomp_mat), replace = FALSE),]
prcomp_mat$ggrepel_labels <- prcomp_mat$Celltype
prcomp_mat$ggrepel_labels[duplicated(prcomp_mat$ggrepel_labels)] <- ""
ggplot(prcomp_mat, aes(x=PC1, y=PC2, color=Celltype)) + geom_point() +
  scale_color_stata() + theme_Nice(angled = FALSE) +
  labs(x=paste0("PC1 (", round(variances_expl2[1]*100, 2), "%)"),
       y=paste0("PC2 (", round(variances_expl2[2]*100, 2), "%)")) +
  geom_label_repel(label=prcomp_mat$ggrepel_labels, show.legend = FALSE, max.overlaps = 100) + 
  coord_fixed()
ggsave(paste0(outdir, "6c.pdf"), width=7, height=7,
       useDingbats=FALSE)

prc2 <- melt(prcomp_mat[,c('PC1', 'PC2', 'NFI/3:NFI', 'MEF2:MADS',
                           'HD/2:homeodomain', 'TBX/4:TBX',
                           'Ebox/CAGCTG:bHLH', 'Ebox/CACCTG:bHLH')], id.vars = c("PC1", "PC2"), variable.name = "Motif family")
prc2 %>% 
  group_split(`Motif family`)  %>%
  map(
    ~ggplot(., aes(PC1, PC2, color = value)) + 
      geom_point(size = 1) +
      scale_colour_gradient2(
        low = "magenta", 
        mid = "black", 
        high = "yellow", 
        midpoint = 0
      ) +
      facet_wrap(~ `Motif family`, labeller = function(x) label_value(x, multi_line = FALSE)) +
      theme_Nice() + theme(legend.position = "right") +
      coord_equal()
  ) %>% 
  plot_grid(plotlist = ., align = 'hv', ncol = 2) 
# (p1 | p2) / (p3 | p4) / (p5 | p6)
ggsave(filename = paste0(outdir, "/6d.pdf"), width = 7,height = 7)



n_bins=14
p1 <- ggplot(data.frame(loo=as.numeric(all_LOO_mat_selection_aggregates_trevino['MEF2:MADS',]),
                        exp=as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'MEF2C',]))),
             aes(x=exp, y=loo)) + theme_Nice() + geom_bin_2d(bins=n_bins) +
  scale_fill_viridis_c(option = "A") + theme(legend.position = "right", 
                                             aspect.ratio=1) +
  labs(x='MEF2C', y='MEF2') +
  ggtitle(paste0("R = ", round(cor(as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'MEF2C',])),
                                   as.numeric(all_LOO_mat_selection_aggregates_trevino['MEF2:MADS',]),
                                   method='spearman'), digits=2)))
p2 <- ggplot(data.frame(loo=as.numeric(all_LOO_mat_selection_aggregates_trevino['Ebox/CAGCTG:bHLH',]),
                        exp=as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'NEUROD2',]))),
             aes(x=exp, y=loo)) + theme_Nice() + geom_bin_2d(bins=n_bins) +
  scale_fill_viridis_c(option = "A") + theme(legend.position = "right", 
                                             aspect.ratio=1) +
  labs(x='NEUROD2', y='Ebox/CAGCTG') +
  ggtitle(paste0("R = ", round(cor(as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'NEUROD2',])),
                                   as.numeric(all_LOO_mat_selection_aggregates_trevino['Ebox/CAGCTG:bHLH',]),
                                   method='spearman'), digits=2)))
p3 <- ggplot(data.frame(loo=as.numeric(all_LOO_mat_selection_aggregates_trevino['HD/2:homeodomain',]),
                        exp=as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'EMX2',]))),
             aes(x=exp, y=loo)) + theme_Nice() + geom_bin_2d(bins=n_bins) +
  scale_fill_viridis_c(option = "A") + theme(legend.position = "right", 
                                             aspect.ratio=1) +
  labs(x='EMX2', y='HD/2') +
  ggtitle(paste0("R = ", round(cor(as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'EMX2',])),
                                   as.numeric(all_LOO_mat_selection_aggregates_trevino['HD/2:homeodomain',]),
                                   method='spearman'), digits=2)))
p4 <- ggplot(data.frame(loo=as.numeric(all_LOO_mat_selection_aggregates_trevino['NFI/3:NFI',]),
                        exp=as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'NFIC',]))),
             aes(x=exp, y=loo)) + theme_Nice() + geom_bin_2d(bins=n_bins) +
  scale_fill_viridis_c(option = "A") + theme(legend.position = "right", 
                                             aspect.ratio=1) +
  labs(x='NFIC', y='NFI/3') +
  ggtitle(paste0("R = ", round(cor(as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'NFIC',])),
                                   as.numeric(all_LOO_mat_selection_aggregates_trevino['NFI/3:NFI',]),
                                   method='spearman'), digits=2)))
p5 <- ggplot(data.frame(loo=as.numeric(all_LOO_mat_selection_aggregates_trevino['TBX/4:TBX',]),
                        exp=as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'TBR1',]))),
             aes(x=exp, y=loo)) + theme_Nice() + geom_bin_2d(bins=n_bins) +
  scale_fill_viridis_c(option = "A") + theme(legend.position = "right", 
                                             aspect.ratio=1) +
  labs(x='TBR1', y='TBX/4') +
  ggtitle(paste0("R = ", round(cor(as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'TBR1',])),
                                   as.numeric(all_LOO_mat_selection_aggregates_trevino['TBX/4:TBX',]),
                                   method='spearman'), digits=2)))
p6 <- ggplot(data.frame(loo=as.numeric(all_LOO_mat_selection_aggregates_trevino['Ebox/CACCTG:bHLH',]),
                        exp=as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'ID4',]))),
             aes(x=exp, y=loo)) + theme_Nice() + geom_bin_2d(bins=n_bins) +
  scale_fill_viridis_c(option = "A") + theme(legend.position = "right", 
                                             aspect.ratio=1) +
  labs(x='ID4', y='Ebox/CACCTG') +
  ggtitle(paste0("R = ", round(cor(as.numeric(logcounts(curr_sce_trevino[rowData(curr_sce_trevino)$Symbol == 'ID4',])),
                                   as.numeric(all_LOO_mat_selection_aggregates_trevino['Ebox/CACCTG:bHLH',]),
                                   method='spearman'), digits=2)))
(p4 | p1) / (p3 | p5) / (p2 | p6)
ggsave(filename=paste0(outdir, '/6e.pdf'), width=7, height=7)


