library(dplyr)
library(parallel)
library(ArchR)
library(RColorBrewer)
library(cowplot)
library(Cairo)


frag_files <- c(
  list.files(
    path = c("fragment/"),
    pattern = ".*_flt_fragment.bed$", full.names = TRUE
  )
)

names(frag_files) <- gsub("(.*)_flt_fragment.bed","\\1",basename(frag_files))

library(pbmcapply)

frag <- pbmclapply(names(frag_files), function(x) {
  df <- read.table(frag_files[[x]])
  df <- df[,c(1,2,3,5,6)]
  df$cell <- x
  colnames(df) <- c("chrom", "start", "end", "mapq", "strand","cell")
  return(df)
}, ignore.interactive = T, mc.cores = 10) %>%
  do.call(rbind, .)


library(gUtils)
frag_gr <- dt2gr(frag)

bins <- read.table('ENSEMBL_2Mb_bin.txt')
colnames(bins) <- c("chrom","start","end")
bins$start<-bins$start+1

bins <- bins[bins$chrom %in% paste0("chr",c(1:22,"X")),]
bins$chrom <- factor(bins$chrom, levels = paste0("chr",c(1:22,"X")))
bins <- bins %>% arrange(chrom, start)
bins_gr <- dt2gr(bins)



frag_gr_ov <- frag_gr %*% bins_gr


frag_cov <- frag_gr_ov %>% as.data.frame() %>% 
  filter(width > 1) %>% 
  dplyr::group_by(cell, subject.id) %>%
  dplyr::summarise(counts = n_distinct(query.id))


frag_cov_bins <- frag_cov %>% 
  dplyr::group_by(cell) %>%
  dplyr::summarise(cov_bins = n_distinct(subject.id))


frag_cov_bins$bin_pct <- frag_cov_bins$cov_bins/nrow(bins)


dim(frag_cov_bins)
nrow(frag_cov_bins[frag_cov_bins$bin_pct>0.7,])
nrow(frag_cov_bins[frag_cov_bins$bin_pct>0.7,])/nrow(frag_cov_bins)

cell_pass <- frag_cov_bins[frag_cov_bins$bin_pct>0.2,]$cell

ggplot(frag_cov_bins, aes(x = bin_pct)) +
  geom_density() + geom_vline(xintercept = 0.7) +
  theme_bw()


library(tidyverse)

frag_cov_mtx <- frag_cov %>%
  reshape2::dcast(cell ~ subject.id, value.var = "counts") %>%
  column_to_rownames("cell")

row_ann <- data.frame(
  row.names = duprm_stats_GM12878$cell,
  protein = duprm_stats_GM12878$protein,
  nFrags = log10(duprm_stats_GM12878$read_counts)
)

col_ann <- data.frame(
  row.names = 1:nrow(bins),
  chrom = bins$chrom
)

library(pheatmap)

unique(bins$chrom)

protein_cols <- setNames(brewer.pal(12,"Set3")[c(10,5)],c('H3K9me3','LaminB1'))
chr_cols <- setNames(c(rep(c("black","grey"),11),"black"), paste0("chr",c(1:22,"X")))

K9_cells <- intersect(rownames(row_ann[row_ann$protein == "H3K9me3",]),cell_pass)
LB1_cells <- intersect(rownames(row_ann[row_ann$protein == "LaminB1",]),cell_pass)

K9_cells <- rownames(slice_max(row_ann[K9_cells,], order_by = nFrags, n=192))
LB1_cells <- rownames(slice_max(row_ann[LB1_cells,], order_by = nFrags, n=192))

row_sums <- rowSums(frag_cov_mtx, na.rm = TRUE)
frag_cov_mtx_lib_norm <- sweep(frag_cov_mtx, 1, row_sums, "/")

H3K9me3_norm_vec <- colSums(frag_cov_mtx[K9_cells,], na.rm = T)/sum(frag_cov_mtx[K9_cells,], na.rm = T)
LaminB1_norm_vec <- colSums(frag_cov_mtx[LB1_cells,], na.rm = T)/sum(frag_cov_mtx[LB1_cells,], na.rm = T)

H3K9me3_frag_cov_mtx_lib_norm_bulk_norm <- t(t(frag_cov_mtx_lib_norm[K9_cells, ])/H3K9me3_norm_vec)
H3K9me3_frag_cov_mtx_lib_norm_bulk_norm[H3K9me3_frag_cov_mtx_lib_norm_bulk_norm >= 0.75 & H3K9me3_frag_cov_mtx_lib_norm_bulk_norm <= 1.25] <- 1
H3K9me3_frag_cov_mtx_lib_norm_bulk_norm[H3K9me3_frag_cov_mtx_lib_norm_bulk_norm > 1.5] <- 2

LaminB1_frag_cov_mtx_lib_norm_bulk_norm <- t(t(frag_cov_mtx_lib_norm[LB1_cells, ])/LaminB1_norm_vec)
LaminB1_frag_cov_mtx_lib_norm_bulk_norm[LaminB1_frag_cov_mtx_lib_norm_bulk_norm >= 0.75 & LaminB1_frag_cov_mtx_lib_norm_bulk_norm <= 1.25] <- 1
LaminB1_frag_cov_mtx_lib_norm_bulk_norm[LaminB1_frag_cov_mtx_lib_norm_bulk_norm > 1.5] <- 2

pheatmap(H3K9me3_frag_cov_mtx_lib_norm_bulk_norm[K9_cells,],
  show_rownames = T, show_colnames = F,
  annotation_row = row_ann,
  annotation_col = col_ann,
  cluster_cols = F,
  annotation_colors = list(chrom = chr_cols, protein = protein_cols),
  clustering_method = "ward.D2",
  breaks = seq(0,2,length.out = 100),
  color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)

pheatmap(LaminB1_frag_cov_mtx_lib_norm_bulk_norm[LB1_cells,],
  show_rownames = T, show_colnames = F,
  annotation_row = row_ann,
  annotation_col = col_ann,
  cluster_cols = F,
  annotation_colors = list(chrom = chr_cols, protein = protein_cols),
  clustering_method = "ward.D2",
  breaks = seq(0,2,length.out = 100),
  color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)
