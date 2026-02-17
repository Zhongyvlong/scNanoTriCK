library(dplyr)
library(parallel)
library(ArchR)
library(RColorBrewer)
library(cowplot)
library(Cairo)


frag_files <- c(
  list.files(
    path = c("/data/csy/zyl_analysis/demultiplex/re_demultiplex_0303_To_0421_data_250725/merged/GM12878/fragment_dup/",
             "/data/csy/zyl_analysis/demultiplex/20250723/GM12878/fragment/"),
    pattern = ".*_flt_fragment.bed$", full.names = TRUE
  )
)

frag_files <- frag_files[grep("O78|O9[5,6]|O88", frag_files)]
names(frag_files) <- gsub("(.*)_flt_fragment.bed","\\1",basename(frag_files))

frag_files




library(pbmcapply)

frag <- pbmclapply(names(frag_files), function(x) {
  df <- read.table(frag_files[[x]])
  df <- df[,c(1,2,3,5,6)]
  df$cell <- x
  colnames(df) <- c("chrom", "start", "end", "mapq", "strand","cell")
  return(df)
}, ignore.interactive = T, mc.cores = 10) %>%
  do.call(rbind, .)


duprm_bam_stats_DuNACHIP250421 <- read.table("/data/csy/zyl_analysis/demultiplex/re_demultiplex_0303_To_0421_data_250725/merged/duprm_bam.stats", header = T)
duprm_bam_stats_DuNACHIP250421 <- duprm_bam_stats_DuNACHIP250421[,c(1,4)]
duprm_bam_stats_DuNACHIP250421$protein <- "protein"
duprm_bam_stats_DuNACHIP250421[grep("^O88",duprm_bam_stats_DuNACHIP250421$cell),"protein"] <- "H3K9me3"

duprm_bam_stats_DuNACHIP250723 <- read.table("/data/csy/zyl_analysis/demultiplex/20250723/GM12878/mapping/duprm_bam.stats", header = T)
duprm_bam_stats_DuNACHIP250723 <- duprm_bam_stats_DuNACHIP250723[,c(1,4)]
duprm_bam_stats_DuNACHIP250723$protein <- "protein"
duprm_bam_stats_DuNACHIP250723[grep("^O78",duprm_bam_stats_DuNACHIP250723$cell),"protein"] <- "H3K9me3"
duprm_bam_stats_DuNACHIP250723[grep("^O95|^O96",duprm_bam_stats_DuNACHIP250723$cell),"protein"] <- "LaminB1"


duprm_stats_GM12878 <- rbind(duprm_bam_stats_DuNACHIP250421, duprm_bam_stats_DuNACHIP250723)
duprm_stats_GM12878<-duprm_stats_GM12878[duprm_stats_GM12878$protein!="protein",]


cell_used <- duprm_stats_GM12878[duprm_stats_GM12878$read_counts>=1000,"cell"]

library(gUtils)
frag_gr <- dt2gr(frag[frag$cell %in% cell_used, ])

frag_gr


bins <- read.table('/data/ass/genome_file/gencode_ensembl/fasta_hg38/ENSEMBL_2Mb_bin.txt')
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



# frag_cov_bins$cell_hap <- paste(frag_cov_bins$cell,frag_cov_bins$haplotype,sep = "_")


cell_pass <- frag_cov_bins[frag_cov_bins$bin_pct>0.2,]$cell

summary(frag_cov_bins$bin_pct)


ggplot(frag_cov_bins, aes(x = bin_pct)) +
  geom_density() + geom_vline(xintercept = 0.7) +
  theme_bw()


library(tidyverse)

# frag_cov_mtx <- frag_cov %>% 
#   mutate(cell_hap = paste(cell, haplotype, sep = "_")) %>%
#   reshape2::dcast(cell_hap ~ subject.id, value.var = "counts") %>%
#   column_to_rownames("cell_hap")

frag_cov_mtx <- frag_cov %>%
  reshape2::dcast(cell ~ subject.id, value.var = "counts") %>%
  column_to_rownames("cell")

row_ann <- data.frame(
  row.names = duprm_stats_GM12878$cell,
  protein = duprm_stats_GM12878$protein,
  nFrags = log10(duprm_stats_GM12878$read_counts)
)


head(row_ann)


col_ann <- data.frame(
  row.names = 1:nrow(bins),
  chrom = bins$chrom
)

library(pheatmap)

unique(bins$chrom)

protein_cols <- setNames(brewer.pal(12,"Set3")[c(10,5)],c('H3K9me3','LaminB1'))
chr_cols <- setNames(c(rep(c("black","grey"),11),"black"), paste0("chr",c(1:22,"X")))
# hap_cols <- setNames(brewer.pal(9,"Set1")[c(2,1)],c('pat','mat'))


K9_cells <- intersect(rownames(row_ann[row_ann$protein == "H3K9me3",]),cell_pass)
LB1_cells <- intersect(rownames(row_ann[row_ann$protein == "LaminB1",]),cell_pass)

K9_cells <- rownames(slice_max(row_ann[K9_cells,], order_by = nFrags, n=192))
LB1_cells <- rownames(slice_max(row_ann[LB1_cells,], order_by = nFrags, n=192))

# options(repr.plot.width = 10, repr.plot.height = 4)
pheatmap(log2(frag_cov_mtx[K9_cells,]),
         show_rownames = F, show_colnames = F,
         annotation_row = row_ann,
         annotation_col = col_ann, 
         cluster_cols = F,
         annotation_colors = list(chrom = chr_cols, protein = protein_cols),
         clustering_method = "ward.D2",
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)


# options(repr.plot.width = 10, repr.plot.height = 4)
pheatmap(log2(frag_cov_mtx[LB1_cells,]),
         show_rownames = F, show_colnames = F,
         annotation_row = row_ann,
         annotation_col = col_ann, 
         cluster_cols = F,
         annotation_colors = list(chrom = chr_cols, protein = protein_cols),
         clustering_method = "ward.D2",
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)


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


p <- pheatmap(H3K9me3_frag_cov_mtx_lib_norm_bulk_norm[K9_cells,],
              show_rownames = T, show_colnames = F,
              annotation_row = row_ann,
              annotation_col = col_ann,
              cluster_cols = F,
              annotation_colors = list(chrom = chr_cols, protein = protein_cols),
              clustering_method = "ward.D2",
              breaks = seq(0,2,length.out = 100),
              color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)


p <- pheatmap(LaminB1_frag_cov_mtx_lib_norm_bulk_norm[LB1_cells,],
              show_rownames = T, show_colnames = F,
              annotation_row = row_ann,
              annotation_col = col_ann,
              cluster_cols = F,
              annotation_colors = list(chrom = chr_cols, protein = protein_cols),
              clustering_method = "ward.D2",
              breaks = seq(0,2,length.out = 100),
              color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)

frag_cov_mtx_norm <- t(scale(t(log2(frag_cov_mtx))))

#H3K9me3
pheatmap(frag_cov_mtx_norm[K9_cells,],
         show_rownames = F, show_colnames = F,
         annotation_row = row_ann,
         annotation_col = col_ann, 
         cluster_cols = F,
         annotation_colors = list(chrom = chr_cols, protein = protein_cols),
         clustering_method = "ward.D2",
         breaks = seq(-4,4,length.out = 100),
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)

#LaminB1
pheatmap(frag_cov_mtx_norm[LB1_cells,],
         show_rownames = F, show_colnames = F,
         annotation_row = row_ann,
         annotation_col = col_ann, 
         cluster_cols = F,
         annotation_colors = list(chrom = chr_cols, protein = protein_cols),
         clustering_method = "ward.D2",
         breaks = seq(-4,4,length.out = 100),
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)



#在log2标准化和中心化基础上，每列减去每列的fragment数均值，相当于校正细胞间测序深度
frag_cov_mtx_norm2 <- rbind(
  sweep(frag_cov_mtx_norm[K9_cells, ], 2, as.numeric(colMeans(frag_cov_mtx_norm[K9_cells, ], na.rm = T)), "-"),
  sweep(frag_cov_mtx_norm[LB1_cells, ], 2, as.numeric(colMeans(frag_cov_mtx_norm[LB1_cells, ], na.rm = T)), "-")
)

# options(repr.plot.width = 10, repr.plot.height = 4)
p <- pheatmap((frag_cov_mtx_norm2[K9_cells, ]),
              show_rownames = T, show_colnames = F,
              annotation_row = row_ann,
              annotation_col = col_ann,
              cluster_cols = F,
              annotation_colors = list(chrom = chr_cols, protein = protein_cols),
              clustering_method = "ward.D2",
              breaks = seq(-3,3,length.out = 100),
              color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)
ggsave("GM12878_H3K9me3_CNV_250911.pdf",
       plot = p,
       device = "pdf",
       height = 10,
       width = 15)


p <- pheatmap::pheatmap((frag_cov_mtx_norm2[LB1_cells, ]),
                        show_rownames = T, show_colnames = F,
                        annotation_row = row_ann,
                        annotation_col = col_ann,
                        cluster_cols = F,
                        annotation_colors = list(chrom = chr_cols, protein = protein_cols),
                        clustering_method = "ward.D2",
                        breaks = seq(-3,3,length.out = 100),
                        color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
)
ggsave("GM12878_LaminB1_CNV_250911.pdf",
       plot = p,
       device = "pdf",
       height = 10,
       width = 15)





library(circlize)
library(RColorBrewer)



#H3K9me3_frag_cov_mtx_lib_norm_bulk_norm[setdiff(K9_cells,c("O88I59","O88I67","O78I67","O78I65","O88I91","O78I91","O88I78","O78I78","O78I85","O88I41","O78I41","O88I64",
#"O88I82","O78I17","O78I63","O88I63","O88I35","O78I35","O88I68","O78I68","O88I11","O78I11")),]
ComplexHeatmap::Heatmap(
  LaminB1_frag_cov_mtx_lib_norm_bulk_norm[setdiff(LB1_cells,c("O95I17","O95I92","O95I36","O95I44","O95I26","O96I28")),][p[["tree_row"]][["order"]],],
  name = "Expression",  # 图例名称
  show_row_names = F,
  show_column_names = FALSE,
  row_names_side = "right",
  
  # 行注释 - 在这里设置行注释的颜色
  #left_annotation = rowAnnotation(nFrags = runif(164), protein = rep("#80B1D3",times=164)),
  #left_annotation = rowAnnotation(nFrags = runif(164), protein = rep("#80B1D3",times=164)),
  
  # 列注释 - 在这里设置列注释的颜色
  #top_annotation = HeatmapAnnotation(
  #  df = col_ann,
  #  col = list(  # 设置列注释的颜色映射
  #    chrom = chr_cols,
  #    protein = protein_cols
  #  )
  #),
  
  # 聚类方法
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  
  # 热图颜色
  col = circlize::colorRamp2(c(0,0.5,1,1.5,2), c("#053061","#A2CAE1","white","#FBE0D3","#940E26")),
  
  
  # 其他参数
  heatmap_legend_param = list(
    title = "Value",
    title_position = "lefttop-rot"
  ),
  
  use_raster = TRUE, 
  raster_quality = 2
)
ggsave("GM12878_LaminB1_CNV_250911_raster.pdf",
       plot = p,
       device = "pdf",
       height = 10,
       width = 15)



ComplexHeatmap::Heatmap(
  H3K9me3_frag_cov_mtx_lib_norm_bulk_norm[setdiff(K9_cells,c("O88I59","O78I59","O88I67","O78I67","O78I65","O88I91","O78I91","O88I78","O78I78","O78I85","O88I41","O78I41","O88I64",
                                                             "O88I82","O78I17","O78I63","O88I63","O88I35","O78I35","O88I68","O78I68","O88I11","O78I11",
                                                             "O78I34","O88I34","O78I71","O88I71","O78I64","O88I64","O78I24","O88I24","O78I82","O88I82",
                                                             "O78I27","O88I27","O78I25","O88I25","O78I85","O88I85","O78I04","O88I04")),][p[["tree_row"]][["order"]],],
  name = "Expression",  # 图例名称
  show_row_names = F,
  show_column_names = FALSE,
  row_names_side = "right",
  
  # 行注释 - 在这里设置行注释的颜色
  #left_annotation = rowAnnotation(nFrags = runif(164), protein = rep("#80B1D3",times=164)),
  #left_annotation = rowAnnotation(nFrags = runif(164), protein = rep("#80B1D3",times=164)),
  
  # 列注释 - 在这里设置列注释的颜色
  #top_annotation = HeatmapAnnotation(
  #  df = col_ann,
  #  col = list(  # 设置列注释的颜色映射
  #    chrom = chr_cols,
  #    protein = protein_cols
  #  )
  #),
  
  # 聚类方法
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  
  # 热图颜色
  col = circlize::colorRamp2(c(0,0.5,1,1.5,2), c("#053061","#A2CAE1","white","#E5B0A8","#940E26")),
  
  
  # 其他参数
  heatmap_legend_param = list(
    title = "Value",
    title_position = "lefttop-rot"
  ),
  
  use_raster = TRUE, 
  raster_quality = 2
)
ggsave("GM12878_H3K9me3_CNV_250911_raster.pdf",
       plot = p,
       device = "pdf",
       height = 10,
       width = 15)


colorRampPalette(c("#053061","#A2CAE1","white","#FAD9C9","#940E26"))(21)




"#053061" "#244E7A" "#436D94" "#638CAD" "#82ABC7" "#A2CAE1" "#B4D4E7" "#C7DFED" "#D9E9F3" "#ECF4F9" "#FFFFFF" "#FEF7F4" "#FDEFE9" "#FBE8DE" "#FBE0D3" "#FAD9C9" "#E5B0A8" "#D18787"
"#BC5F67" "#A83646" "#940E26"












