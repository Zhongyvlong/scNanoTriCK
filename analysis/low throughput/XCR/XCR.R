library(dplyr)
library(ggplot2)
library(Cairo)
library(RColorBrewer)
library(GenomicRanges)


library(ArchR)
set.seed(1)
addArchRThreads(threads = 30)
addArchRGenome("mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)





all_cell_DNA_RNA_hap_count_remove_by_quantile <- all_cell_DNA_RNA_hap_count %>%
  dplyr::filter(DNA_sex=="female" & DNA_tagged_count >= 1000 & DNA_x_tagged_count > 0)

p1 <- ggviolin(all_cell_DNA_RNA_hap_count_remove_by_quantile, 
  x = "celltype", y = "scaled_DNA_chrX_mat",
  fill = "celltype", palette = c("#71C9E5","#F4AF5B","#C0986B","#AD768B"),
  add = c("boxplot"), 
  add.params = list(outlier.shape = NA, width=0.08),
  alpha = 1
) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Cell Type",
    y = "scale X mat",
    title = "female cell"
  ) +
  facet_wrap(~histone) + ylim(c(0,0.15))

all_cell_DNA_RNA_hap_count_remove_by_quantile <- all_cell_DNA_RNA_hap_count %>%
  dplyr::filter(DNA_sex=="female" & DNA_tagged_count >= 1000 & DNA_x_tagged_count > 0)

p2 <- ggviolin(all_cell_DNA_RNA_hap_count_remove_by_quantile, 
  x = "celltype", y = "scaled_DNA_chrX_pat",
  fill = "celltype", palette = c("#71C9E5","#F4AF5B","#C0986B","#AD768B"),
  add = c("boxplot"), 
  add.params = list(outlier.shape = NA, width=0.08),
  alpha = 1
) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Cell Type",
    y = "scale X pat",
    title = "female cell"
  ) +
  facet_wrap(~histone) + ylim(c(0,0.15))

p1+p2
