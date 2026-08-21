library(ggpointdensity)
library(ggplot2)





signal_in_bin_with_different_repetitive_element_density <-
  read.table("Figure2b_ChIP-seq_GM12878_H3K9me3.txt",
              header = F)
ggplot(signal_in_bin_with_different_repetitive_element_density, aes(x = log2(Sum_Signal), y = repeat_percentage)) +
  ggrastr::rasterise(geom_pointdensity(method = 'neighbors', adjust = 0.5), dpi = 300) + # neighbors/kde2d
  scale_color_viridis_c() +
  labs(
    title = "Repeat Percentage vs Log2(Signal)",
    x = "Log2(Signal Value)",
    y = "Repeat Percentage"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )




signal_in_bin_with_different_repetitive_element_density <-
  read.table("Figure2b_scNanoTriCK_GM12878_H3K9me3.txt",
              header = F)
ggplot(signal_in_bin_with_different_repetitive_element_density, aes(x = log2(Sum_Signal), y = repeat_percentage)) +
  ggrastr::rasterise(geom_pointdensity(method = 'neighbors', adjust = 0.5), dpi = 300) + # neighbors/kde2d
  scale_color_viridis_c() +
  labs(
    title = "Repeat Percentage vs Log2(Signal)",
    x = "Log2(Signal Value)",
    y = "Repeat Percentage"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )
















