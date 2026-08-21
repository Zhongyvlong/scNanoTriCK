library(ggpointdensity)
library(ggplot2)


signal_in_bin_with_different_repetitive_element_density <-
  read.table("Ex4b_ChIP-seqs_GM12878_LaminB1.txt",
              header = F)
ggplot(signal_in_bin_with_different_repetitive_element_density, aes(x = log2(Sum_Signal), y = repeat_percentage)) +
  ggrastr::rasterise(geom_pointdensity(method = 'neighbors', adjust = 0.5), dpi = 300) +
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
  read.table("Ex4b_scNanoTriCK_GM12878_LaminB1.txt",
              header = F)
ggplot(signal_in_bin_with_different_repetitive_element_density, aes(x = log2(Sum_Signal), y = repeat_percentage)) +
  ggrastr::rasterise(geom_pointdensity(method = 'neighbors', adjust = 0.5), dpi = 300) +
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
