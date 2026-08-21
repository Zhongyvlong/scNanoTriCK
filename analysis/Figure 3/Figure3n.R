library(dplyr)
library(ggplot2)


df_stats <- read.table("Figure3n_data_table.txt", header = T)

width <- 0.5

ggplot() +
  geom_rect(data = df_stats,
            aes(xmin = as.numeric(factor(celltype)) - width/2,
                xmax = as.numeric(factor(celltype)) + width/2,
                ymin = lower, ymax = upper,
                fill = celltype),
            color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(factor(celltype)) - width/2,
                   xend = as.numeric(factor(celltype)) + width/2,
                   y = middle, yend = middle),
               color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(factor(celltype)),
                   xend = as.numeric(factor(celltype)),
                   y = ymin, yend = lower),
               linetype = "dashed", color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(factor(celltype)),
                   xend = as.numeric(factor(celltype)),
                   y = upper, yend = ymax),
               linetype = "dashed", color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(factor(celltype)) - width/4,
                   xend = as.numeric(factor(celltype)) + width/4,
                   y = ymin, yend = ymin),
               color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(factor(celltype)) - width/4,
                   xend = as.numeric(factor(celltype)) + width/4,
                   y = ymax, yend = ymax),
               color = "black") +
  scale_x_continuous(breaks = seq_along(levels(factor(df_long$celltype))),
                     labels = levels(factor(df_long$celltype))) +
  scale_fill_manual(values = c("#71C9E5","#F4AF5B","#C0986B","#AD768B")) +
  scale_color_manual(values = c("#71C9E5","#F4AF5B","#C0986B","#AD768B")) +
  facet_wrap(~allele, ncol = 2) +
  coord_cartesian(ylim = c(0, 5)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "Celltype", y = "Feature counts",
       title = "Distribution of features across celltypes (Xm vs Xp)")

