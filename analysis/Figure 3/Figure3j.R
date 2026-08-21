library(ggplot2)
library(dplyr)


df_stats <- read.table("Figure3j-1_data_table.txt", header = T)

width <- 0.5

ggplot() +
  geom_rect(data = df_stats,
            aes(xmin = as.numeric(celltype) - width/2,
                xmax = as.numeric(celltype) + width/2,
                ymin = lower,
                ymax = upper,
                fill = celltype),
            color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype) - width/2,
                   xend = as.numeric(celltype) + width/2,
                   y = middle,
                   yend = middle),
               color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype),
                   xend = as.numeric(celltype),
                   y = ymin,
                   yend = lower),
               linetype = "dashed", color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype),
                   xend = as.numeric(celltype),
                   y = upper,
                   yend = ymax),
               linetype = "dashed", color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype) - width/4,
                   xend = as.numeric(celltype) + width/4,
                   y = ymin,
                   yend = ymin),
               color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype) - width/4,
                   xend = as.numeric(celltype) + width/4,
                   y = ymax,
                   yend = ymax),
               color = "black") +
  scale_fill_manual(values = c("#71C9E5","#F4AF5B","#C0986B","#AD768B")) +
  facet_wrap(~histone) +
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
  coord_cartesian(ylim = c(0, 0.15))



df_stats <- read.table("Figure3j-2_data_table.txt", header = T)

width <- 0.5

p2 <- ggplot() +
  geom_rect(data = df_stats,
            aes(xmin = as.numeric(celltype) - width/2,
                xmax = as.numeric(celltype) + width/2,
                ymin = lower,
                ymax = upper,
                fill = celltype),
            color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype) - width/2,
                   xend = as.numeric(celltype) + width/2,
                   y = middle,
                   yend = middle),
               color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype),
                   xend = as.numeric(celltype),
                   y = ymin,
                   yend = lower),
               linetype = "dashed", color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype),
                   xend = as.numeric(celltype),
                   y = upper,
                   yend = ymax),
               linetype = "dashed", color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype) - width/4,
                   xend = as.numeric(celltype) + width/4,
                   y = ymin,
                   yend = ymin),
               color = "black") +
  geom_segment(data = df_stats,
               aes(x = as.numeric(celltype) - width/4,
                   xend = as.numeric(celltype) + width/4,
                   y = ymax,
                   yend = ymax),
               color = "black") +
  scale_fill_manual(values = c("#71C9E5","#F4AF5B","#C0986B","#AD768B")) +
  facet_wrap(~histone) +
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
  coord_cartesian(ylim = c(0, 0.15))
