library(ChIPpeakAnno)
library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(RColorBrewer)



# annotation_results <-
# 	read.table("Figure2d_NGS_GM12878_H3K9me3_peak_annotation_for_repetitive_elements.txt",
# 		header = T)
# annotation_results <-
# 	read.table("Figure2d_scNanoTriCK_GM12878_H3K9me3_peak_annotation_for_repetitive_elements.txt",
# 		header = T)


stats_df <- annotation_results %>%
  count(annotated_element) %>%
  mutate(Percentage = round(n / sum(n) * 100, 1),
         Label = paste0(annotated_element, "\n", n, " (", Percentage, "%)"))

# 按数量排序
stats_df <- stats_df %>% arrange(desc(n))

# 版本1：带标签的饼图
ggplot(stats_df, aes(x = "", y = n, fill = annotated_element)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = Label), 
            position = position_stack(vjust = 0.5),
            size = 4, color = "black", fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Distribution of Peaks in Repetitive Elements",
       subtitle = paste("Total peaks:", sum(stats_df$n)),
       fill = "Element Type") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12))


# 设定顺序
order_levels <- c("LINE", "SINE", "LTR", "Simple_repeat", "Satellite", "OthersRepeat", "unrepeat")

# 把 annotated_element 转成有序因子
stats_df$annotated_element <- factor(stats_df$annotated_element, levels = order_levels)

# 按顺序排序数据框
stats_df <- stats_df[order(stats_df$annotated_element), ]

# 提取数据
counts <- stats_df$n
elements <- stats_df$annotated_element
percent <- round(100 * counts / sum(counts), 2)
labels <- paste0(elements, " (", percent, "%)")

# 画饼图
p <- pie(
  counts,
  labels = labels,
  col = brewer.pal(length(counts), "Set3")
)

