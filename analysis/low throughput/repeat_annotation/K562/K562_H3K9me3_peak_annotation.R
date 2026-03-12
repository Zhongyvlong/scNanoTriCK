library(ChIPpeakAnno)
library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(RColorBrewer)

priority_order <- c(LINE = 1, LTR = 2, SINE = 3, Simple_repeat = 4, Satellite = 5, OthersRepeat = 6)

bed_files <- c(
  LINE = "LINE.bed",
  SINE = "SINE.bed", 
  LTR = "LTR.bed",
  Satellite = "Satellite.bed",
  Simple_repeat = "Simple_repeat.bed",
  OthersRepeat = "OthersRepeat.bed"
)

cat("Reading peak file...\n")
# peaks <- read.table("ENCFF963GZJ_K562_H3K9me3_hg38.bed")
# peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns = F,
#                                   start.field = "V2",seqnames.field = "V1",end.field = "V3")
# peaks <- read.table("GSM7964462_K562_H3K9me3.bed")
# peaks <- makeGRangesFromDataFrame(peaks, keep.extra.columns = F,
#                                   start.field = "V2",seqnames.field = "V1",end.field = "V3")
# peaks <- import("scNanoTriCK_K562_H3K9me3.broadPeak")

cat("Total peaks:", length(peaks), "\n")

cat("Reading repetitive element files...\n")
all_repeats <- GRanges()
for (name in names(bed_files)) {
  cat("  Reading", name, "...\n")
  repeat_gr <- import(bed_files[[name]], format = "BED")
  repeat_gr$type <- name
  all_repeats <- c(all_repeats, repeat_gr)
}
cat("Total repetitive elements:", length(all_repeats), "\n")

annotation_results <- data.frame(
  peak_id = 1:length(peaks),
  chromosome = seqnames(peaks),
  start = start(peaks),
  end = end(peaks),
  peak_length = width(peaks),
  annotated_element = "unrepeat",
  annotation_method = "none",
  all_overlaps = NA_character_,
  overlap_condition = NA_character_,
  stringsAsFactors = FALSE
)

cat("Finding overlaps...\n")
overlaps <- findOverlaps(peaks, all_repeats)
cat("Initial overlaps found:", length(overlaps), "\n")

if (length(overlaps) > 0) {
  cat("Calculating overlap ratios...\n")
  overlap_ranges <- pintersect(peaks[queryHits(overlaps)], all_repeats[subjectHits(overlaps)])
  
  overlap_lengths <- width(overlap_ranges)
  repeat_lengths <- width(all_repeats[subjectHits(overlaps)])
  peak_lengths <- width(peaks[queryHits(overlaps)])
  
  overlap_ratios_repeat <- overlap_lengths / repeat_lengths
  overlap_ratios_peak <- overlap_lengths / peak_lengths
  
  valid_indices <- which(overlap_ratios_repeat >= 2/3 | overlap_ratios_peak >= 2/3)
  cat("Valid overlaps (ratio >= 2/3 for repeat OR peak):", length(valid_indices), "\n")
  
  repeat_only <- sum(overlap_ratios_repeat >= 2/3 & overlap_ratios_peak < 2/3)
  peak_only <- sum(overlap_ratios_repeat < 2/3 & overlap_ratios_peak >= 2/3)
  both_conditions <- sum(overlap_ratios_repeat >= 2/3 & overlap_ratios_peak >= 2/3)
  
  cat("Overlaps satisfying repeat ratio >= 2/3 only:", repeat_only, "\n")
  cat("Overlaps satisfying peak ratio >= 2/3 only:", peak_only, "\n")
  cat("Overlaps satisfying both conditions:", both_conditions, "\n")
  cat("Total valid overlaps:", repeat_only + peak_only + both_conditions, "\n")
  
  if (length(valid_indices) > 0) {
    valid_overlaps_df <- data.frame(
      peak_id = queryHits(overlaps)[valid_indices],
      repeat_type = all_repeats$type[subjectHits(overlaps)[valid_indices]],
      overlap_ratio_repeat = overlap_ratios_repeat[valid_indices],
      overlap_ratio_peak = overlap_ratios_peak[valid_indices],
      condition = ifelse(overlap_ratios_repeat[valid_indices] >= 2/3 & 
                           overlap_ratios_peak[valid_indices] >= 2/3, "both",
                         ifelse(overlap_ratios_repeat[valid_indices] >= 2/3, "repeat_only", "peak_only"))
    )
    
    cat("Grouping overlaps by peak...\n")
    peak_valid_overlaps <- split(valid_overlaps_df, valid_overlaps_df$peak_id)
    
    cat("Annotating peaks by priority...\n")
    
    for (peak_id in names(peak_valid_overlaps)) {
      peak_idx <- as.numeric(peak_id)
      peak_data <- peak_valid_overlaps[[peak_id]]
      
      overlapping_types <- unique(peak_data$repeat_type)
      annotation_results$all_overlaps[peak_idx] <- paste(overlapping_types, collapse = ";")
      
      condition_counts <- table(peak_data$condition)
      annotation_results$overlap_condition[peak_idx] <- paste(
        names(condition_counts), condition_counts, sep = ":", collapse = ";"
      )
      
      valid_types <- overlapping_types[overlapping_types %in% names(priority_order)]
      if (length(valid_types) > 0) {
        priorities <- priority_order[valid_types]
        annotation_results$annotated_element[peak_idx] <- names(which.min(priorities))
        annotation_results$annotation_method[peak_idx] <- "overlap"
      }
    }
  }
}

cat("\n=== Final Overlap Annotation Statistics ===\n")
final_stats <- table(annotation_results$annotated_element)
print(final_stats)

cat("\nAnnotation methods:\n")
print(table(annotation_results$annotation_method))

if (sum(!is.na(annotation_results$overlap_condition)) > 0) {
  cat("\nOverlap condition distribution:\n")
  all_conditions <- unlist(strsplit(annotation_results$overlap_condition[!is.na(annotation_results$overlap_condition)], ";"))
  condition_table <- table(all_conditions)
  print(condition_table)
}

stats_df <- annotation_results %>%
  count(annotated_element) %>%
  mutate(Percentage = round(n / sum(n) * 100, 1),
         Label = paste0(annotated_element, "\n", n, " (", Percentage, "%)"))

stats_df <- stats_df %>% arrange(desc(n))

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


order_levels <- c("LINE", "SINE", "LTR", "Simple_repeat", "Satellite", "OthersRepeat", "unrepeat")

stats_df$annotated_element <- factor(stats_df$annotated_element, levels = order_levels)

stats_df <- stats_df[order(stats_df$annotated_element), ]

counts <- stats_df$n
elements <- stats_df$annotated_element
percent <- round(100 * counts / sum(counts), 2)
labels <- paste0(elements, " (", percent, "%)")

pie(
  counts,
  labels = labels,
  col = brewer.pal(length(counts), "Set3")
)
