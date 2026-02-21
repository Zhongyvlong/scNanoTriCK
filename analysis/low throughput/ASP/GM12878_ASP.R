library(dplyr)
library(parallel)
library(ArchR)
library(RColorBrewer)
library(cowplot)

addArchRGenome("hg38")
set.seed(1)

addArchRThreads(threads = 10)

#haplotyped fragment file directory
frag_files <- c(
  list.files(
    path = c("/data/csy/zyl_analysis/demultiplex/re_demultiplex_0303_To_0421_data_250725/merged/GM12878/phasing_bam/",
             "/data/csy/zyl_analysis/demultiplex/20250723/GM12878/phasing_bam/"),
    pattern = "_duprm_phased_HP.*.bed.gz$", full.names = TRUE
  )
)

frag_files <- frag_files[grep("O8[3,4,9]", frag_files, invert = T)]
names(frag_files) <- gsub("(.*)_duprm_phased_(.*).bed.gz","\\1_\\2",basename(frag_files))

library(pbmcapply)

frag <- pbmclapply(names(frag_files), function(x) {
  df <- read.table(frag_files[[x]])
  colnames(df) <- c("chrom", "start", "end", "mapq", "strand")
  df$file <- x
  df$cell <- gsub("_.*", "", x)
  df$haplotype <- gsub(".*_", "", x)
  return(df)
}, ignore.interactive = T, mc.cores = 10) %>%
  do.call(rbind, .)

cell_used <- unique(stringr::str_split(names(frag_files), pattern = "_", simplify = T)[,1])
frag_flt <- frag[frag$cell %in% cell_used, ]

#peak file directory
peak_files <- list.files(path = "/data/csy/zyl_analysis/demultiplex/20250723/GM12878/callpeak", pattern = ".broadPeak", full.names = T) %>%
  grep(pattern = "CTCF", invert = T, value = T) %>% grep(pattern = "GM12878", value = T)
peak_files <- peak_files[c(2,4)]
names(peak_files) <- c("GM12878_H3K9me3","GM12878_LaminB1")

peak_set <- lapply(peak_files, function(x){
  df <- read.table(x)
  colnames(df) <- c("chrom","start","end","name","fc")
  df <- df[,1:5]
  return(df)
})

#convert peak region dataframe to grange object
library(gUtils)
frag_gr <- dt2gr(frag_flt)

library(BSgenome.Hsapiens.UCSC.hg38)
reset_seqlen <- function(gr, ref = BSgenome.Hsapiens.UCSC.hg38){
  gr@seqinfo@seqlengths <- ref@seqinfo[gr@seqinfo@seqnames]@seqlengths
  return(gr)
}

#set chromosome length
frag_gr <- frag_gr %>% reset_seqlen()

frag_gr_grl <- lapply(names(peak_set), function(x){
  gr <- frag_gr[frag_gr$cell %in% cell_used, ]
})
names(frag_gr_grl) <- names(peak_set)

#extract chromatin modification signal site
frag_flank <- flank(frag_gr, width = 1)

frag_flank_grl <- lapply(names(peak_set), function(x){
  gr <- frag_flank[frag_flank$cell %in% cell_used, ]
})
names(frag_flank_grl) <- names(peak_set)

#extend chromatin modification signal site +/- 75 bp
frag_slop_grl <- lapply(names(peak_set), function(x){
  resize(frag_gr[frag_gr$cell %in% cell_used, ], fix = 'center',width = 150)
}) 
names(frag_slop_grl) <- names(peak_set)

peak_set_gr <- lapply(peak_set, function(df){
  dt2gr(df) %>% reset_seqlen()
})

#extract signal
ov <- lapply(names(frag_flank_grl), function(x){
  frag_flank_grl[[x]] %*% peak_set_gr[[x]]
})
names(ov) <- names(frag_flank_grl)

#create signal dataframe
asp_counts <- lapply(ov, function(gr){
  gr %>%
    as.data.frame() %>%
    dplyr::group_by(subject.id, haplotype) %>%
    dplyr::summarise(
      ncell = n_distinct(cell),
      nfrags = length(cell)
    )
})

asp_counts_df <- lapply(asp_counts, function(df) {
  full_join(df[df$haplotype == "HP1", c(1, 3, 4)],
            df[df$haplotype == "HP2", c(1, 3, 4)],
            suffix = c("_pat", "_mat"), by = "subject.id"
  ) %>% mutate_all(~ coalesce(., 0)) %>%
    mutate(total_frags = nfrags_pat + nfrags_mat) %>%
    mutate(pat_ratio = nfrags_pat / (nfrags_pat + nfrags_mat))
})

#call allel specific peak
asp_counts_df <- lapply(asp_counts_df, function(df) {
  df$p_value <- apply(df, 1, function(x) {
    pat <- x[3]
    mat <- x[5]
    if (!is.na(mat) & !is.na(pat)) {
      res <- binom.test(c(mat, pat), p = 0.5)
      return(res$p.value)
    } else {
      return(NA)
    }
  })
  df$fdr <- p.adjust(df$p_value, method = "BH")
  
  df$significance <- "N.S."
  try(df[!is.na(df$fdr) & df$fdr < 0.05, ]$significance <- "*", silent = T)
  try(df[!is.na(df$fdr) & df$fdr < 0.01, ]$significance <- "**", silent = T)
  try(df[!is.na(df$fdr) & df$fdr < 0.001, ]$significance <- "***", silent = T)
  df$significance <- factor(df$significance, levels = c("N.S.","*","**","***"))
  
  df$type <- "none"
  df[!is.na(df$pat_ratio) & df$pat_ratio > 0.5, ]$type <- "pat"
  df[!is.na(df$pat_ratio) & df$pat_ratio < 0.5, ]$type <- "mat"
  df[!is.na(df$pat_ratio) & df$significance == "N.S.",]$type <- "none"
  return(df)
})

peak_set <- lapply(peak_set, function(df){
  df$idx <- 1:nrow(df)
  return(df)
})

asp_counts_df <- lapply(names(asp_counts_df), function(x){
  left_join(asp_counts_df[[x]],peak_set[[x]][,c("chrom","start","end","idx","fc")], by = c("subject.id" = "idx"))
})
names(asp_counts_df) <- names(peak_set_gr)

ggplot(asp_counts_df$GM12878_LaminB1, aes(x = log10(total_frags), y = pat_ratio)) +
  geom_point(aes(color = fdr)) +
  scale_color_viridis_c(direction = -1) +
  theme_bw()

ggplot(asp_counts_df$GM12878_H3K9me3, aes(x = log10(total_frags), y = pat_ratio)) +
  geom_point(aes(color = fdr)) +
  scale_color_viridis_c(direction = -1) +
  theme_bw()


asp_counts_df$GM12878_LaminB1$type <- factor(asp_counts_df$GM12878_LaminB1$type, levels = c("none","mat","pat"))
asp_counts_df$GM12878_H3K9me3$type <- factor(asp_counts_df$GM12878_H3K9me3$type, levels = c("none","mat","pat"))


library(ggrastr)
p <- ggplot(asp_counts_df$GM12878_LaminB1, aes(x = log10(total_frags), y = pat_ratio)) +
  geom_point_rast(aes(color = type)) +
  scale_color_manual(values = c("darkgrey",brewer.pal(9,"Set1")[c(1,2)])) +
  theme_bw() + labs(title = "GM12878 LaminB1") +
  xlab("Number of fragments") +
  ylab("Paternal/(Paternal + Maternal)") +
  theme(panel.grid = element_blank())
p



p <- ggplot(asp_counts_df$GM12878_H3K9me3, aes(x = log10(total_frags), y = pat_ratio)) +
  geom_point_rast(aes(color = type)) +
  scale_color_manual(values = c("darkgrey",brewer.pal(9,"Set1")[c(1,2)])) +
  theme_bw() + labs(title = "GM12878 H3K9me3") +
  xlab("Number of fragments") +
  ylab("Paternal/(Paternal + Maternal)") +
  theme(panel.grid = element_blank())
p






asp_df <- lapply(asp_counts_df, function(df){
  df[!is.na(df$fdr) & df$fdr < 0.05,] %>% droplevels()
})
names(peak_set)

summary(asp_counts_df$GM12878_LaminB1[asp_counts_df$GM12878_LaminB1$significance != "N.S.",]$total_frags)
summary(asp_counts_df$GM12878_H3K9me3[asp_counts_df$GM12878_H3K9me3$significance != "N.S.",]$total_frags)



asp_chr_sum <- lapply(names(peak_set), function(x) {
  pk_sum <- table(peak_set[[x]]$chrom) %>%
    as.data.frame() %>%
    mutate(Var1 = factor(Var1, levels = paste0("chr", c(1:22, "X")))) %>%
    `colnames<-`(c("chrom", "npeaks"))
  
  asp_sum <- table(asp_df[[x]][, c("chrom", "type")]) %>%
    as.data.frame() %>%
    tidyr::spread(type, Freq) %>%
    dplyr::mutate(nASP = pat + mat) %>%
    dplyr::mutate(chrom = factor(chrom, levels = paste0("chr", c(1:22, "X")))) %>%
    dplyr::arrange(chrom)
  
  asp_sum <- left_join(asp_sum, pk_sum, by = "chrom")
  asp_sum$ASP_ratio <- asp_sum$nASP / asp_sum$npeaks
  asp_sum$pASP_ratio <- asp_sum$pat / asp_sum$npeaks
  asp_sum$mASP_ratio <- asp_sum$mat / asp_sum$npeaks
  
  return(asp_sum)
})
names(asp_chr_sum) <- names(peak_set)

asp_chr_sum$GM12878_LaminB1



asp_chr_sum_m <- lapply(asp_chr_sum, function(df) reshape2::melt(df[,c(1,6:8)]))


asp_chr_sum_m$GM12878_LaminB1$chrom <- factor(asp_chr_sum_m$GM12878_LaminB1$chrom, levels = paste0("chr",c(1:22, "X")))
asp_chr_sum_m$GM12878_H3K9me3$chrom <- factor(asp_chr_sum_m$GM12878_H3K9me3$chrom, levels = paste0("chr",c(1:22, "X")))



asp_chr_sum_m$GM12878_LaminB1$chrom





library(patchwork)
p1 <- ggplot(asp_chr_sum_m$GM12878_LaminB1, aes(x = chrom, y = value * 100)) +
  geom_col(aes(fill = chrom == "chrX")) +
  facet_wrap(vars(variable), ncol = 1) +
  theme_bw() +
  xlab("") + labs(title = "GM12878 LaminB1") +
  ylab("ASP/(total peaks) %") +
  scale_fill_manual("", values = brewer.pal(9, "Set1")[c(9, 1)]) +
  scale_x_discrete(label = c(1:22, "X"), breaks = paste0("chr",c(1:22, "X")), drop = F) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )
p1



p1 <- ggplot(asp_chr_sum_m$GM12878_H3K9me3, aes(x = chrom, y = value * 100)) +
  geom_col(aes(fill = chrom == "chrX")) +
  facet_wrap(vars(variable), ncol = 1) +
  theme_bw() +
  xlab("") + labs(title = "GM12878 H3K9me3") +
  ylab("ASP/(total peaks) %") +
  scale_fill_manual("", values = brewer.pal(9, "Set1")[c(9, 1)]) +
  scale_x_discrete(label = c(1:22, "X")) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )
p1



asp_chr_sum_absolute_count <- lapply(asp_chr_sum, function(df) reshape2::melt(df[,c(1,2:4)]))

p1 <- ggplot(asp_chr_sum_absolute_count$GM12878_LaminB1, aes(x = chrom, y = value)) +
  geom_col(aes(fill = chrom == "chrX")) +
  facet_wrap(vars(variable), ncol = 1) +
  theme_bw() +
  xlab("") + labs(title = "GM12878 LaminB1") +
  ylab("absolute ASP count") +
  scale_fill_manual("", values = brewer.pal(9, "Set1")[c(9, 1)]) +
  scale_x_discrete(label = c(1:22, "X"), breaks = paste0("chr",c(1:22, "X")), drop = F) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )
p1



p1 <- ggplot(asp_chr_sum_absolute_count$GM12878_H3K9me3, aes(x = chrom, y = value)) +
  geom_col(aes(fill = chrom == "chrX")) +
  facet_wrap(vars(variable), ncol = 1) +
  theme_bw() +
  xlab("") + labs(title = "GM12878 H3K9me3") +
  ylab("absolute ASP count") +
  scale_fill_manual("", values = brewer.pal(9, "Set1")[c(9, 1)]) +
  scale_x_discrete(label = c(1:22, "X")) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )
p1





asp_counts_flt <- lapply(asp_counts_df, function(df) {
  df %>%
    filter(type != "none") %>%
    mutate(strand = ".") %>%
    ungroup() %>%
    select(chrom, start, end, type, pat_ratio, strand)
})



write.table(asp_counts_flt$GM12878_LaminB1[asp_counts_flt$GM12878_LaminB1$type == "pat",],'./GM12878_LaminB1_ASP_paternal.bed',row.names = F,col.names = F,sep = '\t',quote = F)
write.table(asp_counts_flt$GM12878_LaminB1[asp_counts_flt$GM12878_LaminB1$type == "mat",],'./GM12878_LaminB1_ASP_maternal.bed',row.names = F,col.names = F,sep = '\t',quote = F)

write.table(asp_counts_flt$GM12878_H3K9me3[asp_counts_flt$GM12878_H3K9me3$type == "pat",],'./GM12878_H3K9me3_ASP_paternal.bed',row.names = F,col.names = F,sep = '\t',quote = F)
write.table(asp_counts_flt$GM12878_H3K9me3[asp_counts_flt$GM12878_H3K9me3$type == "mat",],'./GM12878_H3K9me3_ASP_maternal.bed',row.names = F,col.names = F,sep = '\t',quote = F)


asp_gr <- lapply(asp_counts_flt, function(df){
  lapply(c("pat","mat"), function(x){
    dt2gr(df[df$type == x,])
  }) %>% `names<-`(c("pat","mat"))
})


#==========Functional elements enrichment==========
library("LOLA")

lolatest <- function(lola,mc.cores = 5){
  require(pbmcapply)
  lola <- as.data.frame(lola)
  rownames(lola) <- lola$description
  lola.split <- split(lola,lola$description)
  test.res <- pbmclapply(lola.split,function(x){
    fisher.test(matrix(as.integer(x[,c("support","b","c","d")]),nrow = 2))
  },mc.cores = mc.cores)
  
  res <- lapply(test.res, function(x){
    c(x$p.value,x$estimate,x$conf.int[1],x$conf.int[2])
  }) %>% Reduce(rbind,.) %>%
    as.data.frame() %>%
    `colnames<-`(c("pval","oddsRatio","confIntLow","confIntHigh"))
  res <- cbind(lola[names(lola.split),
                    c("collection","description","support",
                      "b","c","d","filename","size")],res)
  res$LogPadj <- -log10(p.adjust(res$pval,method = "BH"))
  return(res)
}


regionDB <- loadRegionDB("/data/ass/genome_file/LOLA_regionDB/", useCache = T, collections = "hg38")
regionDB$regionAnno

region_sele <- c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only",
                 "intragenic","intergenic","Exon","Intron",'noncoding',
                 "promoter",'promoterUp1000Down500','HCP','ICP','LCP',
                 'enhancer',
                 'CGI','promoterCGI','nonpromoterCGI',
                 "LINE","L1","L2","SINE","Alu","MIR","LTR",'ERV1','ERVK','ERVL','ERVL-MaLR','Low_complexity','Satellite','Simple_repeat','IAP'
)


rgb2hex <- function(rgb) {
  r = rgb[1]
  g = rgb[2]
  b = rgb[3]
  hex_value <- paste0("#", sprintf("%02X%02X%02X", r, g, b))
  return(hex_value)
}

lapply(list(c(255, 0, 0),c(255, 167, 0),c(255, 205, 0),c(255, 170, 170),c(0, 176, 240)), rgb2hex)



region_sele_color <- setNames(c('#FF0000','#FFA700','#FFCD00','#FFAAAA','#00B0F0',
                                brewer.pal(11,"Set3")[c(4,3,7,5,1)],
                                brewer.pal(9,"Spectral")[c(1:5)],
                                brewer.pal(9,"YlOrRd")[5],
                                brewer.pal(9,"PuOr")[c(1:3)],
                                brewer.pal(9,"Blues")[c(9:7)],
                                brewer.pal(9,"Reds")[c(8:6)],
                                brewer.pal(9,"Greens")[c(9:5)],
                                brewer.pal(8,"Accent")[c(8,5,7,6)]),
                              region_sele)


K9_mASP_Func_lola <- runLOLA(
  userSets = asp_gr$GM12878_H3K9me3$mat[asp_gr$GM12878_H3K9me3$mat@seqnames == "chrX"],
  userUniverse = dt2gr(peak_set$GM12878_H3K9me3[peak_set$GM12878_H3K9me3$chr == "chrX",]),
  regionDB, cores = 5
)

K9_mASP_Func_lola2 <- lolatest(K9_mASP_Func_lola)
head(K9_mASP_Func_lola2)

K9_mASP_Func_lola2[K9_mASP_Func_lola2$description == "L1",]


K9_mASP_Func_lola2 <-  K9_mASP_Func_lola2 %>%
  filter(description %in% region_sele) %>%
  mutate(description = factor(description, levels = region_sele))
colnames(K9_mASP_Func_lola2)


ggplot(K9_mASP_Func_lola2, aes(x = description, y = log2(oddsRatio), fill = description)) +
  geom_bar(position = position_dodge(width = 0.9), aes(y = log2(oddsRatio)), stat = "identity") +
  geom_errorbar(aes(ymin = log2(confIntLow), ymax = log2(confIntHigh)), position = position_dodge(width = 0.9), width = 0.5) +
  geom_vline(xintercept = c(5.5,11.5,14.5,15.5,18.5), linetype = "dashed") +
  theme_bw() + labs(title = "LOLA enrichment of GM12878 H3K9me3 mASP") +
  scale_fill_manual(values = region_sele_color) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


ggplot(K9_mASP_Func_lola2, aes(x = description, y = support/(support + c), fill = description)) +
  geom_bar(position = position_dodge(width = 0.9), aes(y = support/(support + c)), stat = "identity") +
  geom_vline(xintercept = c(5.5,11.5,14.5,15.5,18.5), linetype = "dashed") +
  theme_bw() + labs(title = "Percentage of mASP overlapping with genomic elements") +
  ylab("Percentage of peaks") +
  scale_fill_manual(values = region_sele_color) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )




LB1_mASP_Func_lola <- runLOLA(
  userSets = asp_gr$GM12878_LaminB1$mat[asp_gr$GM12878_LaminB1$mat@seqnames == "chrX"],
  userUniverse = dt2gr(peak_set$GM12878_LaminB1[peak_set$GM12878_LaminB1$chr == "chrX",]),
  regionDB, cores = 5
)

LB1_mASP_Func_lola2 <- lolatest(LB1_mASP_Func_lola)
head(LB1_mASP_Func_lola2)

LB1_mASP_Func_lola2[LB1_mASP_Func_lola2$description == "L1",]


LB1_mASP_Func_lola2 <-  LB1_mASP_Func_lola2 %>%
  filter(description %in% region_sele) %>%
  mutate(description = factor(description, levels = region_sele))
colnames(LB1_mASP_Func_lola2)


ggplot(LB1_mASP_Func_lola2, aes(x = description, y = log2(oddsRatio), fill = description)) +
  geom_bar(position = position_dodge(width = 0.9), aes(y = log2(oddsRatio)), stat = "identity") +
  geom_errorbar(aes(ymin = log2(confIntLow), ymax = log2(confIntHigh)), position = position_dodge(width = 0.9), width = 0.5) +
  geom_vline(xintercept = c(5.5,11.5,14.5,15.5,18.5), linetype = "dashed") +
  theme_bw() + labs(title = "LOLA enrichment of GM12878 LaminB1 mASP") +
  scale_fill_manual(values = region_sele_color) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )


ggplot(LB1_mASP_Func_lola2, aes(x = description, y = support/(support + c), fill = description)) +
  geom_bar(position = position_dodge(width = 0.9), aes(y = support/(support + c)), stat = "identity") +
  geom_vline(xintercept = c(5.5,11.5,14.5,15.5,18.5), linetype = "dashed") +
  theme_bw() + labs(title = "Percentage of mASP overlapping with genomic elements") +
  ylab("Percentage of peaks") +
  scale_fill_manual(values = region_sele_color) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
