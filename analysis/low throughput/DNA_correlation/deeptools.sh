#K562 H3K9me3
multiBigwigSummary BED-file \
    --BED het.call.human.K562.hg38.bed \
    -b K562_H3K9me3_merged_scNanoTriCK-TileSize-100-normMethod-nFrags-ArchR.bw \
       GSM7964462_K562_H3K9me3_scNanoSeq_CUT-Tag.bw \
       ENCFF600IBI_K562_H3K9me3_hg38.bigWig \
    --labels K562_H3K9me3_scNanoTriCK_merged scNanoSeq_CUT-TAG_H3K9me3 ChIP-seq_H3K9me3 \
    -out K562_ChIP-seq_H3K9me3_peak_region_correlation.npz \
    -p 60 \
    --outRawCounts K562_ChIP-seq_H3K9me3_peak_region_correlation_cor.tab


#K562 H3K27me3
multiBigwigSummary BED-file \
    --BED Polycomb_domain_bins_10kb.bed \
    -b K562_H3K27me3_merged_scNanoTriCK-TileSize-100-normMethod-nFrags-ArchR.bw \
       GSM7964461_K562_H3K27me3_scNanoSeq_CUT-Tag.bw \
       ENCFF508LLH_K562_H3K27me3_hg38-TileSize-100-normMethod-nFrags-ArchR.bw \
    --labels K562_H3K27me3_scNanoTriCK_merged scNanoSeq_CUT-TAG_H3K27me3 ChIP-seq_H3K27me3 \
    -out K562_ChIP-seq_H3K27me3_peak_region_correlation.npz \
    -p 60 \
    --outRawCounts K562_ChIP-seq_H3K27me3_peak_region_correlation_cor.tab

#GM12878 H3K9me3
multiBigwigSummary BED-file \
    --BED heterochromatin_domain_bins_100kb.bed \
    -b GM12878_H3K9me3_merged_scNanoTriCK-TileSize-100-normMethod-nFrags-ArchR.bw \
       GSM7964476_GM12878_H3K9me3.bw \
       ENCFF349VRN_log2_enrichment_over_Input_5kb.bigwig \
    --labels GM12878_H3K9me3_scNanoTriCK_merged GM12878_H3K9me3_scNanoSeq_CUT-Tag ChIP-seq_H3K9me3 \
    -out GM12878_H3K9me3_heterochromatin_domain_correlation.npz \
    -p 60 \
    --outRawCounts GM12878_H3K9me3_heterochromatin_domain_correlation_cor.tab
















































