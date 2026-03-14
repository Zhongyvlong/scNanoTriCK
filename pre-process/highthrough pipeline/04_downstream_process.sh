#!/bin/bash

##### dual accordance function
accordance() {
    base=$(basename "$1" "_tagged_rc_1_tagged_2.fastq.gz")

    # check files
    if [ ! -f "demultiplex/trim/${base}_tagged_rc_1_tagged_2.fastq.gz" ]; then
        echo "ERROR: File ${base}_tagged_rc_1_tagged_2.fastq.gz not found" >&2
        return 1
    fi

    # fastq information
    zcat "demultiplex/trim/${base}_tagged_rc_1_tagged_2.fastq.gz" | awk -v filter="demultiplex/trim/${base}_tagged_1_tagged_2_discor.fastq" '{{
        if(NR % 4 == 1) {{
            split($1, name, " ")
            split(name[1], parts, "_")
            if (parts[2] == parts[4]) {{
                print # output read
                getline; print
                getline; print
                getline; print
            }} else {{
                print > filter # output filter files
                getline; print > filter
                getline; print > filter
                getline; print > filter
            }}
        }}
    }}' | pigz > "demultiplex/trim/${base}_tagged_1_tagged_2_accor.fastq.gz"

    pigz "demultiplex/trim/${base}_tagged_1_tagged_2_discor.fastq"
}
export -f accordance

# running
find demultiplex/trim/ -name "*_tagged_rc_1_tagged_2.fastq.gz" | parallel -j 35 accordance {}

##### output final files
find demultiplex/trim/ -name "*_tagged_1_tagged_2_accor.fastq.gz" | parallel -j 15 '
        base=$(basename {} "_tagged_1_tagged_2_accor.fastq.gz")
        cat demultiplex/trim/${base}_tagged_1_tagged_2_accor.fastq.gz demultiplex/trim/${base}_tagged_rc_1_filtered_2.fastq.gz demultiplex/trim/${base}_filtered_rc_1_tagged_2.fastq.gz |
        seqkit seq -r -p -j 35 | pigz > demultiplex/trim/${base}_bc_merged.fastq.gz
'

find demultiplex/trim/ -name "*_filtered_rc_1_filtered_2.fastq.gz" | parallel -j 15 '
        base=$(basename {} "_filtered_rc_1_filtered_2.fastq.gz")
        cat demultiplex/trim/${base}_filtered_rc_1_filtered_2.fastq.gz demultiplex/trim/${base}_tagged_1_tagged_2_discor.fastq.gz |
        seqkit seq -r -p -j 35 | pigz > demultiplex/trim/${base}_bc_filtered.fastq.gz
'

for i in `ls demultiplex/raw_data/split/*.fastq.gz`
do base=$(basename $i ".fastq.gz")
seqkit stats -a -T -j 35 demultiplex/raw_data/split/${base}.fastq.gz demultiplex/trim/split/${base}_trimmed.fastq.gz demultiplex/trim/{${base}_tagged_1.fastq.gz,${base}_filtered_1.fastq.gz,${base}_filtered_rc_1_tagged_2.fastq.gz,${base}_filtered_rc_1_filtered_2.fastq.gz,${base}_tagged_1_tagged_2_accor.fastq.gz,${base}_tagged_1_tagged_2_discor.fastq.gz,${base}_bc_merged.fastq.gz} > demultiplex/trim/${base}_fastq_stats.txt
done

echo -e "file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\tQ1\\tQ2\\tQ3\\tsum_gap\\tN50\\tN50_num\\tQ20(%)\\tQ30(%)\\tAvgQual\\tGC(%)\\tsum_n" > demultiplex/trim/O94I10_fastq_stats.txt
find demultiplex/trim/ -name "*_fastq_stats.txt" | parallel -j 35 '
        base=$(basename {} "_fastq_stats.txt")
        cat demultiplex/trim/${base}_fastq_stats.txt | grep -P -v "^file\\tformat" >> demultiplex/trim/O94I10_fastq_stats.txt
'

##### mapping reads to genome
mkdir demultiplex/mapping
find demultiplex/trim/ -name "*_bc_merged.fastq.gz" | parallel -j 35 '
        base=$(basename {} "_bc_merged.fastq.gz")
        minimap2 -t 10 \
            -ax map-ont \
            --secondary=no \
            minimap2_index/mm10/mm10_ont.mmi \
            demultiplex/trim/${base}_bc_merged.fastq.gz 2>demultiplex/logs/${base}_mm2.log |
        samtools view -o demultiplex/mapping/${base}_mm2.bam -
'

add_barcode() {
    base=$(basename "$1" "_mm2.bam")
    if [ ! -f "demultiplex/mapping/${base}_mm2.bam" ]; then
        echo "ERROR: File ${base}_mm2.bam not found" >&2
        return 1
    fi
    cat <(samtools view -H demultiplex/mapping/${base}_mm2.bam) \
            <(samtools view demultiplex/mapping/${base}_mm2.bam |
            awk -vOFS="\t" '{
                    split($1, name, "_")
                    print $0,"CB:Z:"name[2],"RX:Z:"name[3]
            }') |
        samtools view -bS - > demultiplex/mapping/${base}_bc.bam
}
export -f add_barcode
find demultiplex/mapping/ -name "*_mm2.bam" | parallel -j 35 add_barcode {}

map_stats() {
    base=$(basename "$1" "_bc.bam")
    if [ ! -f "demultiplex/mapping/${base}_bc.bam" ]; then
        echo "ERROR: File ${base}_bc.bam not found" >&2
        return 1
    fi
    bedtools bamtobed -i demultiplex/mapping/${base}_bc.bam |
        awk '{split($4, name, "_");print $1,$2,$3,name[2],$5,$6
        }' OFS="\t" |
        awk '{if($5>=20){q20=1}else{q20=0};if($5>30){q30=1}else{q30=0}
        };{print $4,$3-$2,q20,q30
        }' OFS="\t" | sort -k1 |
        bedtools groupby -g 1 -c 1,2,2,2,2,2,3,4 -o count,sum,min,max,mean,median,sum,sum |
        awk '{print $1,$2,$3,$4,$5,$6,$7,$8/$2,$9/$2}' OFS="\t" > demultiplex/mapping/${base}_bc.bam.stats
}
export -f map_stats
find demultiplex/mapping/ -name "*_bc.bam" | parallel -j 35 map_stats {}

Rscript ./NanoCode/merge_bam_stats.r "demultiplex/mapping" "_bc.bam.stats" "demultiplex/mapping/O94I10_bc.bam.stats"

#filter read
# filter alignments, only unique mapped reads with clipping length < 150 were remained.
mkdir tmp
find demultiplex/mapping/ -name "*_bc.bam" | parallel -j 55 '
        base=$(basename {} "_bc.bam")
        samtools view -ShuF 2308 -q 30 -e "sclen < 150 && hclen < 150" {} |
        samtools sort -T tmp/${base} - -o demultiplex/mapping/${base}_bc_q30.bam
'
#summary q30 bam
find demultiplex/mapping/ -name "*_bc_q30.bam" | parallel -j 55 '
        base=$(basename {} "_bc_q30.bam")
        echo {}
        bedtools bamtobed -i {} | awk -v OFS="\t" '\''{{split($4, name, "_")} {print $1,$2,$3,name[2],$5,$6}}'\'' |
        awk -v OFS="\t" '\''{
                if($5 >= 20) {q20=1} else {q20=0}
                if($5 > 30) {q30=1} else {q30=0}
                {print $4,$3-$2,q20,q30}}'\'' |sort -k1 |
        bedtools groupby -g 1 -c 1,2,2,2,2,2,3,4 -o count,sum,min,max,mean,median,sum,sum |
        awk -v OFS="\t" '\''{print $1,$2,$3,$4,$5,$6,$7,$8/$2,$9/$2}'\'' > demultiplex/mapping/${base}_bc_q30.bam.stats
'
Rscript ./NanoCode/merge_bam_stats.r "demultiplex/mapping" "_bc_q30.bam.stats" "demultiplex/mapping/O94I10_bc_q30.bam.stats"
samtools merge -@ 30 -o demultiplex/mapping/O94I10_bc_q30.bam demultiplex/mapping/*_bc_q30.bam

##### deduplicate reads
# dedup with umi_tools, since UMIs tend to be overestimated due to the sequencing errors of Nanopore, and each genomic loci usaually was single copy in scATAC,
# so we only consider the mapped coordinates of alignments when deduplication.
sample="O94I10"
mkdir demultiplex/fragment
samtools index "demultiplex/mapping/${sample}_bc_q30.bam"
umi_tools dedup --extract-umi-method=tag --per-cell --ignore-umi --cell-tag=CB --temp-dir=tmp \
            -I "demultiplex/mapping/${sample}_bc_q30.bam" -L demultiplex/logs/"${sample}_dedup.log" -S "demultiplex/mapping/${sample}_duprm.bam"
bedtools bamtobed -i demultiplex/mapping/${sample}_duprm.bam | awk '{split($4, name, "_");print $1,$2,$3,name[2],$5,$6}' OFS="\t" |awk '!seen[$1,$2,$3,$4]++' > demultiplex/fragment/${sample}_fragment.bed
bgzip demultiplex/fragment/${sample}_fragment.bed
tabix demultiplex/fragment/${sample}_fragment.bed.gz
zcat demultiplex/fragment/${sample}_fragment.bed.gz |
        awk '{if($5>=20){q20=1}else{q20=0};if($5>30){q30=1}else{q30=0}
        };{print $4,$3-$2,q20,q30}' OFS="\t" | sort -k1 |
        bedtools groupby -g 1 -c 1,2,2,2,2,2,3,4 -o count,sum,min,max,mean,median,sum,sum |
        awk '{print $1,$2,$3,$4,$5,$6,$7,$8/$2,$9/$2}' OFS="\t" > demultiplex/mapping/${sample}_dedup.bam.stats
# flank and slop
mkdir demultiplex/bigwig
zcat demultiplex/fragment/${sample}_fragment.bed.gz | grep ^chr | bedtools flank -b 1 -i stdin -g mm10.chrom.size > demultiplex/fragment/${sample}_flank.bed
cat demultiplex/fragment/${sample}_flank.bed | sort -k1,1 -k2,2 | bedtools genomecov -bg -i stdin -g mm10.chrom.size > demultiplex/bigwig/${sample}_flank.bedGraph
bedGraphToBigWig demultiplex/bigwig/${sample}_flank.bedGraph mm10.chrom.size demultiplex/bigwig/${sample}_flank.bw
pigz demultiplex/bigwig/${sample}_flank.bedGraph

bedtools slop -b 75 -i demultiplex/fragment/${sample}_flank.bed -g mm10.chrom.size > demultiplex/fragment/${sample}_slop.bed
cat demultiplex/fragment/${sample}_slop.bed | sort -k1,1 -k2,2 | bedtools genomecov -bg -i stdin -g mm10.chrom.size > demultiplex/bigwig/${sample}_slop.bedGraph
bedGraphToBigWig demultiplex/bigwig/${sample}_slop.bedGraph mm10.chrom.size demultiplex/bigwig/${sample}_slop.bw
pigz demultiplex/bigwig/${sample}_slop.bedGraph

##### output for signac
#!usr/bin/
#pipline
sample="O94I10"
awk -v OFS="\t" -F'\t' '{
    key = $1 "\t" $2 "\t" $3 "\t" $4
    count[key]++
}
END {
    for (key in count) {
        print key "\t" count[key]
    }
}' demultiplex/fragment/${sample}_flank.bed | sort -k 1,1 -k 2,2n |awk '{print$1"\t"$2"\t"$3"\t'$sample'#"$4"\t"$5}' > demultiplex/fragment/${sample}_signac_fragfromflank.bed
bgzip -@ 50 demultiplex/fragment/${sample}_signac_fragfromflank.bed -o demultiplex/fragment/${sample}_signac_fragfromflank.bed.gz
tabix -@ 50 demultiplex/fragment/${sample}_signac_fragfromflank.bed.gz




