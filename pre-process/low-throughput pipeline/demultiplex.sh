conda activate nanopore


minimap2_index="/data/ass/genome_file/gencode_ensembl/minimap2_index/hg38/hg38_ont.mmi"


#generate barcode index file
bash dual_barcode.sh bc_index

mkdir raw_data
nanoplexer -b barcode.fa -t 80 -p raw_data/ pass.fastq.gz
echo "Demultiplex finished! Compressing data ..."
rm raw_data/unclassified.fastq


#cut adaptor
mkdir trim
find raw_data/ -name "*fastq.gz" | parallel -j 55 '
	base=$(basename {} ".fastq.gz")
	pychopper -m edlib \
		-b Tn5_seq.fa \
		-c primer_config.txt \
		-r trim/${base}_trimed_report.pdf \
		-w trim/${base}_rescued.fastq \
		-Q 7 -z 10 -t 1 {} trim/${base}_trimed.fastq 2> trim/${base}_pychopper.log
	cat trim/${base}_rescued.fastq >> trim/${base}_trimed.fastq
	pigz trim/${base}_trimed.fastq
'
rm trim/*fastq


#summary raw fastq information
find merged_raw_data/ -name "*fastq.gz" | parallel -j 55 '
	base=$(basename {} ".fastq.gz")
	seqkit stats -a -T -j 1 {} >> raw_fastq.stats.tmp
'
awk '{if(FNR==1 || FNR%2==0) {print $0}}' raw_fastq.stats.tmp |
	sed 's/raw_data\///' | sed 's/'.'fastq'.'gz//' > raw_fastq.stats
rm raw_fastq.stats.tmp


#summary trim fastq information
find merged_trim/ -name "*fastq.gz" | parallel -j 15 '
	base=$(basename {} ".fastq.gz")
	seqkit stats -a -T -j 1 {} >> trim_fastq.stats.tmp
'
awk '{if(FNR==1 || FNR%2==0) {print $0}}' trim_fastq.stats.tmp |
	sed 's/trim\///' | sed 's/_trimed'.'fastq'.'gz//' > trim_fastq.stats
rm trim_fastq.stats.tmp


#summary bulk raw/trim information
raw_used=`cat raw_fastq.stats | sed 1d | grep -v un | awk '{{sum+=$5}};END{{print sum}}'`
raw_total=`cat raw_fastq.stats | sed 1d | awk '{{sum+=$5}};END{{print sum}}'`
raw_un=`cat raw_fastq.stats | sed 1d | grep un | awk '{{sum+=$5}};END{{print sum}}'`
trim_total=`cat trim_fastq.stats | sed 1d | awk '{{sum+=$5}};END{{print sum}}'`
trim_used=`cat trim_fastq.stats | sed 1d | grep -v un | awk '{{sum+=$5}};END{{print sum}}'`
trim_un=`cat trim_fastq.stats | sed 1d | grep un | awk '{{sum+=$5}};END{{print sum}}'`
echo -e "raw_total\t$raw_total" >> datasize_summary;
echo -e "raw_classified\t$raw_used" >> datasize_summary
echo -e "raw_unclassified\t$raw_un" >> datasize_summary
echo -e "trim_total\t$trim_total" >> datasize_summary
echo -e "trim_classified\t$trim_used" >> datasize_summary
echo -e "trim_unclassified\t$trim_un" >> datasize_summary


#minimap2 mapping
mkdir mapping
find trim/ -name "*fastq.gz" | parallel -j 64 '
	base=$(basename {} ".fastq.gz")
	minimap2 -t 1 \
		-ax map-ont \
		--secondary=no \
		-K 500M \
		$minimap2_index \
		trim/${base}.fastq.gz 2> mapping/${base}_mm2.log | samtools view -o mapping/${base}_mm2.bam
'


#summary map bam
find mapping/ -name "*_mm2.bam" | parallel -j 55 '
	base=$(basename {} "_trimed_mm2.bam")
	echo {}
	bedtools bamtobed -i {} | awk -v OFS="\t" '\''{{split($4, name, "|")} {print $1,$2,$3,name[2],$5,$6}}'\'' |
	awk -v OFS="\t" '\''{
		if($5 >= 20) {q20=1} else {q20=0}
		if($5 > 30) {q30=1} else {q30=0}
		if($6 == "+") {str=1} else {str=0}
		{print "1",$4,$3-$2,q20,q30,str}}'\'' |
	bedtools groupby -g 1 -c 1,2,3,3,3,3,3,4,5,6 \
		-o count,count_distinct,sum,min,max,mean,median,sum,sum,sum |
	awk -v OFS="\t" -v c=${base} '\''{print c,"map",$2,$3,$4,$5,$6,$7,$8,$9/$2,$10/$2,$11/$2}'\'' > mapping/${base}_map.bam.stats
'


#filter read
#soft/hard clipping length < 150
mkdir tmp
find mapping/ -name "*_mm2.bam" | parallel -j 55 '
	base=$(basename {} "_trimed_mm2.bam")
	samtools view -ShuF 2308 -q 30 -e "sclen < 150 && hclen < 150" {} |
	samtools sort -T tmp/${base} -o mapping/${base}_q30.bam
'

#summary q30 bam
find mapping/ -name "*_q30.bam" | parallel -j 55 '
	base=$(basename {} "_q30.bam")
	echo {}
	bedtools bamtobed -i {} | awk -v OFS="\t" '\''{{split($4, name, "|")} {print $1,$2,$3,name[2],$5,$6}}'\'' |
	awk -v OFS="\t" '\''{
		if($5 >= 20) {q20=1} else {q20=0}
		if($5 > 30) {q30=1} else {q30=0}
		if($6 == "+") {str=1} else {str=0}
		{print "1",$4,$3-$2,q20,q30,str}}'\'' |
	bedtools groupby -g 1 -c 1,2,3,3,3,3,3,4,5,6 \
		-o count,count_distinct,sum,min,max,mean,median,sum,sum,sum |
	awk -v OFS="\t" -v c=${base} '\''{print c,"filter",$2,$3,$4,$5,$6,$7,$8,$9/$2,$10/$2,$11/$2}'\'' > mapping/${base}_filter.bam.stats
'


#remove deduplicated reads
find mapping/ -name "*_q30.bam" | parallel -j 55 '
	base=$(basename {} "_q30.bam")
	samtools rmdup -s {} mapping/${base}_duprm.bam
'


#stat duprm bam
find mapping/ -name "*_duprm.bam" | parallel -j 55 '
	base=$(basename {} "_duprm.bam")
	echo {}
	bedtools bamtobed -i {} | awk -v OFS="\t" '\''{{split($4, name, "|")} {print $1,$2,$3,name[2],$5,$6}}'\'' |
	awk -v OFS="\t" '\''{
		if($5 >= 20) {q20=1} else {q20=0}
		if($5 > 30) {q30=1} else {q30=0}
		if($6 == "+") {str=1} else {str=0}
		{print "1",$4,$3-$2,q20,q30,str}}'\'' |
	bedtools groupby -g 1 -c 1,2,3,3,3,3,3,4,5,6 \
		-o count,count_distinct,sum,min,max,mean,median,sum,sum,sum |
	awk -v OFS="\t" -v c=${base} '\''{print c,"duprm",$2,$3,$4,$5,$6,$7,$8,$9/$2,$10/$2,$11/$2}'\'' > mapping/${base}_duprm.bam.stats
'



#merge stats
echo -e "cell\ttype\talign_counts\tread_counts\tdatasize\tmin_len\tmax_len\tmean_len\tmedian_len\tq20\tq30\torientation" > mapping/map_bam.stats
echo -e "cell\ttype\talign_counts\tread_counts\tdatasize\tmin_len\tmax_len\tmean_len\tmedian_len\tq20\tq30\torientation" > mapping/filter_bam.stats
echo -e "cell\ttype\talign_counts\tread_counts\tdatasize\tmin_len\tmax_len\tmean_len\tmedian_len\tq20\tq30\torientation" > mapping/duprm_bam.stats

for i in `ls mapping/*_trimed_mm2.bam`
do
	base=$(basename ${i} "_trimed_mm2.bam")
	cat mapping/${base}_map.bam.stats >> mapping/map_bam.stats
	cat mapping/${base}_filter.bam.stats >> mapping/filter_bam.stats
	cat mapping/${base}_duprm.bam.stats >> mapping/duprm_bam.stats
done


#generate fragment file, cut site flank bed and cut site slop bed
mkdir fragment
find mapping/ -name "*_duprm.bam" | parallel -j 32 '
	base=$(basename {} "_duprm.bam")
	echo {}
	bedtools bamtobed -i {} | grep "chr" | 
		awk -v OFS="\t" -v CB=${base} '\''{print $1,$2,$3,CB,$5,$6}'\'' | sort -k 1,1 -k 2,2n | awk '\''!seen[$1,$2,$3]++'\'' > fragment/${base}_fragment.bed

	uni_read=$(samtools view -c mapping/${base}_q30.bam)
	dedup_read=$(wc -l <fragment/${base}_fragment.bed)
	echo -e "${base}\t${uni_read}\t${dedup_read}\t$(awk -v u=${uni_read} -v d=${dedup_read} '\''BEGIN{printf "%.4f", (u - d) / u}'\'')" >> duplicated_rate.txt

	awk -v OFS="\t" '\''{if($3-$2 > 1000) {print $0}}'\'' fragment/${base}_fragment.bed > fragment/${base}_flt_fragment.bed

	cat fragment/${base}_fragment.bed | grep "^chr" |
		bedtools flank -b 1 -i stdin -g /data/ass/genome_file/gencode_ensembl/fasta_hg38/hg38.chrom.size > fragment/${base}_flank.bed
	bedtools slop -b 75 -i fragment/${base}_flank.bed -g /data/ass/genome_file/gencode_ensembl/fasta_hg38/hg38.chrom.size > fragment/${base}_slop.bed

	cat fragment/${base}_flt_fragment.bed | grep "^chr" |
		bedtools flank -b 1 -i stdin -g /data/ass/genome_file/gencode_ensembl/fasta_hg38/hg38.chrom.size > fragment/${base}_flt_flank.bed
	bedtools slop -b 75 -i fragment/${base}_flt_flank.bed -g /data/ass/genome_file/gencode_ensembl/fasta_hg38/hg38.chrom.size > fragment/${base}_flt_slop.bed
'



#merge single cell
mkdir merged_result

cat fragment/*flank.bed > merged_results/merged_flank.bed
cat fragment/*slop.bed > merged_results/merged_slop.bed
cat fragment/*fragment.bed |
	awk -v n="merged" -v OFS="\t" '{{print $1,$2,$3,n,$5,$6}}' | sort -Vk 1,1 -k 2,2n > merged_results/merged_fragment.bed
bgzip merged_results/merged_fragment.bed
tabix merged_results/merged_fragment.bed.gz

#filter fragment
zcat merged_results/merged_fragment.bed.gz | awk -v OFS="\t" '$3-$2 > 1000 {print $0}' > merged_results/merged_flt_fragment.bed
bgzip merged_results/merged_flt_fragment.bed
tabix merged_results/merged_flt_fragment.bed.gz





































