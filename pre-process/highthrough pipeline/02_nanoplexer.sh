source activate nanopore
# nanoplexer
mkdir demultiplex/raw_data
nanoplexer -b demultiplex/barcode.fa \
        -t 20 -p demultiplex/raw_data/ \
        HUNAchip/pass.fastq.gz

echo "Demultiplex finished! Compressing data ..."
#rm demultiplex/raw_data/unclassified.fastq
pigz -p 35 demultiplex/raw_data/*.fastq
echo "Compress finished!"

#split files
seqkit split2 -s 500000 -t dna -j 40 -1 demultiplex/raw_data/O93I09.fastq.gz -O demultiplex/raw_data/split

# cutadapt
find demultiplex/raw_data/split/ -name "*fastq.gz" | parallel -j 30 '
        base=$(basename {} ".fastq.gz")
        cutadapt -q 7 -e 0.25 -j 40 -m 100 \
                -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT...AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
                -o demultiplex/trim/split/${base}_trimmed.fastq.gz \
                demultiplex/raw_data/split/${base}.fastq.gz \
                > demultiplex/logs/${base}_cutadapt.log
'

