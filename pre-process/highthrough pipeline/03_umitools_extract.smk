### umi-tools extract cell barcode
# input samples and chunks
SAMPLES = ["O94I10"] # ["sample1", "sample2", "sample3"]
CHUNKS = [f"{i:03d}" for i in range(1, 26)]  # 001 to 025
bc_regex = '(?P<cell_1>.{6}){e<=1}(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=6}(?P<cell_2>.{6}){e<=1}(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=7}(?P<cell_3>.{6}){e<=1}(?P<discard_3>AGATGTGTATAAGAGACAG){e<=4}.*'
whitelist = "./NanoCode/sci_whitelist.txt"
rule all:
    input:
        expand("demultiplex/trim/{sample}_{chunk}_tagged_1.fastq.gz",
               sample=SAMPLES,    
               chunk=CHUNKS)      
# demultiplex barcodes and UMIs in 5'
# The 'Utilities.py' of umi_tools should be modified according to https://github.com/CGATOxford/UMI-tools/pull/447 to avoid raising error when ignore umi
rule umitools_run1:
    input:
        fq = "demultiplex/trim/split/{sample}_{chunk}_trimmed.fastq.gz",
        whitelist = whitelist,
    output:
        tagged = temp("demultiplex/tmp_file/{sample}_{chunk}_tagged_1.fastq.gz"),
        filtered = temp("demultiplex/tmp_file/{sample}_{chunk}_filtered_1.fastq.gz"),
    params:
        lambda wildcards: bc_regex,
    log:
        "demultiplex/logs/{sample}_{chunk}_umitools_run1.log",
    shell: """
        # umi_tools version 1.1.5
        set +u; source activate nanopore; set -u
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern='{params}' \
            --whitelist={input.whitelist} \
            --filtered-out {output.filtered} \
            -I {input.fq} \
            -L {log} -S {output.tagged}
    """
# save files
rule preserve_files:
    input:
        run1_tagged = rules.umitools_run1.output.tagged,
        run1_filtered = rules.umitools_run1.output.filtered,
    output:
        out_run1_tagged = "demultiplex/trim/{sample}_{chunk}_tagged_1.fastq.gz",
        out_run1_filtered = "demultiplex/trim/{sample}_{chunk}_filtered_1.fastq.gz",
    shell:
        """
        cp {input.run1_tagged} {output.out_run1_tagged}
        cp {input.run1_filtered} {output.out_run1_filtered}
        """


### reverse sequence
rule all:
    input:
        expand("demultiplex/trim/{sample}_{chunk}_tagged_rc_1.fastq.gz",
               sample=SAMPLES,   
               chunk=CHUNKS)     
# demultiplex barcodes and UMIs in 5'
# The 'Utilities.py' of umi_tools should be modified according to https://github.com/CGATOxford/UMI-tools/pull/447 to avoid raising error when ignore umi
# make reverse-complement of the run1 results,
# Since `umi_tools extract` is very slow when processing 3' barcodes
rule reverse_comp_tagged:
    input:
        "demultiplex/trim/{sample}_{chunk}_tagged_1.fastq.gz",
    output:
        temp("demultiplex/tmp_file/{sample}_{chunk}_tagged_rc_1.fastq.gz"),
    threads: 35,
    shell: """
        set +u; source activate nanopore; set -u
        seqkit seq -r -p -j {threads} {input} | pigz > {output}
    """

rule reverse_comp_filtered:
    input:
        "demultiplex/trim/{sample}_{chunk}_filtered_1.fastq.gz",
    output:
        temp("demultiplex/tmp_file/{sample}_{chunk}_filtered_rc_1.fastq.gz"),
    threads: 35,
    shell: """
        set +u; source activate nanopore; set -u
        seqkit seq -r -p -j {threads} {input} | pigz > {output}
    """

# save files
rule preserve_files:
    input:
        rc_tagged = rules.reverse_comp_tagged.output,
        rc_filtered = rules.reverse_comp_filtered.output,
    output:
        out_rc_tagged = "demultiplex/trim/{sample}_{chunk}_tagged_rc_1.fastq.gz",
        out_rc_filtered = "demultiplex/trim/{sample}_{chunk}_filtered_rc_1.fastq.gz",
    shell:
        """
        cp {input.rc_tagged} {output.out_rc_tagged}
        cp {input.rc_filtered} {output.out_rc_filtered}
        """


### umi-tools extract in reverse seq
rule all:
    input:
        expand("demultiplex/trim/{sample}_{chunk}_tagged_rc_1_tagged_2.fastq.gz",
               sample=SAMPLES,    
               chunk=CHUNKS)      
# demultiplex barcodes and UMIs in 5'
# The 'Utilities.py' of umi_tools should be modified according to https://github.com/CGATOxford/UMI-tools/pull/447 to avoid raising error when ignore umi

# extract the 3' barcodes and UMIs
rule umitools_run2_tagged:
    input:
        fq = "demultiplex/trim/{sample}_{chunk}_tagged_rc_1.fastq.gz",
        whitelist = whitelist,
    output:
        tagged=temp("demultiplex/tmp_file/{sample}_{chunk}_tagged_rc_1_tagged_2.fastq.gz"),
        filtered=temp("demultiplex/tmp_file/{sample}_{chunk}_tagged_rc_1_filtered_2.fastq.gz"),
    params:
        lambda wildcards: bc_regex,
    log:
        "demultiplex/logs/{sample}_{chunk}_umitools_run2_tagged.log",
    shell: """
        # umi_tools version 1.1.5
        set +u; source activate nanopore; set -u
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern='{params}' \
            --whitelist={input.whitelist} \
            --filtered-out {output.filtered} \
            -I {input.fq} \
            -L {log} -S {output.tagged}
    """

rule umitools_run2_filtered:
    input:
        fq = "demultiplex/trim/{sample}_{chunk}_filtered_rc_1.fastq.gz",
        whitelist = whitelist,
    output:
        tagged=temp("demultiplex/tmp_file/{sample}_{chunk}_filtered_rc_1_tagged_2.fastq.gz"),
        filtered=temp("demultiplex/tmp_file/{sample}_{chunk}_filtered_rc_1_filtered_2.fastq.gz"),
    params:
        lambda wildcards: bc_regex,
    log:
        "demultiplex/logs/{sample}_{chunk}_umitools_run2_filtered.log",
    shell: """
        # umi_tools version 1.1.5
        set +u; source activate nanopore; set -u
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern='{params}' \
            --whitelist={input.whitelist} \
            --filtered-out {output.filtered} \
            -I {input.fq} \
            -L {log} -S {output.tagged}
    """
# save files
rule preserve_files:
    input:
        run2_tagged_tagged = rules.umitools_run2_tagged.output.tagged,
        run2_tagged_filtered = rules.umitools_run2_tagged.output.filtered,
        run2_filtered_tagged = rules.umitools_run2_filtered.output.tagged,
        run2_filtered_filtered = rules.umitools_run2_filtered.output.filtered,
    output:
        out_run2_tagged_tagged = "demultiplex/trim/{sample}_{chunk}_tagged_rc_1_tagged_2.fastq.gz",
        out_run2_tagged_filtered = "demultiplex/trim/{sample}_{chunk}_tagged_rc_1_filtered_2.fastq.gz",
        out_run2_filtered_tagged = "demultiplex/trim/{sample}_{chunk}_filtered_rc_1_tagged_2.fastq.gz",
        out_run2_filtered_filtered = "demultiplex/trim/{sample}_{chunk}_filtered_rc_1_filtered_2.fastq.gz",
    shell:
        """
        cp {input.run2_tagged_tagged} {output.out_run2_tagged_tagged}
        cp {input.run2_tagged_filtered} {output.out_run2_tagged_filtered}
        cp {input.run2_filtered_tagged} {output.out_run2_filtered_tagged}
        cp {input.run2_filtered_filtered} {output.out_run2_filtered_filtered}
        """

###









