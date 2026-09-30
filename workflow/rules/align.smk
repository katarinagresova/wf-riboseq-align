# Ribo-seq: map contaminant-free reads to human + yeast transcriptome, then
# split the BAM by species. HUMAN_TRANSCRIPTOME_FA is the only thing `mode`
# changes: the RNA-seq-filtered fasta, or the configured one as is.


# Yeast goes LAST: split_bed_transcriptome relies on that order.
rule combined_transcriptome:
    input:
        human=HUMAN_TRANSCRIPTOME_FA,
        yeast=config["yeast_transcriptome_fa"],
    output:
        f"{RESULTS_DIR}/reference/transcriptome.combined_human_yeast.fa",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.human} {input.yeast} > {output}"


rule star_transcript_index:
    input:
        f"{RESULTS_DIR}/reference/transcriptome.combined_human_yeast.fa",
    output:
        index=directory(f"{RESULTS_DIR}/star_index/transcriptome"),
        chrom_sizes=f"{RESULTS_DIR}/star_index/transcriptome/chrNameLength.txt",
    log:
        f"{LOG_DIR}/star/transcriptome_index.log",
    params:
        log_prefix=f"{LOG_DIR}/star/transcriptome_index.",
    conda:
        "../envs/star.yaml"
    threads: 8
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output.index} --genomeFastaFiles {input} "
        "--genomeSAindexNbases 11 --genomeChrBinNbits 12 "
        "--outFileNamePrefix {params.log_prefix} > {log} 2>&1"


# Multimappers are kept (up to 255 loci); STAR marks one alignment per read
# primary, chosen at random among equally good loci.
rule star_transcript:
    input:
        fastq=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.fastq.gz",
        index=f"{RESULTS_DIR}/star_index/transcriptome",
    output:
        f"{RESULTS_DIR}/star/transcriptome/{{sample}}/{{sample}}.bam",
    params:
        prefix=f"{RESULTS_DIR}/star/transcriptome/{{sample}}/{{sample}}.transcript_",
    conda:
        "../envs/star.yaml"
    threads: 12
    shell:
        r"""
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --outSAMtype BAM Unsorted \
            --outSAMmode NoQS \
            --outSAMattributes NH NM \
            --seedSearchLmax 10 \
            --outFilterMultimapNmax 255 \
            --outFilterMismatchNmax 2 \
            --outFilterMultimapScoreRange 0 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFileNamePrefix {params.prefix} \
            --readFilesIn {input.fastq} \
            --readFilesCommand zcat

        samtools sort -@ {threads} {params.prefix}Aligned.out.bam -o {output}
        rm {params.prefix}Aligned.out.bam
        """


rule split_bed_transcriptome:
    input:
        chrom_sizes=f"{RESULTS_DIR}/star_index/transcriptome/chrNameLength.txt",
        human_fa=HUMAN_TRANSCRIPTOME_FA,
        yeast_fa=config["yeast_transcriptome_fa"],
    output:
        human=f"{RESULTS_DIR}/reference/transcripts.human.bed",
        yeast=f"{RESULTS_DIR}/reference/transcripts.yeast.bed",
    conda:
        "../envs/coreutils.yaml"
    shell:
        r"""
        human_n=$(grep -c '^>' {input.human_fa})
        yeast_n=$(grep -c '^>' {input.yeast_fa})
        total_n=$(wc -l < {input.chrom_sizes})
        if [ $((human_n + yeast_n)) -ne "$total_n" ]; then
            echo "index has $total_n refs, fastas have $human_n human + $yeast_n yeast" >&2
            exit 1
        fi
        awk -v OFS='\t' -v n="$human_n" 'NR <= n {{print $1, 0, $2}}' {input.chrom_sizes} > {output.human}
        awk -v OFS='\t' -v n="$human_n" 'NR > n {{print $1, 0, $2}}' {input.chrom_sizes} > {output.yeast}
        """


rule split_bam_transcriptome:
    input:
        bam=f"{RESULTS_DIR}/star/transcriptome/{{sample}}/{{sample}}.bam",
        bed=f"{RESULTS_DIR}/reference/transcripts.{{species}}.bed",
    output:
        bam=f"{RESULTS_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam.bai",
    wildcard_constraints:
        species="human|yeast",
    conda:
        "../envs/star.yaml"
    shell:
        r"""
        samtools view -b -L {input.bed} {input.bam} > {output.bam}
        samtools index {output.bam}
        """
