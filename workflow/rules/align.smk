# Ribo-seq: map contaminant-free reads to human + spike-in transcriptome, then
# split the BAM by species. HUMAN_TRANSCRIPTOME_FA is the only thing `mode`
# changes: the RNA-seq-filtered fasta, or the configured one as is.


rule combined_transcriptome:
    input:
        human=HUMAN_TRANSCRIPTOME_FA,
        spike_in=config["spike_in_transcriptome_fa"],
    output:
        f"{RESULTS_DIR}/reference/transcriptome.combined_human_spike_in.fa",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.human} {input.spike_in} > {output}"


rule star_transcript_index:
    input:
        f"{RESULTS_DIR}/reference/transcriptome.combined_human_spike_in.fa",
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
        spike_in_fa=config["spike_in_transcriptome_fa"],
    output:
        human=f"{RESULTS_DIR}/reference/transcripts.human.bed",
        spike_in=f"{RESULTS_DIR}/reference/transcripts.spike_in.bed",
    params:
        script=workflow.source_path("../scripts/split_bed_transcriptome.py"),
    log:
        f"{LOG_DIR}/align/split_bed_transcriptome.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.script} {input.chrom_sizes} {input.human_fa} {input.spike_in_fa} "
        "{output.human} {output.spike_in} 2> {log}"


rule split_bam_transcriptome:
    input:
        bam=f"{RESULTS_DIR}/star/transcriptome/{{sample}}/{{sample}}.bam",
        bed=f"{RESULTS_DIR}/reference/transcripts.{{species}}.bed",
    output:
        bam=f"{RESULTS_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam.bai",
    wildcard_constraints:
        species="human|spike_in",
    conda:
        "../envs/star.yaml"
    shell:
        r"""
        samtools view -b -L {input.bed} {input.bam} > {output.bam}
        samtools index {output.bam}
        """
