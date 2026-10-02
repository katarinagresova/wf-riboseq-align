# Ribo-seq: map contaminant-free reads to human + spike-in transcriptome, then
# split the BAM by species. HUMAN_TRANSCRIPTOME_FA is the only thing `mode`
# changes: the RNA-seq-filtered fasta, or the configured one as is.


rule combined_transcriptome:
    input:
        human=HUMAN_TRANSCRIPTOME_FA,
        spike_in=config["spike_in_transcriptome_fa"],
    output:
        f"{EXP_DIR}/reference/transcriptome.combined_human_spike_in.fa",
    log:
        f"{EXP_LOG_DIR}/align/combined_transcriptome.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.human} {input.spike_in} > {output} 2> {log}"


rule star_transcript_index:
    input:
        f"{EXP_DIR}/reference/transcriptome.combined_human_spike_in.fa",
    output:
        index=directory(f"{EXP_DIR}/star_index/transcriptome"),
        chrom_sizes=f"{EXP_DIR}/star_index/transcriptome/chrNameLength.txt",
    log:
        f"{EXP_LOG_DIR}/star/transcriptome_index.log",
    params:
        log_prefix=f"{EXP_LOG_DIR}/star/transcriptome_index.",
    conda:
        "../envs/star.yaml"
    threads: 8
    resources:
        mem_mb=16000,
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output.index} --genomeFastaFiles {input} "
        "--genomeSAindexNbases 11 --genomeChrBinNbits 12 "
        "--outFileNamePrefix {params.log_prefix} > {log} 2>&1"


# Multimappers are kept (up to 255 loci); STAR marks one alignment per read
# primary, chosen at random among equally good loci. CLEAN_FASTQ: Snakefile.
rule star_transcript:
    input:
        fastq=CLEAN_FASTQ,
        index=f"{EXP_DIR}/star_index/transcriptome",
    output:
        f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.bam",
    log:
        f"{EXP_LOG_DIR}/star/transcriptome/{{sample}}.log",
    params:
        prefix=f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.transcript_",
        star_args=RIBO_TRANSCRIPTOME_STAR_ARGS,
    conda:
        "../envs/star.yaml"
    threads: 12
    resources:
        mem_mb=24000,
    shell:
        r"""
        exec > {log} 2>&1
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --outSAMtype BAM Unsorted \
            --outSAMmode NoQS \
            --outSAMattributes NH NM \
            {params.star_args} \
            --outFileNamePrefix {params.prefix} \
            --readFilesIn {input.fastq} \
            --readFilesCommand zcat

        samtools sort -@ {threads} {params.prefix}Aligned.out.bam -o {output}
        rm {params.prefix}Aligned.out.bam
        """


rule split_bed_transcriptome:
    input:
        chrom_sizes=f"{EXP_DIR}/star_index/transcriptome/chrNameLength.txt",
        human_fa=HUMAN_TRANSCRIPTOME_FA,
        spike_in_fa=config["spike_in_transcriptome_fa"],
        script=workflow.source_path("../scripts/split_bed_transcriptome.py"),
    output:
        human=f"{EXP_DIR}/reference/transcripts.human.bed",
        spike_in=f"{EXP_DIR}/reference/transcripts.spike_in.bed",
    log:
        f"{EXP_LOG_DIR}/align/split_bed_transcriptome.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.chrom_sizes} {input.human_fa} {input.spike_in_fa} "
        "{output.human} {output.spike_in} 2> {log}"


rule split_bam_transcriptome:
    input:
        bam=f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.bam",
        bed=f"{EXP_DIR}/reference/transcripts.{{species}}.bed",
    output:
        bam=f"{EXP_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam",
        bai=f"{EXP_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam.bai",
    wildcard_constraints:
        species="human|spike_in",
    log:
        f"{EXP_LOG_DIR}/align/split_bam/{{species}}/{{sample}}.log",
    conda:
        "../envs/star.yaml"
    shell:
        r"""
        exec 2> {log}
        samtools view -b -L {input.bed} {input.bam} > {output.bam}
        samtools index {output.bam}
        """
