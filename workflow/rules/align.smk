# Ribo-seq: map contaminant-free reads to human + spike-in transcriptome, then
# split the BAM by species.
#
# The reference and its index do not depend on the experiment: they live in
# RESULTS_DIR, not EXP_DIR, so a multi-experiment run builds them once.


rule combined_transcriptome:
    input:
        human=config["human_transcriptome_fa"],
        spike_in=config["spike_in_transcriptome_fa"],
    output:
        f"{RESULTS_DIR}/reference/transcriptome.combined_human_spike_in.fa",
    log:
        f"{LOG_DIR}/combined_transcriptome.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.human} {input.spike_in} > {output} 2> {log}"


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
    resources:
        mem_mb=16000,
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output.index} --genomeFastaFiles {input} "
        "--genomeSAindexNbases 11 --genomeChrBinNbits 12 "
        "--outFileNamePrefix {params.log_prefix} > {log} 2>&1"


# Multimappers are kept (up to 255 loci); STAR marks one alignment per read
# primary, chosen at random among equally good loci. CLEAN_FASTQ: common.smk.
rule star_transcript:
    input:
        fastq=CLEAN_FASTQ,
        index=f"{RESULTS_DIR}/star_index/transcriptome",
    output:
        bam=temp(f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.bam"),
        log_final=f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.transcript_Log.final.out",
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

        samtools sort -@ {threads} {params.prefix}Aligned.out.bam -o {output.bam}
        rm {params.prefix}Aligned.out.bam
        """


# The libraries are stranded: drop the antisense alignments and fix NH, MAPQ and
# the primary flag of the reads that had any (see the script).
rule keep_sense:
    input:
        bam=f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.bam",
        script=workflow.source_path("../scripts/filter_alignments.py"),
    output:
        temp(f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.sense.bam"),
    log:
        f"{EXP_LOG_DIR}/align/keep_sense/{{sample}}.log",
    conda:
        "../envs/pysam.yaml"
    shell:
        "python {input.script} sense {input.bam} {output} 2> {log}"


rule split_bed_transcriptome:
    input:
        chrom_sizes=f"{RESULTS_DIR}/star_index/transcriptome/chrNameLength.txt",
        human_fa=config["human_transcriptome_fa"],
        spike_in_fa=config["spike_in_transcriptome_fa"],
        script=workflow.source_path("../scripts/split_bed_transcriptome.py"),
    output:
        human=f"{RESULTS_DIR}/reference/transcripts.human.bed",
        spike_in=f"{RESULTS_DIR}/reference/transcripts.spike_in.bed",
    log:
        f"{LOG_DIR}/split_bed_transcriptome.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.chrom_sizes} {input.human_fa} {input.spike_in_fa} "
        "{output.human} {output.spike_in} 2> {log}"


# One BAM per species. A read aligned to both is in neither: it may come from
# either (see the script). The log counts it as both_species.
rule split_bam_transcriptome:
    input:
        bam=f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.sense.bam",
        bed=f"{RESULTS_DIR}/reference/transcripts.{{species}}.bed",
        script=workflow.source_path("../scripts/filter_alignments.py"),
    output:
        bam=f"{EXP_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam",
        bai=f"{EXP_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam.bai",
    wildcard_constraints:
        species="human|spike_in",
    log:
        f"{EXP_LOG_DIR}/align/split_bam/{{species}}/{{sample}}.log",
    conda:
        "../envs/pysam.yaml"
    shell:
        r"""
        exec 2> {log}
        python {input.script} refs {input.bed} {input.bam} {output.bam}
        python -c 'import sys, pysam; pysam.index(sys.argv[1])' {output.bam}
        """
