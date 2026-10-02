# Remove rRNA / tRNA / other contaminant reads: map to the contaminant set and
# keep what does NOT align, plus the aligned reads that fit the transcriptome
# better (the competitive filter).
#
# The contaminant set and its index do not depend on the experiment: they live
# in RESULTS_DIR, not EXP_DIR, so a multi-experiment run builds them once.


# Not part of the default targets: run explicitly (needs internet access) to
# produce a contaminants_fa from public sources instead of a hand-maintained
# file, then point `contaminants_fa` in config.yaml at the output.
rule build_contaminants:
    input:
        script=workflow.source_path("../scripts/build_contaminants.py"),
    output:
        f"{RESULTS_DIR}/reference/contaminants_built.fa",
    log:
        f"{LOG_DIR}/build_contaminants.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {output} 2> {log}"


rule number_contaminants:
    input:
        fa=config["contaminants_fa"],
        script=workflow.source_path("../scripts/number_contaminants.py"),
    output:
        f"{RESULTS_DIR}/reference/contaminants_numbered.fa",
    log:
        f"{LOG_DIR}/number_contaminants.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.fa} {output} 2> {log}"


rule star_contaminant_index:
    input:
        f"{RESULTS_DIR}/reference/contaminants_numbered.fa",
    output:
        directory(f"{RESULTS_DIR}/star_index/contaminants"),
    log:
        f"{LOG_DIR}/star/contaminant_index.log",
    params:
        log_prefix=f"{LOG_DIR}/star/contaminant_index.",
    conda:
        "../envs/star.yaml"
    threads: 8
    shell:
        "STAR --genomeSAindexNbases 9 --runThreadN {threads} --runMode genomeGenerate "
        "--genomeFastaFiles {input} --genomeDir {output} "
        "--outFileNamePrefix {params.log_prefix} > {log} 2>&1"


# Aligns the trimmed reads to the contaminants once: the unaligned reads are the
# clean fastq, the aligned ones (one alignment each, with AS) are kept for the
# competitive filter below.
rule star_contaminant:
    input:
        fastq=f"{EXP_DIR}/trim_reads/{{sample}}.fastq.gz",
        index=f"{RESULTS_DIR}/star_index/contaminants",
    output:
        fastq=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.fastq.gz",
        log_final=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
        bam=temp(f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Aligned.out.bam"),
    log:
        f"{EXP_LOG_DIR}/star/contaminant/{{sample}}.log",
    params:
        prefix=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_",
        tmp_dir=f"{EXP_DIR}/filter_reads/{{sample}}/_tmpSTAR",
    conda:
        "../envs/star.yaml"
    threads: 8
    resources:
        mem_mb=16000,
    shell:
        r"""
        exec > {log} 2>&1
        rm -rf {params.tmp_dir}
        STAR \
            --genomeDir {input.index} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax 10000 \
            --outSAMmultNmax 1 \
            --alignIntronMax 1 \
            --alignEndsType Extend5pOfRead1 \
            --outTmpDir {params.tmp_dir} \
            --genomeLoad NoSharedMemory \
            --outSAMattributes NH HI AS NM MD \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.prefix} \
            --outReadsUnmapped Fastx \
            --readFilesIn {input.fastq}

        too_many=$(grep "Number of reads mapped to too many loci" {params.prefix}Log.final.out | cut -f2)
        if [ "$too_many" != "0" ]; then
            echo "ERROR: $too_many reads exceeded --outFilterMultimapNmax 10000 in contaminant mapping; raise the limit." >&2
            exit 1
        fi

        gzip -c {params.prefix}Unmapped.out.mate1 > {output.fastq}
        rm -f {params.prefix}Unmapped.out.mate1
        """


# The competitive filter: of the reads star_contaminant aligned, put back the
# ones whose best sense alignment to the transcriptome scores higher than their
# contaminant alignment (scripts/contaminant_compete.py). The transcriptome STAR
# call uses star_transcript's settings (RIBO_TRANSCRIPTOME_STAR_ARGS,
# common.smk), plus the AS attribute.
rule contaminant_compete_align:
    input:
        bam=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Aligned.out.bam",
        transcriptome_index=f"{RESULTS_DIR}/star_index/transcriptome",
    output:
        removed=temp(f"{EXP_DIR}/filter_reads/{{sample}}/compete/removed.fastq.gz"),
        contaminant=temp(f"{EXP_DIR}/filter_reads/{{sample}}/compete/contaminant.tsv.gz"),
        transcriptome=temp(f"{EXP_DIR}/filter_reads/{{sample}}/compete/transcriptome.tsv.gz"),
        log_final=f"{EXP_DIR}/filter_reads/{{sample}}/compete/{{sample}}.transcriptome_Log.final.out",
    # the sample in STAR's file names: MultiQC names a STAR report after its file
    params:
        prefix=f"{EXP_DIR}/filter_reads/{{sample}}/compete/{{sample}}.",
        star_args=RIBO_TRANSCRIPTOME_STAR_ARGS,
    log:
        f"{EXP_LOG_DIR}/contaminant_compete/{{sample}}.align.log",
    conda:
        "../envs/star.yaml"
    threads: 8
    resources:
        mem_mb=24000,
    shell:
        r"""
        exec > {log} 2>&1
        rm -rf {params.prefix}transcriptome_tmpSTAR
        samtools fastq {input.bam} | gzip > {output.removed}
        # read, contaminant record, AS
        samtools view {input.bam} |
            awk -F'\t' -v OFS='\t' '{{for (i = 12; i <= NF; i++) if ($i ~ /^AS:i:/) print $1, $3, substr($i, 6)}}' |
            gzip > {output.contaminant}

        # star_transcript's alignment of those reads
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.transcriptome_index} \
            --outTmpDir {params.prefix}transcriptome_tmpSTAR \
            --outSAMtype BAM Unsorted \
            --outSAMmode NoQS \
            --outSAMattributes NH AS NM \
            {params.star_args} \
            --outFileNamePrefix {params.prefix}transcriptome_ \
            --readFilesIn {output.removed} \
            --readFilesCommand zcat
        # read, flag, AS
        samtools view {params.prefix}transcriptome_Aligned.out.bam |
            awk -F'\t' -v OFS='\t' '{{for (i = 12; i <= NF; i++) if ($i ~ /^AS:i:/) print $1, $2, substr($i, 6)}}' |
            gzip > {output.transcriptome}
        rm -f {params.prefix}transcriptome_Aligned.out.bam
        """


rule contaminant_compete:
    input:
        clean=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.fastq.gz",
        contam_log=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
        removed=f"{EXP_DIR}/filter_reads/{{sample}}/compete/removed.fastq.gz",
        contaminant=f"{EXP_DIR}/filter_reads/{{sample}}/compete/contaminant.tsv.gz",
        transcriptome=f"{EXP_DIR}/filter_reads/{{sample}}/compete/transcriptome.tsv.gz",
        script=workflow.source_path("../scripts/contaminant_compete.py"),
    output:
        fastq=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.competitive.fastq.gz",
        stats=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.compete.tsv",
        records=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.compete_records.tsv",
    log:
        f"{EXP_LOG_DIR}/contaminant_compete/{{sample}}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=16000,
    shell:
        "python {input.script} {input.clean} {input.contam_log} {input.removed} {input.contaminant} "
        "{input.transcriptome} {output.fastq} {output.stats} {output.records} 2> {log}"
