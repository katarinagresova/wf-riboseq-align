# Remove rRNA / tRNA / other contaminant reads: map to the contaminant set and
# keep what does NOT align.


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


rule star_contaminant:
    input:
        fastq=f"{RESULTS_DIR}/trim_reads/{{sample}}.fastq.gz",
        index=f"{RESULTS_DIR}/star_index/contaminants",
    output:
        fastq=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.fastq.gz",
        log_final=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
    params:
        prefix=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.contam_",
        tmp_dir=f"{RESULTS_DIR}/filter_reads/{{sample}}/_tmpSTAR",
    conda:
        "../envs/star.yaml"
    threads: 8
    shell:
        r"""
        rm -rf {params.tmp_dir}
        STAR \
            --genomeDir {input.index} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outSAMunmapped Within \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax 10000 \
            --outSAMmultNmax 1 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outTmpDir {params.tmp_dir} \
            --genomeLoad NoSharedMemory \
            --outSAMattributes NH HI AS NM MD \
            --outSAMtype SAM \
            --outFileNamePrefix {params.prefix} \
            --outReadsUnmapped Fastx \
            --readFilesIn {input.fastq}

        too_many=$(grep "Number of reads mapped to too many loci" {params.prefix}Log.final.out | cut -f2)
        if [ "$too_many" != "0" ]; then
            echo "ERROR: $too_many reads exceeded --outFilterMultimapNmax 10000 in contaminant mapping; raise the limit." >&2
            exit 1
        fi

        gzip -c {params.prefix}Unmapped.out.mate1 > {output.fastq}
        rm -f {params.prefix}Unmapped.out.mate1 {params.prefix}Aligned.out.sam
        """


# contaminant_filter: competitive. star_contaminant discards every read with a
# contaminant hit; these two rules put back the ones whose best sense alignment
# to the transcriptome scores higher than their contaminant alignment
# (scripts/contaminant_compete.py). star_contaminant is left as it is (a change
# to its command would rerun it, and everything after it, in remove mode), so
# the reads are aligned to the contaminants once more, keeping the alignments.
# Both STAR calls must use the same parameters as star_contaminant and
# star_transcript, plus the AS attribute.
rule contaminant_compete_align:
    input:
        fastq=f"{RESULTS_DIR}/trim_reads/{{sample}}.fastq.gz",
        contaminant_index=f"{RESULTS_DIR}/star_index/contaminants",
        transcriptome_index=f"{RESULTS_DIR}/star_index/transcriptome",
    output:
        removed=temp(f"{RESULTS_DIR}/filter_reads/{{sample}}/compete/removed.fastq.gz"),
        contaminant=temp(f"{RESULTS_DIR}/filter_reads/{{sample}}/compete/contaminant.tsv.gz"),
        transcriptome=temp(f"{RESULTS_DIR}/filter_reads/{{sample}}/compete/transcriptome.tsv.gz"),
    params:
        prefix=f"{RESULTS_DIR}/filter_reads/{{sample}}/compete/",
    log:
        f"{LOG_DIR}/contaminant_compete/{{sample}}.align.log",
    conda:
        "../envs/star.yaml"
    threads: 8
    shell:
        r"""
        exec > {log} 2>&1
        rm -rf {params.prefix}contaminant_tmpSTAR {params.prefix}transcriptome_tmpSTAR
        # star_contaminant's alignment; only the aligned reads, one alignment each
        STAR \
            --genomeDir {input.contaminant_index} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax 10000 \
            --outSAMmultNmax 1 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outTmpDir {params.prefix}contaminant_tmpSTAR \
            --genomeLoad NoSharedMemory \
            --outSAMattributes NH HI AS NM MD \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.prefix}contaminant_ \
            --readFilesIn {input.fastq}
        samtools fastq {params.prefix}contaminant_Aligned.out.bam | gzip > {output.removed}
        # read, contaminant record, AS
        samtools view {params.prefix}contaminant_Aligned.out.bam |
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
            --seedSearchLmax 10 \
            --outFilterMultimapNmax 255 \
            --outFilterMismatchNmax 2 \
            --outFilterMultimapScoreRange 0 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFileNamePrefix {params.prefix}transcriptome_ \
            --readFilesIn {output.removed} \
            --readFilesCommand zcat
        # read, flag, AS
        samtools view {params.prefix}transcriptome_Aligned.out.bam |
            awk -F'\t' -v OFS='\t' '{{for (i = 12; i <= NF; i++) if ($i ~ /^AS:i:/) print $1, $2, substr($i, 6)}}' |
            gzip > {output.transcriptome}
        rm -f {params.prefix}contaminant_Aligned.out.bam {params.prefix}transcriptome_Aligned.out.bam
        """


rule contaminant_compete:
    input:
        clean=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.fastq.gz",
        contam_log=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
        removed=f"{RESULTS_DIR}/filter_reads/{{sample}}/compete/removed.fastq.gz",
        contaminant=f"{RESULTS_DIR}/filter_reads/{{sample}}/compete/contaminant.tsv.gz",
        transcriptome=f"{RESULTS_DIR}/filter_reads/{{sample}}/compete/transcriptome.tsv.gz",
        script=workflow.source_path("../scripts/contaminant_compete.py"),
    output:
        fastq=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.competitive.fastq.gz",
        stats=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.compete.tsv",
        records=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.compete_records.tsv",
    log:
        f"{LOG_DIR}/contaminant_compete/{{sample}}.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.clean} {input.contam_log} {input.removed} {input.contaminant} "
        "{input.transcriptome} {output.fastq} {output.stats} {output.records} 2> {log}"
