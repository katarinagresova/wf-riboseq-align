# Matched total RNA-seq (paired-end): map to the UNFILTERED human transcriptome
# and quantify with salmon. Runs in both modes; in mode "filtered" these quants
# also decide which transcripts the ribo reference keeps (autofilter.smk).


rule star_transcript_index_rnaseq:
    input:
        config["human_transcriptome_fa"],
    output:
        index=directory(f"{RESULTS_DIR}/star_index/unfiltered_transcriptome"),
        chrom_sizes=f"{RESULTS_DIR}/star_index/unfiltered_transcriptome/chrNameLength.txt",
    log:
        f"{LOG_DIR}/star/unfiltered_transcriptome_index.log",
    params:
        log_prefix=f"{LOG_DIR}/star/unfiltered_transcriptome_index.",
    conda:
        "../envs/star.yaml"
    threads: 8
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate "
        "--genomeDir {output.index} --genomeFastaFiles {input} "
        "--genomeSAindexNbases 11 --genomeChrBinNbits 12 "
        "--outFileNamePrefix {params.log_prefix} > {log} 2>&1"


# Every alignment of a multimapping fragment is written (no --outSAMmultNmax),
# so salmon distributes it over its transcripts by EM. With
# --outSAMmultNmax 1, as in the eIF pipeline, salmon saw one arbitrary
# alignment per fragment and assigned it there.
rule star_transcript_rnaseq:
    input:
        fastq1=lambda wc: samples.loc[wc.sample, "fastq_1"],
        fastq2=lambda wc: samples.loc[wc.sample, "fastq_2"],
        index=f"{RESULTS_DIR}/star_index/unfiltered_transcriptome",
    output:
        bam=f"{RESULTS_DIR}/star/transcriptome_rnaseq/{{sample}}/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/star/transcriptome_rnaseq/{{sample}}/{{sample}}.bam.bai",
        unsorted_bam=temp(f"{RESULTS_DIR}/star/transcriptome_rnaseq/{{sample}}/{{sample}}.transcript_Aligned.out.bam"),
    params:
        prefix=f"{RESULTS_DIR}/star/transcriptome_rnaseq/{{sample}}/{{sample}}.transcript_",
    conda:
        "../envs/star.yaml"
    threads: 8
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
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFileNamePrefix {params.prefix} \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --readFilesCommand zcat

        samtools sort -@ {threads} {output.unsorted_bam} -o {output.bam}
        samtools index {output.bam}
        """


rule salmon_bam:
    input:
        bam=f"{RESULTS_DIR}/star/transcriptome_rnaseq/{{sample}}/{{sample}}.transcript_Aligned.out.bam",
        fasta=config["human_transcriptome_fa"],
    output:
        f"{RESULTS_DIR}/salmon/{{sample}}/quant.sf",
    params:
        out_dir=f"{RESULTS_DIR}/salmon/{{sample}}",
    log:
        f"{LOG_DIR}/salmon/{{sample}}.log",
    conda:
        "../envs/salmon.yaml"
    threads: 4
    shell:
        r"""
        salmon quant -p {threads} --seqBias -t {input.fasta} -l A -a {input.bam} \
            --output {params.out_dir} > {log} 2>&1

        # salmon only warns when mates are not adjacent in the BAM (e.g. a
        # coordinate-sorted one), and then quantifies them wrongly
        if grep -q "suspicious pair" {log}; then
            echo "ERROR: salmon reported suspicious pairs (see {log}); mates must be adjacent in {input.bam}." >&2
            exit 1
        fi
        """
