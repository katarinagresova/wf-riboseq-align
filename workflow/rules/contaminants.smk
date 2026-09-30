# Remove rRNA / tRNA / other contaminant reads: map to the contaminant set and
# keep what does NOT align.


# Not part of the default targets: run explicitly (needs internet access) to
# produce a contaminants_fa from public sources instead of a hand-maintained
# file, then point `contaminants_fa` in config.yaml at the output.
rule build_contaminants:
    output:
        f"{RESULTS_DIR}/reference/contaminants_built.fa",
    params:
        script=workflow.source_path("../scripts/build_contaminants.py"),
    log:
        f"{LOG_DIR}/build_contaminants.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.script} {output} 2> {log}"


rule number_contaminants:
    input:
        config["contaminants_fa"],
    output:
        f"{RESULTS_DIR}/reference/contaminants_numbered.fa",
    params:
        script=workflow.source_path("../scripts/number_contaminants.py"),
    log:
        f"{LOG_DIR}/number_contaminants.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.script} {input} {output} 2> {log}"


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
        f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.fastq.gz",
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

        gzip -c {params.prefix}Unmapped.out.mate1 > {output}
        rm -f {params.prefix}Unmapped.out.mate1 {params.prefix}Aligned.out.sam
        """
