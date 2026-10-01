# Ribo-seq read pre-processing: adapter trim -> UMI-aware collapse -> strip the
# 4+4 nt randomers. RNA-seq reads go to STAR untrimmed (rnaseq.smk).


rule cutadapt_reads:
    input:
        lambda wc: samples.loc[wc.sample, "fastq_1"],
    output:
        f"{EXP_DIR}/cutadapt_reads/{{sample}}.fastq.gz",
    log:
        f"{EXP_LOG_DIR}/cutadapt/{{sample}}.log",
    # via params, not {config[...]} in the shell string: as a module, the shell
    # would format against the importing workflow's config
    params:
        adapter=config["adapter"],
        minlength=config["minlength"],
        maxlength=config["maxlength"],
        qualcutoff=config["qualcutoff"],
    conda:
        "../envs/cutadapt.yaml"
    threads: 4
    shell:
        "zcat {input} | "
        "cutadapt --cores {threads} -a {params.adapter} "
        "--minimum-length {params.minlength} --maximum-length {params.maxlength} "
        "--quality-cutoff {params.qualcutoff} - 2> {log} | "
        "gzip > {output}"


# Collapse BEFORE stripping the randomers: they are part of the sequence here,
# so identical reads with different randomers stay distinct molecules.
rule collapse_reads:
    input:
        fastq=f"{EXP_DIR}/cutadapt_reads/{{sample}}.fastq.gz",
        script=workflow.source_path("../scripts/collapse_reads.py"),
    output:
        f"{EXP_DIR}/collapse_reads/{{sample}}.fastq.gz",
    log:
        f"{EXP_LOG_DIR}/collapse_reads/{{sample}}.log",
    conda:
        "../envs/python.yaml"
    shell:
        "zcat {input.fastq} | python {input.script} {wildcards.sample} 2> {log} | gzip > {output}"


rule trim_reads:
    input:
        fastq=f"{EXP_DIR}/collapse_reads/{{sample}}.fastq.gz",
        script=workflow.source_path("../scripts/remove_randomers.py"),
    output:
        f"{EXP_DIR}/trim_reads/{{sample}}.fastq.gz",
    log:
        f"{EXP_LOG_DIR}/trim_reads/{{sample}}.log",
    conda:
        "../envs/python.yaml"
    shell:
        "zcat {input.fastq} | python {input.script} 2> {log} | gzip > {output}"
