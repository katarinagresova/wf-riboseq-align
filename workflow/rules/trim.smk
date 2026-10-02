# Ribo-seq read pre-processing: adapter trim -> UMI-aware collapse -> strip the
# 4+4 nt randomers.


rule cutadapt_reads:
    input:
        lambda wc: samples.loc[wc.sample, "fastq_1"],
    output:
        temp(f"{EXP_DIR}/cutadapt_reads/{{sample}}.fastq.gz"),
    log:
        f"{EXP_LOG_DIR}/cutadapt/{{sample}}.log",
    # via params, not {config[...]} in the shell string: as a module, the shell
    # would format against the importing workflow's config
    params:
        adapter=config["adapter"],
        # minlength / maxlength are of the insert: the read still has its
        # 4+4 nt randomers here (trim_reads strips them)
        minlength=config["minlength"] + 8,
        maxlength=config["maxlength"] + 8,
        qualcutoff=config["qualcutoff"],
    conda:
        "../envs/cutadapt.yaml"
    threads: 4
    shell:
        # -O 6: default -O 3 lets a 3 nt chance match to the adapter trim real
        # insert bases off the 3' end; -O 6 is long enough to be a real adapter
        "cutadapt --cores {threads} -a {params.adapter} -O 6 "
        "--minimum-length {params.minlength} --maximum-length {params.maxlength} "
        "--quality-cutoff {params.qualcutoff} -o {output} {input} > {log} 2>&1"


# Collapse BEFORE stripping the randomers: they are part of the sequence here,
# so identical reads with different randomers stay distinct molecules.
rule collapse_reads:
    input:
        fastq=f"{EXP_DIR}/cutadapt_reads/{{sample}}.fastq.gz",
        script=workflow.source_path("../scripts/collapse_reads.py"),
    output:
        temp(f"{EXP_DIR}/collapse_reads/{{sample}}.fastq.gz"),
    log:
        f"{EXP_LOG_DIR}/collapse_reads/{{sample}}.log",
    conda:
        "../envs/pysam.yaml"
    # holds every distinct read in memory: 4.2 GB peak on 41M reads
    resources:
        mem_mb=16000,
    shell:
        "python {input.script} {input.fastq} {wildcards.sample} 2> {log} | gzip > {output}"


rule trim_reads:
    input:
        f"{EXP_DIR}/collapse_reads/{{sample}}.fastq.gz",
    output:
        f"{EXP_DIR}/trim_reads/{{sample}}.fastq.gz",
    log:
        f"{EXP_LOG_DIR}/trim_reads/{{sample}}.log",
    conda:
        "../envs/cutadapt.yaml"
    threads: 4
    shell:
        # Cut the 4 nt randomer off both ends of each read and keep both in
        # the read name, e.g. read "r" starting with ACCG and ending with TATA
        # becomes "r_ACCG:TATA".
        "cutadapt --cores {threads} -u 4 -u -4 "
        "--rename '{{header}}_{{cut_prefix}}:{{cut_suffix}}' -o {output} {input} > {log} 2>&1"
