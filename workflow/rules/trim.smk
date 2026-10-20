# Ribo-seq read pre-processing: adapter trim -> UMI-aware collapse -> strip the
# randomers (config randomer_5p + randomer_3p nt, default 4 + 4).


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
        # randomers here (trim_reads strips them)
        minlength=config["minlength"] + config["randomer_5p"] + config["randomer_3p"],
        maxlength=config["maxlength"] + config["randomer_5p"] + config["randomer_3p"],
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
    # holds every distinct read in memory
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
    params:
        randomer_5p=config["randomer_5p"],
        randomer_3p=config["randomer_3p"],
    conda:
        "../envs/cutadapt.yaml"
    threads: 4
    shell:
        # Cut the randomers off both ends of each read and keep both in the
        # read name, e.g. (4 + 4) read "r" starting with ACCG and ending with
        # TATA becomes "r_ACCG:TATA".
        "cutadapt --cores {threads} -u {params.randomer_5p} -u -{params.randomer_3p} "
        "--rename '{{header}}_{{cut_prefix}}:{{cut_suffix}}' -o {output} {input} > {log} 2>&1"
