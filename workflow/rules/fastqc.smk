# FastQC at two points: raw Ribo-seq, and the trimmed Ribo-seq that
# star_contaminant maps (after cutadapt_reads, collapse_reads and trim_reads,
# so collapsed: its duplication plot is ~flat). Reports only; no rule reads them.
#
# FastQC names its report after the input file, so each job links its fastq
# into its own tmp dir as <name>.fastq.gz and moves the report out. The tmp dir
# also takes FastQC's temporary files (--dir): compute nodes' /tmp can be full.
FASTQC_SHELL = (
    "rm -rf {params.tmp} && mkdir -p {params.tmp} && "
    "ln -s $(readlink -f {input.fastq}) {params.tmp}/{params.name}.fastq.gz && "
    "fastqc --quiet --threads 1 --dir {params.tmp} --outdir {params.tmp} "
    "{params.tmp}/{params.name}.fastq.gz > {log} 2>&1 && "
    "mv {params.tmp}/{params.name}_fastqc.html {output.html} && "
    "mv {params.tmp}/{params.name}_fastqc.zip {output.zip} && "
    "rm -rf {params.tmp}"
)


rule fastqc_ribo_raw:
    input:
        fastq=lambda wc: samples.loc[wc.sample, "fastq_1"],
    output:
        html=f"{EXP_DIR}/fastqc/ribo_raw/{{sample}}_fastqc.html",
        zip=f"{EXP_DIR}/fastqc/ribo_raw/{{sample}}_fastqc.zip",
    params:
        name="{sample}",
        tmp=f"{EXP_DIR}/fastqc/ribo_raw/{{sample}}_tmp",
    log:
        f"{EXP_LOG_DIR}/fastqc/ribo_raw/{{sample}}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        FASTQC_SHELL


rule fastqc_ribo_trimmed:
    input:
        fastq=f"{EXP_DIR}/trim_reads/{{sample}}.fastq.gz",
    output:
        html=f"{EXP_DIR}/fastqc/ribo_trimmed/{{sample}}_fastqc.html",
        zip=f"{EXP_DIR}/fastqc/ribo_trimmed/{{sample}}_fastqc.zip",
    params:
        name="{sample}",
        tmp=f"{EXP_DIR}/fastqc/ribo_trimmed/{{sample}}_tmp",
    log:
        f"{EXP_LOG_DIR}/fastqc/ribo_trimmed/{{sample}}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        FASTQC_SHELL
