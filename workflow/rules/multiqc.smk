# One MultiQC report per experiment: cutadapt, bowtie_align's read classes (its
# assign tsv, as custom content) and FastQC. Sections and sample names:
# workflow/multiqc_config.yaml. qc/summary.tsv stays: MultiQC does not show the
# species split. bowtie's own log is not shown: its "aligned" mixes contaminant
# and transcript hits.
MULTIQC_REPORTS = [
    f"{EXP_LOG_DIR}/cutadapt/{{sample}}.log",
    f"{EXP_DIR}/bowtie/{{sample}}.assign.tsv",
    f"{EXP_DIR}/fastqc/ribo_raw/{{sample}}_fastqc.zip",
    f"{EXP_DIR}/fastqc/ribo_trimmed/{{sample}}_fastqc.zip",
]


rule multiqc:
    input:
        reports=lambda wc: [
            f for pattern in MULTIQC_REPORTS for f in experiment_files(pattern, wc.get("experiment"), "ribo")
        ],
        config=workflow.source_path("../multiqc_config.yaml"),
    output:
        html=f"{EXP_DIR}/multiqc/multiqc_report.html",
        data=directory(f"{EXP_DIR}/multiqc/multiqc_data"),
    log:
        f"{EXP_LOG_DIR}/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --config {input.config} --outdir $(dirname {output.html}) "
        "{input.reports} > {log} 2>&1"
