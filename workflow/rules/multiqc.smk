# One MultiQC report per experiment: cutadapt, the three STAR steps
# (star_contaminant, contaminant_compete_align, star_transcript) and FastQC.
# Sections and sample names: workflow/multiqc_config.yaml. qc/summary.tsv
# stays: MultiQC does not count rescued reads or the species split.
MULTIQC_REPORTS = [
    f"{EXP_LOG_DIR}/cutadapt/{{sample}}.log",
    f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
    f"{EXP_DIR}/filter_reads/{{sample}}/compete/{{sample}}.transcriptome_Log.final.out",
    f"{EXP_DIR}/star/transcriptome/{{sample}}/{{sample}}.transcript_Log.final.out",
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
