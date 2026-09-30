# mode "filtered" only: drop transcripts the matched RNA-seq does not express,
# giving the human half of the ribo reference. In mode "unfiltered" nothing
# requests these outputs, so the rules never run.


rule filter_rnaseq:
    input:
        fa=config["human_transcriptome_fa"],
        quant=expand(f"{RESULTS_DIR}/salmon/{{sample}}/quant.sf", sample=RNA_SAMPLES),
        script=workflow.source_path("../scripts/autofilter_rnaseq.py"),
    output:
        f"{RESULTS_DIR}/reference/rnaseq_filter_blacklist_txid.txt",
    params:
        min_tpm=config["autofilter"]["min_tpm"],
        min_samples=config["autofilter"]["min_samples"],
    log:
        f"{LOG_DIR}/autofilter/filter_rnaseq.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.fa} {output} "
        "{params.min_tpm} {params.min_samples} {input.quant} 2> {log}"


rule make_filtered_resources:
    input:
        fa=config["human_transcriptome_fa"],
        gtf=config["human_transcriptome_gtf"],
        blacklist=f"{RESULTS_DIR}/reference/rnaseq_filter_blacklist_txid.txt",
        script=workflow.source_path("../scripts/filter_transcriptome.py"),
    output:
        fa=f"{RESULTS_DIR}/reference/human_transcriptome.rnaseq_filtered.fa",
        gtf=f"{RESULTS_DIR}/reference/human_transcriptome.rnaseq_filtered.gtf",
    log:
        f"{LOG_DIR}/autofilter/make_filtered_resources.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.fa} {input.gtf} {input.blacklist} "
        "{output.fa} {output.gtf} 2> {log}"
