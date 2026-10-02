# The contaminant set (rRNA, tRNA, snoRNA, ...): bowtie_index (align.smk)
# indexes it with the transcriptome, and bowtie_align removes the reads that fit
# a contaminant at least as well as a transcript.


# Not part of the default targets: run explicitly (needs internet access) to
# rebuild the contaminant set from public sources, then copy the output over
# resources/contaminants_built.fa (git-tracked; the default contaminants_fa).
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
