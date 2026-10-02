# Matched total RNA-seq (paired-end): quantify on the UNFILTERED human
# transcriptome with salmon, the genome as decoys. Runs in both modes; in mode
# "filtered" these quants also decide which transcripts the ribo reference
# keeps (autofilter.smk).
#
# The index does not depend on the experiment: it lives in RESULTS_DIR, not
# EXP_DIR, so a multi-experiment run builds it once.


# salmon's decoy-aware reference: the transcripts, then every genome sequence
# (decoys must come last), and the decoys' names.
rule salmon_gentrome:
    input:
        transcriptome=config["human_transcriptome_fa"],
        genome=config["human_genome_fa"],
    output:
        gentrome=temp(f"{RESULTS_DIR}/salmon_index/gentrome.fa"),
        decoys=f"{RESULTS_DIR}/salmon_index/decoys.txt",
    log:
        f"{LOG_DIR}/salmon_index/gentrome.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        r"""
        exec 2> {log}
        grep '^>' {input.genome} | cut -c2- | cut -d' ' -f1 > {output.decoys}
        cat {input.transcriptome} {input.genome} > {output.gentrome}
        """


# A fragment that aligns better to the genome than to any transcript (intronic,
# intergenic, an unannotated copy) is counted for no transcript. By default
# salmon drops a transcript whose sequence duplicates another's (25 in the
# HCT116 set) from the index and quant.sf, and filter_rnaseq would blacklist it
# as TPM 0: --keepDuplicates.
rule salmon_index:
    input:
        gentrome=f"{RESULTS_DIR}/salmon_index/gentrome.fa",
        decoys=f"{RESULTS_DIR}/salmon_index/decoys.txt",
    output:
        directory(f"{RESULTS_DIR}/salmon_index/transcriptome_genome_decoys"),
    log:
        f"{LOG_DIR}/salmon_index/index.log",
    conda:
        "../envs/salmon.yaml"
    threads: 16
    resources:
        mem_mb=24000,
    shell:
        "salmon index -p {threads} -t {input.gentrome} -d {input.decoys} -k 31 --keepDuplicates "
        "-i {output} > {log} 2>&1"


rule salmon_quant:
    input:
        fastq1=lambda wc: samples.loc[wc.sample, "fastq_1"],
        fastq2=lambda wc: samples.loc[wc.sample, "fastq_2"],
        index=f"{RESULTS_DIR}/salmon_index/transcriptome_genome_decoys",
    output:
        f"{EXP_DIR}/salmon/{{sample}}/quant.sf",
    params:
        out_dir=f"{EXP_DIR}/salmon/{{sample}}",
    log:
        f"{EXP_LOG_DIR}/salmon/{{sample}}.log",
    conda:
        "../envs/salmon.yaml"
    threads: 8
    resources:
        mem_mb=24000,
    shell:
        "salmon quant -p {threads} -i {input.index} -l A -1 {input.fastq1} -2 {input.fastq2} "
        "--seqBias --gcBias --output {params.out_dir} > {log} 2>&1"
