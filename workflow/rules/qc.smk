# Per-library Ribo-seq read counts: the contaminant step (from STAR's
# Log.final.out) and the species split. A read aligned to both species stays in
# both BAMs; it is counted as both_species, not moved. Counts are collapsed
# reads, the unit STAR reports as input reads.


rule qc_library:
    input:
        contam_log=f"{RESULTS_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
        human=f"{RESULTS_DIR}/split_bam/transcriptome/human/{{sample}}.bam",
        spike_in=f"{RESULTS_DIR}/split_bam/transcriptome/spike_in/{{sample}}.bam",
    output:
        f"{RESULTS_DIR}/qc/{{sample}}.stats.tsv",
    conda:
        "../envs/star.yaml"
    shell:
        r"""
        export LC_ALL=C
        star_stat() {{ grep "$1" {input.contam_log} | cut -f2; }}
        input=$(star_stat "Number of input reads")
        unique=$(star_stat "Uniquely mapped reads number")
        multi=$(star_stat "Number of reads mapped to multiple loci")

        # comm columns: 1 = human only, 2 = spike-in only, 3 = both
        comm <(samtools view {input.human} | cut -f1 | sort -u) \
             <(samtools view {input.spike_in} | cut -f1 | sort -u) |
        awk -F'\t' -v OFS='\t' -v sample={wildcards.sample} \
            -v input="$input" -v unique="$unique" -v multi="$multi" '
            $1 != "" {{ h++; next }}
            $2 != "" {{ s++; next }}
            {{ b++ }}
            END {{
                mapped = h + s + b
                print "sample", "input_reads", "contaminant_unique", "contaminant_multi", "clean_reads",
                      "human_only", "spike_in_only", "both_species", "spike_in_fraction"
                print sample, input, unique, multi, input - unique - multi,
                      h + 0, s + 0, b + 0, (mapped ? s / mapped : 0)
            }}' > {output}
        """


rule qc_summary:
    input:
        expand(f"{RESULTS_DIR}/qc/{{sample}}.stats.tsv", sample=RIBO_SAMPLES),
    output:
        f"{RESULTS_DIR}/qc/summary.tsv",
    params:
        script=workflow.source_path("../scripts/qc_summary.py"),
        min_spike_in_fraction=config.get("min_spike_in_fraction", 0.01),
    log:
        f"{LOG_DIR}/qc/qc_summary.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {params.script} {output} {params.min_spike_in_fraction} {input} 2> {log}"
