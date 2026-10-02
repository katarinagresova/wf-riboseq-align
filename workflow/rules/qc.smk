# Per-library Ribo-seq read counts: the contaminant step (from STAR's
# Log.final.out) and the species split. A read aligned to both species is in
# neither BAM; it is counted as both_species (from the split logs, which must
# agree). Counts are collapsed reads, the unit STAR reports as input reads.
# contaminant_unique / _multi are what star_contaminant aligned;
# contaminant_rescued of them are put back (contaminant_compete) and counted in
# clean_reads.


rule qc_library:
    input:
        contam_log=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
        rescue=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.compete.tsv",
        human=f"{EXP_DIR}/split_bam/transcriptome/human/{{sample}}.bam",
        spike_in=f"{EXP_DIR}/split_bam/transcriptome/spike_in/{{sample}}.bam",
        human_log=f"{EXP_LOG_DIR}/align/split_bam/human/{{sample}}.log",
        spike_in_log=f"{EXP_LOG_DIR}/align/split_bam/spike_in/{{sample}}.log",
    output:
        f"{EXP_DIR}/qc/{{sample}}.stats.tsv",
    log:
        f"{EXP_LOG_DIR}/qc/{{sample}}.log",
    conda:
        "../envs/star.yaml"
    shell:
        r"""
        exec 2> {log}
        export LC_ALL=C
        star_stat() {{ grep "$1" {input.contam_log} | cut -f2; }}
        input=$(star_stat "Number of input reads")
        unique=$(star_stat "Uniquely mapped reads number")
        multi=$(star_stat "Number of reads mapped to multiple loci")
        rescued=$(awk -F'\t' '$1 == "rescued" {{print $2}}' {input.rescue})
        split_stat() {{ awk -F'\t' '$1 == "both_species" {{print $2}}' "$1"; }}
        both=$(split_stat {input.human_log})
        if [ -z "$both" ] || [ "$both" != "$(split_stat {input.spike_in_log})" ]; then
            echo "both_species missing from the split logs, or differs between them" >&2
            exit 1
        fi
        # reads = primary records: STAR and keep_sense leave one per read
        human=$(samtools view -c -F 0x100 {input.human})
        spike_in=$(samtools view -c -F 0x100 {input.spike_in})

        awk -v OFS='\t' -v sample={wildcards.sample} -v input="$input" -v unique="$unique" -v multi="$multi" \
            -v rescued="$rescued" -v h="$human" -v s="$spike_in" -v b="$both" 'BEGIN {{
                mapped = h + s + b
                print "sample", "input_reads", "contaminant_unique", "contaminant_multi", "contaminant_rescued",
                      "clean_reads", "human_only", "spike_in_only", "both_species", "spike_in_fraction"
                print sample, input, unique, multi, rescued, input - unique - multi + rescued,
                      h, s, b, (mapped ? s / mapped : 0)
            }}' > {output}
        """


rule qc_summary:
    input:
        tsvs=lambda wc: experiment_files(f"{EXP_DIR}/qc/{{sample}}.stats.tsv", wc.get("experiment"), "ribo"),
        script=workflow.source_path("../scripts/qc_summary.py"),
    output:
        f"{EXP_DIR}/qc/summary.tsv",
    params:
        min_spike_in_fraction=config["min_spike_in_fraction"],
    log:
        f"{EXP_LOG_DIR}/qc/qc_summary.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {output} {params.min_spike_in_fraction} {input.tsvs} 2> {log}"
