# Per-library Ribo-seq read counts: the contaminant step (bowtie_align's assign
# tsv) and the species split. A read aligned to both species is in neither BAM;
# it is counted as both_species (from the split logs, which must agree). Counts
# are collapsed reads, bowtie_align's input reads. contaminant_unique / _multi
# are the reads with a contaminant alignment; contaminant_rescued of them fit
# the transcriptome better (assign_reads.py's rescued) and are counted in
# clean_reads.


rule qc_library:
    input:
        assign=f"{EXP_DIR}/bowtie/{{sample}}.assign.tsv",
        human=f"{EXP_DIR}/split_bam/transcriptome/human/{{sample}}.bam",
        spike_in=f"{EXP_DIR}/split_bam/transcriptome/spike_in/{{sample}}.bam",
        human_log=f"{EXP_LOG_DIR}/align/split_bam/human/{{sample}}.log",
        spike_in_log=f"{EXP_LOG_DIR}/align/split_bam/spike_in/{{sample}}.log",
    output:
        f"{EXP_DIR}/qc/{{sample}}.stats.tsv",
    log:
        f"{EXP_LOG_DIR}/qc/{{sample}}.log",
    conda:
        "../envs/pysam.yaml"
    shell:
        r"""
        exec 2> {log}
        export LC_ALL=C
        assign_stat() {{
            awk -F'\t' -v col="$1" 'NR == 1 {{for (i = 1; i <= NF; i++) if ($i == col) c = i}} NR == 2 && c {{print $c}}' \
                {input.assign}
        }}
        input=$(assign_stat input_reads)
        unique=$(assign_stat contaminant_unique)
        multi=$(assign_stat contaminant_multi)
        rescued=$(assign_stat rescued)
        split_stat() {{ awk -F'\t' '$1 == "both_species" {{print $2}}' "$1"; }}
        both=$(split_stat {input.human_log})
        if [ -z "$both" ] || [ "$both" != "$(split_stat {input.spike_in_log})" ]; then
            echo "both_species missing from the split logs, or differs between them" >&2
            exit 1
        fi
        # reads = primary records: assign_reads.py leaves one per read
        primaries() {{ python -c 'import sys, pysam; print(pysam.view("-c", "-F", "0x100", sys.argv[1]).strip())' "$1"; }}
        human=$(primaries {input.human})
        spike_in=$(primaries {input.spike_in})

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
