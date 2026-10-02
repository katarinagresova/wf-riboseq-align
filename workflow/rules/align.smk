# Ribo-seq: align the trimmed reads to the human + spike-in transcriptome and
# the contaminants in one pass, keep the reads that fit the transcriptome
# better, then split the BAM by species.
#
# The index and the BEDs do not depend on the experiment: they live in
# RESULTS_DIR, not EXP_DIR, so a multi-experiment run builds them once.


# One index of the transcriptome (human + spike-in) and the contaminants:
# bowtie_align aligns each read to all of them at once. The record names must be
# unique across the three fastas (bowtie keeps a header up to the first space),
# and contaminants_fa must hold each sequence once, as build_contaminants writes
# it: a second copy would only split reads between the two (one strand: a
# reverse-complement copy is allowed).
rule bowtie_index:
    input:
        human=config["human_transcriptome_fa"],
        spike_in=config["spike_in_transcriptome_fa"],
        contaminants=config["contaminants_fa"],
    output:
        multiext(f"{RESULTS_DIR}/bowtie_index/ribo", ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt",
                 ".rev.2.ebwt"),
    log:
        f"{LOG_DIR}/bowtie_index.log",
    params:
        prefix=lambda wildcards, output: output[0].removesuffix(".1.ebwt"),
    conda:
        "../envs/bowtie.yaml"
    threads: 8
    resources:
        mem_mb=4000,
    shell:
        r"""
        exec > {log} 2>&1
        awk -v contaminants={input.contaminants} '
            function end_record() {{
                if (file != contaminants) return
                if (seq in first) {{ print "ERROR: " name " has the same sequence as " first[seq]; bad = 1 }}
                else first[seq] = name
            }}
            /^>/ {{
                if (name != "") end_record()
                name = substr($1, 2); seq = ""; file = FILENAME
                if (name in where) {{ print "ERROR: record name " name " in " where[name] " and in " FILENAME; bad = 1 }}
                where[name] = FILENAME
                next
            }}
            FILENAME == contaminants {{ seq = seq toupper($0) }}
            END {{
                if (name != "") end_record()
                if (bad) print "ERROR: record names must be unique across the fastas, and contaminants_fa needs each sequence once, as build_contaminants writes it"
                exit bad
            }}
        ' {input.human} {input.spike_in} {input.contaminants}
        bowtie-build --threads {threads} {input.human},{input.spike_in},{input.contaminants} {params.prefix}
        """


# Aligns the trimmed reads once to the transcriptome and the contaminants
# (BOWTIE_ARGS, common.smk); assign_reads.py keeps the reads that fit the
# transcriptome better than any contaminant, with their best alignments (see
# the script). Fails if bowtie skipped a read or the read counts differ.
rule bowtie_align:
    input:
        fastq=f"{EXP_DIR}/trim_reads/{{sample}}.fastq.gz",
        index=rules.bowtie_index.output,
        contaminants=config["contaminants_fa"],
        script=workflow.source_path("../scripts/assign_reads.py"),
    output:
        bam=temp(f"{EXP_DIR}/bowtie/{{sample}}.bam"),
        stats=f"{EXP_DIR}/bowtie/{{sample}}.assign.tsv",
        records=f"{EXP_DIR}/bowtie/{{sample}}.contaminant_records.tsv",
    log:
        bowtie=f"{EXP_LOG_DIR}/bowtie/{{sample}}.log",
        assign=f"{EXP_LOG_DIR}/bowtie/{{sample}}.assign.log",
    params:
        index=lambda wildcards, input: input.index[0].removesuffix(".1.ebwt"),
        bowtie_args=BOWTIE_ARGS,
        short_length=SHORT_READ_LENGTH,
        short_mismatches=SHORT_READ_MISMATCHES,
    conda:
        "../envs/bowtie.yaml"
    threads: 8
    resources:
        mem_mb=8000,
    shell:
        r"""
        bowtie -p {threads} {params.bowtie_args} --sam -x {params.index} {input.fastq} 2> {log.bowtie} |
            python {input.script} {wildcards.sample} {input.contaminants} {params.short_length} \
                {params.short_mismatches} {output.bam}.unsorted.bam {output.stats} {output.records} 2> {log.assign}
        if grep -q "Exhausted best-first chunk memory" {log.bowtie}; then
            echo "ERROR: bowtie skipped reads (Exhausted best-first chunk memory): raise --chunkmbs" >> {log.assign}
            exit 1
        fi
        processed=$(awk -F': ' '$1 == "# reads processed" {{print $2}}' {log.bowtie})
        assigned=$(awk -F'\t' 'NR == 2 {{print $2}}' {output.stats})
        if [ "$processed" != "$assigned" ]; then
            echo "ERROR: bowtie processed $processed reads, assign_reads.py assigned $assigned" >> {log.assign}
            exit 1
        fi
        python -c 'import sys, pysam; pysam.sort("-@", sys.argv[1], "-m", "512M", "-o", sys.argv[2], sys.argv[3])' \
            {threads} {output.bam} {output.bam}.unsorted.bam
        rm {output.bam}.unsorted.bam
        """


rule split_bed_transcriptome:
    input:
        human_fa=config["human_transcriptome_fa"],
        spike_in_fa=config["spike_in_transcriptome_fa"],
        script=workflow.source_path("../scripts/split_bed_transcriptome.py"),
    output:
        human=f"{RESULTS_DIR}/reference/transcripts.human.bed",
        spike_in=f"{RESULTS_DIR}/reference/transcripts.spike_in.bed",
    log:
        f"{LOG_DIR}/split_bed_transcriptome.log",
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.human_fa} {input.spike_in_fa} {output.human} {output.spike_in} 2> {log}"


# One BAM per species. A read aligned to both is in neither: it may come from
# either (see the script). The log counts it as both_species.
rule split_bam_transcriptome:
    input:
        bam=f"{EXP_DIR}/bowtie/{{sample}}.bam",
        bed=f"{RESULTS_DIR}/reference/transcripts.{{species}}.bed",
        script=workflow.source_path("../scripts/filter_alignments.py"),
    output:
        bam=f"{EXP_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam",
        bai=f"{EXP_DIR}/split_bam/transcriptome/{{species}}/{{sample}}.bam.bai",
    wildcard_constraints:
        species="human|spike_in",
    log:
        f"{EXP_LOG_DIR}/align/split_bam/{{species}}/{{sample}}.log",
    conda:
        "../envs/pysam.yaml"
    shell:
        r"""
        exec 2> {log}
        python {input.script} {input.bed} {input.bam} {output.bam}
        python -c 'import sys, pysam; pysam.index(sys.argv[1])' {output.bam}
        """
