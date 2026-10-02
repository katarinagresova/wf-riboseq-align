# Remove rRNA / tRNA / other contaminant reads: map to the contaminant set and
# keep what does NOT align, plus the aligned reads that fit the transcriptome
# better (the competitive filter).
#
# The contaminant set and its index do not depend on the experiment: they live
# in RESULTS_DIR, not EXP_DIR, so a multi-experiment run builds them once.


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


# contaminants_fa goes to STAR as is: build_contaminants writes it with unique
# names and one record per sequence, and a fasta of the user's own must be so
# too, else an error. STAR needs unique names (it keeps the header up to the
# first space); a sequence twice would only split reads between the copies
# (one strand: a reverse-complement copy is allowed).
#
# --genomeChrBinNbits 8: STAR pads every record to a multiple of 2^bits. The
# default 18 (262 kb) made 1 Mb of short records a 1.9 GB index; STAR's
# recommended min(18, log2(max(length / records, read length))) is ~7-10 for
# contaminant sets, and 256 nt bins exceed any Ribo-seq read.
rule star_contaminant_index:
    input:
        config["contaminants_fa"],
    output:
        directory(f"{RESULTS_DIR}/star_index/contaminants"),
    log:
        f"{LOG_DIR}/star/contaminant_index.log",
    params:
        log_prefix=f"{LOG_DIR}/star/contaminant_index.",
    conda:
        "../envs/star.yaml"
    threads: 8
    shell:
        r"""
        exec > {log} 2>&1
        awk '
            function end_record() {{
                if (seq in first) {{ print "ERROR: " name " has the same sequence as " first[seq]; bad = 1 }}
                else first[seq] = name
            }}
            /^>/ {{
                if (n++) end_record()
                name = substr($1, 2); seq = ""
                if (name in names) {{ print "ERROR: duplicate record name " name; bad = 1 }}
                names[name]
                next
            }}
            {{ seq = seq toupper($0) }}
            END {{
                end_record()
                if (bad) print "ERROR: contaminants_fa needs unique record names and sequences, as build_contaminants writes it"
                exit bad
            }}
        ' {input}
        STAR --genomeSAindexNbases 9 --genomeChrBinNbits 8 --runThreadN {threads} --runMode genomeGenerate \
            --genomeFastaFiles {input} --genomeDir {output} \
            --outFileNamePrefix {params.log_prefix}
        """


# Aligns the trimmed reads to the contaminants once: the unaligned reads are the
# clean fastq, the aligned ones (one alignment each, with AS) are kept for the
# competitive filter below.
rule star_contaminant:
    input:
        fastq=f"{EXP_DIR}/trim_reads/{{sample}}.fastq.gz",
        index=f"{RESULTS_DIR}/star_index/contaminants",
    output:
        fastq=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.fastq.gz",
        log_final=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
        bam=temp(f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Aligned.out.bam"),
    log:
        f"{EXP_LOG_DIR}/star/contaminant/{{sample}}.log",
    params:
        prefix=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_",
        tmp_dir=f"{EXP_DIR}/filter_reads/{{sample}}/_tmpSTAR",
        star_args=CONTAMINANT_STAR_ARGS,
    conda:
        "../envs/star.yaml"
    threads: 8
    resources:
        mem_mb=16000,
    shell:
        r"""
        exec > {log} 2>&1
        rm -rf {params.tmp_dir}
        STAR \
            --genomeDir {input.index} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outMultimapperOrder Random \
            --outSAMmultNmax 1 \
            {params.star_args} \
            --outTmpDir {params.tmp_dir} \
            --genomeLoad NoSharedMemory \
            --outSAMattributes NH HI AS NM MD \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.prefix} \
            --outReadsUnmapped Fastx \
            --readFilesIn {input.fastq}

        too_many=$(grep "Number of reads mapped to too many loci" {params.prefix}Log.final.out | cut -f2)
        if [ "$too_many" != "0" ]; then
            echo "ERROR: $too_many reads exceeded --outFilterMultimapNmax 10000 in contaminant mapping; raise the limit." >&2
            exit 1
        fi

        gzip -c {params.prefix}Unmapped.out.mate1 > {output.fastq}
        rm -f {params.prefix}Unmapped.out.mate1
        """


# The competitive filter: of the reads star_contaminant aligned, put back the
# ones whose best sense alignment to the transcriptome scores higher than their
# best sense contaminant alignment (scripts/contaminant_compete.py). The
# transcriptome STAR call uses star_transcript's settings
# (RIBO_TRANSCRIPTOME_STAR_ARGS, common.smk), plus the AS attribute.
# star_contaminant reported one of a read's best contaminant alignments, in
# either orientation. Where it is antisense, the read is aligned to the
# contaminants again with star_contaminant's settings (CONTAMINANT_STAR_ARGS),
# keeping all its alignments, for the best sense one whatever its score: under
# STAR's default --outFilterMultimapScoreRange of 1, a sense alignment more than
# 1 below the antisense one would not be reported.
rule contaminant_compete_align:
    input:
        bam=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Aligned.out.bam",
        contaminant_index=f"{RESULTS_DIR}/star_index/contaminants",
        transcriptome_index=f"{RESULTS_DIR}/star_index/transcriptome",
    output:
        removed=temp(f"{EXP_DIR}/filter_reads/{{sample}}/compete/removed.fastq.gz"),
        contaminant=temp(f"{EXP_DIR}/filter_reads/{{sample}}/compete/contaminant.tsv.gz"),
        contaminant_sense=temp(f"{EXP_DIR}/filter_reads/{{sample}}/compete/contaminant_sense.tsv.gz"),
        transcriptome=temp(f"{EXP_DIR}/filter_reads/{{sample}}/compete/transcriptome.tsv.gz"),
        log_final=f"{EXP_DIR}/filter_reads/{{sample}}/compete/{{sample}}.transcriptome_Log.final.out",
    # the sample in STAR's file names: MultiQC names a STAR report after its file
    params:
        prefix=f"{EXP_DIR}/filter_reads/{{sample}}/compete/{{sample}}.",
        contaminant_star_args=CONTAMINANT_STAR_ARGS,
        star_args=RIBO_TRANSCRIPTOME_STAR_ARGS,
    log:
        f"{EXP_LOG_DIR}/contaminant_compete/{{sample}}.align.log",
    conda:
        "../envs/star.yaml"
    threads: 8
    resources:
        mem_mb=24000,
    shell:
        r"""
        exec > {log} 2>&1
        rm -rf {params.prefix}contaminant_tmpSTAR {params.prefix}transcriptome_tmpSTAR
        samtools fastq {input.bam} | gzip > {output.removed}
        # read, flag, contaminant record, AS
        samtools view {input.bam} |
            awk -F'\t' -v OFS='\t' '{{for (i = 12; i <= NF; i++) if ($i ~ /^AS:i:/) print $1, $2, $3, substr($i, 6)}}' |
            gzip > {output.contaminant}

        # all contaminant alignments of the reads whose reported one is antisense
        samtools fastq -f 16 {input.bam} | gzip > {params.prefix}antisense.fastq.gz
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.contaminant_index} \
            --outTmpDir {params.prefix}contaminant_tmpSTAR \
            --outSAMtype BAM Unsorted \
            --outSAMmode NoQS \
            --outSAMattributes AS \
            {params.contaminant_star_args} \
            --outSAMmultNmax -1 \
            --outFilterMultimapScoreRange 1000 \
            --outFileNamePrefix {params.prefix}contaminant_ \
            --readFilesIn {params.prefix}antisense.fastq.gz \
            --readFilesCommand zcat
        cat {params.prefix}contaminant_Log.final.out
        antisense=$(samtools view -c -f 16 {input.bam})
        aligned=$(awk -F'|' '/Uniquely mapped reads number|Number of reads mapped to multiple loci/ {{n += $2}}
            END {{print n + 0}}' {params.prefix}contaminant_Log.final.out)
        if [ "$aligned" != "$antisense" ]; then
            echo "ERROR: $aligned of the $antisense reads with an antisense contaminant alignment aligned again" >&2
            exit 1
        fi
        # read, AS of each sense alignment
        samtools view -F 16 {params.prefix}contaminant_Aligned.out.bam |
            awk -F'\t' -v OFS='\t' '{{for (i = 12; i <= NF; i++) if ($i ~ /^AS:i:/) print $1, substr($i, 6)}}' |
            gzip > {output.contaminant_sense}
        # also so that MultiQC does not take its Log.final.out for a second report of the sample (path */compete/*)
        rm -f {params.prefix}antisense.fastq.gz {params.prefix}contaminant_Aligned.out.bam \
            {params.prefix}contaminant_Log.final.out {params.prefix}contaminant_Log.out \
            {params.prefix}contaminant_Log.progress.out {params.prefix}contaminant_SJ.out.tab

        # star_transcript's alignment of those reads
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.transcriptome_index} \
            --outTmpDir {params.prefix}transcriptome_tmpSTAR \
            --outSAMtype BAM Unsorted \
            --outSAMmode NoQS \
            --outSAMattributes NH AS NM \
            {params.star_args} \
            --outFileNamePrefix {params.prefix}transcriptome_ \
            --readFilesIn {output.removed} \
            --readFilesCommand zcat
        # read, flag, AS
        samtools view {params.prefix}transcriptome_Aligned.out.bam |
            awk -F'\t' -v OFS='\t' '{{for (i = 12; i <= NF; i++) if ($i ~ /^AS:i:/) print $1, $2, substr($i, 6)}}' |
            gzip > {output.transcriptome}
        rm -f {params.prefix}transcriptome_Aligned.out.bam
        """


rule contaminant_compete:
    input:
        clean=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.fastq.gz",
        contam_log=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.contam_Log.final.out",
        removed=f"{EXP_DIR}/filter_reads/{{sample}}/compete/removed.fastq.gz",
        contaminant=f"{EXP_DIR}/filter_reads/{{sample}}/compete/contaminant.tsv.gz",
        contaminant_sense=f"{EXP_DIR}/filter_reads/{{sample}}/compete/contaminant_sense.tsv.gz",
        transcriptome=f"{EXP_DIR}/filter_reads/{{sample}}/compete/transcriptome.tsv.gz",
        script=workflow.source_path("../scripts/contaminant_compete.py"),
    output:
        fastq=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.competitive.fastq.gz",
        stats=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.compete.tsv",
        records=f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.compete_records.tsv",
    log:
        f"{EXP_LOG_DIR}/contaminant_compete/{{sample}}.log",
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb=16000,
    shell:
        "python {input.script} {input.clean} {input.contam_log} {input.removed} {input.contaminant} "
        "{input.contaminant_sense} {input.transcriptome} {output.fastq} {output.stats} {output.records} 2> {log}"
