# Config and samples.csv -> the globals and helpers the rules use.

import re

import pandas as pd
from snakemake.utils import validate

# also fills in the defaults (RESULTS_DIR, LOG_DIR, min_spike_in_fraction)
validate(config, "../schemas/config.schema.yaml")

RESULTS_DIR = config["RESULTS_DIR"]
LOG_DIR = config["LOG_DIR"]

samples = pd.read_csv(config["samples"], dtype=str).set_index("sample_id", drop=False)
validate(samples, "../schemas/samples.schema.yaml")
if samples.index.duplicated().any():
    raise ValueError(f"samples.csv: duplicate sample_id {sorted(set(samples.index[samples.index.duplicated()]))}")

# {sample} only ever matches a whole sample_id, so a pattern like fastqc's
# {sample}_R{mate} cannot split one differently
wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples.index),


# Optional `experiment` column: several experiments in one run. Each gets its
# own outputs, incl. its own RNA-seq-filtered reference, under
# {RESULTS_DIR}/{experiment}/ (logs under {LOG_DIR}/{experiment}/). Without the
# column the run is one experiment, directly in RESULTS_DIR / LOG_DIR. The
# experiment-independent outputs (contaminant set and index, unfiltered RNA-seq
# index) are always directly in RESULTS_DIR (reference/, star_index/) and
# LOG_DIR (star/), so they are built once.
if "experiment" in samples.columns:
    if samples["experiment"].isna().any():
        raise ValueError("samples.csv: with an experiment column, every row needs an experiment")
    EXPERIMENTS = sorted(set(samples["experiment"]))
    EXP_DIR = f"{RESULTS_DIR}/{{experiment}}"
    EXP_LOG_DIR = f"{LOG_DIR}/{{experiment}}"

    wildcard_constraints:
        experiment="|".join(re.escape(e) for e in EXPERIMENTS),
else:
    EXPERIMENTS = [None]
    EXP_DIR = RESULTS_DIR
    EXP_LOG_DIR = LOG_DIR


def experiment_samples(experiment, read_type):
    """One experiment's sample_ids of one read type (None: the whole run), in samples.csv order."""
    sub = samples if experiment is None else samples[samples["experiment"] == experiment]
    return list(sub.loc[sub["read_type"] == read_type, "sample_id"])


def experiment_files(pattern, experiment, read_type, **wildcards):
    """`pattern` expanded over one experiment's samples of one read type."""
    if experiment is not None:
        wildcards["experiment"] = experiment
    return expand(pattern, sample=experiment_samples(experiment, read_type), **wildcards)


MODE = config["mode"]
if MODE == "filtered":
    no_rna = [e for e in EXPERIMENTS if not experiment_samples(e, "rna")]
    if no_rna:
        raise ValueError("mode 'filtered' builds the ribo reference from the RNA-seq, "
                         "but samples.csv has no read_type=rna rows"
                         + ("" if no_rna == [None] else f" for experiment(s) {no_rna}"))
    HUMAN_TRANSCRIPTOME_FA = f"{EXP_DIR}/reference/human_transcriptome.rnaseq_filtered.fa"
elif MODE == "unfiltered":
    HUMAN_TRANSCRIPTOME_FA = config["human_transcriptome_fa"]

# What star_transcript maps: star_contaminant's clean fastq plus the reads
# contaminant_compete puts back.
CLEAN_FASTQ = f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.competitive.fastq.gz"


# star_transcript's alignment settings. contaminant_compete_align uses the same
# ones: it compares a read's contaminant AS with the AS star_transcript would
# give it.
RIBO_TRANSCRIPTOME_STAR_ARGS = (
    "--seedSearchLmax 10 "
    "--outFilterMultimapNmax 255 "
    "--outFilterMismatchNmax 2 "
    "--outFilterMultimapScoreRange 0 "
    "--outFilterIntronMotifs RemoveNoncanonical"
)
