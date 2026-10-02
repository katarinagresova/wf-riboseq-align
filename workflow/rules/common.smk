# Config and samples.csv -> the globals and helpers the rules use.

import re

import pandas as pd
from snakemake.utils import validate

# also fills in the defaults (RESULTS_DIR, LOG_DIR, min_spike_in_fraction)
validate(config, "../schemas/config.schema.yaml")
# the set this repo ships; a path the schema cannot default (relative to the
# repo, not to where snakemake runs), and through the source cache as a module
config.setdefault("contaminants_fa", workflow.source_path("../../resources/contaminants_built.fa"))

RESULTS_DIR = config["RESULTS_DIR"]
LOG_DIR = config["LOG_DIR"]

# read_type=rna rows are validated, then skipped (every target is of the ribo
# rows), so an importer can pass the samples.csv it quantifies its RNA-seq from
samples = pd.read_csv(config["samples"], dtype=str).set_index("sample_id", drop=False)
validate(samples, "../schemas/samples.schema.yaml")
if samples.index.duplicated().any():
    raise ValueError(f"samples.csv: duplicate sample_id {sorted(set(samples.index[samples.index.duplicated()]))}")

# {sample} only ever matches a whole sample_id
wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples.index),


# Optional `experiment` column: several experiments in one run. Each gets its
# own outputs under {RESULTS_DIR}/{experiment}/ (logs under
# {LOG_DIR}/{experiment}/). Without the column the run is one experiment,
# directly in RESULTS_DIR / LOG_DIR. The experiment-independent outputs (the
# contaminant set, the Ribo-seq reference and both STAR indexes) are always
# directly in RESULTS_DIR (reference/, star_index/) and LOG_DIR (star/), so
# they are built once.
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


# What star_transcript maps: star_contaminant's clean fastq plus the reads
# contaminant_compete puts back.
CLEAN_FASTQ = f"{EXP_DIR}/filter_reads/{{sample}}/{{sample}}.competitive.fastq.gz"


# star_transcript's alignment settings. contaminant_compete_align uses the same
# ones: it compares a read's contaminant AS with the AS star_transcript would
# give it. --outMultimapperOrder Random makes star_transcript's primary pick
# among equally-good loci random, as its comment claims (STAR's default order
# is fixed, not random). It is not reproducible with > 1 thread, seed or not
# (777 is STAR's default): each thread seeds its RNG with runRNGseed *
# (thread + 1), and which thread maps which reads varies from run to run, so a
# multimapper's primary can change between runs; its alignments do not.
# Accepted (2026-10-02). contaminant_compete_align ignores which alignment is
# primary, so this does not change its output.
RIBO_TRANSCRIPTOME_STAR_ARGS = (
    "--seedSearchLmax 10 "
    "--outFilterMultimapNmax 255 "
    "--outFilterMismatchNmax 2 "
    "--outFilterMultimapScoreRange 0 "
    "--alignIntronMax 1 "
    "--alignEndsType Extend5pOfRead1 "
    "--outMultimapperOrder Random "
    "--runRNGseed 777"
)


# star_contaminant's alignment settings. contaminant_compete_align realigns
# with the same ones the reads whose one reported contaminant alignment
# (star_contaminant: --outSAMmultNmax 1) is antisense, for their best sense
# contaminant AS.
CONTAMINANT_STAR_ARGS = (
    "--outFilterMultimapNmax 10000 "
    "--alignIntronMax 1 "
    "--alignEndsType Extend5pOfRead1"
)
