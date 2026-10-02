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
# contaminant set, the species BEDs and the bowtie index) are always directly
# in RESULTS_DIR (reference/, bowtie_index/) and LOG_DIR, so they are built
# once.
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


# bowtie_align: forward strand only (stranded libraries), end to end with at
# most 2 mismatches, up to 256 alignments per read, best first, in input order
# (the same output for any thread count). 256 in the best stratum = more than
# 255 loci: assign_reads.py drops the read as too_many_loci.
BOWTIE_ARGS = "-v 2 --norc -k 256 --best --reorder"
# A read shorter than SHORT_READ_LENGTH keeps only its alignments with at most
# SHORT_READ_MISMATCHES mismatches. 2 = no cap, until it is chosen
# (PLAN_BOWTIE.md V3).
SHORT_READ_LENGTH = 23
SHORT_READ_MISMATCHES = 2
