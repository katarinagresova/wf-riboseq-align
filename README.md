# wf-riboseq-align

A Snakemake workflow that takes raw Ribo-seq fastq files to **transcriptome-aligned BAM files**, split into the human
reads and the spike-in reads. It also quantifies the matched total RNA-seq with salmon. It is a generalised,
reproducible version of the alignment part of the eIF project pipeline (Gabriel Villamil, Frederick Korbel). It can be
run on its own or imported into another Snakemake workflow.

## What it does

Ribo-seq (single-end, with 4+4 nt randomers around the insert):

1. `cutadapt_reads`: trim the 3' adapter, keeping reads of 18–100 nt
2. `collapse_reads`: collapse identical reads. The randomers are still part of the read at this point, so PCR
   duplicates collapse but distinct molecules stay distinct
3. `trim_reads`: strip the randomers and append them to the read name
4. `star_contaminant`: discard reads that map to rRNA / tRNA / other contaminants
5. `star_transcript`: map to the human + spike-in transcriptome (multimappers kept, up to 255 loci)
6. `split_bam_transcriptome`: split into `human/` and `spike_in/` BAM files, sorted and indexed
7. `qc_library`, `qc_summary`: per-library read counts: contaminant step (input, removed as unique hits, removed as
   multi-copy hits, clean) and species split (human only, spike-in only, both, spike-in fraction = spike-in-only /
   all mapped). Reads aligned to both species stay in both BAM files and are only counted. The run fails if any
   library's spike-in fraction is below `min_spike_in_fraction` (default 0.01; set 0 for libraries without spike-in)

Total RNA-seq (paired-end): map to the human transcriptome with STAR, then quantify with `salmon quant`.

## The one choice: `mode`

The two eIF pipelines differ in a single step: which human transcriptome the Ribo-seq reads are mapped to.

| `mode` | Human reference for Ribo-seq | Used by |
|---|---|---|
| `filtered` | Transcripts the matched RNA-seq expresses: TPM ≥ `min_tpm` in ≥ `min_samples` libraries. Filtered per run. | eIF depletion experiments |
| `unfiltered` | `human_transcriptome_fa` as given | eIF harringtonine experiments |

Both modes append the spike-in transcriptome. For libraries without a spike-in this is harmless: the spike_in BAM is
simply (nearly) empty. `filtered` needs `read_type=rna` rows in `samples.csv`.

## Inputs

`config/samples.csv`, one row per library:

| column | |
|---|---|
| `sample_id` | free label, used in output names. It also becomes the read-name prefix of Ribo-seq reads |
| `read_type` | `ribo` or `rna` |
| `fastq_1` | fastq.gz (R1 for RNA-seq) |
| `fastq_2` | R2 fastq.gz for RNA-seq; empty for Ribo-seq |

`config/config.yaml` contains `mode`, the cutadapt parameters, the contaminant fasta, the human transcriptome
fasta + gtf, the spike-in transcriptome fasta, `min_spike_in_fraction`, and the `autofilter` thresholds. Every key is commented in the file.
Fill in the placeholder paths with your own data and references before running.

## Outputs (under `RESULTS_DIR`, default `results/`)

```
split_bam/transcriptome/human/<sample>.bam(.bai)     Ribo-seq, human transcripts
split_bam/transcriptome/spike_in/<sample>.bam(.bai)  Ribo-seq, spike-in reads
salmon/<sample>/quant.sf                             RNA-seq quantification
qc/<sample>.stats.tsv, qc/summary.tsv                Ribo-seq read counts per library (collapsed reads)
reference/                                           the references these were aligned to
  transcriptome.combined_human_spike_in.fa           (the BAM @SQ lines refer to this)
  human_transcriptome.rnaseq_filtered.{fa,gtf}       (mode filtered)
  rnaseq_filter_blacklist_txid.txt                   (mode filtered: the transcripts dropped)
star/, filter_reads/, ...                            intermediates, incl. STAR Log.final.out
```

Per-step logs are written to `LOG_DIR` (default `logs/`). `logs/collapse_reads/` holds read counts, read-length
distributions and base composition.

## Setup

1. **Conda or mamba**, with the `bioconda` and `conda-forge` channels reachable. Used both to create the env below
   and, via `use-conda: true` in the profiles, to auto-build one isolated env per rule from `workflow/envs/*.yaml`
   on first run.
2. **A `snake` env** with Snakemake and pandas (the `Snakefile` itself does `import pandas`, outside any per-rule
   env):
   ```bash
   mamba create -n snake -c conda-forge -c bioconda "snakemake=9.11" pandas
   ```
   Running on SLURM also needs the executor plugin in the same env:
   ```bash
   mamba install -n snake -c conda-forge -c bioconda snakemake-executor-plugin-slurm
   ```
3. **`config/config.yaml` and `config/samples.csv`** — fill in your own paths (see Inputs above).
4. **The execution profile you'll use**:
   - Local (`workflow/profiles/default/config.yaml`): adjust `cores` to your machine.
   - SLURM (`workflow/profiles/slurm/config.yaml`): set `slurm_account` and `slurm_partition` to your cluster's
     values. The shipped defaults are placeholders and will fail on any other cluster.

## Run it

```bash
./snakemake.sh -n          # dry run: always do this first
./snakemake.sh             # local (workflow/profiles/default)
sbatch slurm_job.sh        # SLURM: one job per rule (workflow/profiles/slurm)
```

Both wrappers activate the `snake` conda env. To use your own config: `./snakemake.sh --configfile my_config.yaml`.

## Use it from another workflow

Import it as a Snakemake module and give it its own config block:

```python
# Snakefile of the importing workflow
configfile: "config/config.yaml"      # contains an `align:` block shaped like this repo's config.yaml

module align:
    snakefile: "/path/to/wf-riboseq-align/workflow/Snakefile"
    config: config["align"]

use rule * from align as align_*
```

Use a different `RESULTS_DIR` for each experiment (for example `results/align/eIF4E_4h`): the `filtered` reference is
built per run. Relative paths in the `align` block are resolved from the importing workflow's directory. Its rules
can then consume, for example, `rules.align_split_bam_transcriptome.output.bam`.

## Reproducibility

- Every rule has a conda env with exact version pins. STAR 2.7.10b, samtools 1.17 and salmon 1.10.2 are the versions
  in the eIF pipeline's container.
- cutadapt is 5.2 on Python 3.13, because no cutadapt build exists for Python 3.14 yet. On 2M eIF4E reads its
  output is byte-identical to 4.4, the container's version, which in turn reproduces the eIF pipeline's trimmed
  reads exactly.
- The helper scripts are Python (3.14) ports of the pipeline's Perl and R scripts. Each was checked against the
  original, or against the eIF pipeline's saved output, on the eIF4E 4 h data. All gave identical output: the
  numbered contaminant fasta, the RNA-seq blacklist (1,673 transcripts), the filtered fasta and gtf, the collapse
  statistics (41.1M → 12.4M reads), the collapsed reads themselves, and the randomer trimming.

## Differences from the eIF pipeline

- **Reads with any contaminant hit are removed.** The eIF pipeline ran the contaminant alignment with
  `--outFilterMultimapNmax 1`. A read hitting more than one contaminant entry then counts as "mapped to too many
  loci", ends up in the unmapped output and passes as clean: 2.26M reads in ribo_07, 34% of its "clean" fastq. Here
  the limit is 10,000, and the rule stops with an error if STAR's `Log.final.out` reports any read mapped to too many
  loci, so a larger copy number cannot leak silently. The limit is this high because a `build_contaminants` fasta
  takes all of Ensembl ncRNA, which likely includes many near-identical snRNA, Y_RNA and 5S copies, and
  `number_contaminants` only removes exact duplicates. `--winAnchorMultimapNmax` stays at its default of 50: raising
  it to 10,000 as well caught only 0.3% more reads, but left 250–470k more reads per library unmapped (and so
  "clean"), including 151k in ribo_07 that map uniquely at the default. Measured on the full eIF4E 4 h data: in
  ribo_07 every other STAR category is identical to the eIF pipeline's, only the leaked reads move from "too many
  loci" to "multiple loci", and the clean fastq drops from 6.58M to 4.32M reads. Across the six libraries, reads in
  the human BAMs drop by 20–38% and in the spike-in BAMs by 55–76%. Every output from `filter_reads` onward
  therefore differs from the eIF pipeline's; only the RNA-seq outputs are unchanged, because they skip the
  contaminant step.
- **Spike-in split in `unfiltered` mode.** The harringtonine pipeline built its Ribo-seq index from the human
  transcriptome only. It then labelled the last 6,612 human transcripts as "yeast", so its `yeast/` BAM files hold
  human reads and the spike-in was never aligned. Here the index is always human + spike-in, and the split checks
  that the index and fasta counts agree.
- **Exact transcript-id matching** when applying the blacklist. The original substring match would over-match
  versioned ids (`ENST…123.1` is contained in `ENST…123.10`). The output is identical on the HCT116 reference.
- **Corrupt contaminant records.** `contaminants.combined_human_yeast.fa` contains stray text in two records (in
  MIR4500HG, and a yeast_HRA1 header fused into tRNA-iMet). The workflow drops the non-nucleotide characters, as R's
  Biostrings did silently, and writes a warning to `logs/number_contaminants.log`.
- **Changed outputs.** The final BAM files are indexed, and `qc/` holds per-library read counts. The
  contaminant-alignment BAM + stats step and everything after the BAM (RiboStan, TPM tables, size factors) are not
  part of this workflow; normalisation stays downstream, as in the eIF pipeline.
- **Deterministic collapsing.** Reads with equal counts are ordered by first appearance. The Perl version used hash
  order, which is random, so read names could differ between runs.
