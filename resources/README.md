# resources/

Scratch space for locally-generated reference files. Everything here is gitignored except
this README and `contaminants_built.fa` (see `.gitignore`). This README documents the
provenance of the files expected here.

## contaminants_built.fa

Produced by the `build_contaminants` rule (`workflow/rules/contaminants.smk`), which runs
`workflow/scripts/build_contaminants.py` in the `python` conda env
(`workflow/envs/python.yaml`: python 3.14.7). The rule is not part of the default targets —
it needs internet access, so it's run explicitly and the output is then copied here. It is
the default `contaminants_fa` (`workflow/rules/common.smk`), which goes to STAR as is.

The script keeps the first record of each sequence (one strand) and prefixes every name with
a running number (`1_human_...`): Ensembl gene symbols repeat (839 `human_Y_RNA`), and STAR
needs unique names. Of the 8,448 downloaded records, 7,398 remain.

Generated 2026-09-30; deduplicated and numbered 2026-10-02 from that download (no re-fetch),
byte-identical to what the former `number_contaminants` rule made of it. Built from:

- Ensembl release 109, `Homo_sapiens.GRCh38.ncrna.fa.gz` and
  `Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz` — all ncRNA biotypes except `lncRNA`
  (rRNA, tRNA, snRNA, snoRNA, scaRNA, miRNA, misc_RNA, ...). Release 109 matches
  `spike_in_transcriptome_fa` elsewhere in the pipeline.
- GtRNAdb `hg38-tRNAs.fa` (eukaryota/Hsapi38) — human nuclear tRNA genes, which Ensembl
  does not annotate.
- NCBI RefSeq, fetched by accession (efetch, db=nuccore): `NR_003286.2` (18S),
  `NR_003287.2` (28S), `NR_003285.2` (5.8S), `NR_023379.1` (5S) — canonical full-length
  rRNAs; Ensembl's human rRNA biotype entries are fragmentary pseudogene copies.
- 12 Illumina TruSeq/RPI small-RNA index-barcode sequences and 2 synthetic RNA
  size-marker oligos, hardcoded in the script (`LITERAL_SEQUENCES`) — wet-lab protocol
  sequences not in any database, copied verbatim from the lab's prior
  `contaminants.combined_human_yeast.fa`.

See the script's docstring and `CLAUDE.md` for further design context.
