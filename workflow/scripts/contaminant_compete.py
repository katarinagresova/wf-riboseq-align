"""Competitive contaminant filter.

Of the reads star_contaminant aligns to a contaminant, this puts one back when its best sense
alignment to the transcriptome scores higher (STAR AS) than its best sense contaminant alignment.
Antisense alignments do not count, on either side: the libraries are stranded, and the
contaminant set is RNA in its own orientation, so a read antisense to a contaminant does not come
from it. star_contaminant reports one of a read's best contaminant alignments; when that one is
antisense, contaminant_compete_align has aligned the read again for its sense alignments. A read
with none is not a contaminant read, and goes back if it has a sense transcriptome alignment.
Measured on eIF4E 4h ribo_07: of the 7.16M reads star_contaminant removed, 31,168 had an antisense
alignment reported, and 2,112 of those now go back that a strand-blind comparison kept discarded
(1,827 ties, 285 contaminant_better), e.g. 127 on SERF2, antisense to MIR1282, and 119 on yeast
CUP1-1/CUP1-2, antisense to RUF5-2.

A tie stays discarded. Ties are reads whose sequence lies in both references, and most are not
footprints: rRNA, snoRNA and 7SL sequence inside transcriptome entries, yeast rRNA inside the
rDNA ORF YLR154C-G. Measured on eIF4E 4h ribo_07 (validation in wf-eIF-deltaTE): 1.68M ties, of
which the few that land in a human CDS pile up (28S rRNA on LSM4 at 775x the rest of its CDS);
the cost is the footprints inside exonic miRNA hairpins (RPS27A).

Classes, per read the contaminant alignment removed:
  no_transcriptome    no alignment to the transcriptome
  antisense_only      antisense alignments only
  contaminant_better  best sense contaminant AS > best sense AS
  tie                 equal: stays discarded
  rescued             best sense AS > best sense contaminant AS, or no sense contaminant
                      alignment: put back

Inputs are contaminant_compete_align's: the removed reads (fastq), read / flag / record / AS of
their reported contaminant alignment, read / AS of the sense contaminant alignments of the reads
where that one is antisense, read / flag / AS of their transcriptome alignments. The read count
must equal what star_contaminant aligned (its Log.final.out).

Writes star_contaminant's clean fastq followed by the rescued reads (a second gzip member), the
class counts (stats.tsv, `class<TAB>reads`, plus a `removed` row and a `rescued_by_strand` row:
the rescued reads whose reported contaminant alignment, antisense, scores at least as high as their
best sense transcriptome one) and the counts per contaminant record (records.tsv; the record is
the reported alignment's, even if antisense). The log also has the reads whose reported
contaminant alignment is antisense (antisense_picks) and how many of them have a sense one
(antisense_picks_with_sense).

Usage: contaminant_compete.py <clean.fastq.gz> <contam_Log.final.out> <removed.fastq.gz>
       <contaminant.tsv.gz> <contaminant_sense.tsv.gz> <transcriptome.tsv.gz> <out.fastq.gz>
       <stats.tsv> <records.tsv>
"""
import gzip
import shutil
import sys
from collections import Counter, defaultdict

CLASSES = ["no_transcriptome", "antisense_only", "contaminant_better", "tie", "rescued"]


def star_removed(log):
    """Reads star_contaminant aligned (unique + multiple loci; it fails on any too-many-loci read)."""
    stats = dict(line.split("|", 1) for line in open(log) if "|" in line)
    stats = {k.strip(): v.strip() for k, v in stats.items()}
    return int(stats["Uniquely mapped reads number"]) + int(stats["Number of reads mapped to multiple loci"])


def main():
    if len(sys.argv) != 10:
        sys.exit(__doc__)
    (clean, contam_log, removed_fq, contaminant_tsv, contaminant_sense_tsv, transcriptome_tsv, out_fq, stats_out,
     records_out) = sys.argv[1:]
    expected = star_removed(contam_log)

    contaminant = {}
    with gzip.open(contaminant_tsv, "rt") as f:
        for line in f:
            name, flag, record, score = line.rstrip("\n").split("\t")
            contaminant[name] = (sys.intern(record), int(score), bool(int(flag) & 16))
    if len(contaminant) != expected:
        sys.exit(f"{len(contaminant)} reads align to the contaminants here, star_contaminant removed {expected}")

    sense_contaminant = {}  # read whose reported contaminant alignment is antisense -> its best sense AS
    with gzip.open(contaminant_sense_tsv, "rt") as f:
        for line in f:
            name, score = line.rstrip("\n").split("\t")
            if not contaminant.get(name, (None, None, False))[2]:
                sys.exit(f"{name} has a sense contaminant alignment in {contaminant_sense_tsv}, but no antisense one "
                         f"reported in {contaminant_tsv}")
            if name not in sense_contaminant or int(score) > sense_contaminant[name]:
                sense_contaminant[name] = int(score)

    sense, antisense = {}, set()
    with gzip.open(transcriptome_tsv, "rt") as f:
        for line in f:
            name, flag, score = line.rstrip("\n").split("\t")
            if int(flag) & 16:
                antisense.add(name)
            elif name not in sense or int(score) > sense[name]:
                sense[name] = int(score)

    per_record = defaultdict(Counter)
    rescued = set()
    rescued_by_strand = antisense_picks = 0
    for name, (record, reported, antisense_pick) in contaminant.items():
        # star_contaminant reported a best alignment: if sense, it is the best sense one
        score = sense_contaminant.get(name) if antisense_pick else reported
        antisense_picks += antisense_pick
        best = sense.get(name)
        if best is None:
            cls = "antisense_only" if name in antisense else "no_transcriptome"
        elif score is None or best > score:
            cls = "rescued"
            rescued.add(name)
            rescued_by_strand += best <= reported
        elif best == score:
            cls = "tie"
        else:
            cls = "contaminant_better"
        per_record[record][cls] += 1

    n = 0
    with open(out_fq, "wb") as out:
        with open(clean, "rb") as f:
            shutil.copyfileobj(f, out)
        with gzip.open(removed_fq, "rt") as f, gzip.GzipFile(fileobj=out, mode="wb", compresslevel=6) as z:
            while True:
                read = [f.readline() for _ in range(4)]
                if not read[0]:
                    break
                n += 1
                if read[0][1:].split()[0] in rescued:
                    z.write("".join(read).encode())
    if n != expected:
        sys.exit(f"{removed_fq} has {n} reads, star_contaminant removed {expected}")

    totals = Counter()
    for c in per_record.values():
        totals.update(c)
    with open(stats_out, "w") as o:
        o.write(f"class\treads\nremoved\t{expected}\n")
        for cls in CLASSES:
            o.write(f"{cls}\t{totals[cls]}\n")
        o.write(f"rescued_by_strand\t{rescued_by_strand}\n")
    with open(records_out, "w") as o:
        o.write("record\treads\t" + "\t".join(CLASSES) + "\n")
        for record, c in sorted(per_record.items(), key=lambda kv: -sum(kv[1].values())):
            o.write(f"{record}\t{sum(c.values())}\t" + "\t".join(str(c[k]) for k in CLASSES) + "\n")
    print(", ".join(f"{cls} {totals[cls]}" for cls in CLASSES), file=sys.stderr)
    print(f"rescued_by_strand\t{rescued_by_strand}", file=sys.stderr)
    print(f"antisense_picks\t{antisense_picks}", file=sys.stderr)
    print(f"antisense_picks_with_sense\t{len(sense_contaminant)}", file=sys.stderr)


if __name__ == "__main__":
    main()
