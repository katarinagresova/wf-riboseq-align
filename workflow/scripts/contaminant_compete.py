"""Competitive contaminant filter (config contaminant_filter: competitive).

star_contaminant discards every read that aligns to a contaminant. This puts a discarded read
back when its best sense alignment to the transcriptome scores higher (STAR AS) than its
contaminant alignment. Antisense alignments do not count: the libraries are stranded.

A tie stays discarded. Ties are reads whose sequence lies in both references, and most are not
footprints: rRNA, snoRNA and 7SL sequence inside transcriptome entries, yeast rRNA inside the
rDNA ORF YLR154C-G. Measured on eIF4E 4h ribo_07 (validation in wf-eIF-deltaTE): 1.68M ties, of
which the few that land in a human CDS pile up (28S rRNA on LSM4 at 775x the rest of its CDS);
the cost is the footprints inside exonic miRNA hairpins (SERF2, RPS27A).

Classes, per read the contaminant alignment removed:
  no_transcriptome    no alignment to the transcriptome
  antisense_only      antisense alignments only
  contaminant_better  contaminant AS > best sense AS
  tie                 equal: stays discarded
  rescued             best sense AS > contaminant AS: put back

Inputs are contaminant_compete_align's: the removed reads (fastq), read / record / AS of their
contaminant alignment, read / flag / AS of their transcriptome alignments. The read count must
equal what star_contaminant removed (its Log.final.out), so the two alignments agree.

Writes star_contaminant's clean fastq followed by the rescued reads (a second gzip member), the
class counts (stats.tsv, `class<TAB>reads`, plus a `removed` row) and the counts per contaminant
record (records.tsv).

Usage: contaminant_compete.py <clean.fastq.gz> <contam_Log.final.out> <removed.fastq.gz>
       <contaminant.tsv.gz> <transcriptome.tsv.gz> <out.fastq.gz> <stats.tsv> <records.tsv>
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
    if len(sys.argv) != 9:
        sys.exit(__doc__)
    clean, contam_log, removed_fq, contaminant_tsv, transcriptome_tsv, out_fq, stats_out, records_out = sys.argv[1:]
    expected = star_removed(contam_log)

    contaminant = {}
    with gzip.open(contaminant_tsv, "rt") as f:
        for line in f:
            name, record, score = line.rstrip("\n").split("\t")
            contaminant[name] = (sys.intern(record), int(score))
    if len(contaminant) != expected:
        sys.exit(f"{len(contaminant)} reads align to the contaminants here, star_contaminant removed {expected}")

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
    for name, (record, score) in contaminant.items():
        best = sense.get(name)
        if best is None:
            cls = "antisense_only" if name in antisense else "no_transcriptome"
        elif best > score:
            cls = "rescued"
            rescued.add(name)
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
    with open(records_out, "w") as o:
        o.write("record\treads\t" + "\t".join(CLASSES) + "\n")
        for record, c in sorted(per_record.items(), key=lambda kv: -sum(kv[1].values())):
            o.write(f"{record}\t{sum(c.values())}\t" + "\t".join(str(c[k]) for k in CLASSES) + "\n")
    print(", ".join(f"{cls} {totals[cls]}" for cls in CLASSES), file=sys.stderr)


if __name__ == "__main__":
    main()
