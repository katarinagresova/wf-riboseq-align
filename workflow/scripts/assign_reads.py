"""Assign each Ribo-seq read from its bowtie alignments (bowtie_align): contaminant, transcriptome or neither.

bowtie aligns every trimmed read once, to the human + spike-in transcripts and the contaminants in one
index, forward strand only (the libraries are stranded, and the contaminants are RNA in their own
orientation), end to end with at most 2 mismatches, up to 256 alignments per read, best (fewest
mismatches) first. This reads its SAM on stdin. A read shorter than <short_length> first loses its
alignments with more than <short_mismatches> mismatches. Then its best stratum (its alignments with the
fewest mismatches left) decides:

  unaligned           no alignment
  over_mismatch_cap   aligned, but only with more mismatches than the cap allows
  too_many_loci       256 alignments in the best stratum (more than 255 loci), none to a contaminant
  contaminant_only    a contaminant in the best stratum, no transcript anywhere
  contaminant_better  a contaminant in the best stratum, a transcript only at a worse one
  tie                 a contaminant and a transcript in the best stratum
  rescued             only transcripts in the best stratum, a contaminant at a worse one
  transcriptome       only transcripts, no contaminant

The first three and the contaminant classes are dropped; a tie stays dropped: its sequence lies in both
references, and most ties are rRNA, snoRNA or 7SL sequence inside transcriptome entries that pile up in
CDSs. A read with 256 alignments may have more that bowtie did not report: its class is from the reported
ones.

A kept read goes to the BAM with its best-stratum alignments only, NH = their number, MAPQ = STAR's
value for that (255 for 1, 3 for 2, 1 for 3-4, 0 for more), and one primary, picked by the CRC32 of
the read name (bowtie marks every alignment primary; the first would favour low reference ids, i.e.
human over spike-in). Tags NH NM MD, no base qualities, header without the contaminants.

Writes the BAM (uncompressed, in input order), one stats row (sample, input_reads, the classes,
contaminant_unique / contaminant_multi = the reads with a contaminant alignment, by the number of
those in their best contaminant stratum) and per contaminant record its reads by class (the record of
the read's first alignment in its best contaminant stratum). Fails unless each read's alignments are
best first.

Usage: assign_reads.py <sample> <contaminants.fa> <short_length> <short_mismatches> <out.bam> <stats.tsv>
       <records.tsv> < bowtie.sam
"""
import sys
import zlib
from collections import Counter, defaultdict

import pysam

CLASSES = ["unaligned", "over_mismatch_cap", "too_many_loci", "contaminant_only", "contaminant_better", "tie",
           "rescued", "transcriptome"]
CONTAMINANT_CLASSES = ["contaminant_only", "contaminant_better", "tie", "rescued"]
MAX_ALIGNMENTS = 256  # bowtie -k


def star_mapq(nh):
    return 255 if nh == 1 else 3 if nh == 2 else 1 if nh <= 4 else 0


def fasta_names(path):
    with open(path) as f:
        return {line[1:].split()[0] for line in f if line.startswith(">")}


def reads(sam):
    """(name, its records), per read: bowtie writes a read's records together."""
    name, group = None, []
    for a in sam:
        if a.query_name != name:
            if group:
                yield name, group
            name, group = a.query_name, []
        group.append(a)
    if group:
        yield name, group


def main():
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    sample, contaminants_fa, short_length, short_mismatches, bam_out, stats_out, records_out = sys.argv[1:]
    short_length, short_mismatches = int(short_length), int(short_mismatches)

    contaminant_names = fasta_names(contaminants_fa)
    totals = Counter()
    contaminant_unique = contaminant_multi = 0
    per_record = defaultdict(Counter)
    with pysam.AlignmentFile("-", "r") as sam:
        names = sam.references
        if not contaminant_names <= set(names):
            sys.exit(f"{len(contaminant_names - set(names))} records of {contaminants_fa} are not in the index")
        is_contaminant = [n in contaminant_names for n in names]
        header = sam.header.to_dict()
        header["SQ"] = [sq for sq in header["SQ"] if sq["SN"] not in contaminant_names]
        out_tid = {}  # input tid -> output tid, transcripts only
        for tid, contaminant in enumerate(is_contaminant):
            if not contaminant:
                out_tid[tid] = len(out_tid)

        with pysam.AlignmentFile(bam_out, "wbu", header=header) as out:
            for name, group in reads(sam):
                if group[0].is_unmapped:
                    totals["unaligned"] += 1
                    continue
                nms = [a.get_tag("NM") for a in group]
                if any(x > y for x, y in zip(nms, nms[1:])):
                    sys.exit(f"{name}: alignments not best first (NM {nms})")
                if group[0].query_length < short_length:
                    n = sum(nm <= short_mismatches for nm in nms)  # best first: the kept ones come first
                    group, nms = group[:n], nms[:n]
                    if not group:
                        totals["over_mismatch_cap"] += 1
                        continue

                n_best = nms.count(nms[0])
                contaminant = [is_contaminant[a.reference_id] for a in group]
                if any(contaminant[:n_best]):
                    cls = ("tie" if not all(contaminant[:n_best]) else
                           "contaminant_better" if not all(contaminant) else "contaminant_only")
                elif n_best == MAX_ALIGNMENTS:
                    cls = "too_many_loci"
                else:
                    cls = "rescued" if any(contaminant) else "transcriptome"
                totals[cls] += 1

                if cls in CONTAMINANT_CLASSES:
                    first = contaminant.index(True)
                    n_contaminant_best = sum(c and nm == nms[first] for c, nm in zip(contaminant, nms))
                    if n_contaminant_best == 1:
                        contaminant_unique += 1
                    else:
                        contaminant_multi += 1
                    per_record[names[group[first].reference_id]][cls] += 1
                if cls not in ("rescued", "transcriptome"):
                    continue

                primary = zlib.crc32(name.encode()) % n_best
                for i, a in enumerate(group[:n_best]):
                    a.reference_id = out_tid[a.reference_id]
                    a.flag = 0 if i == primary else 256
                    a.mapping_quality = star_mapq(n_best)
                    a.set_tags([("NH", n_best), ("NM", nms[i]), ("MD", a.get_tag("MD"))])
                    a.query_qualities = None
                    out.write(a)

    input_reads = sum(totals.values())
    with open(stats_out, "w") as o:
        o.write("\t".join(["sample", "input_reads", *CLASSES, "contaminant_unique", "contaminant_multi"]) + "\n")
        o.write("\t".join(map(str, [sample, input_reads, *(totals[c] for c in CLASSES), contaminant_unique,
                                    contaminant_multi])) + "\n")
    with open(records_out, "w") as o:
        o.write("record\treads\t" + "\t".join(CONTAMINANT_CLASSES) + "\n")
        for record, c in sorted(per_record.items(), key=lambda kv: -sum(kv[1].values())):
            o.write(f"{record}\t{sum(c.values())}\t" + "\t".join(str(c[k]) for k in CONTAMINANT_CLASSES) + "\n")
    print(f"input_reads\t{input_reads}", file=sys.stderr)
    for cls in CLASSES:
        print(f"{cls}\t{totals[cls]}", file=sys.stderr)
    print(f"contaminant_unique\t{contaminant_unique}\ncontaminant_multi\t{contaminant_multi}", file=sys.stderr)


if __name__ == "__main__":
    main()
