"""Collapse identical reads into one record carrying their count.

Port of collapse_reads.pl from the eIF pipeline. fastq(.gz) -> fastq on stdout,
one record per distinct sequence named @<library>_<rank>_x<count>, most
abundant first (ties in order of first appearance), with the quality string of
the first occurrence. Read and composition statistics go to stderr.

Usage: collapse_reads.py <in.fastq.gz> <library> > out.fastq
"""
import sys
from collections import Counter

import pysam

EDGE = 20


def composition(seqs):
    """Per-position A/C/G/T fraction over the first EDGE nt of seqs."""
    # EDGE characters per read (shorter reads padded with "."), so position i
    # of every read is block[i::EDGE]
    block = "".join(s[:EDGE].ljust(EDGE, ".") for s in seqs)
    for i in range(EDGE):
        column = block[i::EDGE]
        # +1 pseudocount per nucleotide, as in the original
        total = len(column) - column.count(".") + 4
        yield [(column.count(nt) + 1) / total for nt in "ACGT"]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    fastq, library = sys.argv[1:]

    counts, qual = Counter(), {}
    for read in pysam.FastxFile(fastq):
        counts[read.sequence] += 1
        qual.setdefault(read.sequence, read.quality)

    # most_common() sorts stably: ties stay in order of first appearance
    for rank, (seq, count) in enumerate(counts.most_common()):
        sys.stdout.write(f"@{library}_{rank}_x{count}\n{seq}\n+\n{qual[seq]}\n")

    raw_len, uniq_len = Counter(), Counter()
    for seq, count in counts.items():
        raw_len[len(seq)] += count
        uniq_len[len(seq)] += 1

    err = sys.stderr
    err.write(f">read statistics\n#property\tcounts_{library}\n"
              f"input\t{counts.total()}\nuniq\t{len(counts)}\n")
    for end, sign, seqs in (("5'", "", counts), ("3'", "-", (s[::-1] for s in counts))):
        err.write(f"\n>nt satistics from {end}end\n## --col_stack=A,C,G,T\n#pos\tA\tC\tG\tT\n")
        for i, freq in enumerate(composition(seqs)):
            err.write(f"{sign}{i}" + "".join(f"\t{f:.3f}" for f in freq) + "\n")
    for title, lengths in (("raw", raw_len), ("unique", uniq_len)):
        err.write(f"\n>read_lengths_{title}\n## --cumsum --steps\n#length\t{library}\n")
        # lengths 15-76 always, others only if seen
        for length in sorted(set(range(15, 77)) | lengths.keys()):
            err.write(f"{length}\t{lengths[length]}\n")


if __name__ == "__main__":
    main()
