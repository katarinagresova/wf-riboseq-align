"""Collapse identical reads into one record carrying their count.

Port of collapse_reads.pl from the eIF pipeline. fastq on stdin -> fastq on
stdout, one record per distinct sequence named @<library>_<rank>_x<count>,
most abundant first (ties in order of first appearance), with the quality
string of the first occurrence. Read and composition statistics go to stderr.

Usage: collapse_reads.py <library> < in.fastq > out.fastq
"""
import sys

import numpy as np

EDGE = 20
NTS = b"ACGT"


def composition(seqs, from_3prime):
    """Per-position A/C/G/T fraction over the first EDGE nt from one end."""
    if from_3prime:
        block = b"".join(s[::-1][:EDGE].ljust(EDGE, b".") for s in seqs)
    else:
        block = b"".join(s[:EDGE].ljust(EDGE, b".") for s in seqs)
    arr = np.frombuffer(block, dtype=np.uint8).reshape(len(seqs), EDGE)
    # +1 pseudocount per nucleotide, as in the original
    total = (arr != ord(".")).sum(axis=0) + len(NTS)
    return {chr(nt): ((arr == nt).sum(axis=0) + 1) / total for nt in NTS}


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    library = sys.argv[1]

    reads = {}
    lines = (line.rstrip(b"\n") for line in sys.stdin.buffer)
    lines = (line for line in lines if line.strip())
    for _name, seq, _plus, qual in zip(lines, lines, lines, lines):
        seq = seq.upper()
        entry = reads.get(seq)
        if entry is None:
            reads[seq] = [1, qual]
        else:
            entry[0] += 1

    ranked = sorted(reads.items(), key=lambda item: -item[1][0])
    out = sys.stdout.buffer
    prefix = f"@{library}_".encode()
    for rank, (seq, (count, qual)) in enumerate(ranked):
        out.write(b"%s%d_x%d\n%s\n+\n%s\n" % (prefix, rank, count, seq, qual))

    raw_len = dict.fromkeys(range(15, 77), 0)
    uniq_len = dict.fromkeys(range(15, 77), 0)
    for seq, (count, _qual) in reads.items():
        raw_len[len(seq)] = raw_len.get(len(seq), 0) + count
        uniq_len[len(seq)] = uniq_len.get(len(seq), 0) + 1

    seqs = list(reads)
    err = sys.stderr
    err.write(f">read statistics\n#property\tcounts_{library}\n"
              f"input\t{sum(raw_len.values())}\nuniq\t{len(reads)}\n\n")
    for title, sign, from_3prime in (("5'", "", False), ("3'", "-", True)):
        if title == "3'":
            err.write("\n")
        freq = composition(seqs, from_3prime)
        err.write(f">nt satistics from {title}end\n## --col_stack=A,C,G,T\n#pos\tA\tC\tG\tT\n")
        for i in range(EDGE):
            err.write(f"{sign}{i}" + "".join("\t%.3f" % freq[nt][i] for nt in "ACGT") + "\n")
    for title, lengths in (("raw", raw_len), ("unique", uniq_len)):
        err.write(f"\n>read_lengths_{title}\n## --cumsum --steps\n#length\t{library}\n")
        for length in sorted(lengths):
            err.write(f"{length}\t{lengths[length]}\n")


if __name__ == "__main__":
    main()
