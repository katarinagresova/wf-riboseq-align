"""Deduplicate contaminant sequences and prefix each name with a running number.

Port of number_contaminants.R from the eIF pipeline: STAR needs unique reference
names, and identical sequences would only split reads between the copies.
Keeps the first occurrence of each sequence; writes 80-column fasta.
Characters outside the IUPAC DNA alphabet are dropped with a warning, which is
what Biostrings did silently (the eIF contaminants fasta has stray text in two
records, so this keeps the index identical to the R version's).

Usage: number_contaminants.py <in.fa> <out.fa>
"""
import sys

WIDTH = 80
DNA_ALPHABET = set("ACGTMRWSYKVHDBN-+.")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    in_fa, out_fa = sys.argv[1:]

    seen = set()
    with open(out_fa, "w") as out:
        for name, seq in read_fasta(in_fa):
            dropped = "".join(c for c in seq.upper() if c not in DNA_ALPHABET)
            if dropped:
                print(f"WARNING {name}: dropped non-DNA characters {dropped!r}", file=sys.stderr)
            seq = "".join(c for c in seq.upper() if c in DNA_ALPHABET)
            if seq in seen:
                continue
            seen.add(seq)
            out.write(f">{len(seen)}_{name}\n")
            for i in range(0, len(seq), WIDTH):
                out.write(seq[i:i + WIDTH] + "\n")


if __name__ == "__main__":
    main()
