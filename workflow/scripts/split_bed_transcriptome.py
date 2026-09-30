"""Split the transcriptome STAR index into human / spike-in BED files by name.

Replaces the count-and-position split in split_bed_transcriptome (align.smk),
which assumed spike-in was appended last to the index. Here each index
reference is looked up by name in the human and spike-in fastas instead, so a
fasta with the right count but wrong content can't pass silently.

Usage: split_bed_transcriptome.py <chrom_sizes> <human_fa> <spike_in_fa> <human_bed> <spike_in_bed>
"""
import sys


def fasta_names(path):
    with open(path) as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    chrom_sizes, human_fa, spike_in_fa, human_bed, spike_in_bed = sys.argv[1:]

    human_set = set(fasta_names(human_fa))
    spike_in_set = set(fasta_names(spike_in_fa))

    dup = human_set & spike_in_set
    if dup:
        sys.exit(f"{len(dup)} transcript name(s) in both {human_fa} and {spike_in_fa}, "
                  f"e.g. {next(iter(dup))}")

    seen = set()
    with open(chrom_sizes) as src, \
            open(human_bed, "w") as human_out, \
            open(spike_in_bed, "w") as spike_in_out:
        for line in src:
            name, length = line.rstrip("\n").split("\t")
            if name in human_set:
                out = human_out
            elif name in spike_in_set:
                out = spike_in_out
            else:
                sys.exit(f"index reference {name!r} is in neither fasta")
            out.write(f"{name}\t0\t{length}\n")
            seen.add(name)

    missing = (human_set | spike_in_set) - seen
    if missing:
        sys.exit(f"{len(missing)} fasta transcript(s) missing from the index, "
                  f"e.g. {next(iter(missing))}")


if __name__ == "__main__":
    main()
