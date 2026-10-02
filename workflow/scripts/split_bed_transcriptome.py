"""Write the human / spike-in BED files: one line per transcript of each fasta, over its whole length.

split_bam_transcriptome keeps the reads on one BED's transcripts. A name in both fastas is an error
(bowtie_index checks that too): its reads could not be told apart.

Usage: split_bed_transcriptome.py <human_fa> <spike_in_fa> <human_bed> <spike_in_bed>
"""
import sys


def fasta_lengths(path):
    lengths = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in lengths:
                    sys.exit(f"{path}: transcript name {name} twice")
                lengths[name] = 0
            else:
                lengths[name] += len(line.strip())
    return lengths


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    human_fa, spike_in_fa, human_bed, spike_in_bed = sys.argv[1:]

    human = fasta_lengths(human_fa)
    spike_in = fasta_lengths(spike_in_fa)
    dup = human.keys() & spike_in.keys()
    if dup:
        sys.exit(f"{len(dup)} transcript name(s) in both {human_fa} and {spike_in_fa}, e.g. {next(iter(dup))}")

    for lengths, bed in ((human, human_bed), (spike_in, spike_in_bed)):
        with open(bed, "w") as out:
            for name, length in lengths.items():
                out.write(f"{name}\t0\t{length}\n")


if __name__ == "__main__":
    main()
