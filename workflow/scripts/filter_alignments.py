"""Keep only the reads whose alignments are all on references in BED (split_bam_transcriptome: one
species' transcripts). A read aligned to both species is in neither species' BAM, as it may come
from either. The kept reads keep all their alignments, so nothing needs fixing: every kept record is
written unchanged.

Two passes over the BAM, so memory holds only the reads with NH > 1 that have an alignment off the
references, and the output keeps the input's order (a sorted BAM stays sorted).

Usage: filter_alignments.py <refs.bed> <in.bam> <out.bam>      (counts to stderr)
"""
import sys

import pysam


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    bed, bam_in, bam_out = sys.argv[1:]

    with open(bed) as f, pysam.AlignmentFile(bam_in) as bam:
        keep_ids = {bam.get_tid(line.split("\t", 1)[0]) for line in f}
    if -1 in keep_ids:
        sys.exit(f"{bed} has a reference that is not in {bam_in}")

    n_unique = 0  # reads with NH 1 whose alignment is off the references
    dropped = set()  # reads with NH > 1 and an alignment off the references
    with pysam.AlignmentFile(bam_in) as f:
        for a in f:
            if a.reference_id not in keep_ids:
                if a.get_tag("NH") == 1:
                    n_unique += 1
                else:
                    dropped.add(a.query_name)

    both = set()  # of them, reads with alignments on the references too
    with pysam.AlignmentFile(bam_in) as f, pysam.AlignmentFile(bam_out, "wb", template=f) as out:
        for a in f:
            if a.reference_id not in keep_ids:
                continue
            if a.query_name in dropped:
                both.add(a.query_name)
                continue
            out.write(a)

    print(f"reads_dropped\t{n_unique + len(dropped)}", file=sys.stderr)  # with an alignment off the references
    print(f"both_species\t{len(both)}", file=sys.stderr)  # of them, with alignments on them too


if __name__ == "__main__":
    main()
