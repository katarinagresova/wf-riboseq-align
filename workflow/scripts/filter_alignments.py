"""Drop alignments from a transcriptome BAM and fix the reads that lose some.

  sense      drop the antisense alignments (keep_sense). The libraries are stranded: a footprint
             has its mRNA's sequence, so it aligns forward to its transcript, and an antisense
             (reverse) alignment is not a footprint of that transcript. STAR cannot align
             single-end reads to one strand only, so they are dropped here.
  refs BED   drop the alignments on references not in BED (split_bam_transcriptome: one species'
             transcripts), so that each species' BAM stands on its own: a read in both has NH =
             its alignments in that BAM, and a primary in each.

A read losing all its alignments is dropped. A read keeping some keeps those, with what the
dropped ones leave stale fixed:
  NH       the number of its kept alignments
  MAPQ     STAR's value for that number: 255 for 1, 3 for 2, 1 for 3-4, 0 for more
  primary  if its primary was dropped, one of the kept alignments, picked by the CRC32 of the
           read name (STAR picks at random among equally good alignments; the first one in the
           BAM would favour low reference ids, i.e. human over spike-in)
Every other record is written unchanged.

Two passes over the BAM, so memory holds only the reads with NH > 1 that lose alignments, and the
output keeps the input's order (a sorted BAM stays sorted). Fails unless each of those reads has
NH minus its dropped alignments kept ones (STAR writes all NH alignments of a read).

Usage: filter_alignments.py sense <in.bam> <out.bam>
       filter_alignments.py refs <refs.bed> <in.bam> <out.bam>      (counts to stderr)
"""
import sys
import zlib

import pysam


def star_mapq(nh):
    return 255 if nh == 1 else 3 if nh == 2 else 1 if nh <= 4 else 0


def main():
    args = sys.argv[1:]
    if args[:1] == ["sense"] and len(args) == 3:
        bed, (bam_in, bam_out) = None, args[1:]
    elif args[:1] == ["refs"] and len(args) == 4:
        bed, bam_in, bam_out = args[1:]
    else:
        sys.exit(__doc__)

    if bed is None:
        def drop(a):
            return a.is_reverse
    else:
        with open(bed) as f, pysam.AlignmentFile(bam_in) as bam:
            keep_ids = {bam.get_tid(line.split("\t", 1)[0]) for line in f}
        if -1 in keep_ids:
            sys.exit(f"{bed} has a reference that is not in {bam_in}")

        def drop(a):
            return a.reference_id not in keep_ids

    n_alignments = 0  # dropped alignments
    n_unique = 0  # reads with NH 1 whose alignment is dropped: nothing to fix
    dropped = {}  # read with NH > 1 -> its number of dropped alignments
    primary_dropped = set()
    with pysam.AlignmentFile(bam_in) as f:
        for a in f:
            if drop(a):
                n_alignments += 1
                if a.get_tag("NH") == 1:
                    n_unique += 1
                    continue
                dropped[a.query_name] = dropped.get(a.query_name, 0) + 1
                if not a.is_secondary:
                    primary_dropped.add(a.query_name)

    kept = {}  # read with dropped alignments -> [its NH after the filter, kept alignments seen]
    with pysam.AlignmentFile(bam_in) as f, pysam.AlignmentFile(bam_out, "wb", template=f) as out:
        for a in f:
            if drop(a):
                continue
            name = a.query_name
            if name in dropped:
                if name not in kept:
                    kept[name] = [a.get_tag("NH") - dropped[name], 0]
                    if kept[name][0] < 1:
                        sys.exit(f"{name}: NH {a.get_tag('NH')}, but {dropped[name]} dropped alignments")
                nh, seen = kept[name]
                kept[name][1] += 1
                a.set_tags([(tag, nh if tag == "NH" else value) for tag, value in a.get_tags()])  # keeps STAR's order
                a.mapping_quality = star_mapq(nh)
                if name in primary_dropped and seen == zlib.crc32(name.encode()) % nh:
                    a.is_secondary = False
            out.write(a)

    for name, (nh, seen) in kept.items():
        if seen != nh:
            sys.exit(f"{name}: {seen} kept alignments, but NH minus its dropped alignments is {nh}")

    print(f"dropped_alignments\t{n_alignments}", file=sys.stderr)
    print(f"reads_with_dropped\t{n_unique + len(dropped)}", file=sys.stderr)
    print(f"reads_dropped\t{n_unique + len(dropped) - len(kept)}", file=sys.stderr)
    print(f"primary_moved\t{len(primary_dropped & kept.keys())}", file=sys.stderr)


if __name__ == "__main__":
    main()
