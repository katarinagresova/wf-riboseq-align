"""Drop the antisense alignments from star_transcript's BAM.

The libraries are stranded: a footprint has its mRNA's sequence, so it aligns forward to its
transcript, and an antisense (reverse) alignment is not a footprint of that transcript. STAR
cannot align single-end reads to one strand only, so they are dropped here. A read with only
antisense alignments is dropped. A read with sense alignments too keeps those, with what the
dropped ones leave stale fixed:
  NH       the number of its sense alignments
  MAPQ     STAR's value for that number: 255 for 1, 3 for 2, 1 for 3-4, 0 for more
  primary  if STAR's primary was antisense, one of the sense alignments, picked by the CRC32 of
           the read name (STAR picks at random among equally good alignments; the first one in
           the BAM would favour low reference ids, i.e. human over spike-in)
Every other record is written unchanged.

Two passes over the coordinate-sorted BAM, so memory holds only the reads with antisense
alignments, and the output keeps the input's order. Fails unless each of those reads has NH minus
its antisense alignments sense alignments (STAR writes all NH alignments of a read).

Usage: keep_sense.py <in.bam> <out.bam>       (counts to stderr)
"""
import sys
import zlib

import pysam


def star_mapq(nh):
    return 255 if nh == 1 else 3 if nh == 2 else 1 if nh <= 4 else 0


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    bam_in, bam_out = sys.argv[1:]

    antisense = {}  # read -> its number of antisense alignments
    antisense_primary = set()
    with pysam.AlignmentFile(bam_in) as f:
        for a in f:
            if a.is_reverse:
                antisense[a.query_name] = antisense.get(a.query_name, 0) + 1
                if not a.is_secondary:
                    antisense_primary.add(a.query_name)

    sense = {}  # read with antisense alignments -> [its NH after the filter, sense alignments seen]
    with pysam.AlignmentFile(bam_in) as f, pysam.AlignmentFile(bam_out, "wb", template=f) as out:
        for a in f:
            if a.is_reverse:
                continue
            name = a.query_name
            if name in antisense:
                if name not in sense:
                    sense[name] = [a.get_tag("NH") - antisense[name], 0]
                    if sense[name][0] < 1:
                        sys.exit(f"{name}: NH {a.get_tag('NH')}, but {antisense[name]} antisense alignments")
                nh, seen = sense[name]
                sense[name][1] += 1
                a.set_tags([(tag, nh if tag == "NH" else value) for tag, value in a.get_tags()])  # keeps STAR's order
                a.mapping_quality = star_mapq(nh)
                if name in antisense_primary and seen == zlib.crc32(name.encode()) % nh:
                    a.is_secondary = False
            out.write(a)

    for name, (nh, seen) in sense.items():
        if seen != nh:
            sys.exit(f"{name}: {seen} sense alignments, but NH minus its antisense alignments is {nh}")

    print(f"antisense_alignments\t{sum(antisense.values())}", file=sys.stderr)
    print(f"reads_with_antisense\t{len(antisense)}", file=sys.stderr)
    print(f"antisense_only_reads\t{len(antisense) - len(sense)}", file=sys.stderr)
    print(f"primary_moved\t{len(antisense_primary & sense.keys())}", file=sys.stderr)


if __name__ == "__main__":
    main()
