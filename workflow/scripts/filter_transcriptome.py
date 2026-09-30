"""Drop blacklisted transcripts from a transcriptome fasta + gtf.

Port of make_resources_autofilter.sh (Frederick Korbel, eIF pipeline). Matches
transcript ids EXACTLY: the original's `grep -F` is a substring match, which
over-matches versioned ids (ENST...123.1 is inside ENST...123.10). GTF lines
without a transcript_id (gene lines, comments) are kept.

Usage: filter_transcriptome.py <fa_in> <gtf_in> <blacklist> <fa_out> <gtf_out>
"""
import re
import sys

TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    fa_in, gtf_in, blacklist, fa_out, gtf_out = sys.argv[1:]

    with open(blacklist) as f:
        remove = {line.strip() for line in f if line.strip()}

    kept = 0
    with open(fa_in) as src, open(fa_out, "w") as out:
        keep = True
        for line in src:
            if line.startswith(">"):
                keep = line[1:].split()[0] not in remove
                kept += keep
            if keep:
                out.write(line)

    with open(gtf_in) as src, open(gtf_out, "w") as out:
        for line in src:
            match = TRANSCRIPT_ID.search(line)
            if match is None or match.group(1) not in remove:
                out.write(line)

    print(f"{len(remove)} transcripts blacklisted, {kept} kept in {fa_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
