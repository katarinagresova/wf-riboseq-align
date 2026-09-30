"""Concatenate the per-library QC tables, then fail if any library's spike-in
fraction (spike-in-only reads / all mapped reads) is below the minimum.

Usage: qc_summary.py <summary_out> <min_spike_in_fraction> <sample.stats.tsv>...
"""
import sys

import pandas as pd


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    summary_out, min_fraction, *stats_files = sys.argv[1:]
    min_fraction = float(min_fraction)

    summary = pd.concat([pd.read_csv(f, sep="\t") for f in stats_files])
    summary.to_csv(summary_out, sep="\t", index=False)

    low = summary[summary["spike_in_fraction"] < min_fraction]
    if len(low):
        sys.exit(f"spike-in fraction below min_spike_in_fraction={min_fraction:g}:\n"
                 + low[["sample", "spike_in_only", "spike_in_fraction"]].to_string(index=False)
                 + "\nSet min_spike_in_fraction: 0 for libraries without a spike-in.")


if __name__ == "__main__":
    main()
