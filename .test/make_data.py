"""Write the synthetic test data in .test/data/ (run from anywhere; seeded, so the output is always the same).

Everything is random sequence, not real data, sized to run the whole workflow in a few minutes:
  human.fa              20 genes, half of them with an exon-skipping second isoform (reads multimap across
                        isoforms), plus HTX021.1, which holds the reverse complement of a HTX011.1 stretch (its
                        reads align to HTX011.1 only: bowtie_align aligns sense only)
  yeast.fa              8 spike-in transcripts; one contains a stretch of a human transcript (its reads align
                        to both species: in neither BAM, counted as both_species)
  contaminants.fa       rRNA / tRNA-like records, plus: a stretch of rRNA that is also inside a human transcript
                        (its reads tie: stay removed),
                        a copy of a human stretch with a mismatch every 12 nt (its reads are rescued), a
                        stretch of rRNA as its own record (its reads align to both: contaminant_multi), the
                        reverse complement of a HTX015.1 stretch (its reads are antisense to it only: plain
                        transcriptome reads) and the reverse complement of the tie stretch (the tie reads hit
                        rRNA_28S only)
  ribo_{a,b}.fastq.gz   single-end 51 nt: 4 nt randomer + insert + 4 nt randomer + adapter, with PCR duplicates,
                        too-short inserts and adapter dimers. At the end: reads antisense to a human transcript
                        (unaligned), reads on the HTX011.1 stretch, reads on the HTX015.1 stretch and 20-22 nt
                        reads with 2 mismatches to a human transcript (over the short-read mismatch cap)
  rna_{a,b}_R{1,2}.fastq.gz  paired-end 2x75 nt, stranded (R1 antisense), none on 4 genes. Only for samples.csv's
                        read_type=rna rows, which the workflow skips; they share the Ribo-seq reads' random stream,
                        so dropping them would change those
"""
import gzip
import os
import random

ADAPTER = "TGGAATTCTCGGGTGCCAAGG"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
rng = random.Random(1)
# the antisense cases (HTX021.1 and the reads at the end of ribo_*) have their own stream, so adding them left
# every other read and record unchanged
antisense_rng = random.Random(2)
# so do the reads on the HTX015.1 stretch, for the competitive filter's strand rule
strand_rng = random.Random(3)
# and the short reads over the mismatch cap
cap_rng = random.Random(4)


def rseq(n, r=rng):
    return "".join(r.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n")
            f.writelines(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))


def write_fastq(path, reads):
    # mtime=0: the same bytes on every run
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as f:
        for name, seq in reads:
            f.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n".encode())


def references():
    rrna = {"rRNA_28S": rseq(3000), "rRNA_18S": rseq(1800), "rRNA_5_8S": rseq(160)}
    trna = {f"tRNA_{i}": rseq(73) for i in range(1, 6)}

    human = {}
    for g in range(1, 21):
        exons = [rseq(rng.randint(120, 400)) for _ in range(rng.randint(3, 5))]
        if g == 1:  # the tie stretch: rRNA sequence inside a human exon
            exons[1] = exons[1][:50] + rrna["rRNA_28S"][1000:1100] + exons[1][50:]
        isoforms = [exons] + ([exons[:1] + exons[2:]] if g % 2 == 0 else [])
        for i, iso in enumerate(isoforms, 1):
            human[f"HTX{g:03d}.{i}"] = "".join(iso)

    # the rescue stretch: a human stretch, as a contaminant with a mismatch every 12 nt
    rescue = list(human["HTX003.1"][200:320])
    for i in range(6, len(rescue), 12):
        rescue[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[rescue[i]]
    contaminants = {**rrna, **trna, "tRNA_1_copy": trna["tRNA_1"], "snoRNA_like_HTX003": "".join(rescue),
                    "rRNA_28S_fragment": rrna["rRNA_28S"][2000:2100]}

    yeast = {f"Y{c}{i:02d}W": rseq(rng.randint(800, 1500)) for i, c in enumerate("ABCDEFGH", 1)}
    # the both-species stretch: a human stretch inside a spike-in transcript
    yeast["YA01W"] = yeast["YA01W"][:300] + human["HTX007.1"][100:200] + yeast["YA01W"][300:]
    return human, yeast, contaminants


def insert_from(seq, lo=26, hi=34, r=rng):
    n = min(r.randint(lo, hi), len(seq))
    p = r.randint(0, len(seq) - n)
    return seq[p:p + n]


def ribo_reads(human, yeast, contaminants, n_molecules):
    tx_names = list(human)
    weights = [rng.lognormvariate(0, 1) for _ in tx_names]
    kinds = ["human", "spike_in", "contaminant", "tie", "rescue", "contaminant_multi", "both_species", "short",
             "dimer", "random"]
    kind_weights = [60, 8, 20, 3, 3, 3, 3, 3, 2, 1]
    contam = [s for name, s in contaminants.items() if name != "snoRNA_like_HTX003"]
    tie = human["HTX001.1"]
    tie = tie[tie.index(contaminants["rRNA_28S"][1000:1100]):][:100]
    rescue = human["HTX003.1"][200:320]
    both_species = human["HTX007.1"][100:200]

    reads = []
    for _ in range(n_molecules):
        kind = rng.choices(kinds, kind_weights)[0]
        if kind == "human":
            insert = insert_from(human[rng.choices(tx_names, weights)[0]])
        elif kind == "spike_in":
            insert = insert_from(rng.choice(list(yeast.values())))
        elif kind == "contaminant":
            insert = insert_from(rng.choice(contam))
        elif kind == "tie":
            insert = insert_from(tie)
        elif kind == "rescue":
            insert = insert_from(rescue)
        elif kind == "contaminant_multi":
            insert = insert_from(contaminants["rRNA_28S_fragment"])
        elif kind == "both_species":
            insert = insert_from(both_species)
        elif kind == "short":
            insert = insert_from(rng.choice(list(human.values())), 4, 9)
        elif kind == "dimer":
            insert = ""
        else:
            insert = rseq(30)
        read = (rseq(4) + insert + rseq(4) + ADAPTER + rseq(20))[:51]
        reads += [read] * rng.choice([1, 1, 1, 1, 2, 2, 3, 5])  # PCR duplicates
    rng.shuffle(reads)
    return [(f"r{i}", r) for i, r in enumerate(reads)]


def antisense_transcript(human):
    # the reverse complement of a HTX011.1 stretch, between random sequence
    r = antisense_rng
    return {"HTX021.1": rseq(150, r) + revcomp(human["HTX011.1"][100:200]) + rseq(150, r)}


def antisense_reads(human, n_molecules, first):
    # alternately: antisense to a human transcript, and on the HTX011.1 stretch that HTX021.1 holds antisense
    r = antisense_rng
    tx_names = list(human)
    reads = []
    for i in range(n_molecules):
        if i % 2:
            insert = revcomp(insert_from(human[r.choice(tx_names)], r=r))
        else:
            insert = insert_from(human["HTX011.1"][100:200], r=r)
        reads.append((f"r{first + i}", (rseq(4, r) + insert + rseq(4, r) + ADAPTER + rseq(20, r))[:51]))
    return reads


def strand_contaminants(human, contaminants):
    # antisense to a human stretch (no sense contaminant alignment), and antisense to the tie stretch (which rRNA_28S
    # holds sense). Not in references(): ribo_reads draws its contaminant reads from that dict
    return {"snoRNA_antisense_HTX015": revcomp(human["HTX015.1"][100:200]),
            "rRNA_28S_antisense": revcomp(contaminants["rRNA_28S"][1000:1100])}


def strand_reads(human, n_molecules, first):
    # on the HTX015.1 stretch that snoRNA_antisense_HTX015 holds antisense
    r = strand_rng
    reads = []
    for i in range(n_molecules):
        insert = insert_from(human["HTX015.1"][100:200], r=r)
        reads.append((f"r{first + i}", (rseq(4, r) + insert + rseq(4, r) + ADAPTER + rseq(20, r))[:51]))
    return reads


def over_cap_reads(human, n_molecules, first):
    # 20-22 nt with 2 mismatches to a human transcript: below SHORT_READ_LENGTH, 1 mismatch at most
    r = cap_rng
    tx_names = list(human)
    reads = []
    for i in range(n_molecules):
        insert = list(insert_from(human[r.choice(tx_names)], 20, 22, r=r))
        for p in (6, 13):
            insert[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[insert[p]]
        insert = "".join(insert)
        reads.append((f"r{first + i}", (rseq(4, r) + insert + rseq(4, r) + ADAPTER + rseq(20, r))[:51]))
    return reads


def rna_reads(human, n_pairs, silent):
    tx_names = [t for t in human if t.split(".")[0] not in silent]
    weights = [rng.lognormvariate(0, 1) for _ in tx_names]
    r1, r2 = [], []
    for i in range(n_pairs):
        seq = human[rng.choices(tx_names, weights)[0]]
        n = rng.randint(150, min(300, len(seq)))
        p = rng.randint(0, len(seq) - n)
        fragment = seq[p:p + n]
        r1.append((f"f{i}/1", revcomp(fragment)[:75]))
        r2.append((f"f{i}/2", fragment[:75]))
    return r1, r2


def main():
    os.makedirs(OUT, exist_ok=True)
    human, yeast, contaminants = references()
    write_fasta(f"{OUT}/human.fa", {**human, **antisense_transcript(human)})
    write_fasta(f"{OUT}/yeast.fa", yeast)
    # without tRNA_1_copy, which only doubles tRNA_1's share of the contaminant reads (it was a duplicate record for
    # number_contaminants to drop; contaminants_fa now goes to bowtie_index as is)
    write_fasta(f"{OUT}/contaminants.fa", {**{n: s for n, s in contaminants.items() if n != "tRNA_1_copy"},
                                          **strand_contaminants(human, contaminants)})

    silent = {"HTX005", "HTX009", "HTX013", "HTX017"}  # single-isoform genes without RNA-seq reads
    for lib in "ab":
        ribo = ribo_reads(human, yeast, contaminants, 1500)
        ribo += antisense_reads(human, 60, len(ribo))
        ribo += strand_reads(human, 30, len(ribo))
        write_fastq(f"{OUT}/ribo_{lib}.fastq.gz", ribo + over_cap_reads(human, 30, len(ribo)))
        r1, r2 = rna_reads(human, 2000, silent)
        write_fastq(f"{OUT}/rna_{lib}_R1.fastq.gz", r1)
        write_fastq(f"{OUT}/rna_{lib}_R2.fastq.gz", r2)


if __name__ == "__main__":
    main()
