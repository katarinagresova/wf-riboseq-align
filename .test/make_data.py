"""Write the synthetic test data in .test/data/ (run from anywhere; seeded, so the output is always the same).

Everything is random sequence, not real data, sized to run the whole workflow in a few minutes:
  human.fa, human.gtf   20 genes, half of them with an exon-skipping second isoform (reads multimap across
                        isoforms); 4 genes get no RNA-seq reads, so mode "filtered" blacklists them
  yeast.fa              8 spike-in transcripts
  contaminants.fa       rRNA / tRNA-like records, plus: an exact duplicate (number_contaminants drops it), a
                        stretch of rRNA that is also inside a human transcript (its reads tie: stay removed), and
                        a copy of a human stretch with a mismatch every 12 nt (its reads are rescued)
  ribo_{a,b}.fastq.gz   single-end 51 nt: 4 nt randomer + insert + 4 nt randomer + adapter, with PCR duplicates,
                        too-short inserts and adapter dimers
  rna_{a,b}_R{1,2}.fastq.gz  paired-end 2x75 nt, stranded (R1 antisense)
"""
import gzip
import os
import random

ADAPTER = "TGGAATTCTCGGGTGCCAAGG"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
rng = random.Random(1)


def rseq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


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

    human, gtf, start = {}, [], 1
    for g in range(1, 21):
        exons = [rseq(rng.randint(120, 400)) for _ in range(rng.randint(3, 5))]
        if g == 1:  # the tie stretch: rRNA sequence inside a human exon
            exons[1] = exons[1][:50] + rrna["rRNA_28S"][1000:1100] + exons[1][50:]
        isoforms = [exons] + ([exons[:1] + exons[2:]] if g % 2 == 0 else [])
        for i, iso in enumerate(isoforms, 1):
            tx, gene = f"HTX{g:03d}.{i}", f"HG{g:03d}"
            human[tx] = "".join(iso)
            attrs = f'gene_id "{gene}"; transcript_id "{tx}";'
            gtf.append(f"chr1\ttest\ttranscript\t{start}\t{start + len(human[tx]) - 1}\t.\t+\t.\t{attrs}\n")
            pos = start
            for exon in iso:
                gtf.append(f"chr1\ttest\texon\t{pos}\t{pos + len(exon) - 1}\t.\t+\t.\t{attrs}\n")
                pos += len(exon)
        start += sum(len(e) for e in exons) + 1000

    # the rescue stretch: a human stretch, as a contaminant with a mismatch every 12 nt
    rescue = list(human["HTX003.1"][200:320])
    for i in range(6, len(rescue), 12):
        rescue[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[rescue[i]]
    contaminants = {**rrna, **trna, "tRNA_1_copy": trna["tRNA_1"], "snoRNA_like_HTX003": "".join(rescue)}

    yeast = {f"Y{c}{i:02d}W": rseq(rng.randint(800, 1500)) for i, c in enumerate("ABCDEFGH", 1)}
    return human, "".join(gtf), yeast, contaminants


def insert_from(seq, lo=26, hi=34):
    n = min(rng.randint(lo, hi), len(seq))
    p = rng.randint(0, len(seq) - n)
    return seq[p:p + n]


def ribo_reads(human, yeast, contaminants, n_molecules):
    tx_names = list(human)
    weights = [rng.lognormvariate(0, 1) for _ in tx_names]
    kinds = ["human", "spike_in", "contaminant", "tie", "rescue", "short", "dimer", "random"]
    kind_weights = [60, 8, 20, 3, 3, 3, 2, 1]
    contam = [s for name, s in contaminants.items() if name != "snoRNA_like_HTX003"]
    tie = human["HTX001.1"]
    tie = tie[tie.index(contaminants["rRNA_28S"][1000:1100]):][:100]
    rescue = human["HTX003.1"][200:320]

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
    human, gtf, yeast, contaminants = references()
    write_fasta(f"{OUT}/human.fa", human)
    with open(f"{OUT}/human.gtf", "w") as f:
        f.write(gtf)
    write_fasta(f"{OUT}/yeast.fa", yeast)
    write_fasta(f"{OUT}/contaminants.fa", contaminants)

    silent = {"HTX005", "HTX009", "HTX013", "HTX017"}  # single-isoform genes without RNA-seq reads
    for lib in "ab":
        write_fastq(f"{OUT}/ribo_{lib}.fastq.gz", ribo_reads(human, yeast, contaminants, 1500))
        r1, r2 = rna_reads(human, 2000, silent)
        write_fastq(f"{OUT}/rna_{lib}_R1.fastq.gz", r1)
        write_fastq(f"{OUT}/rna_{lib}_R2.fastq.gz", r2)


if __name__ == "__main__":
    main()
