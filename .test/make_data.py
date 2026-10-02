"""Write the synthetic test data in .test/data/ (run from anywhere; seeded, so the output is always the same).

Everything is random sequence, not real data, sized to run the whole workflow in a few minutes:
  human.fa, human.gtf   20 genes, half of them with an exon-skipping second isoform (reads multimap across
                        isoforms); 4 genes get no RNA-seq reads, so mode "filtered" blacklists them
  genome.fa             the genome of human.gtf: the genes' exons with introns between them (salmon's decoys)
  yeast.fa              8 spike-in transcripts; one contains a stretch of a human transcript (its reads are in
                        both species' BAMs: both_species)
  contaminants.fa       rRNA / tRNA-like records, plus: an exact duplicate (number_contaminants drops it), a
                        stretch of rRNA that is also inside a human transcript (its reads tie: stay removed),
                        a copy of a human stretch with a mismatch every 12 nt (its reads are rescued), and a
                        stretch of rRNA as its own record (its reads align to both: contaminant_multi)
  ribo_{a,b}.fastq.gz   single-end 51 nt: 4 nt randomer + insert + 4 nt randomer + adapter, with PCR duplicates,
                        too-short inserts and adapter dimers
  rna_{a,b}_R{1,2}.fastq.gz  paired-end 2x75 nt, stranded (R1 antisense), plus pre-mRNA fragments (introns
                        included: salmon assigns most of them to a decoy)
"""
import gzip
import os
import random

ADAPTER = "TGGAATTCTCGGGTGCCAAGG"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
rng = random.Random(1)
# the genome and the pre-mRNA reads: a stream of their own, so the rest is the same as before they were added
genome_rng = random.Random(2)


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

    # genome: 1 kb between genes, 200-600 nt introns. pre_mrna: each gene's unspliced sequence, by transcript stem
    human, gtf, genome, pre_mrna = {}, [], [rseq(1000, genome_rng)], {}
    for g in range(1, 21):
        exons = [rseq(rng.randint(120, 400)) for _ in range(rng.randint(3, 5))]
        if g == 1:  # the tie stretch: rRNA sequence inside a human exon
            exons[1] = exons[1][:50] + rrna["rRNA_28S"][1000:1100] + exons[1][50:]
        gene_start = sum(map(len, genome))
        coords = []  # 1-based, per exon
        for k, exon in enumerate(exons):
            if k:
                genome.append(rseq(genome_rng.randint(200, 600), genome_rng))
            start = sum(map(len, genome)) + 1
            coords.append((start, start + len(exon) - 1))
            genome.append(exon)
        pre_mrna[f"HTX{g:03d}"] = "".join(genome)[gene_start:]
        genome.append(rseq(1000, genome_rng))
        isoforms = [list(range(len(exons)))] + ([[0] + list(range(2, len(exons)))] if g % 2 == 0 else [])
        for i, iso in enumerate(isoforms, 1):
            tx, gene = f"HTX{g:03d}.{i}", f"HG{g:03d}"
            human[tx] = "".join(exons[k] for k in iso)
            attrs = f'gene_id "{gene}"; transcript_id "{tx}";'
            gtf.append(f"chr1\ttest\ttranscript\t{coords[iso[0]][0]}\t{coords[iso[-1]][1]}\t.\t+\t.\t{attrs}\n")
            gtf += [f"chr1\ttest\texon\t{coords[k][0]}\t{coords[k][1]}\t.\t+\t.\t{attrs}\n" for k in iso]

    # the rescue stretch: a human stretch, as a contaminant with a mismatch every 12 nt
    rescue = list(human["HTX003.1"][200:320])
    for i in range(6, len(rescue), 12):
        rescue[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[rescue[i]]
    contaminants = {**rrna, **trna, "tRNA_1_copy": trna["tRNA_1"], "snoRNA_like_HTX003": "".join(rescue),
                    "rRNA_28S_fragment": rrna["rRNA_28S"][2000:2100]}

    yeast = {f"Y{c}{i:02d}W": rseq(rng.randint(800, 1500)) for i, c in enumerate("ABCDEFGH", 1)}
    # the both-species stretch: a human stretch inside a spike-in transcript
    yeast["YA01W"] = yeast["YA01W"][:300] + human["HTX007.1"][100:200] + yeast["YA01W"][300:]
    return human, "".join(gtf), {"chr1": "".join(genome)}, pre_mrna, yeast, contaminants


def insert_from(seq, lo=26, hi=34):
    n = min(rng.randint(lo, hi), len(seq))
    p = rng.randint(0, len(seq) - n)
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


def rna_reads(human, n_pairs, silent):
    tx_names = [t for t in human if t.split(".")[0] not in silent]
    weights = [rng.lognormvariate(0, 1) for _ in tx_names]
    # a generator: each draw happens just before its fragment's, as before pairs() existed
    return pairs((human[rng.choices(tx_names, weights)[0]] for _ in range(n_pairs)), rng)


def pre_mrna_reads(pre_mrna, n_pairs, silent):
    genes = [g for g in pre_mrna if g not in silent]
    return pairs((pre_mrna[genome_rng.choice(genes)] for _ in range(n_pairs)), genome_rng)


def pairs(seqs, r):
    out = []
    for seq in seqs:
        n = r.randint(150, min(300, len(seq)))
        p = r.randint(0, len(seq) - n)
        fragment = seq[p:p + n]
        out.append((revcomp(fragment)[:75], fragment[:75]))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    human, gtf, genome, pre_mrna, yeast, contaminants = references()
    write_fasta(f"{OUT}/human.fa", human)
    with open(f"{OUT}/human.gtf", "w") as f:
        f.write(gtf)
    write_fasta(f"{OUT}/genome.fa", genome)
    write_fasta(f"{OUT}/yeast.fa", yeast)
    write_fasta(f"{OUT}/contaminants.fa", contaminants)

    silent = {"HTX005", "HTX009", "HTX013", "HTX017"}  # single-isoform genes without RNA-seq reads
    for lib in "ab":
        write_fastq(f"{OUT}/ribo_{lib}.fastq.gz", ribo_reads(human, yeast, contaminants, 1500))
        rna = rna_reads(human, 2000, silent) + pre_mrna_reads(pre_mrna, 300, silent)
        write_fastq(f"{OUT}/rna_{lib}_R1.fastq.gz", [(f"f{i}/1", a) for i, (a, _) in enumerate(rna)])
        write_fastq(f"{OUT}/rna_{lib}_R2.fastq.gz", [(f"f{i}/2", b) for i, (_, b) in enumerate(rna)])


if __name__ == "__main__":
    main()
