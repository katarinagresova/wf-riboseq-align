"""Build contaminants.fa from public sources instead of a hand-maintained file.

For human and the yeast spike-in, pulls every ncRNA biotype except lncRNA
(rRNA, tRNA, snRNA, snoRNA, scaRNA, miRNA, misc_RNA, ...) from Ensembl. Ensembl
does not annotate human nuclear tRNA genes, so those come from GtRNAdb
instead; and its human rRNA biotype entries are fragmentary pseudogene copies,
so the four canonical full-length rRNAs are added from NCBI RefSeq. Also
carries over 12 Illumina TruSeq/RPI index-barcode sequences and 2 synthetic
RNA size-marker oligos: wet-lab protocol sequences that aren't in any
database, copied verbatim from the lab's existing contaminants.combined_human_yeast.fa.

Usage: build_contaminants.py <out.fa>
"""
import gzip
import sys
import urllib.request

ENSEMBL_RELEASE = 109  # matches spike_in_transcriptome_fa elsewhere in this pipeline
HUMAN_NCRNA_URL = (
    f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/fasta/homo_sapiens/"
    "ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
)
YEAST_NCRNA_URL = (
    f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/fasta/saccharomyces_cerevisiae/"
    "ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz"
)
GTRNADB_HUMAN_URL = "https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accs}&rettype=fasta&retmode=text"
NCBI_RRNA_ACCESSIONS = ["NR_003286.2", "NR_003287.2", "NR_003285.2", "NR_023379.1"]  # 18S, 28S, 5.8S, 5S

EXCLUDE_BIOTYPES = {"lncRNA"}

# Not from any database: Illumina TruSeq small-RNA multiplexing index primers
# (RPI1-12) and synthetic RNA size-marker oligos, run alongside libraries in
# the wet-lab protocol. Copied verbatim from contaminants.combined_human_yeast.fa.
LITERAL_SEQUENCES = {
    "human_RPI1": "ATCACGA",
    "human_RPI2": "CGATGTA",
    "human_RPI3": "TTAGGCA",
    "human_RPI4": "TGACCAA",
    "human_RPI5": "ACAGTGA",
    "human_RPI6": "GCCAATA",
    "human_RPI7": "CAGATCA",
    "human_RPI8": "ACTTGAA",
    "human_RPI9": "GATCAGA",
    "human_RPI10": "TAGCTTA",
    "human_RPI11": "GGCTACA",
    "human_RPI12": "CTTGTAA",
    "human_marker_27nt": "ATGTACACGGAGTCGAGCTCAACCCGC",
    "human_marker_30nt": "ATGTACACGGAGTCGAGCTCAACCCGCAAC",
}


def fetch(url):
    # GtRNAdb rejects urllib's default User-Agent with 403.
    request = urllib.request.Request(url, headers={"User-Agent": "curl/8.0"})
    with urllib.request.urlopen(request) as resp:
        data = resp.read()
    if url.endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode()


def parse_fasta(text):
    name, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = line[1:], []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def header_tag(header, tag):
    for token in header.split():
        if token.startswith(tag + ":"):
            return token[len(tag) + 1:]
    return None


def ensembl_ncrna(url, prefix):
    for header, seq in parse_fasta(fetch(url)):
        if header_tag(header, "gene_biotype") in EXCLUDE_BIOTYPES:
            continue
        name = header_tag(header, "gene_symbol") or header.split()[0]
        yield f"{prefix}_{name}", seq


def gtrnadb(url, prefix):
    for header, seq in parse_fasta(fetch(url)):
        yield f"{prefix}_{header.split()[0]}", seq


def ncbi_rrna(accessions, prefix):
    url = NCBI_EFETCH_URL.format(accs=",".join(accessions))
    for header, seq in parse_fasta(fetch(url)):
        yield f"{prefix}_{header.split()[0]}", seq


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out_fa, = sys.argv[1:]

    records = []
    records += ensembl_ncrna(HUMAN_NCRNA_URL, "human")
    records += gtrnadb(GTRNADB_HUMAN_URL, "human")
    records += ncbi_rrna(NCBI_RRNA_ACCESSIONS, "human")
    records += ensembl_ncrna(YEAST_NCRNA_URL, "yeast")
    records += LITERAL_SEQUENCES.items()

    with open(out_fa, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n{seq}\n")


if __name__ == "__main__":
    main()
