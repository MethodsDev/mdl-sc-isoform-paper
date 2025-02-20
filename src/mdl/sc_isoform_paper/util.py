import csv
import gzip
import itertools
from collections import defaultdict
from collections.abc import Generator
from pathlib import Path

import pysam
import sparse

import polars as pl
import numpy as np


# standard set of chromosomes, not including alternate contigs
STD_CHRS = frozenset({f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"})

GTF_ROW = tuple[str, ...]
GTF_ATTRS = tuple[tuple[str, str], ...]


# protein-coding categories: including the IG and TR genes
PROTEIN_CODING_CATEGORIES = frozenset(
    {
        "protein_coding",
        "IG_C_gene",
        "IG_D_gene",
        "IG_J_gene",
        "IG_V_gene",
        "TR_C_gene",
        "TR_D_gene",
        "TR_J_gene",
        "TR_V_gene",
    }
)

# most of the IG and TR transcripts have these tags and are discarded
BAD_TAGS = frozenset({"readthrough_transcript", "cds_end_NF", "cds_start_NF"})


def _open(file_path: Path, mode: str = "rt"):
    if file_path.suffix == ".gz":
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def calc_mt_pct(
    m: sparse.GCXS, mt_ix: list[int] | np.ndarray[int]
) -> np.ndarray[float]:
    """calculate percent mitochondrial reads"""
    numis = m.sum(1).todense()
    return m[:, mt_ix].sum(1).todense() / np.maximum(numis, 1)


def read_fasta(fasta_file: Path) -> dict[str, str]:
    """
    read a fasta file and return the sequences for each chromosome, not including
    variants and other non-standard assemblies
    """
    chrs = defaultdict(list)
    with _open(fasta_file) as fh:
        for line in fh:
            if line.startswith(">chr"):
                c = line.strip().split()[0][1:]
            elif c in STD_CHRS:
                chrs[c].append(line.strip())

    chrs = {c: "".join(chrs[c]) for c in chrs}

    return chrs


def read_gtf(gtf_file: Path) -> Generator[tuple[GTF_ROW, GTF_ATTRS]]:
    """
    Read a GTF file and yield each non-comment line as a tuple of two values. The
    first is the full row, the second is a parsed version of the final (attr) column,
    split by semi-colons with quotes removed
    """
    with _open(gtf_file) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if r[0][0] == "#":
                continue
            rd = [v.strip().split() for v in r[8].strip().split(";") if v]
            # strip out malformed entries (seen in isoquant output)
            rd = [v for v in rd if len(v) == 2]
            rd = tuple((k, v.strip('"')) for k, v in rd)
            yield tuple(r), rd


def rd_k(rd: GTF_ATTRS, key: str) -> str:
    """
    Return the first value in the attribute list `rd` that uses the key
    """
    return next(v for k, v in rd if k == key)


def hz_k(rd: GTF_ATTRS, key: str, value: str = None) -> bool:
    """
    Return True if the key is present in this attribute list. Optionally, check that
    a given key-value pair is present
    """
    return any(k == key and (value is None or v == value) for k, v in rd)


def gtf_to_dataframe(
    gtf_file: Path, pc_cats=PROTEIN_CODING_CATEGORIES, bad_tags=BAD_TAGS
) -> pl.DataFrame:
    """
    Reads a GTF file and converts it into a polars DataFrame for efficient processing
    """
    feature_data = defaultdict(list)
    cols = [
        ("feature", pl.Categorical),
        ("transcript_id", str),
        ("gene_type", pl.Categorical),
        ("transcript_type", pl.Categorical),
        ("chromosome", pl.Categorical),
        ("is_fwd", bool),
        ("gencode_basic", bool),
        ("start", int),
        ("end", int),
    ]

    for r, rd in read_gtf(gtf_file):
        if r[2] == "gene":
            continue
        if any(k == "tag" and v in bad_tags for k, v in rd):
            continue

        c = r[0]
        is_fwd = r[6] == "+"

        tid = rd_k(rd, "transcript_id")
        gt = rd_k(rd, "gene_type")
        tt = rd_k(rd, "transcript_type")
        basic = hz_k(rd, "tag", "basic")

        if tt in pc_cats and r[2] in {"transcript", "exon", "CDS", "UTR"}:
            ft_type = r[2]
        elif r[2] == "transcript":
            ft_type = "nc_tx"
        else:
            continue

        for (col, _), v in zip(
            cols, (ft_type, tid, gt, tt, c, is_fwd, basic, int(r[3]), int(r[4]))
        ):
            feature_data[col].append(v)

    return pl.from_dict(feature_data, schema=dict(cols))


def read_bed(bed_file: Path) -> Generator[tuple[str | int | list[int], ...]]:
    with _open(bed_file) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            r[1] = int(r[1])
            r[2] = int(r[2])
            r[3] = r[3].split("|")
            r[9] = int(r[9])
            r[10] = [int(v) for v in r[10].split(",")[:-1]]
            r[11] = [int(v) + r[1] for v in r[11].split(",")[:-1]]
            yield r


def bed_to_dataframe(bed_file: Path) -> pl.DataFrame:
    """
    Reads a BED file and converts it into a polars DataFrame for efficient processing
    """
    exon_data = defaultdict(list)
    cols = [
        ("chromosome", pl.Categorical),
        ("is_fwd", bool),
        ("transcript_id", str),
        ("gene_name", str),
        ("exon_i", int),
        ("n_after", int),
        ("exon_start", int),
        ("exon_end", int),
    ]

    for r in read_bed(bed_file):
        if r[3][1] == "protein_coding":
            assert r[5] in "-+"
            for i, (exon_start, exon_len) in enumerate(zip(r[11], r[10])):
                is_fwd = r[5] == "+"
                n_after = r[9] - i - 1
                exon_end = exon_start + exon_len

                for (col, _), v in zip(
                    cols,
                    [r[0], is_fwd, r[3][0], r[3][2], i, n_after, exon_start, exon_end],
                ):
                    exon_data[col].append(v)

    return pl.from_dict(exon_data, schema=dict(cols))


def filter_reads(
    aligned_bam: pysam.AlignmentFile, min_mapq: int = 60, as_threshold: float = 0.9
):
    """
    Reads an aligned BAM file and yields reads that:
      a) are primary alignments
      b) ...to the standard chromosomes
      c) ...with a high mapping quality or alignment score

    In practice this will be similar to mapq==60 for good sequencing data,
    but it allows a little bit more data to be used
    """
    yield from (
        a
        for a in aligned_bam
        if a.is_mapped
        and a.reference_name in STD_CHRS
        and (a.mapq >= min_mapq or a.get_tag("AS") / a.qlen > as_threshold)
        and not (a.is_secondary or a.is_supplementary)
    )


def filter_gtf(gtf_file: Path, new_gtf: Path, count_file: Path, n: int = 5):
    """
    Filters an IsoQuant GTF to transcripts at or above some minimum count, based on
    the quantification from the transcript_model_count file
    """
    with _open(count_file) as fh:
        tx_counts = {
            r[0]: float(r[1])
            for r in csv.reader(itertools.islice(fh, 1, None), delimiter="\t")
        }

    with _open(gtf_file) as fh, _open(new_gtf, "wt") as out:
        for r in csv.reader(fh, delimiter="\t"):
            if r[0][0] == "#":
                print("\t".join(r), file=out)
            else:
                rd = [v.strip().split() for v in r[8].strip().split(";") if v]
                rd = [v for v in rd if len(v) == 2]
                rd = tuple((k, v.strip('"')) for k, v in rd)

                if not (
                    hz_k(rd, "transcript_id")
                    and tx_counts.get(rd_k(rd, "transcript_id"), 0) < n
                ):
                    print("\t".join(r), file=out)
