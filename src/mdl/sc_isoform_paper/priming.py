import itertools
import logging

from collections import Counter, defaultdict
from enum import IntEnum, auto
from pathlib import Path

import ncls
import numpy as np
import polars as pl
import pysam

from . import isoquant
from . import util
from . import rev_comp

log = logging.getLogger(__name__)


# an enum to represent the various types of priming events we might observe, classified
#  based on their genomic context
#   'motif' is a poly-adenylation signal
#   'GPA' is 'genomic poly(A)', a stretch of adenines in the genome
#   'annotated' or 'PAS' is a cleavage/poly-adenylation site from the reference
class Priming(IntEnum):
    # in transcript, not in CDS
    ANNO_GPA = auto()  # an annotated PAS and motif, with a GPA
    GOOD = auto()  # everything looks good: annotated PAS, has PAM, no GPA
    TX_GPA_PAS = auto()  # no motif but annotated PAS, has GPA
    TX_PAS = auto()  # no motif but annotated PAS
    TX_GPA_MOTIF = auto()  # has motif, no annotated PAS, and has GPA
    TX_MOTIF = auto()  # has motif but no annotated PAS, no GPA
    TX_GPA = auto()  # no motif, has GPA
    TX_UNK = auto()  # in a transcript, nothing else was flagged

    # in transcript but in CDS
    CDS_GPA = auto()  # has a GPA
    CDS_OTHER = auto()  # no priming site identified

    # not in protein-coding transcripts
    AS_TX_GPA_NC = auto()  # antisense to a transcript and in an noncoding RNA, has GPA
    NC_GPA = auto()  # in a noncoding transcript, with GPA
    AS_TX_NC = auto()  # antisense to a transcript and also in a noncoding RNA, no GPA
    NC_OTHER = auto()  # in a noncoding transcript region, no GPA

    AS_TX_GPA = auto()  # antisense to a transcript, with GPA
    INTERGENIC_GPA = auto()  # not in a transcript or noncoding region, but has GPA
    AS_TX = auto()  # antisense to a transcript, no GPA
    INTERGENIC_OTHER = auto()  # not in a transcript or noncoding region, no GPA

    MITO = auto()  # maps to mitochondria

    def __str__(self):
        return self.name


_CONDITIONS = [
    # transcript reads, with and without annotation, motif, genonmic polya
    *(
        dict(anno_PAS=b1, motif=b2, g_polya=b3, in_tx=True, in_cds=False, in_mito=False)
        for b1, b2, b3 in itertools.product((True, False), repeat=3)
    ),
    # cds reads
    *(dict(g_polya=b, in_tx=True, in_cds=True, in_mito=False) for b in (True, False)),
    # non-coding and intergenic reads:
    *(
        dict(in_nc=b1, g_polya=b2, in_tx_AS=b3, in_tx=False, in_mito=False)
        for b1, b2, b3 in itertools.product((True, False), repeat=3)
    ),
    # mito reads
    dict(in_mito=True),
]


_FLAGS = {
    f: i
    for i, f in enumerate(
        (
            "anno_PAS",
            "motif",
            "g_polya",
            "in_tx",
            "in_tx_AS",
            "in_cds",
            "in_nc",
            "in_mito",
        )
    )
}


def _check_flags(t, **kwargs):
    return all(t[_FLAGS[k]] == v for k, v in kwargs.items())


def _make_classes(flags, conditions) -> dict[tuple[bool, ...], Priming]:
    b_to_class = []
    for b in itertools.product((False, True), repeat=len(flags)):
        i = {i for i, cond in enumerate(conditions) if _check_flags(b, **cond)}
        if len(i) == 1:
            b_to_class.append((b, Priming(i.pop() + 1)))
        else:
            raise ValueError(f"{len(i)} classes for condition {b}")

    assert len(set(b_to_class)) == 2 ** len(flags)

    return dict(b_to_class)


B_TO_CLASS = _make_classes(_FLAGS, _CONDITIONS)


def calc_longest_a(seq: str, polya_window: int) -> int:
    best_n = 0
    curr_n = 0
    for i, c in enumerate(seq):
        if c == "A":
            curr_n += 1
        else:
            best_n = max(curr_n, best_n)
            curr_n = 0
            if i > polya_window:
                break

    return max(best_n, curr_n)


class ReferenceRanges:
    MAX_POLYA = 90  # longest stretch of A or T in GRCh38 genome

    def __init__(
        self,
        reference_fasta: Path,
        reference_gtf: Path,
        polya_motif_file: Path,
        polya_annotations: Path,
    ):
        self.chrs = util.read_fasta(reference_fasta)

        with open(polya_motif_file) as fh:
            self.polya_motifs = {line.strip().upper() for line in fh}

        self.reference_gtf = reference_gtf
        self.polya_annotations = polya_annotations
        self.ranges = None

    def build_ranges(self):
        """
        Build NCLS ranges for overlap finding, using a GTF reference. When using
        multiple processes, this should be performed within each subprocess because the
        NCLS objects can't be shared
        """
        if self.ranges is not None:
            log.warning("Overwriting existing feature ranges")

        if (
            self.reference_gtf.suffixes[-2:] == [".gtf", ".gz"]
            or self.reference_gtf.suffix == ".gtf"
        ):
            log.debug(f"Converting {self.reference_gtf} to dataframe")
            feature_df = util.gtf_to_dataframe(self.reference_gtf)
        elif self.reference_gtf.suffix in {".pq", ".parquet"}:
            log.debug(f"Loading parquet file {self.reference_gtf}")
            feature_df = pl.read_parquet(self.reference_gtf)
        else:
            raise ValueError(f"Unrecognized file format {self.reference_gtf}")

        log.debug("Filtering to gencode basic tag")
        feature_df = feature_df.filter(pl.col("gencode_basic"))

        feature_ranges = defaultdict(dict)
        for (feature_type, chrom, is_fwd), df in feature_df.group_by(
            "feature", "chromosome", "is_fwd"
        ):
            # NOTE: need to add 1 to the end of segments because GTF uses inclusive ranges
            feature_ranges[feature_type][chrom, is_fwd] = ncls.NCLS(
                df["start"].to_numpy(), df["end"].to_numpy() + 1, np.arange(df.shape[0])
            )

        pa_extent = defaultdict(list)

        for r, _ in util.read_gtf(self.polya_annotations):
            if r[2] == "polyA_site":
                pa_extent[r[0], r[6] == "+"].append((int(r[3]), int(r[4])))

        for chrom, is_fwd in pa_extent:
            starts, ends = zip(*pa_extent[chrom, is_fwd])
            feature_ranges["polyA_site"][chrom, is_fwd] = ncls.NCLS(
                np.array(starts), np.array(ends) + 1, np.arange(len(starts))
            )

        self.ranges = feature_ranges

    def get_pa_motifs(self, reference_seq):
        return self.polya_motifs.intersection(
            reference_seq[i : i + 6] for i in range(len(reference_seq) - 5)
        )

    def check_feature(self, feature_type, chrom, fwd, start, end) -> bool:
        if (chrom, fwd) in self.ranges[feature_type]:
            return bool(
                next(
                    self.ranges[feature_type][chrom, fwd].find_overlap(start, end),
                    False,
                )
            )
        else:
            return False


class PrimingClassifier:
    GOOD_PRIMING_TAGS = frozenset(
        {Priming.ANNO_GPA, Priming.GOOD, Priming.TX_PAS, Priming.TX_MOTIF}
    )

    def __init__(
        self,
        *,
        feature_pre: int = 5,
        feature_post: int = 5,
        motif_pre: int = 20,
        motif_post: int = 30,
        pas_pre: int = 5,
        pas_post: int = 20,
        polya_window: int = 20,
        polya_max_n: int = 12,
        polya_max_len: int = 6,
    ):
        # overlap with GTF features (transcript, exon, UTR, etc)
        self.feature_pre = feature_pre
        self.feature_post = feature_post

        # overlap this window with polya motif sequences
        self.motif_pre = motif_pre
        self.motif_post = motif_post

        # search for a polya signal from reference
        self.pas_pre = pas_pre
        self.pas_post = pas_post

        self.polya_window = polya_window
        self.polya_max_n = polya_max_n
        self.polya_max_len = polya_max_len

    def read_flags(
        self, reference: ReferenceRanges, chrom: str, start: int, end: int, fwd: bool
    ):
        feature_pre = self.feature_pre
        feature_post = self.feature_post
        pas_pre = self.pas_pre
        pas_post = self.pas_post

        if fwd:
            polya_motif_window = reference.chrs[chrom][
                end - self.motif_pre : end + self.motif_post
            ]

            genomic_polya_window = reference.chrs[chrom][
                end : end + reference.MAX_POLYA + self.polya_window
            ]
        else:
            polya_motif_window = rev_comp(
                reference.chrs[chrom][start - self.motif_post : start + self.motif_pre]
            )
            genomic_polya_window = rev_comp(
                reference.chrs[chrom][
                    start - reference.MAX_POLYA - self.polya_window : start
                ]
            )

            start, end = end, start
            feature_pre, feature_post = feature_post, feature_pre
            pas_pre, pas_post = pas_post, pas_pre

        # read coordinates are 0-indexed, gtf is 1-indexed
        start += 1
        end += 1

        annotated_pas = reference.check_feature(
            "polyA_site", chrom, fwd, end - pas_pre, end + pas_post
        )
        n_motifs = len(reference.get_pa_motifs(polya_motif_window))

        n_genomic_a = genomic_polya_window[: self.polya_window].count("A")
        longest_a = calc_longest_a(genomic_polya_window, self.polya_window)

        in_tx = reference.check_feature(
            "transcript", chrom, fwd, end - feature_pre, end + feature_post
        )
        in_tx_AS = reference.check_feature(
            "transcript", chrom, not fwd, start - feature_post, start + feature_pre
        )
        in_cds = reference.check_feature(
            "CDS", chrom, fwd, end - feature_pre, end + feature_post
        )
        in_nc_tx = reference.check_feature(
            "nc_tx", chrom, fwd, end - feature_pre, end + feature_post
        )

        in_mito = chrom == "chrM"

        return (
            annotated_pas,
            n_motifs,
            n_genomic_a,
            longest_a,
            in_tx,
            in_tx_AS,
            in_cds,
            in_nc_tx,
            in_mito,
        )

    def classify_read(self, ref: ReferenceRanges, a: pysam.AlignedSegment) -> Priming:
        chrom = a.reference_name
        start = a.reference_start
        end = a.reference_end
        is_fwd = a.is_forward

        anno_pas, n_motifs, n_A, max_A, *flags = self.read_flags(
            ref, chrom, start, end, is_fwd
        )
        has_motif = n_motifs > 0
        genomic_polya = n_A > self.polya_max_n or max_A > self.polya_max_len

        return B_TO_CLASS[(anno_pas, has_motif, genomic_polya, *flags)]

    def tag_bam_with_read_stats(
        self,
        bam_file: Path,  # mapped bam file
        anno_bam: Path,  # output bam
        reference: ReferenceRanges,  # reference data
        isoquant_path: Path = None,  # e.g. path/to/OUT
        use_prefix: bool = True,  # filter isoquant files by bam file name
        full_tag: bool = True,  # output the priming tag classification
        simple_tag: bool = False,  # output boolean good/bad priming calls
    ):
        # need to do this inside each subprocess because NCLS objects can't be pickled
        reference.build_ranges()

        if use_prefix:
            # the smrtcell ID for this run
            bam_name = bam_file.name.split(".")[0]
        else:
            bam_name = None

        if isoquant_path is not None:
            log.debug("Loading IsoQuant files")
            rn_to_tx, rn_to_model = isoquant.read_assignment_data(
                next(isoquant_path.glob("*.transcript_models.gtf")),
                next(isoquant_path.glob("*.read_assignments.tsv.gz")),
                next(isoquant_path.glob("*.transcript_model_reads.tsv.gz")),
                bam_name,
            )

        log.debug("Annotating BAM file %s", bam_file)
        with (
            pysam.AlignmentFile(bam_file, "rb", threads=8) as fh,
            pysam.AlignmentFile(anno_bam, "wb", template=fh, threads=8) as out,
        ):
            for a in util.filter_reads(fh):
                # XC: priming classification
                priming_class = self.classify_read(reference, a)
                if full_tag:
                    a.set_tag("XC", str(priming_class))
                # XB: boolean good/bad priming classification
                if simple_tag:
                    a.set_tag("XB", int(priming_class in self.GOOD_PRIMING_TAGS))

                if isoquant_path is not None:
                    # adding various isoquant annotations to the BAM. These may not be
                    # present for a given read (or were ambiguous). Default to "-"
                    #   XI: isoform (transcript id) from read assignment
                    #   XQ: isoquant assignment type
                    #   XS: SQANTI-like match classification (from isoquant)
                    for t, v in zip(
                        ("XI", "XQ", "XS"), rn_to_tx.get(a.query_name, "---")
                    ):
                        a.set_tag(t, str(v), "Z")

                    # YT: transcript id after isoform model construction
                    # YC: transcript type (known, nic, nnic, or novel gene)
                    for t, v in zip(("YT", "YC"), rn_to_model.get(a.query_name, "--")):
                        a.set_tag(t, str(v), "Z")

                out.write(a)
        log.debug("Finished writing annotation to %s", anno_bam)


def count_classes_and_isoquant(key, bam_file):
    class_counter = Counter()
    with pysam.AlignmentFile(bam_file, "rb", threads=2) as fh:
        for a in fh:
            tx = a.get_tag("YT")
            if (iq := a.get_tag("YC")) != "-":
                iq = isoquant.IsoQuantClass[iq]
            p = Priming[a.get_tag("XC")]
            class_counter[tx, iq, p] += 1

    return key, class_counter


def tx_count_breakdown(tx_big, good_classes=PrimingClassifier.GOOD_PRIMING_TAGS):
    tx_class = defaultdict(lambda: defaultdict(set))
    tx_count = defaultdict(Counter)
    tx_good = defaultdict(Counter)
    tx_bad = defaultdict(Counter)

    for k in tx_big:
        for (tx, tx_type, cl), v in tx_big[k].items():
            if tx_type != "-" and tx != "-":
                tx_class[k][tx].add(tx_type)
                tx_count[k][tx] += v
                if cl in good_classes:
                    tx_good[k][tx] += v
                else:
                    tx_bad[k][tx] += v

    assert max(len(v) for k in tx_class for v in tx_class[k].values()) == 1
    tx_class = {k: {tx: tx_class[k][tx].pop() for tx in tx_class[k]} for k in tx_class}

    tx_ratio = defaultdict(dict)

    for k in tx_big:
        for tx in tx_count[k]:
            good_sum = tx_good[k][tx]
            bad_sum = tx_bad[k][tx]

            if good_sum + bad_sum:
                tx_ratio[k][tx] = good_sum / (good_sum + bad_sum)

    return tx_class, tx_count, tx_ratio
