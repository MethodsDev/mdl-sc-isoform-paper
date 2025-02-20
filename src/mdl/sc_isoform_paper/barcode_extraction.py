import itertools
from collections import defaultdict
from pathlib import Path
from typing import Callable

import pysam

import bouncer

from . import rev_comp


BarcodeUMIFunc = Callable[[str, int, int, bool, int], str]


def extract_3p_bc_umi(seq: str, start: int, end: int, is_fwd: bool, buffer: int):
    # 10x 3'
    # (adapterA_TSO, cDNA, polyA, 12bp umi, 16bp bc, adapterB_RC),
    # (adapterB, 16bp bc, 12bp umi, polyT, cDNA, adapterA_TSO_RC),

    # fluent
    # (adapterA_TSO, cDNA, polyA, 12bp umi, 39bp bc, 0-3bp phase, adapterB_RC),
    # (adapterB, 0-3bp phase, 39bp bc, 12bp umi, polyT, cDNA, adapterA_TSO_RC),
    if is_fwd:
        return rev_comp(seq[end - buffer : end])
    else:
        return seq[start : start + buffer]


def extract_5p_bc_umi(seq: str, start: int, end: int, is_fwd: bool, buffer: int):
    # (adapterA_TSO, 16bp bc, 10bp umi, adapterA_TSO_SPLIT_SLS, cDNA, polyA, adapterB_RC)
    # (adapterB, polyT, cDNA, adapterA_TSO_SPLIT_SLS_RC, 10bp umi, 16bp bc, adapterA_TSO_RC),

    if is_fwd:
        return seq[start : start + buffer]
    else:
        return rev_comp(seq[end - buffer : end])


def find_barcode(seq: str, umi_len: int, barcode_set: bouncer.BarcodeSet):
    bc_matches = defaultdict(list)
    bc_dist = set()
    for bc_m, bcc, d in barcode_set.lookup_substrings(seq[:-umi_len]):
        bc_matches[bc_m].append(bcc)
        bc_dist.add(d)

    if len(bc_dist) == 0:
        return None

    assert len(bc_dist) == 1
    bc_dist = bc_dist.pop()

    if len(bc_matches) == 1:
        bc_m, bccs = bc_matches.popitem()
        if bc_dist == 0:
            ix = seq.find(bc_m) + len(bc_m)
            umi = seq[ix : ix + umi_len]
        else:
            # TODO: what is best strategy here?
            ix = min(seq.find(bcc) + len(bcc) for bcc in bccs)
            umi = seq[ix : ix + umi_len]
        return bc_m, umi

    return None


def check_ch_and_barcode(
    read: pysam.AlignedSegment,
    bc_umi_func: BarcodeUMIFunc,
    umi_size: int,
    buffer_size: int,
    barcode_set: bouncer.BarcodeSet,
):
    """
    Given a read, return a tagged read if it is proper and has a valid barcode.
    """
    if read.get_tag("lb") != "Proper":
        return None

    # whether this read is forward (A -> B) or not
    is_fwd = read.get_tag("sf") == 0

    # location of cdna sequence
    cdna_slice = slice(*read.get_tag("cd"))

    # need the coordinates for right inside the adapters
    targets = read.get_tag("ch").split(",")
    inner_start = int(targets[0].split(":")[2])  # end of the adapter on 5' end
    inner_end = int(targets[-1].split(":")[1])  # start of the adapter on 3' end

    bc_umi = find_barcode(
        bc_umi_func(read.query, inner_start, inner_end, is_fwd, buffer_size),
        umi_size,
        barcode_set,
    )
    if bc_umi is not None:
        if is_fwd:
            d = {
                "seq": read.query_sequence[cdna_slice],
                "qual": read.qual[cdna_slice],
                "flag": "4",
            }
        else:
            d = {
                "seq": rev_comp(read.query_sequence[cdna_slice]),
                "qual": read.qual[cdna_slice][::-1],
                "flag": "4",
            }
        return read.to_dict() | d, *bc_umi

    return None


barcode_set = None  # to be set globally before running


def write_barcode_bam(
    bam_file: Path,
    tagged_bam: Path,
    tag_func: BarcodeUMIFunc,
    umi_size: int,
    buffer_size: int,
):
    # TODO: can this be avoided with better PyO3 code?
    if barcode_set is None:
        raise ValueError("global barcode_set has not been initialized")

    with (
        pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=8) as fh,
        pysam.AlignmentFile(tagged_bam, "wb", template=fh, threads=8) as out,
    ):
        map_gen = map(
            check_ch_and_barcode,
            fh,
            itertools.repeat(tag_func),
            itertools.repeat(umi_size),
            itertools.repeat(buffer_size),
            itertools.repeat(barcode_set),
        )

        for i, (a, bc, umi) in enumerate(filter(None, map_gen)):
            a = pysam.AlignedSegment.from_dict(a, None)
            a.set_tag("CB", bc)
            a.set_tag("UB", umi)
            out.write(a)

    return bam_file.name, i
