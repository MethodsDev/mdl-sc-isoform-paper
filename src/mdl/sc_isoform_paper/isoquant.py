import csv
import gzip
import itertools
import re
from collections import Counter
from enum import Enum, auto
from pathlib import Path

from .util import read_gtf, STD_CHRS, rd_k


class IsoQuantClass(Enum):
    KNOWN_KNOWN = auto()  # gene and transcript from the reference
    KNOWN_NIC = auto()  # novel-in-catalog transcript from known gene
    KNOWN_NNIC = auto()  # novel-not-in-catalog transcript from known gene
    NOVEL_NNIC = auto()  # novel transcript from a novel gene

    def __str__(self):
        return self.name


def calc_tx_class(g_id, tx_id):
    if g_id[:4] == "ENSG":
        if tx_id[:4] == "ENST":
            return IsoQuantClass.KNOWN_KNOWN
        elif tx_id.startswith("transcript"):
            if tx_id[-4:] == ".nic":
                return IsoQuantClass.KNOWN_NIC
            elif tx_id[-5:] == ".nnic":
                return IsoQuantClass.KNOWN_NNIC
            else:
                raise ValueError(tx_id)
        else:
            raise ValueError(tx_id)
    elif g_id.startswith("novel_gene"):
        if tx_id[-5:] == ".nnic":
            return IsoQuantClass.NOVEL_NNIC
        else:
            raise ValueError(tx_id)
    else:
        raise ValueError(g_id)


def read_iq_model_file(isoquant_tx_model_file: Path):
    tx_to_gene = dict()
    tx_to_class = dict()

    for r, rd in read_gtf(isoquant_tx_model_file):
        if r[0] not in STD_CHRS:
            continue

        if r[2] == "transcript":
            g_id = rd_k(rd, "gene_id")
            tx_id = rd_k(rd, "transcript_id")
            tx_to_gene[tx_id] = g_id
            tx_to_class[tx_id] = calc_tx_class(g_id, tx_id)

    return tx_to_gene, tx_to_class


def read_assignment_data(
    isoquant_tx_model_file: Path,
    isoquant_read_file: Path,
    isoquant_model_reads: Path,
    prefix: str = None,
):
    # read in the generated GTF
    _, tx_to_class = read_iq_model_file(isoquant_tx_model_file)

    iq_count = Counter()
    rn_to_tx = dict()
    with gzip.open(isoquant_read_file, "rt") as fh:
        for r in csv.reader(itertools.islice(fh, 3, None), delimiter="\t"):
            if (prefix is None or r[0].startswith(prefix)) and r[1] in STD_CHRS:
                iq_count[r[0]] += 1
                m = re.search("Classification=([a-z_]+)", r[8])
                if m:
                    iq = m.group(1)
                else:
                    iq = "-"
                rn_to_tx[r[0]] = (r[3], r[5], iq)

    rn_to_tx = {rn: v for rn, v in rn_to_tx.items() if iq_count[rn] == 1}

    iq_count = Counter()
    rn_to_model = dict()
    with gzip.open(isoquant_model_reads, "rt") as fh:
        for r in csv.reader(itertools.islice(fh, 1, None), delimiter="\t"):
            if prefix is None or r[0].startswith(prefix):
                iq_count[r[0]] += 1
                if r[1] in tx_to_class:
                    rn_to_model[r[0]] = (r[1], tx_to_class[r[1]])

    rn_to_model = {rn: v for rn, v in rn_to_model.items() if iq_count[rn] == 1}

    return rn_to_tx, rn_to_model
