import csv
from collections import Counter, defaultdict
from pathlib import Path


CONFIG_DICT = {
    "TSO_adapters": ["adapterA"],
    "RT_adapters": ["adapterB"],
    "min_rq": 0.99,
    "min_polyA_match": 10,
    "max_err_rate": 0.05,
    "max_err_rate_polya": 0.05,
    "terminal_adapter_search_buffer": 150,
    "terminal_polyA_search_buffer": 150,
    "n_threads": 8,
}


SAMPLE_CONFIG = {
    "PIPseq": {
        "adapterA": "AAGCAGTGGTATCAACGCAGAG",
        "adapterB": "CTACACGACGCTCTTCCGATCT",
        "adapterA_SLS": "TGAATGGG",
    },
    "10x 3'": {
        "adapterA": "AAGCAGTGGTATCAACGCAGAG",
        "adapterB": "CTACACGACGCTCTTCCGATCT",
        "adapterA_SLS": "TACATGGG",
    },
    "10x 5'": {
        "adapterA": "CTACACGACGCTCTTCCGATCT",
        "adapterB": "AAGCAGTGGTATCAACGCAGAG",
        "adapterA_SPLIT_SLS": "TTTCTTATATGGG",
    },
}


def read_sample_reports(run_path: Path):
    """Read in the structure_counts file from the marti reports directory"""
    sample_heads = defaultdict(Counter)

    for fp in run_path.glob("*/reports/structure_counts.tsv"):
        sample = fp.parent.parent.name
        with fp.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sample_heads[sample][row["read_type"]] += int(row["count"])

    return sample_heads


def total_samples(run_path: Path, sample_key: dict[int, str]):
    """Sum up a bunch of structure counts by sample. This assumes that the directories
    are named as *.*.[int] and that these ints are keys into `sample_key`"""

    sample_heads = read_sample_reports(run_path)
    sample_totals = defaultdict(Counter)

    for s in sample_heads:
        sample_totals[sample_key[int(s.split(".")[2])]] += sample_heads[s]

    return sample_totals
