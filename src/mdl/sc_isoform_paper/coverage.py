import itertools

from collections import defaultdict

import numpy as np

import pysam


from .util import read_bed


SPLICE_MATCHES = ("full_splice_match", "incomplete_splice_match")


def read_callback(tx, is_fwd, priming_tag, splice_tag):
    def fn(r):
        return (
            r.get_tag("XI") == tx
            and r.is_forward == is_fwd
            and r.get_tag("XC") == priming_tag
            and r.get_tag("XS") == splice_tag
        )

    return fn


class TranscriptData:
    def __init__(self, gencode_bed):
        self._tx = []
        self.exons = dict()
        self.strand = dict()
        self.length = dict()
        self.loc = dict()

        for r in read_bed(gencode_bed):
            if r[3][1] == "protein_coding":
                assert r[5] in "+-"
                assert not any(e_len == 0 for e_len in r[10])
                self._tx.append(tx := r[3][0])
                exons = np.array(
                    [(e_s, e_s + e_len) for e_s, e_len in zip(r[11], r[10])]
                )

                self.exons[tx] = exons
                self.strand[tx] = r[5]
                self.length[tx] = np.diff(exons, 1).sum()
                self.loc[tx] = (r[0], r[1], r[2])

        self.last_exon_r = self.calc_last_exon_ratio()

    def __iter__(self):
        return iter(self._tx)

    def __len__(self):
        return len(self._tx)

    def calc_exon_index(self, tx):
        tx_range = np.arange(*self.loc[tx][1:])
        return (
            (tx_range >= self.exons[tx][:, 0, None])
            & (tx_range < self.exons[tx][:, 1, None])
        ).any(0)

    def calc_last_exon_ratio(self):
        last_exon_pos = dict()
        last_exon_r = dict()

        for tx, exons in self.exons.items():
            assert list(exons[:, 0]) == sorted(exons[:, 0])
            assert list(exons[:, 1]) == sorted(exons[:, 1])
            exon_lens = np.diff(exons, 1)
            if len(exon_lens) == 1:
                continue
            if self.strand[tx] == "+":
                last_exon_pos[tx] = exon_lens[:-1].sum()
            else:
                last_exon_pos[tx] = exon_lens[1:].sum()
            last_exon_r[tx] = last_exon_pos[tx] / exon_lens.sum()

        return last_exon_r


transcript_data = None  # to be set globally before running


def share_tx_data(tx_data):
    global transcript_data
    transcript_data = tx_data


def calc_cov_from_bam(k, input_bam, tx, priming_tags, splice_tags):
    tx_s = transcript_data.strand[tx]
    exon_ix = transcript_data.calc_exon_index(tx)

    per_tx_depth = dict()
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        for p_tag, s_tag in itertools.product(priming_tags, splice_tags):
            tx_d = np.vstack(
                bam.count_coverage(
                    *transcript_data.loc[tx],
                    read_callback=read_callback(tx, tx_s == "+", p_tag, s_tag),
                ),
                dtype=np.int32,
            ).sum(0)[exon_ix]
            if np.any(tx_d > 0):
                per_tx_depth[k, p_tag, s_tag] = tx_d

    return tx, per_tx_depth


def max_depth_per_tx(keys, tx_d, ps=("good", "bad"), cs=SPLICE_MATCHES):
    p_cl_list = list(itertools.product(ps, cs))

    tx_max = dict()

    for key in keys:
        tx_ids = set.union(set(), *(tx_d.get((key, p, cl), []) for p, cl in p_cl_list))
        for tx in tx_ids:
            # normalize to total depth for this transcript
            tx_max[key, tx] = (
                np.vstack(
                    [
                        tx_d[key, p, cl][tx]
                        for p, cl in p_cl_list
                        if (key, p, cl) in tx_d and tx in tx_d[key, p, cl]
                    ]
                )
                .sum(0)
                .max()
            )
    return tx_max


def bin_tx_by_length(tx_data, all_tx_ids, length_bins):
    binned_tx = defaultdict(list)
    for tx in all_tx_ids:
        ix = np.searchsorted(length_bins, tx_data.length[tx])
        binned_tx[length_bins[ix]].append(tx)

    return dict(binned_tx)


def normalize_depth(
    tx_data, keys, tx_d, binned_tx, n_bins=1000, ps=("good", "bad"), cs=SPLICE_MATCHES
):
    tx_max = max_depth_per_tx(keys, tx_d, ps, cs)

    p_cl_list = list(itertools.product(ps, cs))

    tx_d2 = {
        len_b: {(key, p, cl): dict() for key in keys for p, cl in p_cl_list}
        for len_b in binned_tx
    }

    for len_b, tx_ids in binned_tx.items():
        for key in keys:
            for tx in tx_ids:
                if (txm := tx_max.get((key, tx), 0)) == 0:
                    continue
                for p, cl in p_cl_list:
                    if (key, p, cl) not in tx_d or tx not in tx_d[key, p, cl]:
                        continue
                    y = np.interp(
                        np.linspace(0, tx_data.length[tx], n_bins),
                        np.arange(tx_data.length[tx]),
                        tx_d[key, p, cl][tx]
                        / txm,  # normalize depth to max total depth for this tx
                    )

                    if tx_data.strand[tx] == "-":
                        y = y[::-1]

                    tx_d2[len_b][key, p, cl][tx] = y

    return dict(tx_d2)


def overall_depth(
    tx_data, keys, tx_d, length_bins, n_bins=1000, ps=("good", "bad"), cs=SPLICE_MATCHES
):
    length_bins = np.asarray(length_bins)

    all_tx_ids = set.union(*(set(tx_d[ksc]) for ksc in tx_d))
    binned_tx = bin_tx_by_length(tx_data, all_tx_ids, length_bins)

    normalized_depth = normalize_depth(tx_data, keys, tx_d, binned_tx, n_bins, ps, cs)

    p_cl_list = list(itertools.product(ps, cs))

    tx_d2 = {
        len_b: {
            (key, p, cl): np.zeros(n_bins, dtype=float)
            for key in keys
            for p, cl in p_cl_list
        }
        for len_b in binned_tx
    }

    for len_b in normalized_depth:
        for key, p, cl in normalized_depth[len_b]:
            for y in normalized_depth[len_b][key, p, cl].values():
                tx_d2[len_b][key, p, cl] += y

    return tx_d2, binned_tx
