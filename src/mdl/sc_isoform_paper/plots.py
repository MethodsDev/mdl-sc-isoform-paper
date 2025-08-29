import itertools
import importlib.resources
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib import ticker
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import gaussian_kde, brunnermunzel

from mdl.sc_isoform_paper.pipseq_barcodes import sequence_to_int, barcode_to_int
from pathlib import Path
from typing import Sequence, Optional, Callable, Dict, Any, Tuple
from concurrent.futures import ProcessPoolExecutor
import seaborn as sns
import pandas as pd
import pysam
import re

from .constants import SAMPLE_COLORS
from .coverage import SPLICE_MATCHES
from .isoquant import IsoQuantClass


# load the custom style for our figure panels
# Note: we use the Open Sans font, which might need to be installed from
# https://fonts.google.com/specimen/Open+Sans
with importlib.resources.path(__package__, "figures.mplstyle") as log_config:
    plt.style.use(log_config)


_SKERA_COLORS = [
    (249 / 255, 157 / 255, 65 / 255),
    (255 / 255, 102 / 255, 204 / 255),
    (141 / 255, 109 / 255, 176 / 255),
    (68 / 255, 168 / 255, 223 / 255),
    (0 / 255, 184 / 255, 196 / 255),
    (106 / 255, 191 / 255, 105 / 255),
    (225 / 255, 106 / 255, 44 / 255),
    (223 / 255, 25 / 255, 149 / 255),
    (95 / 255, 36 / 255, 159 / 255),
    (19 / 255, 131 / 255, 198 / 255),
    (0 / 255, 156 / 255, 162 / 255),
    (0 / 255, 157 / 255, 78 / 255),
    (195 / 255, 74 / 255, 33 / 255),
    (0 / 255, 84 / 255, 150 / 255),
    (0 / 255, 91 / 255, 66 / 255),
    (49 / 255, 53 / 255, 65 / 255),
]


_LINE_STYLES = [
    "solid",
    "dotted",
    (0, (1, 1)),
    (0, (1, 10)),
    (5, (10, 3)),
    (0, (5, 10)),
    (0, (5, 5)),
    (0, (5, 1)),
    (0, (3, 10, 1, 10)),
    (0, (3, 5, 1, 5)),
    (0, (3, 1, 1, 1)),
    (0, (3, 5, 1, 5, 1, 5)),
    (0, (3, 10, 1, 10, 1, 10)),
    (0, (3, 1, 1, 1, 1, 1)),
]


def plot_dists(ax, dists, log=False, colors=None, labels=None, title=None):
    if colors is None:
        colors = ["#006DB6"] * len(dists)
    if log:
        dists = [np.log10(d) for d in dists]

    p = ax.violinplot(dists, showmedians=True, widths=0.8)

    for c, b in zip(colors, p["bodies"]):
        b.set_color(c)
        b.set_alpha(1.0)
    for c in ["cmaxes", "cmins", "cmedians"]:
        p[c].set_color("#808080")
    p["cbars"].remove()

    if title is not None:
        ax.set_title(title)
    if labels is not None:
        ax.set_xticks(range(1, len(labels) + 1), labels)


def plot_marker_dotplot(
    labels, markers, marker_nz, marker_counts, z_ix=None, output_file=None
):
    if z_ix is None:
        z_ix = np.arange(marker_nz.shape[1])

    fig, ax = plt.subplots(1, 1, figsize=(6, 8))
    ax.scatter(
        *np.indices(marker_nz.shape),
        s=100 * (marker_nz[:, z_ix].flatten()) ** 2,
        c=marker_counts[:, z_ix].flatten(),
        cmap="Blues",
        norm=colors.LogNorm(),
    )

    ax.set_xticks(np.arange(len(markers)), markers, rotation=90)

    ax.yaxis.remove_overlapping_locs = False
    ax.set_yticks(
        np.arange(len(labels)),
        [labels[i][1] for i in z_ix],
        minor=True,
    )
    lbl_cat = Counter()
    lbl_order = []
    for i in z_ix:
        if (lbl := labels[i][0]) not in lbl_cat:
            lbl_order.append(lbl)
        lbl_cat[lbl] += 1

    ax.set_yticks(
        np.array([lbl_cat[lbl] for lbl in lbl_order]).cumsum()
        - np.array([lbl_cat[lbl] for lbl in lbl_order]) / 2
        - 0.5,
        lbl_order,
        fontsize="medium",
        minor=False,
    )
    ax.tick_params(axis="y", which="major", pad=60, length=0)

    if output_file is not None:
        plt.savefig(output_file)

    plt.show()


# function for plotting ambient contamination: % of cross-species reads
# in the single cells (non-doublets), with human and mouse cells combined into a single panel
def plot_ambient_distributions(
    labels, h_data_list, m_data_list, hi_cutoff=10000, lo_cutoff=1000, output_file=None
):
    fig, axs = plt.subplots(
        1, 2, figsize=(12, 4), sharex=True, gridspec_kw={"bottom": 0.2}
    )

    for lbl, h_data, m_data in zip(labels, h_data_list, m_data_list):
        t_data = h_data + m_data
        ix_h = (h_data > hi_cutoff) & (m_data < lo_cutoff)
        ix_m = (m_data > hi_cutoff) & (h_data < lo_cutoff)
        print(
            lbl,
            np.median(
                np.hstack([m_data[ix_h] / t_data[ix_h], h_data[ix_m] / t_data[ix_m]])
            ),
        )

        for i in (0, 1):
            axs[i].hist(
                np.hstack([m_data[ix_h] / t_data[ix_h], h_data[ix_m] / t_data[ix_m]]),
                bins=np.linspace(0, 0.1, 101),
                log=bool(i),
                density=True,
                histtype="step",
                label=lbl,
                color=SAMPLE_COLORS[lbl],
            )

    for ax in axs:
        ax.set_ylabel("# of Barcodes")
        ax.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
        ax.set_xticks([0, 0.05, 0.1])

    fig.legend(
        *axs[0].get_legend_handles_labels(), ncols=len(labels), loc="lower center"
    )
    fig.suptitle(
        f"Cross-species contamination in cells (A > {hi_cutoff:,}, B < {lo_cutoff:,})"
    )

    if output_file is not None:
        plt.savefig(output_file)
    plt.show()


# function for plotting the barnyard data
def plot_barn(fig, ax, title, h_data, m_data, cutoff=None):
    # version of 'Blues' colormap that is pure white at the bottom
    cmap = colors.LinearSegmentedColormap.from_list(
        "BluesW",
        [
            (0, (1.0, 1.0, 1.0)),
            (0.9, (0 / 255, 109 / 255, 182 / 255)),
            (1.0, (0 / 255, 109 / 255, 182 / 255)),
        ],
    )

    m = ax.hexbin(
        np.maximum(h_data, 1),
        np.maximum(m_data, 1),
        xscale="log",
        yscale="log",
        cmap=cmap,
        norm=colors.LogNorm(vmin=1, vmax=500000, clip=True),
    )

    axins1 = inset_axes(
        ax,
        width="30%",
        height="4%",
        loc="upper left",
        bbox_to_anchor=(0.05, 0.0, 0.9, 1),
        bbox_transform=ax.transAxes,
    )
    axins1.xaxis.set_ticks_position("bottom")
    fig.colorbar(m, cax=axins1, orientation="horizontal")

    if cutoff is not None:
        ax.axvline(cutoff, linestyle=":", color="#63666A")
        ax.axhline(cutoff, linestyle=":", color="#63666A")

    ax.set_xlabel("Human")
    ax.set_ylabel("Mouse")
    ax.set_title(title)


def plot_concat_and_ligations(
    len_to_concat, counts, exclude_nums, arraysize=16, xmax=25000
):
    fig, ax = plt.subplots(
        1,
        4,
        figsize=(16, 5),
        gridspec_kw={"width_ratios": [15, 3, 12, 1], "wspace": 0.05},
    )
    ax[0].set_title("MAS-seq array length and concatentation factor")

    binsize = 250
    bins = np.arange(0, xmax + 1, binsize)

    cat_to_len = defaultdict(list)
    for rl, cat in len_to_concat:
        cat_to_len[cat].append(rl)

    ax[0].hist(
        [cat_to_len[i] for i in sorted(cat_to_len)],
        bins=bins,
        color=_SKERA_COLORS[:arraysize],
        histtype="stepfilled",
        stacked=True,
        label=[f"{i}x" for i in sorted(cat_to_len)],
    )

    ax[0].set_xlim(0, xmax + binsize)
    ax[0].set_xlabel("Read length (kb)", fontsize="x-large")
    ax[0].set_ylabel("Number of Reads", fontsize="x-large")
    ax[0].set_xticks(
        range(0, xmax + 1, binsize * 20),
        labels=range(0, xmax + 1, binsize * 20),
    )

    ax[0].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x // 1000}"))
    y_fmt = ticker.ScalarFormatter(useMathText=True)
    y_fmt.set_powerlimits((-3, 3))
    ax[0].yaxis.set_major_formatter(y_fmt)

    ax[0].tick_params(axis="both", labelsize="x-large")

    # legend
    for idx, c in enumerate(_SKERA_COLORS[:arraysize]):
        bottom = idx + 0.2 + (0.1 * idx)
        cbox = Rectangle((0.2, bottom), 1, 1, lw=0, color=c)
        ax[1].add_patch(cbox)
        ax[1].text(1.4, bottom + 0.5, str(idx + 1) + "x", ha="left", va="center")
    ax[1].set_xlim(0, 2)
    ax[1].set_ylim(0, bottom + 1.2)
    ax[1].tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=False,
        labelleft=False,
        right=False,
        labelright=False,
        top=False,
        labeltop=False,
    )
    ax[1].axis(False)

    edge_len = len(counts)

    m = ax[2].matshow(counts, cmap="Blues")
    ax[2].set_title("MAS ligation heatmap")
    if not exclude_nums:
        for i in range(edge_len):  # row
            for j in range(edge_len):  # column
                ax.text(
                    i,
                    j,
                    str(counts[j][i]),
                    ha="center",
                    va="center",
                    fontsize=6,
                    zorder=10,
                    color="w" if i == j - 1 else "k",
                )

    # heatmap params
    letters = [chr(x + 65) for x in range(edge_len)]
    ax[2].set_xlim(-0.5, edge_len - 0.5)
    ax[2].set_ylim(-0.5, edge_len - 0.5)
    ax[2].set_xticks(np.arange(edge_len))
    ax[2].set_yticks(np.arange(edge_len))
    ax[2].set_xticklabels(letters)
    ax[2].set_yticklabels(letters)
    ax[2].set_xlabel("Adapter 1", fontsize="x-large")
    ax[2].set_ylabel("Adapter 2", fontsize="x-large")
    ax[2].tick_params(
        axis="both",
        which="both",
        bottom=True,
        labelbottom=True,
        left=True,
        labelleft=True,
        right=False,
        labelright=False,
        top=False,
        labeltop=False,
    )

    fig.colorbar(m, cax=ax[3])

    return fig


def plot_stacked_concat(skera_read_lengths, arraysize=16, normalize=False, labels=None):
    fig, ax = plt.subplots(
        1,
        2,
        figsize=(8, 5),
        gridspec_kw={"width_ratios": [15, 3], "wspace": 0.05},
    )
    ax[0].set_title("MAS-seq array length and concatentation factor")

    keys = list(skera_read_lengths)

    cat_hists = dict()
    for k in keys:
        uvals, cat_hists[k] = np.unique(
            [cat for _, cat in skera_read_lengths[k]], return_counts=True
        )
        assert np.array_equal(uvals, np.arange(1, arraysize + 1))

    if normalize:
        for k in keys:
            cat_hists[k] = cat_hists[k] / cat_hists[k].sum()

    x = np.arange(len(keys))
    base_height = np.zeros(len(keys))
    for i in reversed(range(arraysize)):
        h = np.array([cat_hists[k][i] for k in keys])
        ax[0].bar(
            x,
            h,
            width=0.8,
            bottom=base_height,
            color=_SKERA_COLORS[:arraysize][i],
            label=f"{i}x",
        )
        base_height += h

    ax[0].set_xlim(-0.5, len(keys) - 0.5)
    ax[0].set_xticks(x, labels=labels or keys)

    # ax[0].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x // 1000}"))
    if normalize:
        ax[0].set_ylabel("Percentage of Reads", fontsize="x-large")
        ax[0].yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1))
    else:
        ax[0].set_ylabel("Number of Reads", fontsize="x-large")
        y_fmt = ticker.ScalarFormatter(useMathText=True)
        y_fmt.set_powerlimits((-3, 3))
        ax[0].yaxis.set_major_formatter(y_fmt)

    ax[0].tick_params(axis="both", labelsize="x-large")

    # legend
    for idx, c in enumerate(_SKERA_COLORS[:arraysize]):
        bottom = idx + 0.2 + (0.1 * idx)
        cbox = Rectangle((0.2, bottom), 1, 1, lw=0, color=c)
        ax[1].add_patch(cbox)
        ax[1].text(1.4, bottom + 0.5, str(idx + 1) + "x", ha="left", va="center")
    ax[1].set_xlim(0, 2)
    ax[1].set_ylim(0, bottom + 1.2)
    ax[1].tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=False,
        labelleft=False,
        right=False,
        labelright=False,
        top=False,
        labeltop=False,
    )
    ax[1].axis(False)

    return fig


def plot_artifacts(
    sample_list, sample_totals, key_list, *, title=None, output_file=None
):
    _, ax = plt.subplots(1, 1, figsize=(12, 8), gridspec_kw={"bottom": 0.3})

    x = np.arange(len(key_list))
    w = 0.8 / len(sample_list)

    for i, s in enumerate(sample_list):
        s_counts = sample_totals[s]

        tot = s_counts.total()
        vs = np.array([s_counts[k] for k in key_list])
        ax.bar(
            x + 0.1 + i * w, vs / tot, width=w, align="edge", color=SAMPLE_COLORS[s[0]]
        )
        for j, (k, v) in enumerate(zip(key_list, vs)):
            ax.annotate(
                f"{v:,d}\n({(v / tot):.2%})",
                (j + 0.1 + (i + 0.5) * w, v / tot),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                rotation=75,
                fontsize="small",
            )

    ax.xaxis.remove_overlapping_locs = False
    ax.set_xticks(
        (
            x[:, None] + np.linspace(0.1 + 0.5 * w, 0.9 - 0.5 * w, len(sample_list))
        ).flatten(),
        [s[0] for s in sample_list] * len(key_list),
        rotation=45,
        ha="center",
        minor=True,
    )
    ax.set_xticks(x + 0.5, key_list, fontsize="medium", minor=False, rotation=45)
    ax.tick_params(axis="x", which="major", pad=50)
    ax.set_ylim(0, 0.18)
    ax.margins(y=0.2)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))

    if title is not None:
        ax.set_title(title)

    if output_file is not None:
        plt.savefig(output_file)
    plt.show()


def plot_isoform_area(
    sample_order,
    tx_class,
    tx_count,
    tx_ratio,
    min_count=5,
    tx_types=tuple(IsoQuantClass),
    output_path=None,
):
    fig, ax = plt.subplots(
        1,
        len(tx_types) * len(sample_order),
        figsize=(2 * len(sample_order) * len(tx_types), 6),
        sharex=True,
        sharey=True,
        gridspec_kw={"hspace": 0.1, "wspace": 0.05, "bottom": 0.05},
    )

    for i, (txt, k) in enumerate(itertools.product(tx_types, sample_order)):
        tx_list = sorted(
            (
                tx
                for tx in tx_count[k]
                if tx_count[k][tx] >= min_count and tx_class[k][tx] == txt
            ),
            key=tx_ratio[k].get,
            reverse=True,
        )

        y = np.arange(len(tx_list))

        good_x1 = 0.0
        good_x2 = np.array([tx_ratio[k][tx] for tx in tx_list])

        bad_x1 = good_x2
        bad_x2 = 1.0

        ax[i].fill_betweenx(
            y[good_x2 > 0.0],
            x1=good_x1,
            x2=good_x2[good_x2 > 0.0],
            color=SAMPLE_COLORS[k[0]],
            rasterized=True,
            label=" ".join(k) if i < 3 else None,
        )

        ax[i].fill_betweenx(
            y,
            x1=bad_x1,
            x2=bad_x2,
            color=SAMPLE_COLORS[k[0]],
            alpha=0.2,
            rasterized=True,
        )

        p_y = sum(1 for tx in tx_list if tx_ratio[k][tx] > 0.0)
        ax[i].axhline(p_y, xmin=0.05, xmax=0.95, linestyle=":", color="black")
        ax[i].annotate(
            f"{p_y:,}", (0.5, p_y), (0, 5), textcoords="offset points", ha="center"
        )

        ax[i].annotate(
            f"{len(tx_list):,}",
            (0.5, len(tx_list)),
            (0, 5),
            textcoords="offset points",
            ha="center",
        )
        ax[i].axis("off")

    fig.legend(loc="lower center", ncol=len(sample_order))

    for i, txt in enumerate(tx_types):
        ax[3 * i + 1].set_title(txt)

    fig.suptitle("Isoquant Transcript Calls")

    if output_path is not None:
        plt.savefig(output_path)
    plt.show()


def plot_exon_ratio_kde(ax, last_exon_r, tx_bin, n_bins=1000):
    last_exon_ratio = np.array([last_exon_r[tx] for tx in tx_bin])
    g = gaussian_kde(last_exon_ratio)
    ax.fill_between(
        np.linspace(0, n_bins, 1000),
        g.evaluate(np.linspace(0, 1, 1000)),
        step="mid",
        alpha=0.3,
    )
    ax.set_ylim(bottom=0)
    ax.tick_params(right=False, bottom=False, labelright=False, labelbottom=False)
    ax.axis("off")


def plot_total_cov(ax, keys, tx_d, *, ps=("good", "bad"), cs=SPLICE_MATCHES):
    for key in keys:
        # normalize everyone to the highest (mean) depth for this experiment
        m = max(tx_d[key, p, cl].max() for p in ps for cl in cs)
        for i, cl in enumerate(cs):
            for p, s in zip(ps, _LINE_STYLES):
                d = tx_d[key, p, cl] / m

                ax[i].step(
                    np.arange(d.shape[0]),
                    d,
                    alpha=0.8,
                    color=SAMPLE_COLORS[key[0]],
                    where="mid",
                    label=f"{key[0]} - {p}",
                    linestyle=s,
                )

            ax[i].set_xticks(
                np.linspace(0, len(d) - 1, 2),
                [f"{b:.0%}" for b in np.linspace(0, 1, 2)],
            )
            ax[i].set_ylim(bottom=0, top=1.05)
            ax[i].set_yticks([0, 1], ["0", "max"])


def plot_all_covs(
    keys,
    binned_tx_depth,
    binned_tx,
    last_exon_r,
    n_bins=1000,
    ps=("good", "bad"),
    cs=SPLICE_MATCHES,
    output_file=None,
):
    fig, axs = plt.subplots(
        len(binned_tx_depth) * 2,
        len(cs),
        squeeze=False,
        sharey="row",
        figsize=(12, 3 * len(binned_tx_depth) + len(ps) * 0.5),
        gridspec_kw={"hspace": 0.05, "height_ratios": [2, 8] * len(binned_tx_depth)},
        layout="constrained",
    )

    for (b0, b1), (kde_ax, ax) in zip(
        itertools.pairwise(itertools.chain([0], sorted(binned_tx_depth))),
        zip(*(iter(axs),) * 2),
    ):
        plot_total_cov(ax, keys, binned_tx_depth[b1], ps=ps, cs=cs)

        plot_exon_ratio_kde(kde_ax[0], last_exon_r, binned_tx[b1], n_bins)
        plot_exon_ratio_kde(kde_ax[1], last_exon_r, binned_tx[b1], n_bins)

        for i, sm in enumerate(cs):
            if len(binned_tx_depth) > 1:
                kde_ax[i].set_title(f"{sm} {b0} < {b1} ({len(binned_tx[b1])})")
            else:
                kde_ax[i].set_title(f"{sm}")

            ax[i].margins(0.05, 0.05)

    handles, labels = ax[1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="outside lower center",
        ncol=3,
    )

    if output_file is not None:
        plt.savefig(output_file)

    plt.show()


### Methods and helpers added for review analysis

# =========================
# Stats + violin utilities
# =========================

def bh_adjust(p: Sequence[float]) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    out = np.empty_like(q)
    out[order] = q
    return out


def _p_to_stars(p: float) -> str:
    if np.isnan(p): return "n/a"
    if p < 1e-4:   return "****"
    if p < 1e-3:   return "***"
    if p < 1e-2:   return "**"
    if p < 0.05:   return "*"
    return "ns"


def test_vs_first(dists, method="brunnermunzel", alternative="two-sided"):
    """Return p-values for comparisons 2..N vs 1"""
    ref = np.asarray(dists[0], float)
    pvals = []
    for arr in dists[1:]:
        x = np.asarray(arr, float)
        if method == "brunnermunzel":
            stat, p = brunnermunzel(x, ref, alternative=alternative)
        elif method == "mannwhitney":
            # method="exact" only for small sizes & no ties; asymptotic is fine here
            stat, p = mannwhitneyu(x, ref, alternative=alternative, method="asymptotic")
        else:
            raise ValueError("method must be 'brunnermunzel' or 'mannwhitney'")
        pvals.append(p)
    return np.array(pvals)

def add_sig_markers(ax, dists_for_plot, qvals, ref_idx=0, line_color="#444",
                    height_frac=0.03, gap_frac=0.02, lw=1.2, fontsize=12):
    """
    Draw stacked brackets from ref_idx to each other group, labeled by q-values.
    dists_for_plot must be in the same scale as the y-axis (i.e., log10 if you plotted log10).
    """
    N = len(dists_for_plot)
    if N <= 1: return

    # y scale
    y_min = min(np.min(v) for v in dists_for_plot if len(v))
    y_max = max(np.max(v) for v in dists_for_plot if len(v))
    y_range = (y_max - y_min) if y_max > y_min else 1.0

    # starting height just above the tallest violin
    y = y_max + gap_frac * y_range
    top_seen = ax.get_ylim()[1]

    for j in range(1, N):
        x1, x2 = ref_idx + 1, j + 1
        h = height_frac * y_range

        # bracket
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], c=line_color, lw=lw)
        # text
        ax.text((x1 + x2) / 2, y + h, _p_to_stars(float(qvals[j-1])),
                ha="center", va="bottom", fontsize=fontsize, color=line_color)

        top_seen = max(top_seen, y + h)
        # stack next bracket higher
        y += (h + gap_frac * y_range)

    # ensure everything is visible
    ymin, ymax = ax.get_ylim()
    if top_seen > ymax:
        ax.set_ylim(ymin, top_seen + gap_frac * y_range)


def test_all_pairs(
    dists: Sequence[Sequence[float]],
    method: str = "brunnermunzel",
    alternative: str = "two-sided",
) -> tuple[list[tuple[int,int]], np.ndarray]:
    from scipy.stats import brunnermunzel, mannwhitneyu
    pairs = [(i, j) for i in range(len(dists)) for j in range(i + 1, len(dists))]
    pvals = []
    for i, j in pairs:
        x = np.asarray(dists[i], float)
        y = np.asarray(dists[j], float)
        if method == "brunnermunzel":
            _, p = brunnermunzel(x, y, alternative=alternative)
        elif method == "mannwhitney":
            _, p = mannwhitneyu(x, y, alternative=alternative, method="asymptotic")
        else:
            raise ValueError("method must be 'brunnermunzel' or 'mannwhitney'")
        pvals.append(p)
    return pairs, np.asarray(pvals, float)


def add_sig_for_pairs(
    ax: plt.Axes,
    dists_for_plot: Sequence[Sequence[float]],
    pairs: Sequence[tuple[int,int]],
    qvals: Sequence[float],
    *,
    line_color: str = "#444",
    height_frac: float = 0.03,
    gap_frac: float = 0.02,
    lw: float = 1.2,
    fontsize: float = 12,
) -> None:
    valid = [np.asarray(v) for v in dists_for_plot if len(v)]
    if not valid:
        return
    y_min = min(np.min(v) for v in valid)
    y_max = max(np.max(v) for v in valid)
    y_range = max(y_max - y_min, 1.0)

    y = y_max + gap_frac * y_range
    levels = []
    for _ in pairs:
        levels.append(y)
        y += (height_frac * y_range + gap_frac * y_range)

    for (i, j), q, y0 in zip(pairs, qvals, levels):
        x1, x2 = i + 1, j + 1
        h = height_frac * y_range
        ax.plot([x1, x1, x2, x2], [y0, y0 + h, y0 + h, y0], c=line_color, lw=lw)
        ax.text((x1 + x2) / 2, y0 + h, _p_to_stars(float(q)),
                ha="center", va="bottom", fontsize=fontsize, color=line_color)

    ymin, ymax = ax.get_ylim()
    if y > ymax:
        ax.set_ylim(ymin, y + gap_frac * y_range)


def plot_dists2(
    ax: plt.Axes,
    dists: Sequence[Sequence[float]],
    *,
    log: bool = False,
    colors: Sequence[str] | None = None,
    labels: Sequence[str] | None = None,
    title: str | None = None,
    clip_percentile: float = 99.9,
):
    if colors is None:
        colors = ["#006DB6"] * len(dists)

    cleaned, keep_idx = [], []
    for i, arr in enumerate(dists):
        a = np.asarray(arr, float)
        a = a[np.isfinite(a)]
        if log:
            a = a[a > 0]
            if a.size:
                a = np.log10(a)
        if a.size == 0:
            continue
        if a.size == 1:
            a = a + np.array([0.0, 1e-9])
        cleaned.append(a); keep_idx.append(i)

    if not cleaned:
        if title: ax.set_title(title)
        ax.text(0.5, 0.5, "No data > 0", ha="center", va="center", transform=ax.transAxes)
        return None

    p = ax.violinplot(cleaned, showmedians=True, widths=0.8)
    for c, body in zip([colors[i] for i in keep_idx], p["bodies"]):
        body.set_color(c); body.set_alpha(1.0)
    for part in ("cmaxes", "cmins", "cmedians"):
        p[part].set_color("#808080")
    p["cbars"].remove()

    if title: ax.set_title(title)
    if labels: ax.set_xticks(range(1, len(cleaned) + 1), [labels[i] for i in keep_idx])

    if log:
        all_logs = np.hstack(cleaned)
        e_min = int(np.floor(np.min(all_logs)))
        top = float(np.percentile(all_logs, clip_percentile))
        e_max_tick = int(np.floor(top))
        majors = np.arange(e_min, e_max_tick + 1)
        if majors.size:
            ax.set_yticks(majors); ax.set_yticklabels([f"$10^{e}$" for e in majors])
            minor = np.log10([v * 10**e for e in majors for v in range(2, 10)])
            minor = minor[minor < top]
            ax.set_yticks(minor, minor=True)
            ax.set_ylim(bottom=e_min, top=top * 1.02)
    return p


def plot_dists_vs_first(
    dists,
    labels,
    *,
    title="UMIs/cell",
    y_label="#UMIs",
    method="brunnermunzel",         # or "mannwhitney"
    alternative="two-sided",
    figsize=(16, 5),
    savepath: Path | None = None,
):
    """
    Log10-violin plot of multiple groups, with significance markers for (2..N) vs first.
    Returns (fig, ax). Saves if `savepath` is provided (slashes in filename sanitized).
    """
    # --- single pass clean so plotting & stats stay aligned ---
    kept_raw, kept_log, kept_labels = [], [], []
    for arr, lab in zip(dists, labels):
        a = np.asarray(arr, float)
        a = a[np.isfinite(a)]
        if a.size == 0:
            continue
        a_pos = a[a > 0]                 # for plotting on log scale
        if a_pos.size == 0:
            continue
        kept_raw.append(a)               # stats on raw
        kept_log.append(np.log10(a_pos)) # brackets on same scale as axis
        kept_labels.append(lab)

    if len(kept_raw) < 2:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.text(0.5, 0.5, "Not enough data to compare", ha="center", va="center", transform=ax.transAxes)
        if savepath is not None:
            fig.savefig(_sanitize_filename(Path(savepath)))
        return fig, ax

    # --- draw violins (directly, to avoid a second round of filtering) ---
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    vp = ax.violinplot(kept_log, showmedians=True, widths=0.8)
    for b in vp["bodies"]:
        b.set_alpha(1.0)
    for k in ("cmaxes", "cmins", "cmedians"):
        vp[k].set_color("#808080")
    vp["cbars"].remove()

    ax.set_title(title)
    ax.set_xticks(range(1, len(kept_labels) + 1), kept_labels)

    # --- nice log decade ticks & headroom for one row of stars ---
    all_logs = np.hstack(kept_log)
    e_min = int(np.floor(all_logs.min()))
    e_max_data = float(np.percentile(all_logs, 99.9))
    majors = np.arange(e_min, int(np.floor(e_max_data)) + 1)
    ax.set_yticks(majors)
    ax.set_yticklabels([f"$10^{e}$" for e in majors])
    minors = np.log10([v * 10**e for e in majors for v in range(2, 10)])
    minors = minors[minors < e_max_data]
    ax.set_yticks(minors, minor=True)
    ax.set_ylim(e_min, e_max_data * 1.02)   # initial; add_sig_markers may expand

    # --- stats: 2..N vs first (on raw), FDR, then brackets (on log) ---
    pvals = test_vs_first(kept_raw, method=method, alternative=alternative)
    qvals = bh_adjust(pvals)
    add_sig_markers(ax, kept_log, qvals, ref_idx=0)

    ax.set_ylabel(y_label)

    # NOTE: avoid fig.tight_layout() when using constrained layout styles
    if savepath is not None:
        fig.savefig(_sanitize_name(Path(savepath)), bbox_inches="tight")
    return fig, ax


# =========================
# Triptych (Fresh vs Fixed)
# =========================

def plot_fixfresh_triptych(
    hto7_c4_numis: Sequence[float],
    hto8_c4_numis: Sequence[float],
    sample_reads: dict,
    sample_read_lens: dict,
    *,
    figure_path=None,
    filename: str = "fig1h_fixfresh_comparison.svg",
):
    fig, ax = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: SR nUMIs (log) + sig
    d0_raw = [hto7_c4_numis, hto8_c4_numis]
    d0_plot = [np.log10(np.asarray(x)[np.asarray(x) > 0]) for x in d0_raw]
    plot_dists(ax[0], d0_raw, log=True, labels=["Fresh", "Fixed"], title="nUMIs per cell (Illumina)")
    pairs0, p0 = test_all_pairs(d0_raw, method="brunnermunzel", alternative="two-sided")
    q0 = bh_adjust(p0); add_sig_for_pairs(ax[0], d0_plot, pairs0, q0)
    ax[0].set_yticks([3, 4], ["$10^3$", "$10^4$"])
    ax[0].set_yticks(np.log10([v * 10**i for i in range(3, 5) for v in range(2, 10)][:-4]), minor=True)

    # Panel 2: PB UMI counts (log) + sig
    d1_raw = [[len(v) for v in sample_reads[k].values()] for k in sample_reads]
    d1_plot = [np.log10(np.asarray(x)[np.asarray(x) > 0]) for x in d1_raw]
    plot_dists(ax[1], d1_raw, log=True, labels=["Fresh", "Fixed"], title="nUMIs per cell (PB)")
    pairs1, p1 = test_all_pairs(d1_raw, method="brunnermunzel", alternative="two-sided")
    q1 = bh_adjust(p1); add_sig_for_pairs(ax[1], d1_plot, pairs1, q1)
    ax[1].set_yticks([2, 3, 4], ["$10^2$", "$10^3$", "$10^4$"])
    ax[1].set_yticks(np.log10([v * 10**i for i in range(1, 5) for v in range(2, 10)][1:-6]), minor=True)

    # Panel 3: PB read lengths (linear)
    plot_dists(
        ax[2],
        [[v for vs in sample_read_lens[k].values() for v in vs] for k in ("fresh_cd14", "fixed_cd14")],
        labels=["Fresh", "Fixed"],
        title="Read length (PB)",
    )
    ax[2].set_ylim(0, 800)

    fig.suptitle("CD14+ cluster stats")
    if figure_path is not None:
        fig.savefig(figure_path / filename)
    plt.show()
    return fig, ax

# =====================================
# Cluster-data prep for split violins
# =====================================

def rn_to_readumi_bcset(tagged_bam: str, bc_set: set[int]) -> dict[int, set[tuple[str, str]]]:
    bc_reads = defaultdict(set)
    with pysam.AlignmentFile(tagged_bam, "rb", check_sq=False, threads=2) as fh:
        for a in fh:
            if not a.has_tag("CB") or not a.has_tag("UB"):
                continue
            cb = a.get_tag("CB")
            if cb is None:
                continue
            if (bc := sequence_to_int(cb)) in bc_set:
                bc_reads[bc].add((a.query, a.get_tag("UB")))
    return bc_reads


def _rn_worker(bam_path: str, bc_set: set[int]) -> dict[int, set[tuple[str, str]]]:
    return rn_to_readumi_bcset(bam_path, bc_set)


def _row_sums_1d(mat):
    s = mat.sum(1)
    if hasattr(s, "todense"):   # sparse.numba_backend
        v = s.todense()
    elif hasattr(s, "A1"):      # SciPy CSR/CSC
        return s.A1
    else:
        v = np.asarray(s)
    return np.asarray(v).ravel()


def prepare_cluster_data(
    tagged_bams: Sequence[str],
    m_hto7, m_hto8,          # scipy sparse matrices (cells x genes), rows align with barcodes
    hto7_sc, hto8_sc,        # Scanpy objects with obs['leiden']
    hto7_barcodes: Sequence[int],
    hto8_barcodes: Sequence[int],
    MAPPING_PAIRS,           # dict{name:(fresh_cid|None,fixed_cid|None)} OR list[(f,x)]
    *,
    n_workers: int = 8,
) -> dict:
    sr_h7 = _row_sums_1d(m_hto7)
    sr_h8 = _row_sums_1d(m_hto8)

    idx7 = {int(bc): i for i, bc in enumerate(hto7_barcodes)}
    idx8 = {int(bc): i for i, bc in enumerate(hto8_barcodes)}

    def _bc_to_cluster_int(adata):
        lab = adata.obs['leiden'].astype(str).astype(int).to_numpy()
        bcs = adata.obs_names.to_numpy()
        d = {}
        for s, cid in zip(bcs, lab):
            if cid == -1:
                continue
            d[barcode_to_int(s.split('-')[0])] = int(cid)
        return d

    fresh_map = _bc_to_cluster_int(hto7_sc)
    fixed_map = _bc_to_cluster_int(hto8_sc)

    fresh_bc_by_cluster = defaultdict(set)
    for bc in hto7_barcodes:
        cid = fresh_map.get(int(bc))
        if cid is not None:
            fresh_bc_by_cluster[cid].add(int(bc))

    fixed_bc_by_cluster = defaultdict(set)
    for bc in hto8_barcodes:
        cid = fixed_map.get(int(bc))
        if cid is not None:
            fixed_bc_by_cluster[cid].add(int(bc))

    bc_union = set(hto7_barcodes) | set(hto8_barcodes)
    pb_reads = defaultdict(set)
    with ProcessPoolExecutor(max_workers=n_workers) as exc:
        for r in exc.map(_rn_worker, tagged_bams, itertools.repeat(bc_union)):
            for bc, pairs in r.items():
                pb_reads[bc].update(pairs)

    def _gather_sr(sr_vec, idx_map, bcs):
        if not bcs:
            return np.array([], dtype=float)
        ii = [idx_map[bc] for bc in bcs if bc in idx_map]
        return sr_vec[ii] if ii else np.array([], dtype=float)

    def _pb_umi_counts(bcs):
        return [len(pb_reads.get(bc, ())) for bc in bcs]

    def _pb_read_lengths(bcs):
        out = []
        for bc in bcs:
            for q, _ in pb_reads.get(bc, ()):
                if q is not None:
                    out.append(len(q))
        return out

    if isinstance(MAPPING_PAIRS, dict):
        mapping_items = list(MAPPING_PAIRS.items())
    else:
        def _name_for(f, x):
            return f"F{f if f is not None else '–'}<->X{x if x is not None else '–'}"
        mapping_items = [(_name_for(f, x), (f, x)) for (f, x) in MAPPING_PAIRS]

    out = {}
    for name, (f, x) in mapping_items:
        b7 = fresh_bc_by_cluster.get(int(f), set()) if f is not None else set()
        b8 = fixed_bc_by_cluster.get(int(x), set()) if x is not None else set()
        out[name] = {
            "illu": (_gather_sr(sr_h7, idx7, b7), _gather_sr(sr_h8, idx8, b8)),
            "pb_umi_counts": {"fresh": _pb_umi_counts(b7), "fixed": _pb_umi_counts(b8)},
            "pb_read_lens":  {"fresh": _pb_read_lengths(b7), "fixed": _pb_read_lengths(b8)},
            "meta": {"fresh_cid": f, "fixed_cid": x, "complete": (f is not None and x is not None)},
        }
    return out

def _sanitize_name(x) -> str:
    return re.sub(r"[\\/]+", "_", str(x))


def plot_cluster_triptychs(
    cluster_data: Dict[Any, Dict[str, Any]],
    *,
    outdir=None,
    readlen_ylim: Tuple[float, float] = (0, 2000),
):
    """
    For each entry in `cluster_data`, make a 1x3 triptych:
      [Illumina nUMIs (log, +sig), PB UMI counts (log, +sig), PB read length (linear)].
    Returns: dict[key] = (fig, ax_array)
    Saves to SVG if `outdir` is provided (filename sanitized).
    """
    results = {}

    def _sig_panel(ax, raw_groups):
        # log10 for plotting (drop non-positive)
        logged = [np.log10(np.asarray(g, float)[np.asarray(g, float) > 0]) for g in raw_groups]
        # violin
        plot_dists2(ax, raw_groups, log=True, labels=["Fresh", "Fixed"])
        # significance (two groups → one pair)
        try:
            pairs, p = test_all_pairs(raw_groups, method="brunnermunzel", alternative="two-sided")
            q = bh_adjust(p)
            add_sig_for_pairs(ax, logged, pairs, q)
        except Exception:
            pass  # if empty/degenerate, just skip the sig overlay

    for key, dat in cluster_data.items():
        # keep old behavior: skip sentinel -1 if keys are ints
        if isinstance(key, int) and key == -1:
            continue

        fig, ax = plt.subplots(1, 3, figsize=(20, 6))

        # Panel 1 — Illumina UMIs (log) + sig
        illu_fresh, illu_fixed = dat["illu"]
        _sig_panel(ax[0], [illu_fresh, illu_fixed])

        # Panel 2 — PB UMI counts (log) + sig
        pb_umis_f = dat["pb_umi_counts"]["fresh"]
        pb_umis_x = dat["pb_umi_counts"]["fixed"]
        _sig_panel(ax[1], [pb_umis_f, pb_umis_x])

        # Panel 3 — PB read lengths (linear)
        plot_dists2(
            ax[2],
            [dat["pb_read_lens"]["fresh"], dat["pb_read_lens"]["fixed"]],
            labels=["Fresh", "Fixed"],
        )
        ax[2].set_ylim(*readlen_ylim)

        fig.suptitle(f"Cluster {key} stats")
        fig.tight_layout()

        # optional save
        if outdir is not None:
            safe = _sanitize_name(key)
            fig.savefig(outdir / f"sens5_manualpairs_violins_cluster{safe}.svg")

        results[key] = (fig, ax)

    return results

# ==========================================
# Split-violin (seaborn) figures with stats
# ==========================================

def _complete_pairs(cluster_data: dict) -> list:
    out = []
    for k, v in cluster_data.items():
        meta = v.get("meta", {})
        if meta.get("complete") is True:
            out.append(k)
    if out:
        return out
    return [k for k in cluster_data.keys() if "–" not in k]


def _ensure_violin_samples(vals, log_transform=False, eps=1e-9):
    a = np.asarray(vals, dtype=float)
    a = a[np.isfinite(a)]
    if log_transform:
        a = a[a > 0]
        if a.size:
            a = np.log10(a)
    if a.size == 0:
        return None
    if a.size == 1:
        a = a + np.array([0.0, eps])
    return a


def _build_df(cluster_data: dict, metric_key: str, log_transform: bool):
    order = _complete_pairs(cluster_data)
    rows, kept_order = [], []
    for pair in order:
        dat = cluster_data[pair]
        if metric_key == "illu":
            fresh_vals, fixed_vals = dat["illu"]
        elif metric_key == "pb_umi_counts":
            fresh_vals = dat["pb_umi_counts"]["fresh"]; fixed_vals = dat["pb_umi_counts"]["fixed"]
        elif metric_key == "pb_read_lens":
            fresh_vals = dat["pb_read_lens"]["fresh"];  fixed_vals = dat["pb_read_lens"]["fixed"]
        else:
            raise ValueError(metric_key)

        f = _ensure_violin_samples(fresh_vals, log_transform)
        x = _ensure_violin_samples(fixed_vals, log_transform)
        if f is None or x is None:
            continue
        kept_order.append(pair)
        rows.extend({"pair": pair, "phase": "Fresh", "value": v} for v in f)
        rows.extend({"pair": pair, "phase": "Fixed", "value": v} for v in x)
    return pd.DataFrame(rows), kept_order


def _max_per_pair_percentile(df: pd.DataFrame, order: list, q: float):
    tops = []
    for p in order:
        v = df.loc[df["pair"] == p, "value"].to_numpy()
        if v.size:
            tops.append(np.percentile(v, q))
    return max(tops) if tops else None


def _pair_arrays_for_test(df: pd.DataFrame, pair: str, log_transform: bool):
    a = df.loc[(df["pair"] == pair) & (df["phase"] == "Fresh"), "value"].to_numpy()
    b = df.loc[(df["pair"] == pair) & (df["phase"] == "Fixed"), "value"].to_numpy()
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    return a, b


def _compute_pair_pvalues(df: pd.DataFrame, order: list, log_transform: bool):
    from scipy import stats
    pvals = {}
    for pair in order:
        a, b = _pair_arrays_for_test(df, pair, log_transform)
        if a.size == 0 or b.size == 0:
            pvals[pair] = np.nan; continue
        try:
            pvals[pair] = brunnermunzel(a, b, alternative="two-sided").pvalue
        except Exception:
            pvals[pair] = np.nan
    return pvals


def _format_log_axis(ax: plt.Axes, all_vals: np.ndarray, *, top_log=None, headroom_frac=0.15, tailroom_abs=0.15):
    if all_vals.size == 0:
        return
    min_log = float(np.min(all_vals))
    max_log = float(np.max(all_vals)) if top_log is None else float(top_log)
    e_min_tick = int(np.floor(min_log))
    e_max_tick = int(np.floor(max_log))
    majors = np.arange(e_min_tick, e_max_tick + 1)
    ax.set_yticks(majors); ax.set_yticklabels([f"$10^{e}$" for e in majors])
    minors = np.log10([v * 10**e for e in majors for v in range(2, 10)])
    ax.set_yticks(minors, minor=True)
    y_bottom = e_min_tick - tailroom_abs
    y_top_data = max_log
    y_range = max(y_top_data - y_bottom, 1e-6)
    y_headroom = headroom_frac * y_range
    ax.set_ylim(y_bottom, y_top_data + y_headroom)


def _annotate_sig_stars_row(ax: plt.Axes, df: pd.DataFrame, order: list, *, log_transform: bool, color="black",
                            bar_halfwidth=0.22, lw=1.3, fontsize=10):
    pvals = _compute_pair_pvalues(df, order, log_transform)
    xs = ax.get_xticks()
    y0, y1 = ax.get_ylim()
    y = y1 - 0.10 * (y1 - y0)
    dy = 0.02 * (y1 - y0)
    for x, pair in zip(xs, order):
        ax.plot([x - bar_halfwidth, x + bar_halfwidth], [y, y], color=color, lw=lw, solid_capstyle="butt")
        ax.text(x, y + dy, _p_to_stars(pvals[pair]), ha="center", va="bottom", fontsize=fontsize, color=color)


def _annotate_sig_stars_at_top(ax: plt.Axes, df: pd.DataFrame, order: list, *, log_transform: bool, color="black",
                               bar_halfwidth=0.22, lw=1.3, fontsize=11, y_frac=0.97, dy_frac=0.02):
    pvals = _compute_pair_pvalues(df, order, log_transform)
    xs = ax.get_xticks()
    trans = ax.get_xaxis_transform()
    for x, pair in zip(xs, order):
        ax.plot([x - bar_halfwidth, x + bar_halfwidth], [y_frac, y_frac],
                color=color, lw=lw, solid_capstyle="butt",
                transform=trans, clip_on=False, zorder=5)
        ax.text(x, y_frac + dy_frac, _p_to_stars(pvals[pair]),
                ha="center", va="bottom", fontsize=fontsize, color=color,
                transform=trans, clip_on=False, zorder=6)


def make_fig_short_reads(cluster_data: dict, title: str = "nUMIs per cell (Illumina)"):
    df, order = _build_df(cluster_data, metric_key="illu", log_transform=True)
    if df.empty:
        raise ValueError("No complete pairs with Illumina data.")
    fig, ax = plt.subplots(figsize=(max(8, 0.8 * len(order) + 2), 6))
    sns.violinplot(
        data=df, x="pair", y="value", hue="phase",
        split=True, inner="quartile", order=order, ax=ax, cut=0.2
    )
    ax.set_title(title); ax.set_xlabel("Cluster pair (Fresh<->Fixed)"); ax.set_ylabel(""); ax.legend(title="Phase")
    _format_log_axis(ax, df["value"].to_numpy(), top_log=df["value"].max(),
                     headroom_frac=0.15, tailroom_abs=0.15)
    _annotate_sig_stars_row(ax, df, order, log_transform=True, color="black")
    for tick in ax.get_xticklabels():
        tick.set_rotation(45); tick.set_ha("right")
    fig.tight_layout()
    return fig, ax


def make_fig_long_reads(
    cluster_data: dict,
    title_top: str = "UMIs per cell (PacBio)",
    title_bottom: str = "Read length (PacBio)",
):
    df_um, order_um = _build_df(cluster_data, metric_key="pb_umi_counts", log_transform=True)
    df_len, order_len = _build_df(cluster_data, metric_key="pb_read_lens",  log_transform=False)
    if df_um.empty and df_len.empty:
        raise ValueError("No complete pairs with PacBio data.")
    order = order_um if len(order_um) >= len(order_len) else order_len

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(10, 1.0 * len(order) + 4), 10))

    if not df_um.empty:
        sns.violinplot(
            data=df_um, x="pair", y="value", hue="phase",
            split=True, inner="quartile", order=order, ax=ax1, cut=0.2
        )
        ax1.set_title(title_top); ax1.set_xlabel(""); ax1.set_ylabel(""); ax1.legend(title="Phase")
        _format_log_axis(ax1, df_um["value"].to_numpy(), top_log=df_um["value"].max(),
                         headroom_frac=0.15, tailroom_abs=0.15)
        _annotate_sig_stars_row(ax1, df_um, order, log_transform=True, color="black")
        for t in ax1.get_xticklabels(): t.set_rotation(45); t.set_ha("right")
    else:
        ax1.axis("off"); ax1.set_title(title_top + " (no data)")

    if not df_len.empty:
        sns.violinplot(
            data=df_len, x="pair", y="value", hue="phase",
            split=True, inner="quartile", order=order, ax=ax2
        )
        ax2.set_title(title_bottom); ax2.set_xlabel("Cluster pair (Fresh<->Fixed)"); ax2.set_ylabel(""); ax2.legend(title="Phase")
        top_lin = _max_per_pair_percentile(df_len, order, q=99.0)
        if top_lin is not None:
            ax2.set_ylim(bottom=0, top=top_lin * 1.02)
        _annotate_sig_stars_at_top(ax2, df_len, order, log_transform=False, color="black")
        for t in ax2.get_xticklabels(): t.set_rotation(45); t.set_ha("right")
    else:
        ax2.axis("off"); ax2.set_title(title_bottom + " (no data)")

    fig.suptitle("PacBio metrics by mapped cluster pairs", y=0.98)
    fig.tight_layout()
    return fig, (ax1, ax2)




