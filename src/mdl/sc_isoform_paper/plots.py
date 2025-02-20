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
from scipy.stats import gaussian_kde

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
