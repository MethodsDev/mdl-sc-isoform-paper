from typing import Dict, Tuple, Any, Optional, Iterable
import pysam
import pandas as pd
from collections import defaultdict
import pyranges as pr
import logging

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages


__all__ = [
    "compare_BAM_annot_PAS",
    "plot_stacked_barplot_pas_site_good_matches",
    "plot_stacked_barplot_pas_site_non_good_matches",
    "save_overview_good_vs_non_pdf",
    "save_stacked_bars_per_dataset_pdf",
    "save_all_stackedbars_pas_pdf_reports",
]

logger = logging.getLogger(__name__)

# --- prepare 20‑unit bins ---
bins       = [0, 20, 40, 60, 80, 100]
bin_labels = pd.IntervalIndex.from_breaks(bins, closed='right')

# --- define good priming categories ---
good_cats = ["ANNO_GPA", "GOOD", "TX_PAS", "TX_MOTIF"]

# List all categories (in the desired order)
categories = [
    "ANNO_GPA", "GOOD", "TX_GPA_PAS", "TX_PAS", "TX_GPA_MOTIF", "TX_MOTIF",
    "TX_GPA", "TX_UNK",
    "CDS_GPA", "CDS_OTHER",
    "AS_TX_GPA_NC", "NC_GPA", "AS_TX_NC", "NC_OTHER",
    "AS_TX_GPA", "INTERGENIC_GPA", "AS_TX", "INTERGENIC_OTHER",
    "MITO"
]

# Manual color palette mapping.
color_map = {
    "ANNO_GPA":         "#1F618D",
    "GOOD":             "#1C9099",
    "TX_GPA_PAS":       "#74add1",
    "TX_PAS":           "#53A454",
    "TX_GPA_MOTIF":     "#85929E",
    "TX_MOTIF":         "#9B59B6",
    "TX_GPA":           "#b8e186",
    "TX_UNK":           "#fed976",

    "CDS_GPA":          "#e7298a",
    "CDS_OTHER":        "#f4a582",

    "AS_TX_GPA_NC":     "#F6222E",
    "NC_GPA":           "#FF6B6B",
    "AS_TX_NC":         "#FF6347",
    "NC_OTHER":         "#ED4C67",

    "AS_TX_GPA":        "#A93226",
    "INTERGENIC_GPA":   "#C44D58",
    "AS_TX":            "#DD7788",

    "INTERGENIC_OTHER": "#D980FA",
    "MITO":             "#85929E"
}

# --- mapping for legend labels ---
legend_names = {
    "TE": "TE — terminal exon",
    "EX": "EX — exonic",
    "IN": "IN — intronic",
    "AL": "AL — alternative (exonic/intronic)",
    "DI": "DI — downstream intergenic",
    "UI": "UI — upstream intergenic",
    "NA": "NA — not available (ambiguous)",
    "nan": "nan"
}


# Manual color for each genomic class
class_color_map = {
    "TE":  "#1F618D",
    "EX":  "#53A454",
    "IN":  "#E67E22",
    "AL":  "#8E44AD",
    "DI":  "#47E3FF",
    "UI":  "#F6222E",
    "NA":  "#7F8C8D",
    "nan": "#BDC3C7"
}

# --- shade color according to bin ---
def shade_color(color, frac):
    r, g, b = color
    return (
        frac * r + (1 - frac),
        frac * g + (1 - frac),
        frac * b + (1 - frac),
    )

# Convert hex → RGB tuples
rgb_palette = {cat: to_rgb(hexc) for cat, hexc in color_map.items()}



# Load PAS sites into two PyRanges (one per strand)
column_names = [
    'Chromosome','Start','End','pas_id','expression',
    'strand','tissue_support','protocol_count','stringency',
    'genomic_class','polyA_signal'
]





# Chunked BAM processor
def compare_BAM_annot_PAS(bam_path, pas_plus, pas_minus, margin=5, chunk_size=200_000):
    sam = pysam.AlignmentFile(bam_path, "rb")

    # overall tallies
    category_counts  = defaultdict(lambda: {'matches':0,'non_matches':0})
    # per-category × bin → matched count
    bin_match_counts = defaultdict(int)
    # per-category × bin × genomic_class → matched count
    class_counts     = defaultdict(int)

    failed_chunks = 0
    buffer = []

    def flush_buffer(buf):
        nonlocal failed_chunks
        if not buf:
            return

        df_reads = pd.DataFrame(buf)
        chroms   = df_reads['Chromosome'].unique().tolist()
        logger.info(f"Flushing chunk size={len(buf)} contigs={chroms}")

        try:
            for strand, pas in (('+', pas_plus), ('-', pas_minus)):
                sub = df_reads[df_reads.strand == strand]
                if sub.empty:
                    continue

                gr     = pr.PyRanges(sub[[
                    'Chromosome','Start','End','XC_tag','read_idx'
                ]])
                joined = gr.join(pas)

                # if no overlaps → all are non-matches
                if joined.df.empty:
                    tot = sub['XC_tag'].value_counts()
                    for tag, t in tot.items():
                        category_counts[tag]['non_matches'] += t
                    continue

                # dedupe so each read counts once
                dfj = joined.df.drop_duplicates(subset=['read_idx'])

                # overall matches / non-matches
                matched = dfj['XC_tag'].value_counts()
                total   = sub['XC_tag'].value_counts()
                for tag, t in total.items():
                    m = int(matched.get(tag, 0))
                    category_counts[tag]['matches']     += m
                    category_counts[tag]['non_matches'] += (t - m)

                # now bin & class breakdown
                # joined.df already has 'bin' & 'genomic_class'
                for _, row in dfj.iterrows():
                    tag = row['XC_tag']
                    b   = row['bin']
                    gc  = row['genomic_class']
                    bin_match_counts[(tag,b)]      += 1
                    class_counts[(tag,b,gc)]       += 1

        except Exception as e:
            failed_chunks += 1
            logger.error(f"Error in flush_buffer for contigs {chroms}: {e}")
        finally:
            buf.clear()

    # stream through the BAM, tagging each read
    for i, read in enumerate(sam.fetch(), start=1):
        if not read.has_tag('XC'):
            continue
        xc    = read.get_tag('XC')
        end   = read.reference_start if read.is_reverse else read.reference_end
        chrom = sam.get_reference_name(read.reference_id)

        buffer.append({
            'Chromosome': chrom,
            'Start':       max(0, end - margin),
            'End':         end + margin,
            'XC_tag':      xc,
            'strand':      '-' if read.is_reverse else '+',
            'read_idx':    i,
        })

        if i % chunk_size == 0:
            logger.info(f"Processed {i:,} reads — flushing")
            flush_buffer(buffer)

    # final flush
    flush_buffer(buffer)
    logger.info(f"Done: {i:,} reads, {failed_chunks} chunk errors")

    return category_counts, bin_match_counts, class_counts



# stacked barplots
## internal helpers
def _collect_classes(class_bin_counts: Dict[Tuple[str, pd.Interval, Any], int],
                     cats_to_plot: Iterable[str]) -> List[Tuple[Any, str]]:
    """Return sorted list of (orig_gc, label) pairs found in the selected categories."""
    raw_gcs = {
        gc for (cat, _, gc) in class_bin_counts.keys() if cat in set(cats_to_plot)
    }
    classes = []
    for gc in raw_gcs:
        label = "nan" if pd.isna(gc) else str(gc)
        classes.append((gc, label))
    classes.sort(key=lambda pair: pair[1])
    return classes

def _make_base_colors(classes: List[Tuple[Any, str]], cmap_name: Optional[str]) -> Dict[str, Tuple[float,float,float]]:
    """Pick a base RGB for each class label, preferring class_color_map then falling back to a cmap."""
    if cmap_name is not None:
        cmap = plt.get_cmap(cmap_name)
    else:
        cmap = plt.get_cmap("tab10" if len(classes) <= 10 else "hsv")

    fallback = {label: to_rgb(cmap(i % cmap.N)) for i, (_, label) in enumerate(classes)}
    base = {
        label: to_rgb(class_color_map[label]) if label in class_color_map else fallback[label]
        for _, label in classes
    }
    return base

def _plot_pas_core(
    *,
    results: Dict[str, Dict[str, Any]],
    dataset_key: str,
    cats_to_plot: List[str],
    title: str,
    ax: Optional[plt.Axes],
    figsize: Optional[Tuple[float, float]],
    width: float,
    legend_outside: bool,
    non_match_color: str,
    cmap_name: Optional[str],
) -> Tuple[plt.Figure, plt.Axes]:

    if dataset_key not in results:
        raise KeyError(f"'{dataset_key}' not in results. Available: {list(results.keys())}")

    data = results[dataset_key]
    if "category_counts" not in data or "class_bin_counts" not in data:
        raise KeyError("results[dataset_key] must include 'category_counts' and 'class_bin_counts'.")

    category_counts: Dict[str, Dict[str, Any]] = data["category_counts"]
    class_bin_counts: Dict[Tuple[str, pd.Interval, Any], int] = data["class_bin_counts"]

    # classes present in the selected categories
    classes = _collect_classes(class_bin_counts, cats_to_plot)
    base_colors = _make_base_colors(classes, cmap_name)

    # figure/axes
    if ax is None:
        if figsize is None:
            # dynamic width: ~0.6 per bar, min 8
            figsize = (max(8.0, 0.6 * max(1, len(cats_to_plot))), 6.0)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    # bars
    x = range(len(cats_to_plot))
    for xi, cat in zip(x, cats_to_plot):
        bottom = 0
        # matches by class & bin
        for orig_gc, label in classes:
            base = base_colors[label]
            for idx, bin_iv in enumerate(bin_labels):
                cnt = class_bin_counts.get((cat, bin_iv, orig_gc), 0)
                if cnt == 0:
                    continue
                frac = (idx + 1) / len(bin_labels)
                col = shade_color(base, frac)
                ax.bar(xi, cnt, width, bottom=bottom, color=col, edgecolor="none")
                bottom += cnt

        # non-matches
        nm = category_counts.get(cat, {}).get("non_matches", 0)
        if nm:
            ax.bar(xi, nm, width, bottom=bottom, color=non_match_color, edgecolor="none")

    # cosmetics
    ax.set_xticks(list(x))
    ax.set_xticklabels(cats_to_plot, rotation=45, ha="right")
    ax.set_ylabel("Read count")
    ax.set_title(title)
    fig.tight_layout()

    # legend
    handles = [Patch(facecolor=base_colors[label], label=legend_names.get(label, label))
               for _, label in classes]
    handles.append(Patch(facecolor=non_match_color, label="Non-match"))

    if legend_outside:
        ax.legend(handles=handles, title="Genomic class",
                  bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
        fig.tight_layout()
    else:
        ax.legend(handles=handles, title="Genomic class")

    return fig, ax



def plot_stacked_barplot_pas_site_good_matches(
    results: Dict[str, Dict[str, Any]],
    dataset_key: str,
    *,
    ax: Optional[plt.Axes] = None,
    figsize: Optional[Tuple[float, float]] = (8, 6),
    title: str = "PAS-site matches for good-priming categories",
    width: float = 0.6,
    legend_outside: bool = True,
    non_match_color: str = "grey",
    cmap_name: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    GOOD categories version (kept for backward compatibility).
    """
    return _plot_pas_core(
        results=results,
        dataset_key=dataset_key,
        cats_to_plot=list(good_cats),
        title=title,
        ax=ax,
        figsize=figsize,
        width=width,
        legend_outside=legend_outside,
        non_match_color=non_match_color,
        cmap_name=cmap_name,
    )

def plot_stacked_barplot_pas_site_non_good_matches(
    results: Dict[str, Dict[str, Any]],
    dataset_key: str,
    *,
    ax: Optional[plt.Axes] = None,
    figsize: Optional[Tuple[float, float]] = None,  # dynamic default
    title: str = "PAS-site matches for non-good priming categories",
    width: float = 0.6,
    legend_outside: bool = True,
    non_match_color: str = "grey",
    cmap_name: Optional[str] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    NON-GOOD categories version. Figuresize defaults to dynamic width.
    """
    data = results[dataset_key]
    category_counts: Dict[str, Dict[str, Any]] = data["category_counts"]

    # keep preferred ordering if `categories` global is present, else use keys()
    all_present = [c for c in (globals().get("categories") or category_counts.keys())
                   if c in category_counts]
    non_good = [c for c in all_present if c not in set(good_cats)]

    return _plot_pas_core(
        results=results,
        dataset_key=dataset_key,
        cats_to_plot=non_good,
        title=title,
        ax=ax,
        figsize=figsize,  # None → dynamic inside core
        width=width,
        legend_outside=legend_outside,
        non_match_color=non_match_color,
        cmap_name=cmap_name,
    )

# ---------- shared helpers for PDFs ----------

def _non_good_categories_for_dataset(category_counts: Dict[str, Dict[str, Any]]) -> List[str]:
    """Return non-good categories present in dataset, preserving your global 'categories' order if available."""
    present = list(category_counts.keys())
    ordered_all = [c for c in (globals().get("categories") or present) if c in present]
    return [c for c in ordered_all if c not in set(good_cats)]

def _compute_global_maxima(results: Dict[str, Dict[str, Any]]) -> Tuple[int, int]:
    """
    Compute shared Y-axis maxima across datasets for good and non-good category panels.
    Uses class_bin_counts + non_matches to get total per category.
    """
    global_max_good = 0
    global_max_non  = 0

    for data in results.values():
        cc  = data["category_counts"]
        cbc = data["class_bin_counts"]

        # GOOD cats
        for cat in good_cats:
            if cat not in cc:
                continue
            # sum all matched counts for this category across bins/classes
            matched = sum(cnt for (c, _iv, _gc), cnt in cbc.items() if c == cat)
            total   = matched + int(cc[cat].get("non_matches", 0))
            if total > global_max_good:
                global_max_good = total

        # NON-GOOD cats (present in this dataset)
        for cat in _non_good_categories_for_dataset(cc):
            matched = sum(cnt for (c, _iv, _gc), cnt in cbc.items() if c == cat)
            total   = matched + int(cc[cat].get("non_matches", 0))
            if total > global_max_non:
                global_max_non = total

    return global_max_good, global_max_non

def _non_good_categories_for_dataset(category_counts: Dict[str, Dict[str, Any]]) -> list:
    present = list(category_counts.keys())
    ordered_all = [c for c in (globals().get("categories") or present) if c in present]
    return [c for c in ordered_all if c not in set(good_cats)]

def _classes_for_cats(cbc: Dict[Tuple[str, pd.Interval, Any], int], cats: Iterable[str]) -> Tuple[list, list]:
    """Return (classes, class_labels) for the selected categories."""
    raw = {gc for (c, _iv, gc) in cbc.keys() if c in set(cats)}
    classes = sorted(raw, key=lambda x: str(x))
    labels  = ["nan" if pd.isna(gc) else str(gc) for gc in classes]
    return classes, labels

# main plotting

def save_overview_good_vs_non_pdf(
    results: Dict[str, Dict[str, Any]],
    output_path: str,
    *,
    dataset_order: Optional[List[str]] = None,
    figsize: Tuple[float, float] = (8, 6),
    colors: Optional[Dict[str, str]] = None,
    rotation: int = 45,
) -> None:
    """
    Save a single-page PDF with side-by-side bars (Good vs Non-good) for each dataset.
    """
    if colors is None:
        colors = {"Good": "#1f77b4", "Non-good": "#ff7f0e"}

    datasets = dataset_order or list(results.keys())

    good_totals: List[int] = []
    non_totals:  List[int] = []

    for name in datasets:
        data = results[name]
        cc   = data["category_counts"]

        good_sum = 0
        for cat in good_cats:
            d = cc.get(cat)
            if d is not None:
                good_sum += int(d.get("matches", 0)) + int(d.get("non_matches", 0))
        good_totals.append(good_sum)

        non_sum = 0
        for cat in _non_good_categories_for_dataset(cc):
            d = cc.get(cat)
            if d is not None:
                non_sum += int(d.get("matches", 0)) + int(d.get("non_matches", 0))
        non_totals.append(non_sum)

    x     = list(range(len(datasets)))
    width = 0.35

    with PdfPages(output_path) as pdf:
        fig, ax = plt.subplots(figsize=figsize)
        ax.bar([xi - width/2 for xi in x], good_totals, width, label="Good",     color=colors["Good"])
        ax.bar([xi + width/2 for xi in x], non_totals,  width, label="Non-good", color=colors["Non-good"])

        ax.set_xticks(x)
        ax.set_xticklabels(datasets, rotation=rotation, ha="right")
        ax.set_ylabel("Read count")
        ax.set_title("Good vs Non-good by dataset")
        ax.legend()
        fig.tight_layout()

        pdf.savefig(fig)
        plt.close(fig)


def save_stacked_bars_per_dataset_pdf(
    results: Dict[str, Dict[str, Any]],
    output_path: str,
    *,
    include_independent: bool = True,
    include_shared: bool = True,
    shared_pad: float = 1.05,  # multiply global max by this
    width: float = 0.6,
    legend_outside: bool = True,
    non_match_color: str = "grey",
    cmap_name: Optional[str] = None,
) -> None:
    """
    For each dataset, save up to four pages:
      1) Good (independent y)
      2) Non-good (independent y)
      3) Good (shared y across datasets)
      4) Non-good (shared y across datasets)
    """

    if not (include_independent or include_shared):
        raise ValueError("At least one of include_independent/include_shared must be True.")

    gmax, nmax = _compute_global_maxima(results)

    with PdfPages(output_path) as pdf:
        for name, data in results.items():
            cc = data["category_counts"]
            non_good = _non_good_categories_for_dataset(cc)

            # --- independent ---
            if include_independent:
                # Good
                fig, ax = _plot_pas_core(
                    results=results,
                    dataset_key=name,
                    cats_to_plot=list(good_cats),
                    title=f"{name} — Good (independent)",
                    ax=None,
                    figsize=None,           # dynamic width: ~0.6 per bar, min 8
                    width=width,
                    legend_outside=legend_outside,
                    non_match_color=non_match_color,
                    cmap_name=cmap_name,
                )
                pdf.savefig(fig); plt.close(fig)

                # Non-good
                fig, ax = _plot_pas_core(
                    results=results,
                    dataset_key=name,
                    cats_to_plot=non_good,
                    title=f"{name} — Non-good (independent)",
                    ax=None,
                    figsize=None,  # dynamic width by number of categories
                    width=width,
                    legend_outside=legend_outside,
                    non_match_color=non_match_color,
                    cmap_name=cmap_name,
                )
                pdf.savefig(fig); plt.close(fig)

            # --- shared y-lims across datasets ---
            if include_shared:
                # Good (shared)
                fig, ax = _plot_pas_core(
                    results=results,
                    dataset_key=name,
                    cats_to_plot=list(good_cats),
                    title=f"{name} — Good (shared)",
                    ax=None,
                    figsize=None,
                    width=width,
                    legend_outside=legend_outside,
                    non_match_color=non_match_color,
                    cmap_name=cmap_name,
                )
                ax.set_ylim(0, gmax * shared_pad)
                fig.tight_layout()
                pdf.savefig(fig); plt.close(fig)

                # Non-good (shared)
                fig, ax = _plot_pas_core(
                    results=results,
                    dataset_key=name,
                    cats_to_plot=non_good,
                    title=f"{name} — Non-good (shared)",
                    ax=None,
                    figsize=None,
                    width=width,
                    legend_outside=legend_outside,
                    non_match_color=non_match_color,
                    cmap_name=cmap_name,
                )
                ax.set_ylim(0, nmax * shared_pad)
                fig.tight_layout()
                pdf.savefig(fig); plt.close(fig)


def save_all_stackedbars_pas_pdf_reports(
    results: Dict[str, Dict[str, Any]],
    overview_pdf_path: str,
    stacked_pdf_path: str,
    *,
    dataset_order: Optional[List[str]] = None,
) -> None:
    """
    Generate both:
      - Overview (Good vs Non-good by dataset) PDF
      - Per-dataset stacked bars PDF (independent + shared y)
    """
    save_overview_good_vs_non_pdf(results, overview_pdf_path, dataset_order=dataset_order)
    save_stacked_bars_per_dataset_pdf(results, stacked_pdf_path)


def save_pie_global_nested_pdf(
    results: Dict[str, Dict[str, Any]],
    output_path: str,
    *,
    dataset_order: Optional[list] = None,
    figsize: Tuple[float, float] = (10, 6),
    good_label: str = "Good",
    non_good_label: str = "Non-good",
    inner_color_good: Optional[Tuple[float, float, float]] = None,
    inner_color_non_good: Optional[Tuple[float, float, float]] = None,
) -> None:
    """
    For each dataset, draw a *single* nested donut:
      - inner: Good vs Non-good totals
      - outer: categories (good first, then non-good), colored by `rgb_palette`
    """
    # sensible defaults matching your notebook
    if inner_color_good is None:
        inner_color_good = rgb_palette.get("GOOD", to_rgb("#1C9099"))
    if inner_color_non_good is None:
        inner_color_non_good = rgb_palette.get("TX_GPA_MOTIF", to_rgb("#85929E"))

    datasets = dataset_order or list(results.keys())

    with PdfPages(output_path) as pdf:
        for name in datasets:
            data = results[name]
            cc   = data["category_counts"]

            # categories in desired order (good first, then non-good)
            non_good = _non_good_categories_for_dataset(cc)
            all_ordered = list(good_cats) + non_good

            # inner sizes/colors (Good vs Non-good)
            inner_good = sum(int(cc.get(c, {}).get("matches", 0)) + int(cc.get(c, {}).get("non_matches", 0))
                             for c in good_cats)
            inner_non  = sum(int(cc.get(c, {}).get("matches", 0)) + int(cc.get(c, {}).get("non_matches", 0))
                             for c in non_good)
            inner_sizes  = [inner_good, inner_non]
            inner_colors = [inner_color_good, inner_color_non_good]

            # outer sizes/colors (per category)
            outer_sizes  = [int(cc.get(c, {}).get("matches", 0)) + int(cc.get(c, {}).get("non_matches", 0))
                            for c in all_ordered]
            outer_colors = [rgb_palette.get(c, to_rgb("#BBBBBB")) for c in all_ordered]

            fig, ax = plt.subplots(figsize=figsize)
            # adjust margins: legend on right, title a bit higher
            fig.subplots_adjust(left=0.05, right=0.75, top=0.8)

            # inner ring
            ax.pie(
                inner_sizes,
                radius=1.0,
                colors=inner_colors,
                labels=[good_label, non_good_label],
                labeldistance=0.7,
                textprops={'fontsize': 9},
                wedgeprops=dict(width=0.3, edgecolor='w'),
            )

            # outer ring
            ax.pie(
                outer_sizes,
                radius=1.3,
                colors=outer_colors,
                labels=all_ordered,
                labeldistance=1.05,
                textprops={'fontsize': 7},
                wedgeprops=dict(width=0.3, edgecolor='w'),
            )

            ax.set_title(f"{name}: Global nested categories", y=1.05, fontsize=12)

            # legend outside on the right
            handles = [Patch(facecolor=rgb_palette.get(c, to_rgb("#BBBBBB")), label=c) for c in all_ordered]
            ax.legend(
                handles=handles,
                title="Category",
                bbox_to_anchor=(1.01, 0.5),
                loc="center left",
                fontsize=7,
                title_fontsize=8,
            )

            pdf.savefig(fig)
            plt.close(fig)


def save_pie_nested_by_class_pdf(
    results: Dict[str, Dict[str, Any]],
    output_path: str,
    *,
    dataset_order: Optional[list] = None,
    figsize: Tuple[float, float] = (8, 6),
) -> None:
    """
    For each dataset, write TWO pages:
      1) Good categories (inner: categories in color)
      2) Non-good categories (inner: categories in color)
    Outer ring = class × bin slices (shaded), plus Non-match per category.
    """
    def _plot_nested(data: Dict[str, Any], cats: list, title: str, pdf: PdfPages) -> None:
        cc  = data["category_counts"]
        cbc = data["class_bin_counts"]

        # inner ring totals/colors
        inner_sizes  = [int(cc.get(c, {}).get("matches", 0)) + int(cc.get(c, {}).get("non_matches", 0)) for c in cats]
        inner_colors = [rgb_palette.get(c, to_rgb("#BBBBBB")) for c in cats]

        # classes (order by label)
        classes, class_labels = _classes_for_cats(cbc, cats)
        class_base = {lab: to_rgb(class_color_map.get(lab, "#888888")) for lab in class_labels}

        # outer ring: class×bin (shaded) + Non-match
        outer_sizes, outer_colors = [], []
        for c in cats:
            for orig_gc, lab in zip(classes, class_labels):
                base = class_base[lab]
                for idx, iv in enumerate(bin_labels):
                    cnt = int(cbc.get((c, iv, orig_gc), 0))
                    if not cnt:
                        continue
                    outer_sizes.append(cnt)
                    outer_colors.append(shade_color(base, (idx + 1) / len(bin_labels)))
            # Non-match slice per category
            nm = int(cc.get(c, {}).get("non_matches", 0))
            if nm:
                outer_sizes.append(nm)
                outer_colors.append((0.5, 0.5, 0.5))

        fig, ax = plt.subplots(figsize=figsize)
        ax.pie(
            inner_sizes,
            radius=1.0,
            colors=inner_colors,
            labels=cats,
            labeldistance=0.7,
            textprops={'fontsize': 8},
            wedgeprops=dict(width=0.3, edgecolor='w'),
        )
        ax.pie(
            outer_sizes,
            radius=1.3,
            colors=outer_colors,
            labels=[None] * len(outer_sizes),
            wedgeprops=dict(width=0.3, edgecolor='none'),
        )
        ax.set_title(title, y=1.02, fontsize=10)

        # legend
        handles = [Patch(facecolor=class_base[lab], label=legend_names.get(lab, lab)) for lab in class_labels]
        handles.append(Patch(facecolor=(0.5, 0.5, 0.5), label="Non-match"))
        ax.legend(
            handles=handles,
            title="Genomic class",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            fontsize=7,
            title_fontsize=8,
        )
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    datasets = dataset_order or list(results.keys())
    with PdfPages(output_path) as pdf:
        for name in datasets:
            data = results[name]
            cc   = data["category_counts"]
            non_good = _non_good_categories_for_dataset(cc)

            _plot_nested(data, list(good_cats), f"{name} — Good categories", pdf)
            _plot_nested(data, non_good,       f"{name} — Non-good categories", pdf)


def save_pie_nested_by_class_greyscale_pdf(
    results: Dict[str, Dict[str, Any]],
    output_path: str,
    *,
    dataset_order: Optional[list] = None,
    figsize: Tuple[float, float] = (8, 6),
    grey_min: float = 0.3,
    grey_max: float = 0.7,
) -> None:
    """
    For each dataset, write TWO pages:
      1) Good categories (inner = greyscale ramp)
      2) Non-good categories (inner = greyscale ramp)
    Outer ring = class × bin (shaded) + Non-match.
    """
    def _plot_nested_grey(data: Dict[str, Any], cats: list, title: str, pdf: PdfPages) -> None:
        cc  = data["category_counts"]
        cbc = data["class_bin_counts"]

        inner_sizes = [int(cc.get(c, {}).get("matches", 0)) + int(cc.get(c, {}).get("non_matches", 0)) for c in cats]
        values = np.linspace(grey_max, grey_min, max(1, len(cats)))
        inner_colors = [(v, v, v) for v in values]

        classes, class_labels = _classes_for_cats(cbc, cats)
        class_base = {lab: to_rgb(class_color_map.get(lab, "#888888")) for lab in class_labels}

        outer_sizes, outer_colors = [], []
        for c in cats:
            for orig_gc, lab in zip(classes, class_labels):
                base = class_base[lab]
                for idx, iv in enumerate(bin_labels):
                    cnt = int(cbc.get((c, iv, orig_gc), 0))
                    if not cnt:
                        continue
                    outer_sizes.append(cnt)
                    outer_colors.append(shade_color(base, (idx + 1) / len(bin_labels)))
            nm = int(cc.get(c, {}).get("non_matches", 0))
            if nm:
                outer_sizes.append(nm)
                outer_colors.append((0.5, 0.5, 0.5))

        fig, ax = plt.subplots(figsize=figsize)
        ax.pie(
            inner_sizes,
            radius=1.0,
            colors=inner_colors,
            labels=cats,
            labeldistance=0.7,
            textprops={'fontsize': 8},
            wedgeprops=dict(width=0.3, edgecolor='white'),
        )
        ax.pie(
            outer_sizes,
            radius=1.3,
            colors=outer_colors,
            labels=[None] * len(outer_sizes),
            wedgeprops=dict(width=0.3, edgecolor='none'),
        )
        ax.set_title(title, y=1.02, fontsize=10)

        handles = [Patch(facecolor=class_base[lab], label=legend_names.get(lab, lab)) for lab in class_labels]
        handles.append(Patch(facecolor=(0.5, 0.5, 0.5), label="Non-match"))
        ax.legend(
            handles=handles,
            title="Genomic class",
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            fontsize=7,
            title_fontsize=8,
        )
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    datasets = dataset_order or list(results.keys())
    with PdfPages(output_path) as pdf:
        for name in datasets:
            data = results[name]
            cc   = data["category_counts"]
            non_good = _non_good_categories_for_dataset(cc)

            _plot_nested_grey(data, list(good_cats), f"{name} — Good categories", pdf)
            _plot_nested_grey(data, non_good,       f"{name} — Non-good categories", pdf)



