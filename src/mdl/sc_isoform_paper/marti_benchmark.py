from matplotlib import pyplot as plt, cm
from matplotlib.patches import Rectangle
import numpy as np
import gzip
import re
from pathlib import Path
from collections import defaultdict, Counter
import pandas as pd
from upsetplot import UpSet, from_memberships
import seaborn as sns


ARTEFACT_RE = re.compile(r'(?:^|\s)artefact=([^\s]+)')

# Timing logs parsing
# ---------- paths ----------
def marti_bam_path(base_name, marti_root):
    p = marti_root / base_name / f"{base_name}.classified.bam"
    return p

def restrander_paths(base_name, restrander_root):
    p_ok  = restrander_root / f"{base_name}.restrander_PB.fastq.gz"

    if "." in base_name:
        head, tail = base_name.split(".", 1)
        unk_base = f"{head}-unknowns.{tail}"
    else:
        unk_base = f"{base_name}-unknowns"
    p_unk = restrander_root / f"{unk_base}.restrander_PB.fastq.gz"
    return p_ok, p_unk

# ---------- parsing ----------
def parse_marti_bam_labels(bam_path):
    """Return dict[read_id] = raw_label from BAM 'lb' tag; missing -> None."""
    import pysam
    out = {}
    if not bam_path.exists():
        return out
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False, threads=2) as bam:
        for aln in bam:
            rid = aln.query_name
            try:
                out[rid] = aln.get_tag("lb")
            except KeyError:
                out[rid] = None
    return out

def _iter_fastq_names(fq_path):
    """Yield header lines' identifiers as (read_id, header_line)."""
    if not fq_path.exists():
        return
    with gzip.open(fq_path, "rt", encoding="utf-8", errors="replace") as fh:
        it = iter(fh)
        for header in it:
            if not header.startswith("@"):
                # fast path: try to realign to next record
                continue
            read_id = header[1:].split()[0]
            yield read_id, header.rstrip("\n")
            _ = next(it, None)  # seq
            _ = next(it, None)  # +
            _ = next(it, None)  # qual

def parse_restrander_labels(ok_fq, unk_fq):
    """Return dict[read_id] = raw_label where ok -> 'Proper'; unknowns -> artefact or 'unknown'."""
    out = {}
    # proper file
    if ok_fq.exists():
        for rid, _hdr in _iter_fastq_names(ok_fq):
            out[rid] = "Proper"
    # unknowns file
    if unk_fq.exists():
        for rid, hdr in _iter_fastq_names(unk_fq):
            m = ARTEFACT_RE.search(hdr)
            out[rid] = m.group(1) if m else "unknown"
    return out

# ---------- categorization ----------
def cat_marti(raw):
    return raw
    # if raw in {"Proper", "TsoTso", "RtRt", "Unk"}:
    #     return raw
    # return "Other"

def cat_restrander(raw):
    if raw == "Proper":
        return "Proper"
    if raw == "TSO-TSO":
        return "TsoTso"
    if raw == "RTP-RTP":
        return "RtRt"
    return "Unk"

# ---------- per-input summary ----------
def summarize_input(base_name, marti_root, restrander_root):
    mbam = marti_bam_path(base_name, marti_root)
    ok_fq, unk_fq = restrander_paths(base_name, restrander_root)

    marti_labels = parse_marti_bam_labels(mbam)
    restr_labels = parse_restrander_labels(ok_fq, unk_fq)

    # unified keyed store: read_id -> (marti_raw, restr_raw)
    keys = set(marti_labels) | set(restr_labels)
    rows = []
    for rid in keys:
        m_raw = marti_labels.get(rid)
        r_raw = restr_labels.get(rid)
        m_cat = cat_marti(m_raw)
        r_cat = cat_restrander(r_raw)
        rows.append((rid, m_raw, r_raw, m_cat, r_cat))

    df = pd.DataFrame(rows, columns=["read_id", "marti_raw", "restr_raw", "marti_cat", "restr_cat"])

    ct = pd.crosstab(df["marti_cat"], df["restr_cat"]).astype(int)
    ct.index.name = "Marti"
    ct.columns.name = "Restrander"

    return df, ct

# ---------- multi-input summary ----------
def summarize_many(base_names, marti_root, restrander_root):
    per_input = []
    agg = Counter()
    for b in base_names:
        _df, ct = summarize_input(b, marti_root, restrander_root)
        # store per-input table flattened
        t = ct.stack().reset_index()
        t.columns = ["marti_cat", "restr_cat", "count"]
        t.insert(0, "input", b)
        per_input.append(t)
        # aggregate
        for (m, r), c in ct.stack().items():
            agg[(m, r)] += int(c)
    per_input_df = pd.concat(per_input, ignore_index=True) if per_input else pd.DataFrame(columns=["input","marti_cat","restr_cat","count"])
    agg_df = (
        pd.DataFrame([{"marti_cat": m, "restr_cat": r, "count": c} for (m, r), c in agg.items()])
        .pivot(index="marti_cat", columns="restr_cat", values="count")
        .fillna(0)
        .astype(int)
    )
    agg_df.index.name = "Marti"
    agg_df.columns.name = "Restrander"
    return per_input_df, agg_df


# --- collapse helper ---
def collapse_marti_ct(ct, keep=("Proper","TsoTso","RtRt","Unk"), other="Other"):
    """
    Return a copy of ct where Marti categories not in `keep` are summed into `other`.
    ct: crosstab with Marti as index and Restrander as columns.
    """
    idx = pd.Series(ct.index, index=ct.index)
    idx = idx.where(idx.isin(keep), other)
    out = ct.groupby(idx).sum()
    # stable row order: keep first, then 'other' if present
    order = [*keep] + ([other] if other in out.index and other not in keep else [])
    out = out.reindex(order, fill_value=0)
    out.index.name = ct.index.name
    out.columns.name = ct.columns.name
    return out


# Plotting
# ----- absolute times per input, grouped by tool -----
def plot_absolute(df_long, metric="cpu_time", title=None, savepath=None):
    d = df_long[df_long["metric"] == metric].copy()
    inputs = sorted(d["input"].unique())
    tools  = sorted(d["tool"].unique())
    ix = {k:i for i,k in enumerate(inputs)}
    tx = {k:i for i,k in enumerate(tools)}
    W = 0.8
    bar_w = W / max(1, len(tools))

    x = np.arange(len(inputs))
    fig, ax = plt.subplots(figsize=(max(8, len(inputs)*0.5), 8))
    for j, tool in enumerate(tools):
        y = np.full(len(inputs), np.nan)
        sub = d[d["tool"] == tool]
        for _, r in sub.iterrows():
            y[ix[r["input"]]] = r["seconds"]
        ax.bar(x + (j - (len(tools)-1)/2)*bar_w, y, width=bar_w, label=tool)

    ax.set_xticks(x)
    ax.set_xticklabels(inputs, rotation=60, ha="right")
    ax.set_ylabel("Seconds")
    ax.set_title(title or f"{metric.replace('_',' ').title()} by input")
    ax.legend(ncol=min(3, len(tools)))
    ax.axhline(0, lw=0.5, color="k")
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=300)
    plt.show()

# ----- heatmap with a direction for stats -----
def plot_heatmap(ct, normalize=None, title=None, savepath=None):
    # counts in display orientation: rows=Restrander, cols=Marti
    B = ct.values.astype(float).T
    ann = ct.values.astype(int).T  # counts for annotations

    if normalize == 'all':
        denom = B.sum()
        D = B / denom if denom else B
        cbar_label = "fraction of all"
    elif normalize == 'row':
        denom = B.sum(axis=1, keepdims=True)          # per displayed row
        D = np.divide(B, denom, out=np.zeros_like(B), where=denom != 0)
        cbar_label = "row fraction"
    elif normalize == 'col':
        denom = B.sum(axis=0, keepdims=True)          # per displayed column
        D = np.divide(B, denom, out=np.zeros_like(B), where=denom != 0)
        cbar_label = "column fraction"
    else:
        D = B
        cbar_label = "count"

    # figure sized for cols=Marti, rows=Restrander
    fig, ax = plt.subplots(figsize=(1.2 * ct.shape[0] + 2, 1.0 * ct.shape[1] + 2))
    im = ax.imshow(D, aspect="auto")

    # axes: x=Marti (columns of D), y=Restrander (rows of D)
    ax.set_xticks(range(ct.shape[0])); ax.set_xticklabels(ct.index, rotation=45, ha='right')
    ax.set_yticks(range(ct.shape[1])); ax.set_yticklabels(ct.columns)

    # annotations: counts and optional %
    for i in range(ct.shape[1]):          # rows in D / Restrander
        for j in range(ct.shape[0]):      # cols in D / Marti
            txt = str(ann[i, j])
            if normalize:
                txt += f"\n{D[i, j]*100:.1f}%"
            ax.text(j, i, txt, ha='center', va='center', fontsize=9)

    ax.set_xlabel("Marti categories", fontsize=12, fontweight="bold", labelpad=10)
    ax.set_ylabel("Restrander categories", fontsize=12, fontweight="bold", labelpad=10)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(cbar_label)

    ax.set_title(title or "Marti vs Restrander")
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath)  # .svg/.png/.pdf by extension
    plt.show()


# ----- Upset plot -----
def plot_upset_from_ct(
    ct,
    title="UpSet: Marti vs Restrander",
    savepath=None,
    top=None,
    exclude_pair=None,
    annotate=True
):
    s = ct.stack().astype(int)
    s = s[s > 0].sort_values(ascending=False)
    if exclude_pair and exclude_pair in s.index:
        s = s.drop(exclude_pair)
    if top:
        s = s.head(top)

    memberships = [(f"Marti:{m}", f"Restrander:{r}") for (m, r) in s.index]
    data = from_memberships(memberships, data=s.values)

    up = UpSet(
        data,
        subset_size="sum",
        show_percentages=not annotate,   # toggle built-in labels to replace with custom
        sort_by="cardinality"
    )
    fig = plt.figure(figsize=(10, 6))
    up.plot(fig=fig)
    plt.suptitle(title)

    bar_ax = max(fig.axes, key=lambda ax: len(ax.patches))
    if annotate:
        total = float(s.values.sum()) if s.size else 1.0
        for patch in bar_ax.patches:
            h = patch.get_height()
            if h <= 0:
                continue
            x = patch.get_x() + patch.get_width()/2
            bar_ax.text(x, h, f"{int(round(h))} ({h/total*100:.1f}%)",
                        ha="left", va="bottom", fontsize=8, rotation=30)
        fig.tight_layout()
        
    # avoid scientific notation like "1e7"
    # bar_ax.ticklabel_format(style="plain", axis="y")
    # or just move it to the left
    off = bar_ax.yaxis.get_offset_text()
    off.set_x(-0.05)      # nudge left; adjust until it clears the bar area
    off.set_ha('right')   # align to the right so it stays tight to axis
    
    if savepath:
        plt.savefig(savepath)
    plt.show()


def plot_upset_from_per_input_df(per_input_df, savepath=None):

    KEEP = ("Proper","TsoTso","RtRt","Unk")
    OTHER = "Other"

    exclude_pair = ("Proper","Proper")
    top = None

    df = per_input_df.copy()
    df = df[df["count"] > 0]  # drop zeros so they don't pollute stacks
    df["marti_collapsed"] = df["marti_cat"].where(df["marti_cat"].isin(KEEP), OTHER)
    df["by"] = df.apply(
        lambda r: r["marti_cat"] if r["marti_collapsed"] == OTHER else r["marti_collapsed"],
        axis=1
    )

    # choose intersections (by collapsed Marti × Restrander)
    totals = (df.groupby(["marti_collapsed","restr_cat"])["count"].sum()
                .sort_values(ascending=False))
    if exclude_pair and exclude_pair in totals.index:
        totals = totals.drop(exclude_pair)
    if top:
        totals = totals.head(top)
    shown = set(totals.index)

    # filter rows to shown intersections; build memberships + data
    sub = df[df.set_index(["marti_collapsed","restr_cat"]).index.isin(shown)].copy()
    sub["memberships"] = list(
        zip("Marti:" + sub["marti_collapsed"], "Restrander:" + sub["restr_cat"])
    )
    upset_input = from_memberships(sub["memberships"].tolist(),
                                   data=sub[["by","count"]])
    by_labels = list(dict.fromkeys(sub["by"]))  # preserve encounter order
    colors = sns.color_palette("tab20", len(by_labels))

    upset = UpSet(
        upset_input, 
        sort_by="cardinality", 
        sum_over="count", 
        subset_size="sum", 
        intersection_plot_elements=0  # disable the default bar chart
    )

    upset.add_stacked_bars(
        by="by", 
        sum_over="count",
        colors=colors,
        title="Marti categories",
        elements=10
    )

    fig = plt.figure(figsize=(12, 7))
    upset.plot(fig=fig)

    # Move the legend to the right of the figure
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg:
            leg.set_bbox_to_anchor((1.05, 1))   # move outside to the right
            leg.set_loc("upper left")           # anchor from upper left

    # annotate totals on stacked bars
    total_sum = float(sub["count"].sum()) or 1.0
    # pick the axis with the most rectangle patches (the stacked-bars axis)
    stack_ax = max(fig.axes, key=lambda a: sum(isinstance(p, Rectangle) for p in a.patches))
    # sum stacks per bar by x-center
    per_bar = {}
    top_y = {}
    for p in stack_ax.patches:
        xc = round(p.get_x() + p.get_width()/2, 6)
        per_bar[xc] = per_bar.get(xc, 0.0) + p.get_height()
        top_y[xc] = max(top_y.get(xc, 0.0), p.get_y() + p.get_height())
    for xc, h in per_bar.items():
        stack_ax.text(
            xc, top_y[xc], f"{int(round(h))} ({h/total_sum*100:.1f}%)",
            ha="left", va="bottom", fontsize=8, rotation=30, zorder=10, clip_on=False
        )
    
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # leave space for the legend
    if savepath:
        plt.savefig(savepath)
    plt.show()