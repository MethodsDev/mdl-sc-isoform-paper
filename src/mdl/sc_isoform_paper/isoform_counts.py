from __future__ import annotations
from pathlib import Path
from collections import defaultdict
from typing import Mapping, Set, Tuple, Literal
import os
import gzip
import shutil

import pandas as pd
import pysam
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite




def build_cell_by_isoform_matrix(
    bam_file: str | os.PathLike,
    cell_to_cluster: Mapping[str, str | int],
    good_priming_categories: Set[str],
    *,
    isoform_tag: Literal["XI", "YT"] = "XI",
    threads: int = 0,
) -> tuple[dict, dict]:
    """Return two nested dicts {cell: {isoform: set(UMIs)}}: (all, good-priming)."""
    cell_isoform_umis = defaultdict(lambda: defaultdict(set))
    cell_goodprim_isoform_umis = defaultdict(lambda: defaultdict(set))

    with pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=threads) as fh:
        for a in fh:
            try:
                cell = a.get_tag("CB")
                umi = a.get_tag("UB")
                isoform = a.get_tag(isoform_tag)
                priming_category = a.get_tag("XC")
            except KeyError:
                continue
            if cell not in cell_to_cluster: 
                continue
            if isoform is None:
                continue

            cell_isoform_umis[cell][isoform].add(umi)
            if priming_category in good_priming_categories:
                cell_goodprim_isoform_umis[cell][isoform].add(umi)

    return cell_isoform_umis, cell_goodprim_isoform_umis


def cell_isoform_set_to_df(cell_isoform_umis: Mapping[str, Mapping[str, Set[str]]]) -> pd.DataFrame:
    rows = []
    for cell, isoforms in cell_isoform_umis.items():
        for isoform, umis in isoforms.items():
            rows.append((cell, isoform, len(umis)))
    df = pd.DataFrame(rows, columns=["cell", "isoform", "count"])
    return (
        df.pivot_table(index="cell", columns="isoform", values="count", fill_value=0)
          .astype("int32")
    )


def make_anndata_from_counts(counts_df: pd.DataFrame, cell_to_cluster: Mapping[str, str | int]) -> sc.AnnData:
    adata = sc.AnnData(counts_df)
    adata.obs["cluster"] = pd.Series(cell_to_cluster).reindex(counts_df.index)
    return adata


def library_prefix(library: tuple[str, ...]) -> str:
    """Reproduce notebook naming exactly."""
    left = library[0].replace(" ", "_") if library else ""
    right = (library[1] + "_") if len(library) == 2 else ""
    return "_".join([left, right])


def write_h5ad(adata: sc.AnnData, path: str | os.PathLike) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(path)


def export_10x_from_h5ad(
    input_h5ad: str | os.PathLike,
    output_dir: str | os.PathLike,
    *,
    bad_isoforms: Set[str] = frozenset({"NA", "."}),
) -> None:
    """
    Export one .h5ad to a single 10X-style folder:
      barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz, and cell_clusters.csv
    """
    input_h5ad = Path(input_h5ad)
    output_dir = Path(output_dir)
    adata = sc.read_h5ad(input_h5ad)
    adata.var_names_make_unique()

    good_columns = ~adata.var_names.isin(bad_isoforms)
    adata = adata[:, good_columns]

    matrix = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)

    output_dir.mkdir(parents=True, exist_ok=True)

    # barcodes
    pd.Series(adata.obs_names).to_csv(output_dir / "barcodes.tsv", index=False, header=False)

    # features
    features = pd.DataFrame(
        {"gene_id": adata.var_names, "gene_name": adata.var_names, "feature_type": "Gene Expression"}
    )
    features.to_csv(output_dir / "features.tsv", sep="\t", index=False, header=False)

    # matrix.mtx (genes x cells)
    mmwrite(str(output_dir / "matrix.mtx"), matrix.T)

    # gzip three files
    for fname in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
        src = output_dir / fname
        dst = output_dir / f"{fname}.gz"
        with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    # clusters
    adata.obs[["cluster"]].to_csv(output_dir / "cell_clusters.csv", index=True, header=True)