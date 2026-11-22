import os
import shutil
import json
import warnings
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
import jinja2
import qiime2

from ._methods import QC_PREFIX, index_columns

ASSET_DIR = os.path.join(os.path.dirname(__file__), "assets", "web")

COLOR_MAP = {
    "Red": "#E41A1C",
    "Orange": "#FF7F00",
    "Yellow": "#FFFF33",
    "Green": "#4DAF4A",
    "Purple": "#88419D",
    "Blue": "#377EB8",
    "Other": "#999999"
}

FDR_ALPHA = 0.05


def _json_for_script(obj):
    """Serializes obj to JSON that is safe to embed inside a <script> tag."""
    return json.dumps(obj).replace("</", "<\\/")


def _finite_or_none(value):
    value = float(value)
    return value if np.isfinite(value) else None


def benjamini_hochberg(pvals):
    """Benjamini-Hochberg FDR-adjusted q-values (same order as input)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order] * n / np.arange(1, n + 1)
    # Enforce monotonicity from the largest p-value downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qvals = np.empty(n)
    qvals[order] = np.minimum(ranked, 1.0)
    return qvals


def analysis_units(indices, metadata, group_column, subject_column=None):
    """
    Index values joined to the group labels, one row per analysis unit:
    samples, or subjects within a group (mean of their samples) when
    `subject_column` is given.
    """
    columns = metadata.columns
    if group_column not in columns:
        raise ValueError("Metadata column %r not found." % group_column)
    if columns[group_column].type != "categorical":
        raise ValueError("group_column %r must be categorical." % group_column)
    if subject_column is not None and subject_column not in columns:
        raise ValueError("Metadata column %r not found." % subject_column)

    meta_df = metadata.to_dataframe()
    keep = [group_column] + ([subject_column] if subject_column else [])
    design = meta_df[keep].dropna()
    design.columns = (["__group__"]
                      + (["__subject__"] if subject_column else []))

    metrics = index_columns(indices)
    merged = indices[metrics].join(design, how="inner")
    if subject_column:
        merged["__subject__"] = merged["__subject__"].astype(str)
        merged = (merged.groupby(["__group__", "__subject__"])[metrics]
                  .mean().reset_index())
    return merged, metrics


def group_comparisons(units, metrics):
    """
    Mann-Whitney U (two groups) or Kruskal-Wallis (more groups) test of every
    index, with effect sizes and Benjamini-Hochberg q-values across indices.

    Effect sizes: for two groups, AUC = P(X_B > X_A) + 0.5 P(X_B = X_A),
    with A and B the groups in sorted order (0.5 = no separation); for more
    groups, epsilon-squared = H / (n - 1).
    """
    results = []
    group_names = sorted(units["__group__"].astype(str).unique())
    for metric in metrics:
        data = {g: units.loc[units["__group__"].astype(str) == g, metric]
                .dropna().to_numpy(dtype=float) for g in group_names}
        data = {g: v for g, v in data.items() if len(v) >= 2}
        if len(data) < 2:
            continue
        names = list(data)
        values = [data[g] for g in names]
        try:
            if len(values) == 2:
                stat, pval = mannwhitneyu(values[1], values[0],
                                          alternative="two-sided")
                method = "Mann-Whitney U"
                effect_name = "AUC (%s > %s)" % (names[1], names[0])
                effect = stat / (len(values[0]) * len(values[1]))
            else:
                stat, pval = kruskal(*values)
                method = "Kruskal-Wallis"
                n_total = sum(len(v) for v in values)
                effect_name = "epsilon-squared"
                effect = stat / (n_total - 1)
        except ValueError as e:
            warnings.warn("Skipping test of '%s': %s" % (metric, e),
                          UserWarning)
            continue
        if np.isnan(pval):
            warnings.warn("Skipping test of '%s': p-value is NaN (e.g. all "
                          "values identical)." % metric, UserWarning)
            continue
        results.append({
            "metric": metric,
            "method": method,
            "groups": ", ".join("%s (n=%d, median=%.3g)"
                                % (g, len(data[g]), np.median(data[g]))
                                for g in names),
            "statistic": float(stat),
            "effect_name": effect_name,
            "effect": float(effect),
            "p_value": float(pval)
        })

    qvals = benjamini_hochberg([r["p_value"] for r in results])
    for result, q in zip(results, qvals):
        result["q_value"] = float(q)
        result["significant"] = bool(q < FDR_ALPHA)
    return results


def _qc_summary(indices):
    rows = []
    labels = {
        QC_PREFIX + "depth": "Reads per sample",
        QC_PREFIX + "assigned_fraction": "Fraction of reads in a complex",
        QC_PREFIX + "species_fraction": "Fraction of reads named to species",
    }
    for col, label in labels.items():
        if col not in indices.columns:
            continue
        values = indices[col].dropna()
        if values.empty:
            continue
        rows.append({
            "label": label,
            "min": float(values.min()),
            "median": float(values.median()),
            "max": float(values.max()),
        })
    return rows


def summarize(
    output_dir: str,
    complexes: pd.DataFrame,
    indices: pd.DataFrame,
    mapping: pd.DataFrame = None,
    metadata: qiime2.Metadata = None,
    group_column: str = None,
    subject_column: str = None
) -> None:
    if (group_column or subject_column) and metadata is None:
        raise ValueError("group_column and subject_column require metadata.")
    if subject_column and not group_column:
        raise ValueError("subject_column requires group_column.")

    data_dir = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    complexes.to_csv(os.path.join(data_dir, "complexes.tsv"), sep="\t",
                     index_label="sample-id")
    indices.to_csv(os.path.join(data_dir, "indices.tsv"), sep="\t",
                   index_label="sample-id")

    # Vendored JS assets (no CDN: the .qzv must work offline)
    for sub in ("css", "js"):
        shutil.copytree(os.path.join(ASSET_DIR, sub),
                        os.path.join(output_dir, sub), dirs_exist_ok=True)

    total_features = len(mapping) if mapping is not None else 0
    resolution_counts = {}
    if mapping is not None and "Resolution" in mapping.columns:
        resolution_counts = mapping["Resolution"].value_counts().to_dict()

    groups = {}
    stats_results = []
    unit_label = "samples"
    if group_column:
        units, metrics = analysis_units(indices, metadata, group_column,
                                        subject_column)
        stats_results = group_comparisons(units, metrics)
        if subject_column:
            unit_label = "subjects (mean of their samples per group)"
        group_series = metadata.to_dataframe()[group_column].dropna()
        groups = {str(k): str(v) for k, v in group_series.items()}

    complexes_long = [
        {"sample_id": str(sample_id), "complex": str(name),
         "abundance": float(val), "group": groups.get(str(sample_id), "")}
        for sample_id, row in complexes.iterrows()
        for name, val in row.items()
    ]
    indices_long = [
        {"sample_id": str(sample_id), "metric": str(metric),
         "value": _finite_or_none(row[metric]),
         "group": groups.get(str(sample_id), "all samples")}
        for sample_id, row in indices.iterrows()
        for metric in index_columns(indices)
    ]

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(ASSET_DIR),
        autoescape=jinja2.select_autoescape(["html"])
    )
    template = env.get_template("index.html")
    html_content = template.render(
        total_samples=len(complexes),
        total_features=total_features,
        resolution_counts=resolution_counts,
        qc_rows=_qc_summary(indices),
        group_column=group_column,
        subject_column=subject_column,
        unit_label=unit_label,
        stats_results=stats_results,
        fdr_alpha=FDR_ALPHA,
        indices_table_html=indices.head(20).to_html(
            classes="table table-striped", border=0, na_rep="NaN"),
        complexes_json=_json_for_script(complexes_long),
        indices_json=_json_for_script(indices_long),
        color_map_json=_json_for_script(COLOR_MAP)
    )
    with open(os.path.join(output_dir, "index.html"), "w",
              encoding="utf-8") as fh:
        fh.write(html_content)
