"""Simulated known-truth validation of collapse_complexes.

    python -I simulation.py <results_dir>

Generates synthetic subgingival communities with known complex fractions and
SMDI values, runs the plugin, and compares recovered values with truth. The
MVP uses two scenarios (clean and 20% of assigned taxa unclassified) at
three sequencing depths.

Writes:

    <results_dir>/tables/simulation.tsv
    <results_dir>/figures/fig6_simulation.png
    <results_dir>/figures/fig6_simulation.svg
"""
import argparse
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
import biom
from scipy import stats as sps

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from q2_periodontitis._methods import (  # noqa: E402
    COMPLEXES,
    SMDI_SPECIES,
    clr_microbiome,
    collapse_complexes,
    load_complex_reference,
    load_signature,
    load_synonyms,
    taxon_key,
    _canonical,
)

DEPTHS = [1000, 5000, 20000]
N_HEALTH = 50
N_PERIODONTITIS = 50
SEED = 20261104
PSEUDOCOUNT = 1e-6
UNCLASSIFIED = "d__Bacteria; g__uncultured"


def build_taxa():
    """Reference and SMDI taxa, with their true complex and taxonomy string."""
    synonyms = load_synonyms()
    reference = load_complex_reference()
    signature = load_signature(SMDI_SPECIES)
    direction = {
        taxon_key(_canonical(taxon, synonyms)): value
        for taxon, value in zip(signature["taxon"], signature["direction"])
    }

    rows = []
    seen = set()

    def add(name, complex_name, level):
        key = taxon_key(_canonical(name, synonyms))
        if key in seen:
            return
        seen.add(key)
        if level == "genus":
            taxonomy = "d__Bacteria; g__%s" % name.replace(" ", "_")
        else:
            genus = name.split()[0]
            taxonomy = "d__Bacteria; g__%s; s__%s" % (
                genus, name.replace(" ", "_"))
        rows.append({
            "taxon": name,
            "key": key,
            "complex": complex_name,
            "level": level,
            "taxonomy": taxonomy,
        })

    for row in reference.itertuples():
        add(row.taxon, row.complex, row.level)
    for name in signature["taxon"]:
        add(name, "Other", "species")
    return pd.DataFrame(rows), direction


def concentration(taxa, direction):
    """Group-specific Dirichlet concentrations for health and periodontitis."""
    health = np.ones(len(taxa))
    periodontitis = np.ones(len(taxa))
    for index, row in taxa.iterrows():
        value = direction.get(row["key"])
        if value == "dysbiotic":
            health[index], periodontitis[index] = 0.35, 2.5
        elif value == "normobiotic":
            health[index], periodontitis[index] = 2.0, 0.6
        elif row["complex"] in ("Red", "Orange"):
            health[index], periodontitis[index] = 0.5, 2.0
        elif row["complex"] in ("Yellow", "Green", "Purple", "Blue"):
            health[index], periodontitis[index] = 1.8, 0.6
    return health, periodontitis


def generate_communities(taxa, direction, n_health=N_HEALTH,
                         n_periodontitis=N_PERIODONTITIS, seed=SEED):
    """Latent relative abundances and group labels."""
    rng = np.random.default_rng(seed)
    health_alpha, periodontitis_alpha = concentration(taxa, direction)
    health = rng.dirichlet(health_alpha * 5, size=n_health)
    periodontitis = rng.dirichlet(periodontitis_alpha * 5,
                                  size=n_periodontitis)
    relative = np.vstack([health, periodontitis])
    groups = np.array(["health"] * n_health
                      + ["periodontitis"] * n_periodontitis)
    return relative, groups


def simulate_counts(relative, depth, seed):
    """Multinomial sequencing counts for each sample."""
    rng = np.random.default_rng(seed)
    return np.vstack([rng.multinomial(depth, row) for row in relative])


def truth_from_relative(relative, taxa, direction):
    """True complex fractions and SMDI from the latent communities."""
    sample_ids = ["S%d" % i for i in range(len(relative))]
    complexes = pd.DataFrame(index=sample_ids)
    true_complex = taxa["complex"].to_numpy()
    for name in COMPLEXES:
        complexes[name] = relative[:, true_complex == name].sum(axis=1)

    frame = pd.DataFrame(relative, index=sample_ids,
                         columns=taxa["key"].to_numpy())
    clr = clr_microbiome(frame, pseudocount=PSEUDOCOUNT)
    dysbiotic = [key for key, value in direction.items()
                 if value == "dysbiotic"]
    normobiotic = [key for key, value in direction.items()
                   if value == "normobiotic"]
    smdi = clr[dysbiotic].mean(axis=1) - clr[normobiotic].mean(axis=1)
    return complexes, smdi


def choose_unassigned(taxa, fraction=0.2, seed=SEED + 1):
    """Indices of assigned taxa to mask in the stress scenario."""
    candidates = np.flatnonzero(taxa["complex"].to_numpy() != "Other")
    count = max(1, int(round(fraction * len(candidates))))
    rng = np.random.default_rng(seed)
    return rng.choice(candidates, size=count, replace=False)


def run_scenario(counts, taxonomy, table, feature_ids, taxa, true_complex,
                 true_smdi):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        complexes, indices, mapping = collapse_complexes(
            table, taxonomy, smdi_pseudocount=PSEUDOCOUNT)

    recovered = complexes[COMPLEXES].to_numpy()
    expected = true_complex[COMPLEXES].to_numpy()
    complex_mae = float(np.mean(np.abs(recovered - expected)))

    recovered_smdi = indices["SMDI"].to_numpy(dtype=float)
    expected_smdi = true_smdi.to_numpy(dtype=float)
    smdi_mae = float(np.mean(np.abs(recovered_smdi - expected_smdi)))
    smdi_pearson = float(np.corrcoef(expected_smdi, recovered_smdi)[0, 1])
    smdi_spearman = float(sps.spearmanr(expected_smdi, recovered_smdi)[0])

    recovered_complex = mapping["Complex"].reindex(feature_ids).to_numpy()
    correct = recovered_complex == taxa["complex"].to_numpy()
    assignment_accuracy = float(counts[:, correct].sum() / counts.sum())
    return {
        "complex_mae": complex_mae,
        "smdi_mae": smdi_mae,
        "smdi_pearson": smdi_pearson,
        "smdi_spearman": smdi_spearman,
        "assignment_accuracy": assignment_accuracy,
    }, indices["SMDI"]


def run_simulation():
    taxa, direction = build_taxa()
    relative, groups = generate_communities(taxa, direction)
    true_complex, true_smdi = truth_from_relative(relative, taxa, direction)
    feature_ids = ["F%d" % i for i in range(len(taxa))]
    sample_ids = ["S%d" % i for i in range(len(relative))]
    unassigned = choose_unassigned(taxa)

    rows = []
    detail = None
    for depth in DEPTHS:
        counts = simulate_counts(relative, depth, SEED + depth)
        table = biom.Table(counts.T, feature_ids, sample_ids)
        for scenario in ("clean", "unassigned"):
            taxonomy = taxa["taxonomy"].tolist()
            if scenario == "unassigned":
                for index in unassigned:
                    taxonomy[index] = UNCLASSIFIED
            frame = pd.DataFrame(
                {"Taxon": taxonomy},
                index=pd.Index(feature_ids, name="Feature ID"))
            metrics, recovered_smdi = run_scenario(
                counts, frame, table, feature_ids, taxa,
                true_complex, true_smdi)
            unassigned_fraction = (
                float(counts[:, unassigned].sum() / counts.sum())
                if scenario == "unassigned" else 0.0)
            rows.append({
                "scenario": scenario,
                "depth": depth,
                "seed": SEED,
                "n_samples": len(relative),
                "n_taxa": len(taxa),
                "unassigned_read_fraction": unassigned_fraction,
                **metrics,
            })
            if scenario == "clean" and depth == 5000:
                detail = {
                    "true_smdi": true_smdi,
                    "recovered_smdi": recovered_smdi,
                    "groups": groups,
                }
    return pd.DataFrame(rows), detail


def plot_simulation(frame, detail, figdir):
    """Writes fig6_simulation.png and .svg."""
    if frame.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame = frame.copy()
    for column in ("depth", "complex_mae", "smdi_spearman"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["depth", "complex_mae", "smdi_spearman"])
    if frame.empty:
        return

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8,
        "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
        "xtick.color": "#52514e", "ytick.color": "#52514e",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.axisbelow": True,
        "grid.color": "#e4e3df", "grid.linewidth": 0.6,
        "legend.frameon": False, "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))

    ax = axes[0]
    if detail is not None:
        true = np.asarray(detail["true_smdi"], dtype=float)
        recovered = np.asarray(detail["recovered_smdi"], dtype=float)
        groups = np.asarray(detail["groups"])
        for group, color in (("health", "#2a78d6"),
                             ("periodontitis", "#eb6834")):
            mask = groups == group
            ax.scatter(true[mask], recovered[mask], s=12, color=color,
                       alpha=0.8, label=group)
        low = min(true.min(), recovered.min())
        high = max(true.max(), recovered.max())
        ax.plot([low, high], [low, high], color="#52514e", lw=0.8, ls="--")
        ax.legend(fontsize=6, loc="upper left")
    ax.set_xlabel("True SMDI (latent communities)")
    ax.set_ylabel("Recovered SMDI")
    ax.set_title("a  SMDI recovery (clean, 5k reads)", fontsize=8,
                 loc="left", color="#0b0b0b")

    ax = axes[1]
    style = {
        "clean": ("#2a78d6", "o"),
        "unassigned": ("#eb6834", "s"),
    }
    for scenario, (color, marker) in style.items():
        group = frame[frame["scenario"] == scenario].sort_values("depth")
        ax.plot(group["depth"], group["complex_mae"], color=color,
                marker=marker, ms=4, lw=1.2, label=scenario)
    ax.set_xscale("log")
    ax.set_xticks(DEPTHS)
    ax.set_xticklabels(["1k", "5k", "20k"])
    ax.set_xlabel("Sequencing depth (reads)")
    ax.set_ylabel("Complex-fraction MAE")
    ax.set_title("b  Complex recovery", fontsize=8, loc="left",
                 color="#0b0b0b")
    ax.legend(fontsize=6)

    os.makedirs(figdir, exist_ok=True)
    fig.tight_layout()
    for extension in ("png", "svg"):
        fig.savefig(os.path.join(figdir, "fig6_simulation.%s" % extension),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", help="results directory")
    args = parser.parse_args()
    results = os.path.abspath(args.results)
    tables = os.path.join(results, "tables")
    figures = os.path.join(results, "figures")
    os.makedirs(tables, exist_ok=True)

    frame, detail = run_simulation()
    frame.to_csv(os.path.join(tables, "simulation.tsv"), sep="\t",
                 index=False, float_format="%.6g")
    plot_simulation(frame, detail, figures)
    print("simulation written to", os.path.join(tables, "simulation.tsv"))


if __name__ == "__main__":
    main()
