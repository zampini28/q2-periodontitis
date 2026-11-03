"""Runtime and memory benchmark for collapse_complexes.

    python -I performance.py <results_dir>
    python -I performance.py <results_dir> --quick
    python -I performance.py <results_dir> --figure

The benchmark uses synthetic sparse feature tables, not biological data.
Each configuration runs in an isolated child process so that peak RSS is
measured per table size. Results are written to:

    <results_dir>/tables/performance.tsv
    <results_dir>/performance_environment.json

With --figure, the figure is written to <results_dir>/figures/.
"""
import argparse
import json
import os
import platform
import resource
import subprocess
import sys
import time
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy import sparse
import biom

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from q2_periodontitis import __version__  # noqa: E402
from q2_periodontitis._methods import (  # noqa: E402
    collapse_complexes,
    load_complex_reference,
)

FEATURES = [100, 1000, 5000, 10000, 25000]
SAMPLES = [10, 100, 500]
QUICK_FEATURES = [50, 100]
QUICK_SAMPLES = [5, 10]
DENSITY = 0.05
SEED = 20261103
REPEATS = 3
_UNASSIGNED = [
    "d__Bacteria; g__uncultured",
    "d__Bacteria; g__metagenome",
]


def taxonomy_templates():
    """QIIME-style taxonomy strings drawn from the built-in reference."""
    reference = load_complex_reference()
    templates = []
    species = reference.loc[reference["level"] == "species", "taxon"]
    for name in species:
        parts = name.split()
        if len(parts) >= 2:
            genus = parts[0]
            templates.append("d__Bacteria; g__%s; s__%s"
                             % (genus, name.replace(" ", "_")))
    genera = reference.loc[reference["level"] == "genus", "taxon"]
    for name in genera:
        templates.append("d__Bacteria; g__%s" % name.replace(" ", "_"))
    return templates + _UNASSIGNED


def make_table(features, samples, density=DENSITY, seed=SEED):
    """Synthetic counts (features x samples) and their taxonomy."""
    if features < 1 or samples < 1:
        raise ValueError("features and samples must be positive")
    rng = np.random.default_rng(seed)
    templates = taxonomy_templates()
    taxa = rng.choice(templates, size=features, replace=True)
    feature_ids = ["F%d" % i for i in range(features)]
    sample_ids = ["S%d" % j for j in range(samples)]

    nnz = max(1, int(round(features * samples * density)))
    rows = rng.integers(0, features, nnz)
    cols = rng.integers(0, samples, nnz)
    values = rng.integers(1, 1000, nnz)
    covered = {int(value) for value in cols}
    for sample in range(samples):
        if sample not in covered:
            rows = np.append(rows, 0)
            cols = np.append(cols, sample)
            values = np.append(values, 1)

    matrix = sparse.coo_matrix(
        (values, (rows, cols)), shape=(features, samples)).tocsr()
    taxonomy = pd.DataFrame(
        {"Taxon": taxa},
        index=pd.Index(feature_ids, name="Feature ID"))
    return biom.Table(matrix, feature_ids, sample_ids), taxonomy


def measure(features, samples, density=DENSITY, seed=SEED):
    """Runs collapse_complexes once and returns its time and peak RSS."""
    table, taxonomy = make_table(features, samples, density, seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        start = time.perf_counter()
        collapse_complexes(table, taxonomy)
        seconds = time.perf_counter() - start

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mb = peak / (1024 * 1024) if sys.platform == "darwin" else peak / 1024
    return {
        "features": features,
        "samples": samples,
        "density": density,
        "nonzeros": int(table.matrix_data.nnz),
        "seconds": seconds,
        "peak_rss_mb": peak_mb,
        "seed": seed,
    }


def environment_metadata():
    """Software versions and hardware used for the benchmark."""
    info = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "biom": getattr(biom, "__version__", "unknown"),
        "q2_periodontitis": __version__,
    }
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    info["total_memory_gb"] = round(kb / 1024 / 1024, 2)
                    break
    except OSError:
        pass
    return info


def run_one(features, samples, density, seed):
    """Runs one configuration in a child process and parses its JSON."""
    command = [
        sys.executable, os.path.abspath(__file__), "--single",
        "--features", str(features),
        "--samples", str(samples),
        "--density", str(density),
        "--seed", str(seed),
    ]
    result = subprocess.run(
        command, check=True, capture_output=True, text=True)
    return json.loads(result.stdout.strip().splitlines()[-1])


def plot_performance(perf, figdir):
    """Writes fig5_performance.png and .svg."""
    if perf.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    perf = perf.copy()
    for column in ("features", "samples", "seconds", "peak_rss_mb"):
        perf[column] = pd.to_numeric(perf[column], errors="coerce")
    perf = perf.dropna(subset=["features", "samples", "seconds",
                               "peak_rss_mb"])
    if perf.empty:
        return

    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#88419d"]
    markers = ["o", "s", "^", "D"]
    sample_counts = sorted(perf["samples"].unique())
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
    for index, sample_count in enumerate(sample_counts):
        group = perf[perf["samples"] == sample_count].sort_values("features")
        style = {
            "color": colors[index % len(colors)],
            "marker": markers[index % len(markers)],
            "ms": 4, "lw": 1.2,
            "label": "%d samples" % int(sample_count),
        }
        axes[0].plot(group["features"], group["seconds"], **style)
        axes[1].plot(group["features"], group["peak_rss_mb"], **style)

    for ax, ylabel, title in (
            (axes[0], "Runtime (s)", "a  Runtime"),
            (axes[1], "Peak RSS (MB)", "b  Peak memory")):
        ax.set_xscale("log")
        ax.set_xlabel("Features")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=8, loc="left", color="#0b0b0b")
        ticks = sorted(perf["features"].unique())
        ax.set_xticks(ticks)
        ax.set_xticklabels([
            "%g" % value if value < 1000 else "%gk" % (value / 1000)
            for value in ticks])
        ax.legend(fontsize=6, title="Samples", title_fontsize=6)
    axes[0].set_yscale("log")

    os.makedirs(figdir, exist_ok=True)
    fig.tight_layout()
    for extension in ("png", "svg"):
        fig.savefig(os.path.join(figdir, "fig5_performance.%s" % extension),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="?", help="results directory")
    parser.add_argument("--quick", action="store_true",
                        help="run a small smoke-test grid")
    parser.add_argument("--figure", action="store_true",
                        help="write fig5_performance")
    parser.add_argument("--single", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("--features", type=int)
    parser.add_argument("--samples", type=int)
    parser.add_argument("--density", type=float, default=DENSITY)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--repeats", type=int, default=REPEATS,
                        help="runs per configuration (default %d)" % REPEATS)
    args = parser.parse_args()

    if args.single:
        if args.features is None or args.samples is None:
            parser.error("--single requires --features and --samples")
        print(json.dumps(measure(
            args.features, args.samples, args.density, args.seed)))
        return

    if not args.results:
        parser.error("the results directory is required")
    results = os.path.abspath(args.results)
    tables = os.path.join(results, "tables")
    os.makedirs(tables, exist_ok=True)
    features = QUICK_FEATURES if args.quick else FEATURES
    samples = QUICK_SAMPLES if args.quick else SAMPLES
    repeats = 1 if args.quick else args.repeats
    if repeats < 1:
        parser.error("--repeats must be positive")

    rows = []
    for sample_count in samples:
        for feature_count in features:
            print("benchmark %d features x %d samples (%d repeats)"
                  % (feature_count, sample_count, repeats), flush=True)
            for _ in range(repeats):
                rows.append(run_one(
                    feature_count, sample_count, args.density, args.seed))

    raw = pd.DataFrame(rows)
    frame = raw.groupby(
        ["features", "samples", "density", "seed"], as_index=False).agg(
            nonzeros=("nonzeros", "first"),
            repeats=("seconds", "size"),
            seconds=("seconds", "median"),
            seconds_min=("seconds", "min"),
            seconds_max=("seconds", "max"),
            peak_rss_mb=("peak_rss_mb", "max"))
    frame.to_csv(
        os.path.join(tables, "performance.tsv"), sep="\t", index=False,
        float_format="%.6g")
    metadata = environment_metadata()
    metadata["repeats"] = repeats
    with open(os.path.join(results, "performance_environment.json"),
              "w") as fh:
        json.dump(metadata, fh, indent=2)
        fh.write("\n")
    if args.figure:
        plot_performance(frame, os.path.join(results, "figures"))
    print("performance written to", os.path.join(tables, "performance.tsv"))


if __name__ == "__main__":
    main()
