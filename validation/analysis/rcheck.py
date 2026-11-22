"""Cross-check and benchmark of the SMDI against R microbiome 1.6.0.

    python -I rcheck.py <data_dir> <results_dir> --rscript <Rscript> [--quick]

Cross-check: the SMDI and SMDI_simplified of collapse_complexes (default
parameters) are compared, sample by sample, with smdi_microbiome.R, which
applies microbiome::transform(x, "clr") of microbiome 1.6.0, the
implementation used to define the index. Both get the same counts and the
same feature -> taxon grouping (feature_units), so the check covers the
aggregation, pseudocount, CLR and index; taxonomy parsing is the plugin's.
Tables checked: the simulated communities (simulation.py), the synthetic
benchmark tables (performance.py) and every processed cohort and classifier
in <data_dir>.

Benchmark: on the synthetic tables of performance.py, the 49-species SMDI
computed by the plugin (compute_smdi, including taxonomy parsing) and by
smdi_microbiome.R (taxon grouping given). Each run is a separate process;
runtime is measured inside the process from loaded tables to index, and peak
RSS is that of the whole process, read by the process itself from VmHWM in
/proc/self/status (Linux). os.wait4 is not used for memory: its ru_maxrss
includes the parent's memory that the child carried until exec.

Writes:

    <results_dir>/tables/r_crosscheck.tsv
    <results_dir>/tables/r_benchmark.tsv
    <results_dir>/r_environment.json
    <results_dir>/figures/fig7_r_benchmark.png and .svg
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time
import warnings

import numpy as np
import pandas as pd
import biom

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import performance  # noqa: E402
import simulation  # noqa: E402
from common import CLASSIFIERS, load_cohorts  # noqa: E402
from q2_periodontitis._methods import (  # noqa: E402
    SMDI_SIMPLIFIED,
    SMDI_SPECIES,
    _canonical,
    collapse_complexes,
    compute_smdi,
    feature_units,
    load_signature,
    load_synonyms,
    taxon_key,
)

R_SCRIPT = os.path.join(HERE, "smdi_microbiome.R")
INDICES = {"SMDI": SMDI_SPECIES, "SMDI_simplified": SMDI_SIMPLIFIED}
TOLERANCE = 1e-9
REPEATS = performance.REPEATS


def signature_units(path, synonyms):
    """Signature taxa as the unit keys of feature_units, with direction."""
    sig = load_signature(path)
    level = sig["level"].iloc[0]
    if level == "species":
        keys = [taxon_key(_canonical(t, synonyms)) for t in sig["taxon"]]
    else:
        keys = [_canonical(t, synonyms).lower() for t in sig["taxon"]]
    return level, pd.DataFrame({"unit": keys, "direction": sig["direction"]})


def write_inputs(table, taxonomy, index, synonyms, workdir):
    """Writes the three input files of smdi_microbiome.R."""
    level, sig = signature_units(INDICES[index], synonyms)
    units = feature_units(table, taxonomy, synonyms, level)
    paths = [os.path.join(workdir, n)
             for n in ("counts.tsv", "units.tsv", "signature.tsv")]
    counts = pd.DataFrame(table.matrix_data.toarray(),
                          index=pd.Index(table.ids(axis="observation"),
                                         name="feature"),
                          columns=table.ids(axis="sample"))
    counts.to_csv(paths[0], sep="\t")
    pd.DataFrame({"feature": counts.index, "unit": units}).to_csv(
        paths[1], sep="\t", index=False)
    sig.to_csv(paths[2], sep="\t", index=False)
    return paths


def run_r(rscript, table, taxonomy, index, synonyms):
    """SMDI from R (Series), R's compute seconds and the peak RSS in MB."""
    with tempfile.TemporaryDirectory() as workdir:
        inputs = write_inputs(table, taxonomy, index, synonyms, workdir)
        out = os.path.join(workdir, "smdi.tsv")
        seconds, peak = timing(
            run_child([rscript, R_SCRIPT] + inputs + [out]))
        smdi = pd.read_csv(out, sep="\t", index_col=0)["SMDI"]
    smdi.index = smdi.index.astype(str)
    return smdi, seconds, peak


def timing(stdout):
    """Seconds and peak RSS (MB) from a child's last line '<s>\\t<MB>'."""
    seconds, peak = stdout.strip().splitlines()[-1].split("\t")
    return float(seconds), float(peak)


def own_peak_mb():
    """Peak RSS of this process since exec (VmHWM), in MB."""
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) / 1024
    raise RuntimeError("VmHWM not found in /proc/self/status")


def run_child(command):
    """Runs a command and returns its stdout."""
    # The QIIME 2 environment sets R_LIBS_USER to its own R 4.3 library,
    # which would make the R 3.6 of microbiome 1.6.0 load the wrong packages
    env = {k: v for k, v in os.environ.items()
           if not k.startswith("R_") and k != "RSTUDIO_WHICH_R"}
    # Output goes to files, so a large output can never block the child
    with tempfile.TemporaryFile("w+") as out, \
            tempfile.TemporaryFile("w+") as err:
        proc = subprocess.run(command, stdout=out, stderr=err, text=True,
                              env=env)
        if proc.returncode:
            err.seek(0)
            raise RuntimeError("%s failed:\n%s" % (command[0], err.read()))
        out.seek(0)
        return out.read()


def compare(plugin, r):
    """Max absolute difference; NaN must match NaN."""
    r = r.reindex(plugin.index)
    nan_mismatch = int((plugin.isna() != r.isna()).sum())
    both = plugin.notna() & r.notna()
    diff = (plugin[both] - r[both]).abs()
    return {"samples": len(plugin), "nan_mismatch": nan_mismatch,
            "max_abs_diff": float(diff.max()) if both.any() else np.nan,
            "pass": bool(nan_mismatch == 0 and (not both.any()
                                                or diff.max() < TOLERANCE))}


def check_table(rscript, source, classifier, table, taxonomy, synonyms):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        _, indices, _ = collapse_complexes(table, taxonomy)
    rows = []
    for index in INDICES:
        r_smdi, _, _ = run_r(rscript, table, taxonomy, index, synonyms)
        rows.append({"source": source, "classifier": classifier,
                     "index": index, "features": table.shape[0],
                     **compare(indices[index], r_smdi)})
        print(source, classifier, index, rows[-1]["max_abs_diff"],
              flush=True)
    return rows


def simulated_tables(quick):
    taxa, direction = simulation.build_taxa()
    relative, _ = simulation.generate_communities(taxa, direction)
    feature_ids = ["F%d" % i for i in range(len(taxa))]
    sample_ids = ["S%d" % i for i in range(len(relative))]
    taxonomy = pd.DataFrame({"Taxon": taxa["taxonomy"].tolist()},
                            index=pd.Index(feature_ids, name="Feature ID"))
    for depth in simulation.DEPTHS[:1] if quick else simulation.DEPTHS:
        counts = simulation.simulate_counts(
            relative, depth, simulation.SEED + depth)
        yield ("simulated %d reads" % depth,
               biom.Table(counts.T, feature_ids, sample_ids), taxonomy)


def cohort_tables(data_dir):
    import qiime2
    for c in load_cohorts().itertuples():
        out = os.path.join(data_dir, c.cohort, "out")
        if not os.path.exists(os.path.join(out, "table.qza")):
            continue
        table = qiime2.Artifact.load(os.path.join(out, "table.qza")).view(
            biom.Table)
        for clf in CLASSIFIERS:
            path = os.path.join(out, "taxonomy-%s.qza" % clf)
            if os.path.exists(path):
                yield c.cohort, clf, table, qiime2.Artifact.load(path).view(
                    pd.DataFrame)


def crosscheck(rscript, data_dir, quick):
    synonyms = load_synonyms()
    rows = []
    for source, table, taxonomy in simulated_tables(quick):
        rows += check_table(rscript, source, "", table, taxonomy, synonyms)
    grid = [(100, 10)] if quick else [(1000, 100), (10000, 500)]
    for features, samples in grid:
        table, taxonomy = performance.make_table(features, samples)
        rows += check_table(rscript, "synthetic %d x %d"
                            % (features, samples), "", table, taxonomy,
                            synonyms)
    if data_dir and not quick:
        for cohort, clf, table, taxonomy in cohort_tables(data_dir):
            rows += check_table(rscript, cohort, clf, table, taxonomy,
                                synonyms)
    return pd.DataFrame(rows)


def plugin_seconds(features, samples):
    """Child-process mode: seconds of compute_smdi on one synthetic table."""
    table, taxonomy = performance.make_table(features, samples)
    synonyms = load_synonyms()
    signature = load_signature(SMDI_SPECIES)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        start = time.perf_counter()
        compute_smdi(table, taxonomy, signature, synonyms)
        return time.perf_counter() - start


def benchmark(rscript, quick):
    synonyms = load_synonyms()
    features = performance.QUICK_FEATURES if quick else performance.FEATURES
    samples = performance.QUICK_SAMPLES if quick else performance.SAMPLES
    repeats = 1 if quick else REPEATS
    rows = []
    for n_samples in samples:
        for n_features in features:
            print("R benchmark %d features x %d samples"
                  % (n_features, n_samples), flush=True)
            table, taxonomy = performance.make_table(n_features, n_samples)
            for _ in range(repeats):
                seconds, peak = timing(run_child([
                    sys.executable, "-I", os.path.abspath(__file__),
                    "--plugin-seconds", str(n_features), str(n_samples)]))
                rows.append(("q2-periodontitis", n_features, n_samples,
                             table.matrix_data.nnz, seconds, peak))
                _, seconds, peak = run_r(rscript, table, taxonomy, "SMDI",
                                         synonyms)
                rows.append(("R microbiome 1.6.0", n_features, n_samples,
                             table.matrix_data.nnz, seconds, peak))
    raw = pd.DataFrame(rows, columns=["implementation", "features",
                                      "samples", "nonzeros", "seconds",
                                      "peak_rss_mb"])
    return raw.groupby(["implementation", "features", "samples"],
                       as_index=False).agg(
        nonzeros=("nonzeros", "first"), repeats=("seconds", "size"),
        seconds=("seconds", "median"), seconds_min=("seconds", "min"),
        seconds_max=("seconds", "max"), peak_rss_mb=("peak_rss_mb", "max"))


def plot_benchmark(bench, figdir):
    """Writes fig7_r_benchmark.png and .svg."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"q2-periodontitis": "#2a78d6", "R microbiome 1.6.0": "#eb6834"}
    markers = ["o", "s", "^", "D"]
    sample_counts = sorted(bench["samples"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for impl, color in colors.items():
        for i, n_samples in enumerate(sample_counts):
            g = bench[(bench["implementation"] == impl)
                      & (bench["samples"] == n_samples)].sort_values(
                          "features")
            style = {"color": color, "marker": markers[i % len(markers)],
                     "ms": 4, "lw": 1.2,
                     "label": "%s, %d samples" % (impl, n_samples)}
            axes[0].plot(g["features"], g["seconds"], **style)
            axes[1].plot(g["features"], g["peak_rss_mb"], **style)
    for ax, ylabel, title in ((axes[0], "Runtime (s)", "a  Runtime"),
                              (axes[1], "Peak RSS (MB)", "b  Peak memory")):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Features")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=8, loc="left")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(True, color="#e4e3df", linewidth=0.6)
    axes[1].legend(fontsize=6, frameon=False)
    os.makedirs(figdir, exist_ok=True)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(figdir, "fig7_r_benchmark.%s" % ext),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "--plugin-seconds":
        seconds = plugin_seconds(int(sys.argv[2]), int(sys.argv[3]))
        print("%.6f\t%.3f" % (seconds, own_peak_mb()))
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", help="processed cohorts ('' for none)")
    parser.add_argument("results")
    parser.add_argument("--rscript", required=True,
                        help="Rscript of the microbiome 1.6.0 environment")
    parser.add_argument("--quick", action="store_true",
                        help="small smoke-test tables only")
    args = parser.parse_args()

    tables = os.path.join(args.results, "tables")
    os.makedirs(tables, exist_ok=True)
    stdout = run_child([args.rscript, R_SCRIPT, "--versions"])
    with open(os.path.join(args.results, "r_environment.json"), "w") as fh:
        fh.write(stdout)

    check = crosscheck(args.rscript, args.data_dir, args.quick)
    check.to_csv(os.path.join(tables, "r_crosscheck.tsv"), sep="\t",
                 index=False, float_format="%.3g")
    bench = benchmark(args.rscript, args.quick)
    bench.to_csv(os.path.join(tables, "r_benchmark.tsv"), sep="\t",
                 index=False, float_format="%.6g")
    plot_benchmark(bench, os.path.join(args.results, "figures"))
    failed = check[~check["pass"]]
    if len(failed):
        print(failed.to_string())
        sys.exit("R cross-check FAILED for %d table(s)" % len(failed))
    print("R cross-check passed for %d tables" % len(check))


if __name__ == "__main__":
    main()
