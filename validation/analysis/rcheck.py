"""Cross-check and benchmark of the SMDI against R microbiome 1.6.0.

    python -I rcheck.py <data_dir> <results_dir> --rscript <Rscript> [--quick]

Cross-check: the SMDI, SMDI_simplified and SMDI_custom of
collapse_complexes are compared, sample by sample, with smdi_microbiome.R, which
applies microbiome::transform(x, "clr") of microbiome 1.6.0, the
implementation used to define the index. Both get the same counts, raw
taxonomy strings, synonyms and signature files; smdi_microbiome.R parses the
taxonomy itself, so the check covers taxonomy parsing, aggregation,
pseudocount, CLR and index. The feature -> taxon units of both sides are
also compared feature by feature (unit_mismatch).
SMDI_custom uses custom_signature.tsv, a test signature (not a biological
one) written with historical names (Bacteroides forsythus, Streptococcus
sanguis, ...) and another HMT spelling ("oral taxon 360"). Each table is
checked with the default parameters, with min_confidence = 0.7 and with a
fixed smdi_pseudocount of 1e-5 (SETTINGS). Taxonomies without a Confidence
column (simulated, synthetic and VSEARCH tables) get confidences drawn from
U(0.5, 1) with a fixed seed, so the exclusion rule is exercised; the
excluded column counts the features it excluded.
Tables checked: the simulated communities (simulation.py), the synthetic
benchmark tables (performance.py) and every processed cohort and classifier
in <data_dir>.

Benchmark: on the synthetic tables of performance.py, the 49-species SMDI
computed by the plugin (compute_smdi, including taxonomy parsing) and by
smdi_microbiome.R (taxonomy parsing not timed). Each run is a separate process;
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
    _confidence_series,
    SYNONYMS,
    _taxon_series,
    collapse_complexes,
    compute_smdi,
    feature_units,
    load_signature,
    load_synonyms,
)

R_SCRIPT = os.path.join(HERE, "smdi_microbiome.R")
CUSTOM_SIGNATURE = os.path.join(HERE, "custom_signature.tsv")
INDICES = {"SMDI": SMDI_SPECIES, "SMDI_simplified": SMDI_SIMPLIFIED,
           "SMDI_custom": CUSTOM_SIGNATURE}
# collapse_complexes parameters -> smdi_microbiome.R options
SETTINGS = [("default", {}),
            ("min_confidence=0.7", {"min_confidence": 0.7}),
            ("smdi_pseudocount=1e-05", {"smdi_pseudocount": 1e-5})]
R_OPTIONS = {"min_confidence": "--min-confidence=%r",
             "smdi_pseudocount": "--pseudocount=%r"}
CONFIDENCE_SEED = 20261122
TOLERANCE = 1e-9
REPEATS = performance.REPEATS


def write_taxonomy(taxonomy, path):
    """Feature id -> raw taxonomy string (and Confidence, when present), as
    smdi_microbiome.R reads it."""
    taxa = _taxon_series(taxonomy)
    out = pd.DataFrame({"feature": taxa.index.astype(str),
                        "Taxon": taxa.values})
    if isinstance(taxonomy, pd.DataFrame) and "Confidence" in taxonomy:
        out["Confidence"] = taxonomy["Confidence"].values
    out.to_csv(path, sep="\t", index=False)


def with_confidence(taxonomy):
    """The taxonomy with a Confidence column: its own, or U(0.5, 1) values
    drawn with a fixed seed when it has none."""
    if _confidence_series(taxonomy) is not None:
        return taxonomy
    taxonomy = taxonomy.copy()
    rng = np.random.default_rng(CONFIDENCE_SEED)
    taxonomy["Confidence"] = rng.uniform(0.5, 1.0, len(taxonomy))
    return taxonomy


def write_inputs(table, taxonomy, index, workdir):
    """Input files of smdi_microbiome.R: counts, raw taxonomy, synonyms and
    signature (the last two are the plugin's data files, unchanged)."""
    paths = [os.path.join(workdir, n) for n in ("counts.tsv", "taxonomy.tsv")]
    counts = pd.DataFrame(table.matrix_data.toarray(),
                          index=pd.Index(table.ids(axis="observation"),
                                         name="feature"),
                          columns=table.ids(axis="sample"))
    counts.to_csv(paths[0], sep="\t")
    write_taxonomy(taxonomy, paths[1])
    return paths + [SYNONYMS, INDICES[index]]


def run_r(rscript, table, taxonomy, index, params=None):
    """SMDI from R (Series), R's units (Series by feature), R's compute
    seconds and the peak RSS in MB. `params` are collapse_complexes
    parameters (SETTINGS), passed as the matching R options."""
    options = [R_OPTIONS[k] % v for k, v in (params or {}).items()]
    with tempfile.TemporaryDirectory() as workdir:
        inputs = write_inputs(table, taxonomy, index, workdir)
        out = os.path.join(workdir, "smdi.tsv")
        units_out = os.path.join(workdir, "units.tsv")
        seconds, peak = timing(run_child(
            [rscript, R_SCRIPT] + inputs + [out, units_out] + options))
        smdi = pd.read_csv(out, sep="\t", index_col=0)["SMDI"]
        units = read_units(units_out)
    smdi.index = smdi.index.astype(str)
    return smdi, units, seconds, peak


def read_units(path):
    """Feature -> unit table written by smdi_microbiome.R."""
    units = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                        quoting=3)
    return pd.Series(units["unit"].values, index=units["feature"].values)


def unit_mismatches(table, taxonomy, index, synonyms, r_units, exclude):
    """Features whose unit (taxon key) differs between plugin and R;
    `exclude` are the plugin's low-confidence features."""
    level = load_signature(INDICES[index])["level"].iloc[0].strip().lower()
    plugin = pd.Series(feature_units(table, taxonomy, synonyms, level,
                                     exclude),
                       index=list(table.ids(axis="observation")))
    return int((plugin != r_units.reindex(plugin.index)).sum())


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


def compare(plugin, r, unit_mismatch=0):
    """Max absolute difference; NaN must match NaN and every feature must
    have the same unit on both sides."""
    r = r.reindex(plugin.index)
    nan_mismatch = int((plugin.isna() != r.isna()).sum())
    both = plugin.notna() & r.notna()
    diff = (plugin[both] - r[both]).abs()
    return {"samples": len(plugin),
            "nan_samples": int(plugin.isna().sum()),
            "nan_mismatch": nan_mismatch,
            "unit_mismatch": unit_mismatch,
            "max_abs_diff": float(diff.max()) if both.any() else np.nan,
            "pass": bool(nan_mismatch == 0 and unit_mismatch == 0
                         and (not both.any() or diff.max() < TOLERANCE))}


def check_table(rscript, source, classifier, table, taxonomy, synonyms):
    taxonomy = with_confidence(taxonomy)
    custom = load_signature(CUSTOM_SIGNATURE)
    rows = []
    for setting, params in SETTINGS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            _, indices, mapping = collapse_complexes(
                table, taxonomy, smdi_signature=custom, **params)
        exclude = set(mapping.index[mapping["Resolution"]
                                    == "low-confidence"])
        for index in INDICES:
            r_smdi, r_units, _, _ = run_r(rscript, table, taxonomy, index,
                                          params)
            mismatch = unit_mismatches(table, taxonomy, index, synonyms,
                                       r_units, exclude)
            rows.append({"source": source, "classifier": classifier,
                         "setting": setting, "index": index,
                         "features": table.shape[0],
                         "excluded": len(exclude),
                         **compare(indices[index], r_smdi, mismatch)})
            print(source, classifier, setting, index,
                  rows[-1]["max_abs_diff"], flush=True)
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
                _, _, seconds, peak = run_r(rscript, table, taxonomy,
                                            "SMDI")
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
