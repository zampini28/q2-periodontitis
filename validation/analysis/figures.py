"""Figures of the validation, from the result tables.

    python -I figures.py <results_dir>

Writes PNG (300 dpi) and SVG to <results_dir>/figures/:

    fig1_forest     per-cohort and pooled subject-level AUC, SMDI and Red
    fig2_roc        ROC curves per case-control cohort, SMDI and Red
    fig3_robustness classifier, pseudocount and rarefaction
    fig4_therapy    change of SMDI and Red after therapy
    fig5_performance runtime and peak memory on synthetic tables
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import (COMPARATOR, PRIMARY, PRIMARY_CLASSIFIER,  # noqa: E402
                    load_cohorts, load_samples, subject_case_control)
from performance import plot_performance  # noqa: E402

# Reference categorical palette (validated: scripts/validate_palette.js)
COLOR = {"SMDI": "#2a78d6", "complex_Red": "#eb6834",
         "SMDI_simplified": "#1baf7a"}
MARKER = {"SMDI": "o", "complex_Red": "s", "SMDI_simplified": "^"}
LABEL = {"SMDI": "SMDI", "complex_Red": "Red complex",
         "SMDI_simplified": "SMDI simplified"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "svg.fonttype": "none"})

results = sys.argv[1]
tables = os.path.join(results, "tables")
figdir = os.path.join(results, "figures")
os.makedirs(figdir, exist_ok=True)
cohorts = load_cohorts().set_index("cohort")


def read(name):
    path = os.path.join(tables, name)
    try:
        return pd.read_csv(path, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame()


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(figdir, "%s.%s" % (name, ext)), dpi=300,
                    bbox_inches="tight")
    plt.close(fig)


def label(cohort):
    c = cohorts.loc[cohort]
    return "%s (%s)" % (cohort, c["region"])


# ------------------------------------------------------------ figure 1
disc = read("discrimination.tsv")
meta = read("meta_auc.tsv")
clf = PRIMARY_CLASSIFIER
sets = [("independent", "primary", "Independent cohorts"),
        ("non-independent", "non-independent",
         "Cohorts used to derive or test the SMDI")]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharex=True)
for ax, (set_name, val, title) in zip(axes, sets):
    d = disc[(disc["classifier"] == clf) & (disc["meta"] == val)]
    names = list(dict.fromkeys(d["cohort"]))
    y = np.arange(len(names))[::-1] + 2
    for k, m in enumerate([PRIMARY, COMPARATOR]):
        off = 0.15 - 0.3 * k
        g = d[d["measure"] == m].set_index("cohort").reindex(names)
        ax.errorbar(g["auc"], y + off,
                    xerr=[g["auc"] - g["ci_low"], g["ci_high"] - g["auc"]],
                    fmt=MARKER[m], color=COLOR[m], ms=4, lw=1, capsize=0,
                    label=LABEL[m])
        pooled = meta[(meta["set"] == set_name) & (meta["classifier"] == clf)
                      & (meta["measure"] == m)
                      & (meta["analysis"] == "REML+HKSJ")] \
            if not meta.empty else pd.DataFrame()
        if not pooled.empty:
            p = pooled.iloc[0]
            yy = 0.6 + off
            ax.fill([p["ci_low"], p["estimate"], p["ci_high"], p["estimate"]],
                    [yy, yy + 0.18, yy, yy - 0.18], color=COLOR[m], lw=0)
            ax.plot([p["pi_low"], p["pi_high"]], [yy, yy], color=COLOR[m],
                    lw=0.8, ls=":")
    ax.set_yticks(list(y) + ([0.6] if not meta.empty else []))
    ax.set_yticklabels([label(n) for n in names]
                       + (["Pooled (REML, HKSJ)"] if not meta.empty else []))
    ax.axvline(0.5, color=INK2, lw=0.8, ls="--")
    ax.set_xlim(0.3, 1.0)
    ax.set_title(title, fontsize=8, color=INK, loc="left")
axes[0].legend(loc="lower left", fontsize=7)
fig.supxlabel("Subject-level AUC, health vs periodontitis (95%% CI); "
              "classifier: %s" % clf, fontsize=8)
fig.tight_layout()
save(fig, "fig1_forest")

# ------------------------------------------------------------ figure 2


def roc(pos, neg):
    """Empirical ROC: higher index = periodontitis."""
    thr = np.unique(np.concatenate([pos, neg]))[::-1]
    tpr = [0.0] + [(pos >= t).mean() for t in thr]
    fpr = [0.0] + [(neg >= t).mean() for t in thr]
    return np.array(fpr), np.array(tpr)


cc = [c for c in cohorts.index if cohorts.loc[c, "case_control"] == "yes"
      and load_samples(results, c, clf) is not None]
ncol = 4
nrow = int(np.ceil(len(cc) / ncol)) or 1
fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 1.9 * nrow),
                         sharex=True, sharey=True, squeeze=False)
for ax, c in zip(axes.flat, cc):
    subj = subject_case_control(load_samples(results, c, clf))
    for m, ls in [(PRIMARY, "-"), (COMPARATOR, "--")]:
        s = subj.dropna(subset=[m])
        fpr, tpr = roc(s.loc[s["status"] == "periodontitis", m].to_numpy(),
                       s.loc[s["status"] == "health", m].to_numpy())
        a = disc[(disc["cohort"] == c) & (disc["classifier"] == clf)
                 & (disc["measure"] == m)]["auc"]
        ax.step(fpr, tpr, where="post", color=COLOR[m], lw=1.4, ls=ls,
                label="%s %.2f" % (LABEL[m], a.iloc[0]) if len(a) else m)
    ax.plot([0, 1], [0, 1], color=GRID, lw=0.8)
    indep = "" if cohorts.loc[c, "independent"] == "yes" else " *"
    ax.set_title(label(c) + indep, fontsize=7, color=INK)
    ax.legend(loc="lower right", fontsize=6, handlelength=1.5)
    ax.set_aspect("equal")
for ax in axes.flat[len(cc):]:
    ax.axis("off")
for ax in axes[-1]:
    ax.set_xlabel("1 - specificity")
for ax in axes[:, 0]:
    ax.set_ylabel("Sensitivity")
fig.tight_layout()
save(fig, "fig2_roc")

# ------------------------------------------------------------ figure 3
eps = read("robustness_epsilon.tsv")
rar = read("robustness_rarefaction.tsv")
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8))
d = disc[(disc["measure"] == PRIMARY)
         & disc["meta"].isin(["primary", "non-independent", "region-only"])]
piv = d.pivot_table(index="cohort", columns="classifier", values="auc")
piv = piv.dropna()
ax = axes[0]
for i, (c, r) in enumerate(piv.iterrows()):
    ax.plot([r["nb"], r["vsearch"]], [i, i], color=GRID, lw=2)
ax.scatter(piv["nb"], range(len(piv)), color=COLOR["SMDI"], marker="o", s=14,
           label="naive Bayes", zorder=3)
ax.scatter(piv["vsearch"], range(len(piv)), facecolor="white",
           edgecolor=COLOR["SMDI"], marker="o", s=14, label="VSEARCH",
           zorder=3)
ax.set_yticks(range(len(piv)))
ax.set_yticklabels(piv.index, fontsize=6)
ax.set_xlabel("SMDI AUC")
ax.set_title("a  Classifier", fontsize=8, loc="left", color=INK)
ax.legend(fontsize=6, loc="lower left")
ax = axes[1]
if not eps.empty:
    e = eps[(eps["classifier"] == clf) & (eps["measure"] == PRIMARY)]
    order = ["original", "1e-06", "1e-05", "0.0001"]
    e = e.assign(epsilon=e["epsilon"].astype(str))
    for c, g in e.groupby("cohort"):
        g = g.set_index("epsilon").reindex(order)
        ax.plot(range(len(order)), g["auc"], color=COLOR["SMDI"], lw=0.8,
                alpha=0.6, marker="o", ms=2.5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["table\n(original)", "1e-6", "1e-5", "1e-4"])
    ax.set_xlabel("SMDI pseudocount")
    ax.set_ylabel("SMDI AUC")
ax.set_title("b  Pseudocount", fontsize=8, loc="left", color=INK)
ax = axes[2]
if not rar.empty:
    r = rar[(rar["classifier"] == clf)]
    for m in [PRIMARY, "SMDI_simplified", COMPARATOR]:
        g = r[r["measure"] == m]
        mean = g.groupby(["cohort", "depth"])["auc"].mean().reset_index()
        delta = mean.merge(g.groupby("cohort")["auc_unrarefied"].first(),
                           on="cohort")
        delta["d"] = delta["auc"] - delta["auc_unrarefied"]
        med = delta.groupby("depth")["d"].median()
        ax.plot(med.index, med.values, color=COLOR[m], marker=MARKER[m],
                ms=4, lw=1.4, label=LABEL[m])
    ax.axhline(0, color=INK2, lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xticks([1000, 5000, 10000])
    ax.set_xticklabels(["1k", "5k", "10k"])
    ax.set_xlabel("Rarefaction depth (reads)")
    ax.set_ylabel("Median change in AUC vs unrarefied")
    ax.legend(fontsize=6)
ax.set_title("c  Rarefaction", fontsize=8, loc="left", color=INK)
fig.tight_layout()
save(fig, "fig3_robustness")

# ------------------------------------------------------------ figure 4
th = read("therapy.tsv")
if not th.empty:
    t = th[(th["classifier"] == clf) & th["measure"].isin([PRIMARY,
                                                           COMPARATOR])]
    t = t.assign(key=t["cohort"] + "  " + t["follow_up"]
                 + np.where(t["arm"] != "all", "  " + t["arm"].astype(str),
                            ""))
    keys = list(dict.fromkeys(t["key"]))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 0.28 * len(keys) + 1.2),
                             sharey=True)
    for ax, m in zip(axes, [PRIMARY, COMPARATOR]):
        g = t[t["measure"] == m].set_index("key").reindex(keys)
        y = np.arange(len(keys))[::-1]
        ax.barh(y, g["rank_biserial"], color=COLOR[m], height=0.6)
        for yy, (_, r) in zip(y, g.iterrows()):
            if pd.notna(r["q"]):
                ax.text(1.02, yy, "q=%.2g" % r["q"], va="center", fontsize=6,
                        color=INK2, transform=ax.get_yaxis_transform())
        ax.axvline(0, color=INK2, lw=0.8)
        ax.set_xlim(-1, 1)
        ax.set_yticks(y)
        ax.set_yticklabels(keys, fontsize=6)
        ax.set_xlabel("Rank-biserial r of paired change\n(< 0: decreased)")
        ax.set_title(LABEL[m], fontsize=8, loc="left", color=INK)
    fig.tight_layout()
    save(fig, "fig4_therapy")

# ------------------------------------------------------------ figure 5
perf = read("performance.tsv")
if not perf.empty:
    plot_performance(perf, figdir)
print("figures written to", figdir)
