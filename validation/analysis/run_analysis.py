"""Runs the pre-specified analyses (docs/ANALYSIS_PLAN.md, sections 6-7).

    python -I run_analysis.py <results_dir>

Reads <results_dir>/samples/*.tsv (collect.py) and writes one TSV per
analysis to <results_dir>/tables/:

    discrimination.tsv     Q1  per-cohort subject-level AUC
    meta_auc.tsv           Q1  random-effects pooled AUCs (+ sensitivity)
    meta_auc_loo.tsv       Q1  leave-one-cohort-out, primary measure
    delong.tsv             Q2  paired DeLong SMDI (and SMDI_simplified) vs Red
    meta_delong.tsv        Q2  pooled AUC differences
    cutoffs.tsv            Q3  sensitivity/specificity at published cut-offs
    cutoffs_pooled.tsv     Q3  pooled sensitivity/specificity
    therapy.tsv            Q4  paired change after therapy
    therapy_arms.tsv       Q4  change, treated vs untreated arm
    prognosis.tsv          Q5  baseline AUC for non-response
    clinical.tsv           Q6  Spearman with clinical measures and stage
    classifier_agreement.tsv   nb vs vsearch
    regions.tsv            same samples, two 16S regions

Benjamini-Hochberg q-values are computed within each family (one question,
one classifier). Nothing in this script chooses cohorts, cut-offs or
measures from the data.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stats as st  # noqa: E402
from common import (CLASSIFIERS, CLINICAL, COMPARATOR, CUTOFFS,  # noqa: E402
                    MEASURES, N_BOOT, PRIMARY, SEED, load_cohorts,
                    load_samples, subject_case_control)

results = sys.argv[1]
tables = os.path.join(results, "tables")
os.makedirs(tables, exist_ok=True)
cohorts = load_cohorts()
qc = pd.read_csv(os.path.join(results, "qc.tsv"), sep="\t")
available = {(r.cohort, r.classifier) for r in qc.itertuples()
             if r.processed == "yes"}


def save(df, name):
    df.to_csv(os.path.join(tables, name), sep="\t", index=False,
              float_format="%.6g")
    return df


def add_q(df, family_cols=("classifier",), p="p"):
    df = df.copy()
    df["q"] = np.nan
    if df.empty:
        return df
    for _, idx in df.groupby(list(family_cols)).groups.items():
        df.loc[idx, "q"] = st.bh(df.loc[idx, p].to_numpy())
    return df


# ----------------------------------------------------------------- Q1
disc_rows, subjects = [], {}
for c in cohorts[cohorts["case_control"] == "yes"].itertuples():
    for clf in CLASSIFIERS:
        if (c.cohort, clf) not in available:
            continue
        subj = subject_case_control(load_samples(results, c.cohort, clf))
        subjects[(c.cohort, clf)] = subj
        for m in MEASURES:
            pos = subj.loc[subj["status"] == "periodontitis", m].dropna()
            neg = subj.loc[subj["status"] == "health", m].dropna()
            if len(pos) < 2 or len(neg) < 2:
                continue
            a = st.auc(pos, neg)
            lo, hi, bse = st.bootstrap_auc(pos, neg, N_BOOT, SEED)
            disc_rows.append({
                "cohort": c.cohort, "classifier": clf, "meta": c.meta,
                "independent": c.independent, "region": c.region,
                "measure": m, "n_health": len(neg),
                "n_periodontitis": len(pos),
                "median_health": neg.median(),
                "median_periodontitis": pos.median(), "auc": a,
                "ci_low": lo, "ci_high": hi, "se_bootstrap": bse,
                "se_hanley_mcneil": st.hanley_mcneil_se(a, len(pos), len(neg)),
                "p": sps.mannwhitneyu(pos, neg,
                                      alternative="two-sided").pvalue})
disc = save(pd.DataFrame(disc_rows), "discrimination.tsv")

flagged = set(qc.loc[qc["qc_flags"].fillna("") != "", "cohort"])
meta_rows, loo_rows = [], []
SETS = {"independent": disc["meta"] == "primary",
        "non-independent": disc["meta"] == "non-independent"}
for set_name, mask in SETS.items():
    for (clf, m), g in disc[mask].groupby(["classifier", "measure"]):
        variants = [("REML+HKSJ", "se_hanley_mcneil", "REML", True, g),
                    ("DL", "se_hanley_mcneil", "DL", False, g),
                    ("REML+HKSJ, bootstrap SE", "se_bootstrap", "REML",
                     True, g)]
        if set_name == "independent":
            variants.append(("REML+HKSJ, QC-flagged excluded",
                             "se_hanley_mcneil", "REML", True,
                             g[~g["cohort"].isin(flagged)]))
        for label, se_col, method, hk, gg in variants:
            if len(gg) < 3:   # plan, section 8: no pooled estimate
                continue
            res = st.random_effects_auc(gg["auc"], gg[se_col], method, hk)
            meta_rows.append({"set": set_name, "classifier": clf,
                              "measure": m, "analysis": label,
                              "cohorts": ";".join(gg["cohort"]), **res})
        if set_name == "independent" and m == PRIMARY and len(g) >= 4:
            for left in g["cohort"]:
                gg = g[g["cohort"] != left]
                res = st.random_effects_auc(gg["auc"],
                                            gg["se_hanley_mcneil"])
                loo_rows.append({"classifier": clf, "left_out": left, **res})
save(pd.DataFrame(meta_rows), "meta_auc.tsv")
save(pd.DataFrame(loo_rows), "meta_auc_loo.tsv")

# ----------------------------------------------------------------- Q2
dl_rows = []
for (cohort, clf), subj in subjects.items():
    meta_set = cohorts.set_index("cohort").loc[cohort, "meta"]
    for m in [PRIMARY, "SMDI_simplified"]:
        d = subj.dropna(subset=[m, COMPARATOR])
        pos = d[d["status"] == "periodontitis"]
        neg = d[d["status"] == "health"]
        diff, se, z, p = st.delong_paired(pos[m], neg[m], pos[COMPARATOR],
                                          neg[COMPARATOR])
        dl_rows.append({"cohort": cohort, "classifier": clf,
                        "meta": meta_set, "measure": m,
                        "comparator": COMPARATOR,
                        "auc_measure": st.auc(pos[m], neg[m]),
                        "auc_comparator": st.auc(pos[COMPARATOR],
                                                 neg[COMPARATOR]),
                        "difference": diff, "se": se, "z": z, "p": p})
delong = add_q(pd.DataFrame(dl_rows))
save(delong, "delong.tsv")
md_rows = []
for set_name, val in [("independent", "primary"),
                      ("non-independent", "non-independent")]:
    for (clf, m), g in delong[delong["meta"] == val].groupby(
            ["classifier", "measure"]):
        g = g[g["se"] > 0]
        if len(g) < 3:
            continue
        res = st.random_effects(g["difference"], g["se"] ** 2)
        md_rows.append({"set": set_name, "classifier": clf, "measure": m,
                        "comparator": COMPARATOR,
                        "cohorts": ";".join(g["cohort"]), **res})
save(pd.DataFrame(md_rows), "meta_delong.tsv")

# ----------------------------------------------------------------- Q3
cut_rows = []
for (cohort, clf), subj in subjects.items():
    meta_set = cohorts.set_index("cohort").loc[cohort, "meta"]
    for m, cut in CUTOFFS:
        d = subj.dropna(subset=[m])
        truth = (d["status"] == "periodontitis").to_numpy()
        pred = (d[m] > cut).to_numpy()
        tp, fn = int((pred & truth).sum()), int((~pred & truth).sum())
        tn, fp = int((~pred & ~truth).sum()), int((pred & ~truth).sum())
        se_lo, se_hi = st.wilson(tp, tp + fn)
        sp_lo, sp_hi = st.wilson(tn, tn + fp)
        cut_rows.append({"cohort": cohort, "classifier": clf,
                         "meta": meta_set, "measure": m, "cutoff": cut,
                         "tp": tp, "fn": fn, "tn": tn, "fp": fp,
                         "sensitivity": tp / (tp + fn), "sens_low": se_lo,
                         "sens_high": se_hi, "specificity": tn / (tn + fp),
                         "spec_low": sp_lo, "spec_high": sp_hi,
                         "accuracy": (tp + tn) / len(d)})
cuts = save(pd.DataFrame(cut_rows), "cutoffs.tsv")
cp_rows = []
for set_name, val in [("independent", "primary"),
                      ("non-independent", "non-independent")]:
    for (clf, m, cut), g in cuts[cuts["meta"] == val].groupby(
            ["classifier", "measure", "cutoff"]):
        if len(g) < 3:
            continue
        sens = st.pooled_proportion(g["tp"], g["tp"] + g["fn"])
        spec = st.pooled_proportion(g["tn"], g["tn"] + g["fp"])
        cp_rows.append({"set": set_name, "classifier": clf, "measure": m,
                        "cutoff": cut, "k": len(g),
                        "sensitivity": sens["estimate"],
                        "sens_low": sens["ci_low"],
                        "sens_high": sens["ci_high"], "sens_i2": sens["i2"],
                        "specificity": spec["estimate"],
                        "spec_low": spec["ci_low"],
                        "spec_high": spec["ci_high"], "spec_i2": spec["i2"]})
save(pd.DataFrame(cp_rows), "cutoffs_pooled.tsv")


# ----------------------------------------------------------------- Q4
def paired_units(df, follow, arm=None):
    """Baseline and follow-up values per paired unit (samples of one unit
    at one time averaged)."""
    unit = df["unit"].where(df["unit"].fillna("") != "", df["subject"])
    d = df.assign(unit=unit)
    keys = ["unit"] + ([arm] if arm else [])
    per = d.groupby(keys + ["time"])[MEASURES].mean()
    times = per.index.get_level_values("time")
    if "baseline" not in times or follow not in times:
        return None
    b = per.xs("baseline", level="time")
    f = per.xs(follow, level="time")
    both = b.join(f, lsuffix="_b", rsuffix="_f", how="inner")
    return both.reset_index()


th_rows, arm_rows = [], []
for c in cohorts[cohorts["therapy_followups"] != ""].itertuples():
    arm = c.therapy_arm_column or None
    for clf in CLASSIFIERS:
        if (c.cohort, clf) not in available:
            continue
        df = load_samples(results, c.cohort, clf)
        if c.case_control == "yes":
            df = df[df["status"] == "periodontitis"]
        for follow in c.therapy_followups.split(","):
            both = paired_units(df, follow, arm)
            if both is None:
                print("WARNING: %s %s: no baseline or %r samples; follow-up "
                      "skipped" % (c.cohort, clf, follow), file=sys.stderr)
                continue
            strata = ([("all", both)] if arm is None
                      else list(both.groupby(arm)))
            for stratum, g in strata:
                for m in MEASURES:
                    d = (g[m + "_f"] - g[m + "_b"]).dropna()
                    if len(d) < 3:
                        continue
                    p = (sps.wilcoxon(d).pvalue if (d != 0).any()
                         else np.nan)
                    th_rows.append({
                        "cohort": c.cohort, "classifier": clf,
                        "independent": c.independent, "follow_up": follow,
                        "arm": stratum, "measure": m, "n": len(d),
                        "median_baseline": g[m + "_b"].median(),
                        "median_follow_up": g[m + "_f"].median(),
                        "median_change": d.median(),
                        "fraction_decreased": (d < 0).mean(),
                        "rank_biserial": st.rank_biserial_paired(d), "p": p})
            if arm is not None:
                labels = sorted(both[arm].unique())
                treated = [x for x in labels
                           if x in ("intervention", "antibiotics")]
                if len(labels) == 2 and treated:
                    other = [x for x in labels if x != treated[0]][0]
                    for m in MEASURES:
                        dt = (both.loc[both[arm] == treated[0], m + "_f"]
                              - both.loc[both[arm] == treated[0], m + "_b"])
                        do = (both.loc[both[arm] == other, m + "_f"]
                              - both.loc[both[arm] == other, m + "_b"])
                        dt, do = dt.dropna(), do.dropna()
                        if len(dt) < 3 or len(do) < 3:
                            continue
                        arm_rows.append({
                            "cohort": c.cohort, "classifier": clf,
                            "follow_up": follow, "treated": treated[0],
                            "reference": other, "measure": m,
                            "n_treated": len(dt), "n_reference": len(do),
                            "median_change_treated": dt.median(),
                            "median_change_reference": do.median(),
                            # AUC = P(change treated < change reference)
                            "auc_larger_decrease": st.auc(-dt, -do),
                            "p": sps.mannwhitneyu(
                                dt, do, alternative="two-sided").pvalue})
therapy = pd.DataFrame(th_rows)
therapy = add_q(therapy, ("classifier", "independent"))
save(therapy, "therapy.tsv")
save(add_q(pd.DataFrame(arm_rows)), "therapy_arms.tsv")

# ----------------------------------------------------------------- Q5
pr_rows = []
for c in cohorts[cohorts["prognosis"] != ""].itertuples():
    for clf in CLASSIFIERS:
        if (c.cohort, clf) not in available:
            continue
        df = load_samples(results, c.cohort, clf)
        df = df[df["time"] == "baseline"]
        unit = df["unit"].where(df["unit"].fillna("") != "", df["subject"])
        per = df.assign(unit=unit).groupby(["unit", "response"])[
            MEASURES].mean().reset_index()
        positive = [r for r in per["response"].unique()
                    if r in ("NR", "non-responder")]
        if not positive:
            continue
        for m in MEASURES:
            pos = per.loc[per["response"] == positive[0], m].dropna()
            neg = per.loc[per["response"] != positive[0], m].dropna()
            if len(pos) < 2 or len(neg) < 2:
                continue
            lo, hi, _ = st.bootstrap_auc(pos, neg, N_BOOT, SEED)
            pr_rows.append({"cohort": c.cohort, "classifier": clf,
                            "independent": c.independent, "measure": m,
                            "n_nonresponders": len(pos),
                            "n_responders": len(neg),
                            "auc_nonresponse": st.auc(pos, neg),
                            "ci_low": lo, "ci_high": hi,
                            "p": sps.mannwhitneyu(
                                pos, neg, alternative="two-sided").pvalue})
save(add_q(pd.DataFrame(pr_rows)), "prognosis.tsv")

# ----------------------------------------------------------------- Q6
cl_rows = []
for c in cohorts[cohorts["clinical"] != "no"].itertuples():
    for clf in CLASSIFIERS:
        if (c.cohort, clf) not in available:
            continue
        df = load_samples(results, c.cohort, clf)
        df = df[(df["time"] == "baseline") & (df["status"] != "gingivitis")
                & (df["status"] != "other")]
        unit = df["unit"].where(df["unit"].fillna("") != "", df["subject"])
        per = df.assign(unit=unit).groupby("unit")[
            MEASURES + CLINICAL].mean()
        for var in CLINICAL:
            if per[var].notna().sum() < 10:
                continue
            for m in MEASURES:
                d = per[[m, var]].dropna()
                rho, p = sps.spearmanr(d[m], d[var])
                cl_rows.append({"cohort": c.cohort, "classifier": clf,
                                "independent": c.independent,
                                "clinical": var, "measure": m, "n": len(d),
                                "rho": rho, "p": p})
save(add_q(pd.DataFrame(cl_rows)), "clinical.tsv")

# --------------------------------------------------------- robustness
ag_rows = []
for c in cohorts.itertuples():
    if not all((c.cohort, k) in available for k in CLASSIFIERS):
        continue
    a = load_samples(results, c.cohort, "nb").set_index("sample-id")
    b = load_samples(results, c.cohort, "vsearch").set_index("sample-id")
    for m in MEASURES:
        d = pd.concat([a[m], b[m]], axis=1, keys=["nb", "vsearch"]).dropna()
        diff = d["vsearch"] - d["nb"]
        ag_rows.append({"cohort": c.cohort, "measure": m, "n": len(d),
                        "spearman": sps.spearmanr(d["nb"],
                                                  d["vsearch"])[0],
                        "median_difference": diff.median(),
                        "iqr_low": diff.quantile(0.25),
                        "iqr_high": diff.quantile(0.75)})
save(pd.DataFrame(ag_rows), "classifier_agreement.tsv")

REGION_PAIRS = [("PRJNA773195", "PRJNA773202", "sample_code"),
                ("SRP009299_V12", "SRP009299_V4", "library")]
rg_rows = []
for first, second, key in REGION_PAIRS:
    for clf in CLASSIFIERS:
        if (first, clf) not in available or (second, clf) not in available:
            continue
        per = []
        for co in (first, second):
            df = load_samples(results, co, clf)
            df = df[df["status"].isin(["health", "periodontitis"])]
            per.append(df.groupby([key, "status", "subject"])[
                MEASURES].mean().reset_index().set_index(key))
        both = per[0].join(per[1][MEASURES], lsuffix="_1", rsuffix="_2",
                           how="inner")
        for m in [PRIMARY, "SMDI_simplified", COMPARATOR]:
            d = both.dropna(subset=[m + "_1", m + "_2"])
            shift = d[m + "_2"] - d[m + "_1"]
            subj = d.groupby(["subject", "status"])[
                [m + "_1", m + "_2"]].mean()
            subj = subj.reset_index()
            aucs = [st.auc(subj.loc[subj["status"] == "periodontitis", col],
                           subj.loc[subj["status"] == "health", col])
                    for col in (m + "_1", m + "_2")]
            rg_rows.append({
                "region_1": first, "region_2": second, "classifier": clf,
                "measure": m, "n_samples": len(d),
                "median_1": d[m + "_1"].median(),
                "median_2": d[m + "_2"].median(),
                "median_shift_2_minus_1": shift.median(),
                "wilcoxon_p": sps.wilcoxon(shift).pvalue,
                "spearman": sps.spearmanr(d[m + "_1"], d[m + "_2"])[0],
                "auc_1": aucs[0], "auc_2": aucs[1]})
save(pd.DataFrame(rg_rows), "regions.tsv")
print("tables written to", tables)
