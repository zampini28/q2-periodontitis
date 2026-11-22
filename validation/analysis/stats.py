"""Statistics of the multi-cohort validation (docs/ANALYSIS_PLAN.md, sec. 6).

All functions take plain NumPy arrays. "pos" are the cases (periodontitis,
non-responders) and "neg" the controls; AUC = P(pos > neg) + 0.5 P(tie).
"""
import numpy as np
from scipy import stats


# ------------------------------------------------------------------ AUC

def auc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return u / (len(pos) * len(neg))


def bootstrap_auc(pos, neg, n_boot=2000, seed=0):
    """Stratified bootstrap: cases and controls resampled separately.
    Returns (2.5 %, 97.5 % percentiles, SD of the bootstrap AUCs)."""
    rng = np.random.default_rng(seed)
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    boot = np.array([auc(rng.choice(pos, len(pos)), rng.choice(neg, len(neg)))
                     for _ in range(n_boot)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return lo, hi, boot.std(ddof=1)


def hanley_mcneil_se(a, n_pos, n_neg):
    """SE of an AUC (Hanley & McNeil 1982, eq. 1-3)."""
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    var = (a * (1 - a) + (n_pos - 1) * (q1 - a * a)
           + (n_neg - 1) * (q2 - a * a)) / (n_pos * n_neg)
    return np.sqrt(var)


def _midrank(x):
    """Mid-ranks (1-based) of x, ties averaged."""
    return stats.rankdata(x, method="average")


def delong(pos, neg):
    """AUCs and their covariance for k measures on the same subjects
    (DeLong et al. 1988, fast algorithm of Sun & Xu 2014).

    pos: k x m array (cases), neg: k x n array (controls).
    Returns (aucs [k], covariance [k x k]).
    """
    pos = np.atleast_2d(np.asarray(pos, float))
    neg = np.atleast_2d(np.asarray(neg, float))
    k, m = pos.shape
    n = neg.shape[1]
    tx = np.array([_midrank(r) for r in pos])
    ty = np.array([_midrank(r) for r in neg])
    tz = np.array([_midrank(r) for r in np.hstack([pos, neg])])
    aucs = (tz[:, :m].sum(axis=1) / m - (m + 1) / 2) / n
    v01 = (tz[:, :m] - tx) / n          # structural components of cases
    v10 = 1 - (tz[:, m:] - ty) / m      # of controls
    sx = np.atleast_2d(np.cov(v01))
    sy = np.atleast_2d(np.cov(v10))
    return aucs, sx / m + sy / n


def delong_paired(pos_a, neg_a, pos_b, neg_b):
    """Paired DeLong test of AUC(a) - AUC(b) (same subjects, same order).
    Returns (difference, SE, z, two-sided p)."""
    aucs, cov = delong(np.vstack([pos_a, pos_b]), np.vstack([neg_a, neg_b]))
    diff = aucs[0] - aucs[1]
    se = np.sqrt(max(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1], 0.0))
    if se == 0:
        return diff, 0.0, np.nan, np.nan
    z = diff / se
    return diff, se, z, 2 * stats.norm.sf(abs(z))


# ------------------------------------------------------- meta-analysis

def _q_statistic(y, v):
    w = 1 / v
    mu = np.sum(w * y) / np.sum(w)
    return np.sum(w * (y - mu) ** 2)


def i_squared(y, v):
    k = len(y)
    q = _q_statistic(y, v)
    return max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0


def tau2_dl(y, v):
    """DerSimonian-Laird estimate of the between-study variance."""
    w = 1 / v
    q = _q_statistic(y, v)
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    return max(0.0, (q - (len(y) - 1)) / c)


def tau2_reml(y, v, tol=1e-12, max_iter=1000):
    """REML estimate of the between-study variance, bounded at 0, by the
    fixed-point iteration of Viechtbauer (2005, eq. 12), started at the
    DerSimonian-Laird estimate."""
    tau2 = tau2_dl(y, v)
    for _ in range(max_iter):
        w = 1 / (v + tau2)
        mu = np.sum(w * y) / np.sum(w)
        new = max(0.0, np.sum(w ** 2 * ((y - mu) ** 2 - v)) / np.sum(w ** 2)
                  + 1 / np.sum(w))
        if abs(new - tau2) < tol:
            return new
        tau2 = new
    raise RuntimeError("REML did not converge")


def random_effects(y, v, method="REML", hksj=True):
    """Random-effects pooled estimate.

    method: "REML" or "DL" (between-study variance).
    hksj: Hartung-Knapp-Sidik-Jonkman CI (t with k-1 df) instead of Wald.
    Returns a dict with estimate, se, ci, tau2, i2, Q, p_Q, k and the 95 %
    prediction interval (t with k-2 df, Higgins et al. 2009).
    """
    y, v = np.asarray(y, float), np.asarray(v, float)
    k = len(y)
    tau2 = tau2_reml(y, v) if method == "REML" else tau2_dl(y, v)
    w = 1 / (v + tau2)
    mu = np.sum(w * y) / np.sum(w)
    se_wald = np.sqrt(1 / np.sum(w))
    if hksj and k > 1:
        q_hk = np.sum(w * (y - mu) ** 2) / (k - 1)
        se = np.sqrt(q_hk / np.sum(w))
        crit = stats.t.ppf(0.975, k - 1)
        p = 2 * stats.t.sf(abs(mu / se), k - 1) if se > 0 else np.nan
    else:
        se = se_wald
        crit = stats.norm.ppf(0.975)
        p = 2 * stats.norm.sf(abs(mu / se))
    q = _q_statistic(y, v)
    if k > 2:
        half = stats.t.ppf(0.975, k - 2) * np.sqrt(tau2 + se ** 2)
        pi = (mu - half, mu + half)
    else:
        pi = (np.nan, np.nan)
    return {"estimate": mu, "se": se, "ci_low": mu - crit * se,
            "ci_high": mu + crit * se, "p": p, "tau2": tau2,
            "i2": i_squared(y, v), "Q": q,
            "p_Q": stats.chi2.sf(q, k - 1) if k > 1 else np.nan,
            "k": k, "pi_low": pi[0], "pi_high": pi[1],
            "method": method + ("+HKSJ" if hksj else "")}


def _inv_logit(x):
    return 1 / (1 + np.exp(-np.asarray(x, float)))


def random_effects_auc(auc, se, method="REML", hksj=True):
    """Random-effects pooled AUC on the logit scale.

    Each AUC enters as logit(AUC) with the delta-method variance
    se**2 / (AUC (1 - AUC))**2. The estimate, CI and prediction interval
    are back-transformed, so they stay inside (0, 1); tau2, I2 and Q refer
    to the logit scale. AUCs of exactly 0 or 1 have no logit and raise.
    """
    a, se = np.asarray(auc, float), np.asarray(se, float)
    if np.any((a <= 0) | (a >= 1)):
        raise ValueError("AUC must be strictly between 0 and 1 for logit "
                         "pooling")
    y = np.log(a / (1 - a))
    v = (se / (a * (1 - a))) ** 2
    res = random_effects(y, v, method, hksj)
    for key in ("estimate", "ci_low", "ci_high", "pi_low", "pi_high"):
        res[key] = float(_inv_logit(res[key]))
    res["se"] = np.nan   # an SE on the logit scale is not an AUC SE
    res["scale"] = "logit"
    return res


# ------------------------------------------------------- proportions

def wilson(k, n, z=1.959963984540054):
    """Wilson score 95 % interval for k successes out of n."""
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return centre - half, centre + half


def pooled_proportion(ks, ns):
    """Random-effects (DerSimonian-Laird, Wald) pooled proportion on the
    logit scale; 0.5 is added to cells of studies with 0 % or 100 %."""
    ks, ns = np.asarray(ks, float), np.asarray(ns, float)
    edge = (ks == 0) | (ks == ns)
    ks = np.where(edge, ks + 0.5, ks)
    ns = np.where(edge, ns + 1.0, ns)
    y = np.log(ks / (ns - ks))
    v = 1 / ks + 1 / (ns - ks)
    res = random_effects(y, v, method="DL", hksj=False)
    inv = lambda x: 1 / (1 + np.exp(-x))  # noqa: E731
    return {"estimate": inv(res["estimate"]), "ci_low": inv(res["ci_low"]),
            "ci_high": inv(res["ci_high"]), "i2": res["i2"], "k": res["k"]}


# --------------------------------------------------------- paired tests

def rank_biserial_paired(d):
    """Matched-pairs rank-biserial correlation of the differences d
    (zeros dropped, as in the Wilcoxon signed-rank test): positive when
    the values increased."""
    d = np.asarray(d, float)
    d = d[d != 0]
    if len(d) == 0:
        return np.nan
    r = stats.rankdata(np.abs(d))
    return (r[d > 0].sum() - r[d < 0].sum()) / r.sum()


def bh(p):
    """Benjamini-Hochberg q-values (NaN kept as NaN)."""
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q
