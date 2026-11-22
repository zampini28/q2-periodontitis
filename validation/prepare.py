"""Builds the QIIME 2 manifest and the standardised metadata of a cohort.

    python -I prepare.py <cohort> <data_dir>

<data_dir> holds <accession>.runs.tsv (ENA file report, see download.sh) and
<accession>/fastq/. Writes <data_dir>/<cohort>/manifest.tsv and
<data_dir>/<cohort>/metadata.tsv. BioSample records are fetched from the
EBI BioSamples API once and cached in <data_dir>/<accession>/biosamples.json.

Every metadata table has the columns below (empty when not applicable);
cohort-specific columns are kept after them.

    sample-id  run  subject  status  time  unit  arm  response  stage
    probing_depth  clinical_attachment_loss  bleeding_on_probing
    plaque_index

* status: health | gingivitis | periodontitis | other (excluded from
  every contrast, e.g. healthy sites of periodontitis subjects)
* time: baseline | the follow-up label (therapy cohorts); baseline for
  case-control cohorts
* unit: the paired unit of a therapy cohort (subject, or site)
* stage: 0 (health) to 4 (stage IV), when published
"""
import gzip
import json
import os
import re
import sys
import urllib.request
from collections import Counter

import pandas as pd

STANDARD = ["sample-id", "run", "subject", "status", "time", "unit", "arm",
            "response", "stage", "probing_depth", "clinical_attachment_loss",
            "bleeding_on_probing", "plaque_index"]

cohort, data = sys.argv[1:3]
# Cohorts split from one accession (e.g. SRP009299_V4) share its runs
ACCESSION = {"SRP009299_V12": "SRP009299", "SRP009299_V4": "SRP009299"}
acc = ACCESSION.get(cohort, cohort)
runs = pd.read_csv(os.path.join(data, acc + ".runs.tsv"), sep="\t",
                   dtype=str).fillna("")
fastq = os.path.abspath(os.path.join(data, acc, "fastq"))
cache_path = os.path.join(data, acc, "biosamples.json")
cache = {}
if os.path.exists(cache_path):
    with open(cache_path) as fh:
        cache = json.load(fh)


def biosample(sample_accession):
    """BioSample characteristics as {lower-case name: text}."""
    if sample_accession not in cache:
        url = "https://www.ebi.ac.uk/biosamples/samples/" + sample_accession
        with urllib.request.urlopen(url, timeout=60) as fh:
            record = json.load(fh)
        cache[sample_accession] = {
            k: v[0]["text"] for k, v in record.get("characteristics",
                                                   {}).items()}
    return {re.sub(r"[_ ]+", " ", k.lower()): v
            for k, v in cache[sample_accession].items()}


def save_cache():
    if cache:
        with open(cache_path, "w") as fh:
            json.dump(cache, fh, indent=1, sort_keys=True)


def attributes():
    """One row of BioSample attributes per run."""
    rows = []
    for _, r in runs.iterrows():
        rows.append({"run": r["run_accession"],
                     **biosample(r["sample_accession"])})
    save_cache()
    return pd.DataFrame(rows)


def files(run):
    """FASTQ files of a run, as published (1 = single, 2 = paired)."""
    return sorted(f for f in os.listdir(fastq) if f.startswith(run + ".")
                  or f.startswith(run + "_"))


def manifest(df, paired):
    if paired:
        return pd.DataFrame({
            "sample-id": df["sample-id"],
            "forward-absolute-filepath": [
                os.path.join(fastq, r + "_1.fastq.gz") for r in df["run"]],
            "reverse-absolute-filepath": [
                os.path.join(fastq, r + "_2.fastq.gz") for r in df["run"]],
        })
    paths = []
    for r in df["run"]:
        single = [f for f in files(r) if f == r + ".fastq.gz"]
        # Paired runs published with R1 only, or R1 used alone
        paths.append(os.path.join(fastq, single[0] if single
                                  else r + "_1.fastq.gz"))
    return pd.DataFrame({"sample-id": df["sample-id"],
                         "absolute-filepath": paths})


# ------------------------------------------------------------- cohorts

def chen2018():
    """SRP075100: sample titles GROUP_SUBJECT_SITE. HP healthy, D1P
    periodontitis before and D2P after non-surgical therapy; SITE is the
    tooth number with a suffix."""
    df = runs[runs["sample_title"].str.match(r"^(HP|D1P|D2P)_")].copy()
    parts = df["sample_title"].str.split("_", expand=True)
    tooth = parts[2].str.extract(r"^(\d+)")[0]
    out = pd.DataFrame({
        "sample-id": df["sample_title"].to_numpy(),
        "run": df["run_accession"].to_numpy(),
        "group": parts[0].to_numpy(),
    })
    out["status"] = out["group"].map(lambda g: "health" if g == "HP"
                                     else "periodontitis")
    out["time"] = out["group"].map(lambda g: "after-therapy" if g == "D2P"
                                   else "baseline")
    out["subject"] = (parts[0].str[0] + "_" + parts[1]).to_numpy()
    out["unit"] = (out["subject"] + "_" + tooth.to_numpy())
    return manifest(out, paired=False), out


def pei2020():
    """SRP226726: library names N<id> (healthy) and P<id> (periodontitis),
    one subgingival sample per subject."""
    out = pd.DataFrame({"run": runs["run_accession"],
                        "sample-id": runs["library_name"]})
    out["status"] = out["sample-id"].str[0].map(
        {"N": "health", "P": "periodontitis"})
    out["subject"] = out["sample-id"]
    return manifest(out, paired=True), out


def read_prefix(run, n=200):
    """Most common first 8 nt of the first reads of a run."""
    path = os.path.join(fastq, run + ".fastq.gz")
    starts = Counter()
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                starts[line[:8]] += 1
            if i >= 4 * n:
                break
    return starts.most_common(1)[0][0] if starts else ""


def griffen2012(region):
    """SRP009299: 454 reads of V1-V2 (27F remnant, start GTTTGATC) and V4
    (primer-free, start ACTGGG) from the same libraries. Libraries
    PK<id>H / PK<id>D are healthy / diseased sites. As in Chen et al.
    (2022), healthy sites of periodontitis subjects are not used as
    controls (status "other")."""
    rows = []
    for run, sample, alias in zip(runs["run_accession"],
                                  runs["sample_accession"],
                                  runs["library_name"]):
        start = read_prefix(run)
        run_region = "V12" if start.startswith("GTTTGATC") else "V4"
        if run_region != region:
            continue
        a = biosample(sample)
        rows.append({"run": run, "sample-id": "%s-%s" % (alias, run),
                     "library": alias,
                     "subject": a.get("individual", alias[:-1]),
                     "site_type": "diseased" if alias.endswith("D")
                     else "healthy",
                     "condition": a.get("condition", "").lower(),
                     "first_bases": start})
    save_cache()
    out = pd.DataFrame(rows)
    perio = set(out.loc[out["site_type"] == "diseased", "subject"])
    out["status"] = [
        "health" if subj not in perio
        else ("periodontitis" if site == "diseased" else "other")
        for subj, site in zip(out["subject"], out["site_type"])]
    return manifest(out, paired=False), out


def kirst2015():
    """SRP051201: library names Lib_CP.<id>.<A|B> (chronic periodontitis)
    and Lib_H.<id>.<A|B> (healthy); A and B are two sites of one subject."""
    out = pd.DataFrame({"run": runs["run_accession"],
                        "library": runs["library_name"]})
    out["sample-id"] = out["library"].str.replace("Lib_", "", regex=False)
    group = out["sample-id"].str.split(".").str[0]
    out["status"] = group.map({"CP": "periodontitis", "H": "health"})
    out["subject"] = out["sample-id"].str.rsplit(".", n=1).str[0]
    return manifest(out, paired=False), out


STAGE = {"health": 0, "stage i periodontitis": 1, "stage ii periodontitis": 2,
         "stage iii periodontitis": 3, "stage iv periodontitis": 4}


def spanish_stages():
    """PRJNA863881: one subgingival sample per subject; BioSample attributes
    give the periodontal status (health, gingivitis, stage I-IV) and the
    subject's mean probing depth, attachment loss, bleeding on probing and
    plaque index."""
    a = attributes()
    out = pd.DataFrame({
        "run": a["run"],
        "sample-id": runs.set_index("run_accession").loc[
            a["run"], "library_name"].to_numpy(),
        "status_detail": a["periodontal health status"],
        "probing_depth": a.get("probing depth"),
        "clinical_attachment_loss": a.get("clinical attachment loss"),
        "bleeding_on_probing": a.get("bleeding on probing"),
        "plaque_index": a.get("plaque index"),
        "age": a.get("age"), "tobacco": a.get("tobacco"),
    })
    detail = out["status_detail"].str.lower().str.replace("_", " ")
    out["status"] = detail.map(lambda s: s if s in ("health", "gingivitis")
                               else "periodontitis")
    out["stage"] = detail.map(STAGE)
    out["subject"] = out["sample-id"]
    return manifest(out, paired=True), out


def byrne2022():
    """PRJNA786436: paired non-responding (NR) and good-response (GR) sites
    of the same patients, at baseline and 3 months after non-surgical
    therapy, with site-level clinical measures."""
    a = attributes()
    out = pd.DataFrame({
        "run": a["run"],
        "sample-id": "S" + a["sample id"],
        "patient": a["patient"], "tooth_site": a["site"],
        "time": a["time(0=0m, 1=3m)"].map({"0": "baseline", "1": "3m"}),
        "response": a["treatment group(0=nr, 1=gr)"].map(
            {"0": "NR", "1": "GR"}),
        "probing_depth": a["probing depth"],
        "clinical_attachment_loss": a["clinical attachment level"],
        "bleeding_on_probing": a["bleeding on probing"],
        "plaque_index": a["plaque index"],
    })
    out["unit"] = out["patient"] + "_" + out["tooth_site"].str.replace(" ", "")
    out["status"] = "periodontitis"
    out["subject"] = out["patient"]
    return manifest(out, paired=False), out


def johnston2021():
    """PRJNA736618: subgingival plaque of patients at day 0 (before) and
    days 1, 7 and 90 after full-mouth mechanical debridement; patient and
    day are in the experiment titles ("Patient 26 at day 90")."""
    rows = []
    for run, title in zip(runs["run_accession"], runs["experiment_title"]):
        m = re.search(r"Patient (\d+) at day (\d+)", title)
        rows.append({"run": run,
                     "sample-id": "P%s_D%s" % (m.group(1), m.group(2)),
                     "subject": "P" + m.group(1),
                     "time": "baseline" if m.group(2) == "0"
                     else "D" + m.group(2)})
    out = pd.DataFrame(rows)
    out["status"] = "periodontitis"
    out["unit"] = out["subject"]
    return manifest(out, paired=True), out


def houston():
    """PRJNA726254: BioSample host_disease 'healthy periodontium',
    'chronic periodontitis' (before) and 'chronic periodontitis-treated'
    (after non-surgical therapy) of the same subjects (host_subject_id)."""
    a = attributes()
    disease = a["host disease"].str.lower()
    out = pd.DataFrame({
        "run": a["run"],
        "subject": a["host subject id"],
        "status": disease.map(lambda s: "health" if "healthy" in s
                              else "periodontitis"),
        "time": disease.map(lambda s: "after-therapy" if "treated" in s
                            else "baseline"),
        "host_disease": a["host disease"],
        "age": a.get("host age"), "smoker": a.get("smoker"),
    })
    out["sample-id"] = out["subject"] + "_" + out["time"]
    out["unit"] = out["subject"]
    return manifest(out, paired=False), out


def elderly_nspt():
    """PRJEB106816: elderly periodontitis patients, subgingival and
    supragingival plaque at baseline and follow-up 2, in a therapy
    (Intervention) and an untreated Control arm. The BioSample description
    reads 'Participant <n>, <site> Plaque, <time>, <arm>'. Only subgingival
    samples are used."""
    a = attributes()
    parts = a["description"].str.extract(
        r"Participant (\d+), (\w+) Plaque, ([^,]+), (\w+)")
    out = pd.DataFrame({"run": a["run"], "subject": "P" + parts[0],
                        "plaque": parts[1], "time_label": parts[2],
                        "arm": parts[3].str.lower()})
    out = out[out["plaque"].str.lower().str.startswith("subgingival")].copy()
    out["time"] = out["time_label"].map(
        lambda t: "baseline" if t == "Baseline" else "follow-up-2")
    out["status"] = "periodontitis"
    out["unit"] = out["subject"]
    out["sample-id"] = out["subject"] + "_" + out["time"]
    return manifest(out, paired=True), out


def thai_stages():
    """PRJNA1294584: BioSample title control / Stage2 / Stage3 / Stage4,
    one sample per subject (host_subject_id). The controls all have
    gingivitis (Sawangpanyangkura et al. 2026, J Clin Periodontol
    53:1332): there is no health group."""
    a = attributes()
    label = a["title"] if "title" in a else runs["sample_title"]
    stage = label.str.extract(r"Stage(\d)")[0]
    out = pd.DataFrame({"run": a["run"], "subject": a["host subject id"],
                        "group_label": label})
    out["status"] = ["gingivitis" if str(x).lower() == "control"
                     else "periodontitis" for x in label]
    out["stage"] = [int(st) if isinstance(st, str) else None
                    for st in stage]
    out["sample-id"] = out["subject"]
    return manifest(out, paired=True), out


def xian_stages():
    """PRJNA1470369: BioSample host_disease healthy / Stage I-IV
    periodontitis, one sample per subject."""
    a = attributes()
    disease = a["host disease"].str.lower()
    out = pd.DataFrame({"run": a["run"], "host_disease": a["host disease"]})
    out["status"] = disease.map(lambda s: "health" if s == "healthy"
                                else "periodontitis")
    out["stage"] = disease.map(lambda s: 0 if s == "healthy" else STAGE[s])
    lib = runs.set_index("run_accession")["library_name"]
    out["sample-id"] = lib.loc[out["run"]].to_numpy()
    out["subject"] = out["sample-id"]
    return manifest(out, paired=True), out


def peridysseq():
    """PRJEB83770: BioSample 'nose/mouth/teeth/throat disorder' Healthy,
    Gingivitis or Periodontitis; one subgingival sample per subject."""
    a = attributes()
    status = a["nose/mouth/teeth/throat disorder"].str.lower().map(
        {"healthy": "health", "gingivitis": "gingivitis",
         "periodontitis": "periodontitis"})
    alias = runs.set_index("run_accession")["sample_alias"]
    out = pd.DataFrame({"run": a["run"], "status": status,
                        "sample-id": alias.loc[a["run"]].to_numpy()})
    out["subject"] = out["sample-id"]
    return manifest(out, paired=True), out


def santiago(region):
    """PRJNA773195 (V1-V3) and PRJNA773202 (V3-V4): the same subgingival
    samples of systemically healthy subjects with or without periodontitis
    (BioSample Health_condition). Mock communities are dropped. The
    OSRG<n> code that ends every sample alias identifies the subject and
    links the two regions (Code_sample_CRG repeats across runs)."""
    a = attributes()
    a = a[a["health condition"].isin(["Healthy", "Periodontitis"])].copy()
    out = pd.DataFrame({
        "run": a["run"],
        "status": a["health condition"].map(
            {"Healthy": "health", "Periodontitis": "periodontitis"}),
        "severity": a.get("periodontitis severity"),
        "smoking": a.get("smoking"), "country": a.get("geo loc name"),
        "code_sample_crg": a.get("code sample crg"),
    })
    alias = runs.set_index("run_accession")["sample_alias"]
    out["sample-id"] = alias.loc[out["run"]].to_numpy()
    out["sample_code"] = out["sample-id"].str.extract(r"(OSRG\d+)$")[0]
    assert out["sample_code"].notna().all() and out["sample_code"].is_unique
    out["subject"] = out["sample_code"]
    out["region"] = region
    return manifest(out, paired=True), out


def bizzarro2016():
    """SRP062274: periodontitis patients treated with scaling and root
    planing with or without systemic amoxicillin + metronidazole, sampled
    at baseline, 3, 6 and 12 months; BioSample attributes give the time
    point, the antibiotic arm and the responder label."""
    a = attributes()
    alias = runs.set_index("run_accession")["sample_alias"]
    out = pd.DataFrame({"run": a["run"],
                        "sample-id": alias.loc[a["run"]].to_numpy(),
                        "timepoint": a["timepoint"],
                        "arm": a["antibiotics"].map(
                            {"0": "placebo", "1": "antibiotics"}),
                        "response": a["responder"].map(
                            {"0": "non-responder", "1": "responder"}),
                        "smoker": a.get("smoker")})
    out["subject"] = out["sample-id"].str.split(".").str[0]
    out["time"] = out["timepoint"].map(
        lambda t: "baseline" if t == "bsl" else t)
    out["status"] = "periodontitis"
    out["unit"] = out["subject"]
    return manifest(out, paired=False), out


COHORTS = {
    "SRP075100": chen2018,
    "PRJNA863881": spanish_stages,
    "PRJNA786436": byrne2022,
    "PRJNA736618": johnston2021,
    "PRJNA726254": houston,
    "PRJEB106816": elderly_nspt,
    "PRJNA1294584": thai_stages,
    "PRJNA1470369": xian_stages,
    "PRJEB83770": peridysseq,
    "PRJNA773202": lambda: santiago("V3-V4"),
    "PRJNA773195": lambda: santiago("V1-V3"),
    "SRP226726": pei2020,
    "SRP009299_V12": lambda: griffen2012("V12"),
    "SRP009299_V4": lambda: griffen2012("V4"),
    "SRP051201": kirst2015,
    "SRP062274": bizzarro2016,
}

if cohort not in COHORTS:
    raise SystemExit("unknown cohort " + cohort)
man, meta = COHORTS[cohort]()
meta = meta.reset_index(drop=True)
for col in STANDARD:
    if col not in meta.columns:
        meta[col] = ""
if (meta["time"] == "").all():
    meta["time"] = "baseline"
meta = meta[STANDARD + [c for c in meta.columns if c not in STANDARD]]

assert meta["sample-id"].is_unique, "duplicated sample IDs"
assert meta["status"].isin(["health", "gingivitis", "periodontitis",
                            "other"]).all(), "unmapped status"
missing = [p for col in man.columns[1:] for p in man[col]
           if not os.path.exists(p)]
assert not missing, "missing FASTQ, e.g. %s" % missing[:3]
out = os.path.join(data, cohort)
os.makedirs(out, exist_ok=True)
man.to_csv(os.path.join(out, "manifest.tsv"), sep="\t", index=False)
meta.to_csv(os.path.join(out, "metadata.tsv"), sep="\t", index=False)
print(meta.groupby(["status", "time"]).agg(
    samples=("sample-id", "size"), subjects=("subject", "nunique")))
