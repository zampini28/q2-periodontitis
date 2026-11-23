"""Taxonomy parsing of the plugin against smdi_microbiome.R on every label a
classifier trained on a reference database can output.

    python -I parser_check.py <ref_dir> <results_dir> --rscript <Rscript>

For each reference taxonomy (eHOMD, SILVA, GTDB, Greengenes2) the labels
are every distinct taxonomy string and every truncation of it at a rank
boundary (what naive Bayes and VSEARCH consensus report when they stop at a
higher rank), each with ';' and '; ' separators, plus 'Unassigned'. Both
parsers assign every label to a species and a genus unit; a database passes
when no label differs. For each index, the signature taxa that no label of
the database reaches are listed: such taxa can only enter the index as
zeros with that database.

Downloads the references into <ref_dir> when absent. Writes:

    <results_dir>/tables/parser_crosscheck.tsv
    <results_dir>/tables/parser_unreachable.tsv
    <results_dir>/parser_references.json   (URL and sha256 of each file)
"""
import argparse
import gzip
import hashlib
import io
import json
import os
import sys
import tempfile
import urllib.request
import zipfile

import biom
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rcheck  # noqa: E402
from q2_periodontitis._methods import (  # noqa: E402
    _canonical, feature_units, load_signature, load_synonyms, taxon_key)

# name, URL, format
REFERENCES = [
    ("eHOMD V15.23",
     "https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.23/"
     "HOMD_16S_rRNA_RefSeq_V15.23.qiime.taxonomy", "tsv"),
    ("SILVA 138 (99%)",
     "https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza", "qza"),
    ("GTDB R226 (bac120)",
     "https://data.gtdb.ecogenomic.org/releases/release226/226.0/"
     "bac120_taxonomy_r226.tsv.gz", "tsv.gz"),
    ("Greengenes2 2022.10 (backbone)",
     "http://ftp.microbio.me/greengenes_release/2022.10/"
     "2022.10.backbone.tax.qza", "qza"),
]


def fetch(url, ref_dir):
    """Path of the downloaded file (downloaded when absent) and its
    sha256."""
    path = os.path.join(ref_dir, url.rsplit("/", 1)[1])
    if not os.path.exists(path):
        os.makedirs(ref_dir, exist_ok=True)
        tmp = path + ".part"
        urllib.request.urlretrieve(url, tmp)
        os.replace(tmp, path)
    with open(path, "rb") as fh:
        return path, hashlib.sha256(fh.read()).hexdigest()


def read_taxa(path, fmt):
    """The taxonomy strings of a reference file."""
    if fmt == "qza":
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist()
                        if n.endswith("/data/taxonomy.tsv"))
            text = z.read(name).decode()
        df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str,
                         keep_default_na=False, quoting=3)
        return df.iloc[:, 1]
    opener = gzip.open if fmt == "tsv.gz" else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, dtype=str,
                         keep_default_na=False, quoting=3)
    return df.iloc[:, 1]


def labels(taxa):
    """Every distinct label and rank truncation, with both separators."""
    out = {"Unassigned"}
    for taxon in set(taxa):
        ranks = [r.strip() for r in taxon.split(";")]
        for k in range(1, len(ranks) + 1):
            out.add(";".join(ranks[:k]))
            out.add("; ".join(ranks[:k]))
    return sorted(out)


def signature_keys(path, synonyms):
    sig = load_signature(path)
    level = sig["level"].iloc[0].strip().lower()
    keys = [taxon_key(_canonical(t.strip(), synonyms)) if level == "species"
            else _canonical(t.strip(), synonyms).lower()
            for t in sig["taxon"]]
    return level, list(sig["taxon"]), keys


def check(name, label_list, rscript, synonyms):
    """Rows of parser_crosscheck.tsv and parser_unreachable.tsv for one
    database."""
    ids = ["L%d" % i for i in range(len(label_list))]
    taxonomy = pd.DataFrame({"Taxon": label_list}, index=pd.Index(ids))
    table = biom.Table(np.ones((len(ids), 1)), ids, ["S1"])
    rows, unreachable, resolved_units = [], [], {}
    with tempfile.TemporaryDirectory() as workdir:
        tax_path = os.path.join(workdir, "labels.tsv")
        rcheck.write_taxonomy(taxonomy, tax_path)
        for level in ("species", "genus"):
            out = os.path.join(workdir, "units-%s.tsv" % level)
            rcheck.run_child([rscript, rcheck.R_SCRIPT, "--units", tax_path,
                              rcheck.SYNONYMS, level, out])
            r = rcheck.read_units(out).reindex(ids)
            plugin = pd.Series(feature_units(table, taxonomy, synonyms,
                                             level), index=ids)
            bad = plugin != r
            resolved = ~plugin.str.startswith("unresolved: ")
            resolved_units[level] = set(plugin[resolved])
            rows.append({"database": name, "level": level,
                         "labels": len(ids),
                         "resolved": int(resolved.sum()),
                         "units": len(resolved_units[level]),
                         "mismatch": int(bad.sum()),
                         "pass": bool(not bad.any())})
            for i in plugin.index[bad][:10]:
                print("  MISMATCH %s | plugin: %s | R: %s"
                      % (taxonomy.loc[i, "Taxon"], plugin[i], r[i]))
    for index, path in rcheck.INDICES.items():
        level, taxa, keys = signature_keys(path, synonyms)
        missing = [t for t, k in zip(taxa, keys)
                   if k not in resolved_units[level]]
        unreachable += [{"database": name, "index": index, "taxon": t}
                        for t in missing]
        rows[-1 if level == "genus" else 0].setdefault(
            "signature_reached", {})[index] = "%d/%d" % (
                len(taxa) - len(missing), len(taxa))
    return rows, unreachable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ref_dir")
    parser.add_argument("results")
    parser.add_argument("--rscript", required=True)
    args = parser.parse_args()

    synonyms = load_synonyms()
    rows, unreachable, provenance = [], [], []
    for name, url, fmt in REFERENCES:
        path, sha = fetch(url, args.ref_dir)
        taxa = read_taxa(path, fmt)
        label_list = labels(taxa)
        print("%s: %d records, %d distinct labels, %d labels checked"
              % (name, len(taxa), taxa.nunique(), len(label_list)),
              flush=True)
        provenance.append({"database": name, "url": url, "sha256": sha,
                           "records": int(len(taxa)),
                           "distinct_labels": int(taxa.nunique())})
        db_rows, db_unreachable = check(name, label_list, args.rscript,
                                        synonyms)
        rows += db_rows
        unreachable += db_unreachable

    for row in rows:
        reached = row.pop("signature_reached", {})
        row["signature_reached"] = "; ".join(
            "%s %s" % kv for kv in reached.items())
    tables = os.path.join(args.results, "tables")
    os.makedirs(tables, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(os.path.join(tables, "parser_crosscheck.tsv"), sep="\t",
                  index=False)
    pd.DataFrame(unreachable, columns=["database", "index", "taxon"]).to_csv(
        os.path.join(tables, "parser_unreachable.tsv"), sep="\t",
        index=False)
    with open(os.path.join(args.results, "parser_references.json"),
              "w") as fh:
        json.dump(provenance, fh, indent=2)
    print(result.to_string())
    if not result["pass"].all():
        sys.exit("Parser cross-check FAILED")
    print("Parser cross-check passed for %d databases" % len(REFERENCES))


if __name__ == "__main__":
    main()
