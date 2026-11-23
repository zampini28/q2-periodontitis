"""Input tables for the SMDI web calculator of Chen, Marsh & Al-Hebshi
(https://bioinformatics.forsyth.org/smdi/, Calculator V1.2), and the values
q2-periodontitis gives for them.

    python -I make_tables.py <out_dir>

The calculator is run by hand (it asks for a typed verification code); see
README.md. Each table has the 49 SMDI species, named with the labels of the
paper's Appendix Data Set 1 (`original_label` in smdi_chen2022_species.tsv,
the naming of the calculator's own sample input), plus 30 other oral species,
in 8 samples (fixed seed):

    dense.txt    every species in every sample: no zero, so no pseudocount
                 under either convention
    zeros.txt    dysbiotic species at 0-2 reads, the others at 50-499 with
                 about 30% zeros, every SMDI species present in at least
                 one sample: the two pseudocount conventions differ by
                 0.1-0.25 per sample
    absent.txt   as zeros.txt without 5 SMDI species (run with zero-fill)

expected.tsv gives, per table and sample, the plugin's SMDI (microbiome
1.6.0 behaviour: when the table has a zero, half the smallest non-zero
relative abundance is added to every entry) and the SMDI with that value
added to zero entries only (the convention the paper's text and the
calculator's guide describe).
"""
import os
import sys
import warnings

import biom
import numpy as np
import pandas as pd

from q2_periodontitis._methods import (
    SMDI_SPECIES, compute_smdi, load_signature, load_synonyms)

SEED = 20261123
SAMPLES = ["S%d" % i for i in range(1, 9)]
ABSENT = ["Treponema medium", "Rothia aeria", "Fretibacterium sp. HMT 359",
          "Mycoplasma faucium", "Prevotella loescheii"]
# Oral species that are not SMDI taxa, named as in the calculator's sample
# input (smdi/demo/upload_s.txt)
OTHERS = [
    "Prevotella denticola", "Prevotella intermedia", "Prevotella nigrescens",
    "Prevotella melaninogenica", "Prevotella oris", "Campylobacter rectus",
    "Campylobacter showae", "Campylobacter gracilis", "Parvimonas micra",
    "Eikenella corrodens", "Capnocytophaga gingivalis",
    "Capnocytophaga ochracea", "Capnocytophaga sputigena",
    "Streptococcus gordonii", "Streptococcus oralis", "Streptococcus mutans",
    "Leptotrichia buccalis", "Selenomonas sputigena", "Dialister invisus",
    "Haemophilus parainfluenzae", "Neisseria flavescens",
    "Lautropia mirabilis", "Kingella oralis", "Actinomyces israelii",
    "Actinomyces gerencseriae", "Veillonella dispar",
    "Peptostreptococcus stomatis", "Porphyromonas catoniae",
    "Tannerella sp._oral_taxon_286", "Leptotrichia hofstadii"]


def taxonomy_string(name):
    return "k__Bacteria; g__%s; s__%s" % (name.split(" ")[0], name)


def smdi_values(counts, plugin_names, signature, synonyms):
    """Plugin SMDI (Series) and the zeros-only variant (ndarray)."""
    ids = ["F%d" % i for i in range(len(plugin_names))]
    table = biom.Table(counts, ids, SAMPLES)
    taxonomy = pd.DataFrame({"Taxon": [taxonomy_string(n)
                                       for n in plugin_names]},
                            index=pd.Index(ids, name="Feature ID"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        smdi, missing = compute_smdi(table, taxonomy, signature, synonyms)
    # zeros-only variant on the same taxon-level table (absent taxa as zeros)
    rel = counts / counts.sum(axis=0)
    eps = rel[rel > 0].min() / 2
    keys = list(plugin_names)
    dys = [keys.index(t) for t, d in zip(signature["taxon"],
                                         signature["direction"])
           if d == "dysbiotic" and t in keys]
    nor = [keys.index(t) for t, d in zip(signature["taxon"],
                                         signature["direction"])
           if d == "normobiotic" and t in keys]
    n_absent_dys = sum(d == "dysbiotic" and t not in keys for t, d in
                       zip(signature["taxon"], signature["direction"]))
    n_absent_nor = sum(d == "normobiotic" and t not in keys for t, d in
                       zip(signature["taxon"], signature["direction"]))
    full = np.vstack([np.where(rel == 0, eps, rel),
                      np.full((n_absent_dys + n_absent_nor, rel.shape[1]),
                              eps)])
    logs = np.log(full)
    clr = logs - logs.mean(axis=0)
    dys += list(range(len(keys), len(keys) + n_absent_dys))
    nor += list(range(len(keys) + n_absent_dys,
                      len(keys) + n_absent_dys + n_absent_nor))
    zeros_only = clr[dys].mean(axis=0) - clr[nor].mean(axis=0)
    return smdi.reindex(SAMPLES), zeros_only, missing


def write_calculator_table(path, labels, counts):
    with open(path, "w") as fh:
        fh.write("#SampleID\t" + "\t".join(SAMPLES) + "\n")
        for label, row in zip(labels, counts):
            fh.write(label + "\t" + "\t".join(str(int(v)) for v in row)
                     + "\n")


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    signature = load_signature(SMDI_SPECIES)
    synonyms = load_synonyms()
    plugin_names = list(signature["taxon"]) + OTHERS
    labels = list(signature["original_label"]) + OTHERS
    rng = np.random.default_rng(SEED)
    dense = rng.integers(1, 2000, (len(labels), len(SAMPLES)))
    # Dysbiotic species at 0-2 reads, all others abundant with ~30% zeros:
    # adding the pseudocount to every entry then raises the rare dysbiotic
    # counts markedly, so the two conventions are easy to tell apart.
    dysbiotic = np.array([d == "dysbiotic" for d in signature["direction"]]
                         + [False] * len(OTHERS))
    zeros = rng.integers(50, 500, dense.shape) * (
        rng.random(dense.shape) > 0.3)
    zeros[dysbiotic] = rng.integers(0, 3, (dysbiotic.sum(), len(SAMPLES)))
    # every SMDI species present in at least one sample
    for i in range(len(signature)):
        if not zeros[i].any():
            zeros[i, rng.integers(len(SAMPLES))] = 1
    keep = [i for i, n in enumerate(plugin_names) if n not in ABSENT]
    tables = {"dense": (plugin_names, labels, dense),
              "zeros": (plugin_names, labels, zeros),
              "absent": ([plugin_names[i] for i in keep],
                         [labels[i] for i in keep], zeros[keep])}
    rows = []
    for name, (names, labs, counts) in tables.items():
        write_calculator_table(os.path.join(out, name + ".txt"), labs, counts)
        plugin, zeros_only, missing = smdi_values(counts, names, signature,
                                                  synonyms)
        has_zero = bool((counts == 0).any()) or bool(missing)
        for j, sample in enumerate(SAMPLES):
            rows.append({"table": name, "sample": sample,
                         "has_zero": has_zero,
                         "absent_smdi_taxa": len(missing),
                         "plugin_microbiome_1.6.0": plugin[sample],
                         "zeros_only_pseudocount": zeros_only[j],
                         "calculator": ""})
    pd.DataFrame(rows).to_csv(os.path.join(out, "expected.tsv"), sep="\t",
                              index=False, float_format="%.6f")


if __name__ == "__main__":
    main()
