import os
import re
import warnings

import numpy as np
import pandas as pd
from scipy import sparse
import biom
import qiime2

# Socransky et al. (1998) complexes plus the Actinomyces ("Blue") group of
# Socransky & Haffajee (2002, Fig. 1). "Other" holds every read that is not
# assigned to a complex.
COMPLEXES = [
    "Red",
    "Orange",
    "Yellow",
    "Green",
    "Purple",
    "Blue",
    "Other"
]

LEVELS = ["species", "genus"]
DIRECTIONS = ["dysbiotic", "normobiotic"]

# Columns of the indices table that describe data quality rather than
# dysbiosis; they are carried along but never differenced or tested.
QC_PREFIX = "qc_"

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
COMPLEX_REFERENCE = os.path.join(DATA_DIR, "socransky_complexes.tsv")
SYNONYMS = os.path.join(DATA_DIR, "synonyms_map.tsv")
PC2014_REFERENCE = os.path.join(DATA_DIR, "perezchaparro2014_moderate.tsv")
SMDI_SIMPLIFIED = os.path.join(DATA_DIR, "smdi_simplified_genus.tsv")
SMDI_SPECIES = os.path.join(DATA_DIR, "smdi_chen2022_species.tsv")


def load_complex_reference():
    """Built-in taxon -> complex reference (columns taxon, complex, level,
    source, evidence)."""
    return pd.read_csv(COMPLEX_REFERENCE, sep="\t", dtype=str)


def load_synonyms():
    """Historical name -> current name table, keyed in lower case."""
    df = pd.read_csv(SYNONYMS, sep="\t", dtype=str)
    return {old.strip().lower(): new.strip()
            for old, new in zip(df["old_name"], df["canonical_name"])}


def load_signature(path):
    """Reads a dysbiosis signature (taxon, direction, level) from a TSV."""
    return pd.read_csv(path, sep="\t", dtype=str)


RANK_PREFIX_RE = re.compile(r"^[a-z]__")

# Database suffixes appended to genus names or epithets: SILVA
# (Prevotella_7, Treponema_2) and GTDB / Greengenes2 (Fusobacterium_A,
# Prevotella_D, nucleatum_D).
SUFFIX_TOKEN_RE = re.compile(r"\d+|[A-Z]{1,2}")
# Tokens after which a number is an identifier, not a suffix ("oral taxon 361")
ID_PREFIX_TOKENS = {"taxon", "hmt", "hot", "ot", "str.", "strain"}
# Human Microbial Taxon numbers: "HMT 360", "HMT-360", "HOT 360",
# "oral taxon 360"
HMT_RE = re.compile(r"(?:\bHMT|\bHOT|oral taxon)[\s_-]*(\d{1,4})\b",
                    re.IGNORECASE)

CLADE_RE = re.compile(r"[\s_]+clade[\s_]+\d+", re.IGNORECASE)

UNINFORMATIVE_PREFIXES = (
    "uncultured", "unidentified", "unclassified", "unknown", "unassigned",
    "ambiguous", "metagenome", "environmental", "incertae"
)


def _is_uninformative(name):
    """True for placeholder labels such as 'uncultured_bacterium',
    'metagenome', 'unidentified' or 'Ambiguous_taxa'."""
    lower = name.strip().lower()
    return (not lower
            or lower.startswith(UNINFORMATIVE_PREFIXES)
            or "metagenome" in lower
            or lower in {"bacterium", "sp", "sp."})


def _clean_name(raw):
    """
    Normalizes a genus or species label: splits on spaces/underscores,
    drops database suffix tokens (Prevotella_7 -> Prevotella,
    Fusobacterium_A nucleatum_D -> Fusobacterium nucleatum) and the brackets
    that mark misclassified genera ([Eubacterium] nodatum).
    """
    # eHOMD 16S clades inside a species ("infantis_clade_431") are not
    # taxa of their own: drop them.
    raw = CLADE_RE.sub("", raw)
    tokens = [t for t in re.split(r"[\s_]+", raw.strip()) if t]
    kept = []
    for tok in tokens:
        if (kept and SUFFIX_TOKEN_RE.fullmatch(tok)
                and kept[-1].lower() not in ID_PREFIX_TOKENS):
            continue
        kept.append(tok)
    name = " ".join(kept)
    # "[Eubacterium] nodatum" -> "Eubacterium nodatum"; "[G-1]" placeholders
    # inside genus labels are kept, since they are part of the eHOMD name.
    return re.sub(r"^\[([A-Za-z]+)\]", r"\1", name)


def hmt_number(name):
    """Human Microbial Taxon number in a label, as an int, or None."""
    if not name:
        return None
    match = HMT_RE.search(name)
    return int(match.group(1)) if match else None


def taxon_key(name):
    """Matching key of a taxon name: "hmt N" for Human Microbial Taxon
    phylotypes (so "Fretibacterium sp. oral taxon 360" and "Fretibacterium
    sp. HMT 360" match, whatever placeholder genus a database version
    uses), the lower-case name otherwise."""
    number = hmt_number(name)
    return "hmt %d" % number if number is not None else name.lower()


def _canonical(name, synonyms):
    if name is None:
        return None
    return synonyms.get(name.lower(), name)


def parse_taxonomy_string(tax_str, synonyms=None):
    """
    Parses QIIME 2 / eHOMD / Greengenes / Greengenes2 / GTDB / SILVA taxonomy
    strings into (genus, species). Database suffixes (SILVA "_7", GTDB "_A")
    are removed, placeholder names ("uncultured", "metagenome", ...) ignored
    and historical names resolved through `synonyms`.

    Multi-species labels ("Streptococcus mitis/oralis") are returned as a
    "/"-separated list of binomials; see `species_candidates`.
    """
    if synonyms is None:
        synonyms = {}

    if not isinstance(tax_str, str) or not tax_str.strip():
        return None, None

    ranks = [r.strip() for r in tax_str.split(";")]
    has_prefixes = any(RANK_PREFIX_RE.match(r) for r in ranks)
    genus_raw = None
    species_raw = None

    for pos, r in enumerate(ranks):
        clean_r = RANK_PREFIX_RE.sub("", r).strip()
        if _is_uninformative(clean_r):
            continue
        if has_prefixes:
            if r.startswith("g__"):
                genus_raw = clean_r
            elif r.startswith("s__"):
                species_raw = clean_r
        elif pos == 5:
            genus_raw = clean_r
        elif pos == 6:
            species_raw = clean_r

    genus = None
    if genus_raw:
        tokens = _clean_name(genus_raw).split(" ")
        # Keep eHOMD placeholder genera whole ("Peptostreptococcaceae [G-6]")
        keep = 2 if len(tokens) > 1 and tokens[1].startswith("[") else 1
        genus = _canonical(" ".join(tokens[:keep]), synonyms)

    species_full = None
    if species_raw:
        species_clean = _clean_name(species_raw)
        if species_clean[:1].isupper():
            # Binomial already present (s__Porphyromonas_gingivalis)
            species_full = species_clean
        elif genus:
            # Epithet only (eHOMD, Greengenes 13_8: g__Porphyromonas;
            # s__gingivalis)
            species_full = f"{genus} {species_clean}"
        if species_full and len(species_full.split(" ")) < 2:
            species_full = None

    if species_full and "/" in species_full:
        species_full = "/".join(
            _canonical(s, synonyms) for s in species_candidates(species_full))
    else:
        species_full = _canonical(species_full, synonyms)

    return genus, species_full


def species_candidates(species):
    """Expands "Streptococcus mitis/oralis" into
    ["Streptococcus mitis", "Streptococcus oralis"]; a single name is
    returned as a one-element list."""
    if not species:
        return []
    if "/" not in species:
        return [species]
    tokens = species.split(" ")
    genus = tokens[0]
    rest = " ".join(tokens[1:])
    candidates = []
    for part in rest.split("/"):
        part = part.strip()
        if not part:
            continue
        candidates.append(part if part[:1].isupper() else f"{genus} {part}")
    return candidates


def build_reference(reference_df, synonyms):
    """
    Validates a taxon -> complex table and returns two case-insensitive
    lookups: species (and subspecies) names, and genus names. Rows without a
    `level` column are species when the taxon has two or more words, genus
    otherwise. Raises ValueError on unknown complexes or levels.
    """
    df = reference_df.copy()
    missing = [c for c in ("taxon", "complex") if c not in df.columns]
    if missing:
        raise ValueError(
            "The reference must have the columns 'taxon' and 'complex'; "
            "missing: %s." % ", ".join(missing))
    df["taxon"] = df["taxon"].astype(str).str.strip()
    df["complex"] = df["complex"].astype(str).str.strip()
    if "level" not in df.columns or df["level"].isna().all():
        df["level"] = np.where(df["taxon"].str.contains(" "),
                               "species", "genus")
    df["level"] = df["level"].fillna("").astype(str).str.strip()
    blank = df["level"] == ""
    df.loc[blank, "level"] = np.where(
        df.loc[blank, "taxon"].str.contains(" "), "species", "genus")

    invalid = sorted(set(df["complex"]) - set(COMPLEXES))
    if invalid:
        offenders = sorted(df.loc[df["complex"].isin(invalid), "taxon"])
        raise ValueError(
            "Invalid complex name(s) %s (taxa: %s). Valid complexes are: %s."
            % (", ".join(repr(c) for c in invalid),
               ", ".join(offenders[:10])
               + (" ..." if len(offenders) > 10 else ""),
               ", ".join(COMPLEXES)))
    bad_levels = sorted(set(df["level"]) - set(LEVELS))
    if bad_levels:
        raise ValueError("Invalid level(s) %s; use 'species' or 'genus'."
                         % ", ".join(repr(v) for v in bad_levels))

    species, genera = {}, {}
    for taxon, complex_name, level in zip(df["taxon"], df["complex"],
                                          df["level"]):
        name = _canonical(taxon, synonyms)
        target = species if level == "species" else genera
        target[name.lower()] = complex_name
    return species, genera


def merge_references(base_df, custom_df):
    """Custom entries are added to, and override, the built-in ones."""
    if custom_df is None or custom_df.empty:
        return base_df
    missing = [c for c in ("taxon", "complex") if c not in custom_df.columns]
    if missing:
        raise ValueError(
            "custom_mapping must have the columns 'taxon' and 'complex'; "
            "missing: %s." % ", ".join(missing))
    custom = custom_df.copy()
    custom["taxon"] = custom["taxon"].astype(str).str.strip()
    keep = ~base_df["taxon"].str.lower().isin(custom["taxon"].str.lower())
    return pd.concat([base_df[keep], custom], ignore_index=True, sort=False)


def assign_feature_to_complex(tax_str, species_lookup, genus_lookup,
                              synonyms):
    """
    Assigns one taxonomic string to a complex. Returns (complex, resolution)
    with resolution one of:

    - "species": the species (or subspecies) is listed in the reference;
    - "genus": the species is unknown or unlisted, and the genus itself is
      listed as a whole (only genera with an explicit genus-level entry);
    - "ambiguous": a multi-species label whose candidates fall in different
      complexes (or only some of them are listed);
    - "unassigned": nothing matched.
    """
    genus, species = parse_taxonomy_string(tax_str, synonyms)

    candidates = species_candidates(species)
    if len(candidates) > 1:
        hits = {species_lookup.get(c.lower()) for c in candidates}
        if len(hits) == 1 and None not in hits:
            return hits.pop(), "species"
        if hits != {None}:
            return "Other", "ambiguous"
    elif species and species.lower() in species_lookup:
        return species_lookup[species.lower()], "species"

    if genus and genus.lower() in genus_lookup:
        return genus_lookup[genus.lower()], "genus"

    return "Other", "unassigned"


def _signature_lookup(path=PC2014_REFERENCE, synonyms=None):
    """Names (lower case) and HMT numbers of the Pérez-Chaparro et al. (2014)
    moderate-evidence taxa."""
    synonyms = synonyms or {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    names = {_canonical(t.strip(), synonyms).lower() for t in df["taxon"]}
    hmts = {int(h) for h in df["hmt"].dropna()}
    return names, hmts


def is_pc2014_taxon(species, names, hmts):
    """True when the parsed species (a single, non-ambiguous name) is one of
    the moderate-evidence taxa, matched by name or by HMT number."""
    if not species or "/" in species:
        return False
    if species.lower() in names:
        return True
    number = hmt_number(species)
    return number is not None and number in hmts


def _confidence_series(taxonomy):
    if "Confidence" not in taxonomy.columns:
        return None
    # Missing / non-numeric confidences carry no information and are never
    # filtered out.
    return pd.to_numeric(taxonomy["Confidence"], errors="coerce")


def _taxon_series(taxonomy):
    if isinstance(taxonomy, pd.Series):
        return taxonomy
    if "Taxon" in taxonomy.columns:
        return taxonomy["Taxon"]
    return taxonomy.iloc[:, 0]


def validate_signature(signature):
    """Checks a dysbiosis signature table and returns it normalized."""
    missing = [c for c in ("taxon", "direction", "level")
               if c not in signature.columns]
    if missing:
        raise ValueError("The signature must have the columns 'taxon', "
                         "'direction' and 'level'; missing: %s."
                         % ", ".join(missing))
    sig = signature.copy()
    for col in ("taxon", "direction", "level"):
        sig[col] = sig[col].astype(str).str.strip()
    sig["direction"] = sig["direction"].str.lower()
    sig["level"] = sig["level"].str.lower()
    bad = sorted(set(sig["direction"]) - set(DIRECTIONS))
    if bad:
        raise ValueError("Invalid direction(s) %s; use 'dysbiotic' or "
                         "'normobiotic'." % ", ".join(repr(b) for b in bad))
    levels = set(sig["level"])
    if len(levels) != 1 or not levels <= set(LEVELS):
        raise ValueError("All signature rows must share one level, "
                         "'species' or 'genus' (found: %s)."
                         % ", ".join(sorted(levels)))
    for direction in DIRECTIONS:
        if not (sig["direction"] == direction).any():
            raise ValueError("The signature needs at least one %s taxon."
                             % direction)
    if sig["taxon"].map(taxon_key).duplicated().any():
        raise ValueError("The signature lists a taxon more than once.")
    return sig


def aggregate_counts(table, taxonomy, synonyms, level, exclude=None):
    """
    Sums the feature counts by species or genus name (samples x taxa).
    Features without a name at that level, and features in `exclude` (e.g.
    low-confidence ones), are kept as their own units, keyed by their full
    taxonomy string, so that the total of every sample is unchanged.
    """
    exclude = exclude or set()
    taxa = _taxon_series(taxonomy)
    feature_ids = list(table.ids(axis="observation"))
    units = []
    for f_id in feature_ids:
        tax_str = taxa.get(f_id, "")
        if f_id in exclude:
            units.append("low-confidence: %s" % f_id)
            continue
        genus, species = parse_taxonomy_string(tax_str, synonyms)
        name = species if level == "species" else genus
        if name and "/" not in name:
            units.append(taxon_key(name) if level == "species"
                         else name.lower())
        else:
            units.append("unresolved: %s" % tax_str)
    codes, uniques = pd.factorize(pd.Series(units))
    indicator = sparse.csr_matrix(
        (np.ones(len(codes)), (codes, np.arange(len(codes)))),
        shape=(len(uniques), len(codes)))
    summed = indicator @ table.matrix_data  # taxa x samples, still sparse
    return pd.DataFrame(summed.toarray().T,
                        index=list(table.ids(axis="sample")),
                        columns=list(uniques))


def clr_microbiome(counts, pseudocount=None):
    """
    Centered log-ratio transform as implemented by the R package microbiome
    (v1.6.0, transform(x, "clr")), the implementation used to define the
    SMDI: counts are closed to relative abundances, a single pseudocount
    equal to half the smallest non-zero relative abundance of the whole
    table is added to every entry, and each sample is centred on the mean
    of its natural logs.

    With `pseudocount`, that fixed value is added to every relative
    abundance instead, so that a sample's CLR no longer depends on the
    other samples of the table.
    """
    depth = counts.sum(axis=1)
    rel = counts.div(depth.where(depth > 0), axis=0)
    values = rel.to_numpy(dtype=float)
    if pseudocount is not None:
        if not pseudocount > 0:
            raise ValueError("The pseudocount must be greater than 0.")
        values = values + pseudocount
    else:
        finite = values[np.isfinite(values)]
        if (finite == 0).any():
            values = values + finite[finite > 0].min() / 2.0
    logs = np.log(values)
    clr = logs - logs.mean(axis=1, keepdims=True)
    return pd.DataFrame(clr, index=counts.index, columns=counts.columns)


def compute_smdi(table, taxonomy, signature, synonyms, exclude=None,
                 pseudocount=None):
    """
    Subgingival Microbial Dysbiosis Index (Chen, Marsh & Al-Hebshi 2022):
    the mean CLR abundance of the dysbiotic taxa minus the mean CLR abundance
    of the normobiotic taxa. The CLR is taken over every taxon of the table
    at the signature's level. Signature taxa absent from the whole table are
    added as all-zero columns (that is, they receive the pseudocount), as
    they would be in a table merged with samples that contain them; a
    warning lists them. Features in `exclude` never count as signature taxa,
    but their reads stay in the table. `pseudocount` is passed to
    clr_microbiome (None: the original, table-dependent pseudocount).

    Returns (smdi, missing_taxa). Samples without reads get NaN.
    """
    sig = validate_signature(signature)
    level = sig["level"].iloc[0]
    counts = aggregate_counts(table, taxonomy, synonyms, level, exclude)

    if level == "species":
        names = [taxon_key(_canonical(t, synonyms)) for t in sig["taxon"]]
    else:
        names = [_canonical(t, synonyms).lower() for t in sig["taxon"]]
    missing = [n for n in names if n not in counts.columns]
    for name in missing:
        counts[name] = 0.0

    clr = clr_microbiome(counts, pseudocount)
    dysbiotic = [n for n, d in zip(names, sig["direction"])
                 if d == "dysbiotic"]
    normobiotic = [n for n, d in zip(names, sig["direction"])
                   if d == "normobiotic"]
    smdi = clr[dysbiotic].mean(axis=1) - clr[normobiotic].mean(axis=1)
    smdi[counts.sum(axis=1) == 0] = np.nan
    return smdi, missing


def collapse_complexes(
    table: biom.Table,
    taxonomy: pd.DataFrame,
    custom_mapping: pd.DataFrame = None,
    smdi_signature: pd.DataFrame = None,
    min_confidence: float = 0.0,
    smdi_pseudocount: float = None
) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    """
    Collapses a feature table into Socransky complexes and computes the SMDI
    (see compute_smdi), PC2014_moderate and the quality-control columns.
    """
    if not 0.0 <= min_confidence <= 1.0:
        raise ValueError("min_confidence must be within [0, 1].")
    if smdi_pseudocount is not None and not 0.0 < smdi_pseudocount < 1.0:
        raise ValueError("smdi_pseudocount must be within (0, 1).")

    synonyms = load_synonyms()
    reference = merge_references(load_complex_reference(), custom_mapping)
    species_lookup, genus_lookup = build_reference(reference, synonyms)
    pc_names, pc_hmts = _signature_lookup(synonyms=synonyms)

    taxa = _taxon_series(taxonomy)
    confidence = _confidence_series(taxonomy)

    sample_ids = list(table.ids(axis="sample"))
    feature_ids = list(table.ids(axis="observation"))
    missing_tax = [f for f in feature_ids if f not in taxa.index]
    if missing_tax:
        raise ValueError(
            "%d feature(s) of the table have no taxonomy (e.g. %s)."
            % (len(missing_tax), ", ".join(missing_tax[:5])))

    # Each feature belongs to exactly one complex, so that every sample's
    # complex profile sums to 1. Features whose classifier confidence is
    # below min_confidence go to "Other".
    feature_assignments = []
    w_matrix = np.zeros((len(feature_ids), len(COMPLEXES)))
    pc_vector = np.zeros(len(feature_ids))
    species_vector = np.zeros(len(feature_ids))
    low_conf_features = set()

    for idx, f_id in enumerate(feature_ids):
        tax_str = str(taxa[f_id])
        conf = np.nan
        if confidence is not None:
            conf = float(confidence.get(f_id, np.nan))
        genus, species = parse_taxonomy_string(tax_str, synonyms)

        low_conf = not np.isnan(conf) and conf < min_confidence
        if low_conf:
            low_conf_features.add(f_id)
            complex_name, resolution = "Other", "low-confidence"
        else:
            complex_name, resolution = assign_feature_to_complex(
                tax_str, species_lookup, genus_lookup, synonyms)
        w_matrix[idx, COMPLEXES.index(complex_name)] = 1.0

        pc_member = (not low_conf
                     and is_pc2014_taxon(species, pc_names, pc_hmts))
        pc_vector[idx] = float(pc_member)
        species_vector[idx] = float(
            bool(species) and "/" not in species and not low_conf)

        feature_assignments.append({
            "Feature ID": f_id,
            "Taxonomy": tax_str,
            "Genus": genus or "",
            "Species": species or "",
            "Complex": complex_name,
            "Resolution": resolution,
            "PC2014": "yes" if pc_member else "no",
            "Confidence": conf
        })

    counts = table.matrix_data  # sparse (M x N)
    depth = np.asarray(counts.sum(axis=0)).ravel()
    safe_depth = np.where(depth > 0, depth, 1.0)
    df_complexes = pd.DataFrame(
        np.asarray(counts.T @ w_matrix) / safe_depth[:, None],
        index=sample_ids, columns=COMPLEXES)

    pc_fraction = np.asarray(counts.T @ pc_vector).ravel() / safe_depth
    species_fraction = (np.asarray(counts.T @ species_vector).ravel()
                        / safe_depth)

    df_indices = pd.DataFrame(index=sample_ids)
    signatures = [("SMDI", load_signature(SMDI_SPECIES)),
                  ("SMDI_simplified", load_signature(SMDI_SIMPLIFIED))]
    if smdi_signature is not None and not smdi_signature.empty:
        signatures.append(("SMDI_custom", smdi_signature))
    for column, signature in signatures:
        smdi, missing = compute_smdi(table, taxonomy, signature, synonyms,
                                     low_conf_features, smdi_pseudocount)
        if missing:
            warnings.warn(
                "%s: %d signature taxa are absent from the table and enter "
                "the index as zeros: %s"
                % (column, len(missing), ", ".join(missing)), UserWarning)
        df_indices[column] = smdi.reindex(sample_ids)

    df_indices["PC2014_moderate"] = pc_fraction
    empty = depth == 0
    df_indices.loc[empty, "PC2014_moderate"] = np.nan
    df_indices[QC_PREFIX + "depth"] = depth
    df_indices[QC_PREFIX + "assigned_fraction"] = np.where(
        empty, np.nan, 1.0 - df_complexes["Other"].to_numpy())
    df_indices[QC_PREFIX + "species_fraction"] = np.where(
        empty, np.nan, species_fraction)

    mapping_columns = [
        "Feature ID", "Taxonomy", "Genus", "Species",
        "Complex", "Resolution", "PC2014", "Confidence"
    ]
    df_mapping = pd.DataFrame(
        feature_assignments, columns=mapping_columns
    ).set_index("Feature ID")
    return df_complexes, df_indices, df_mapping


def _normalize_label(value):
    """String form of a metadata value; integral floats lose their '.0'
    so numeric time points (e.g. 0.0) match labels such as '0'."""
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    return str(value)


def index_columns(indices):
    """Dysbiosis index columns, without the quality-control ones."""
    return [c for c in indices.columns if not str(c).startswith(QC_PREFIX)]


def delta_indices(
    indices: pd.DataFrame,
    metadata: qiime2.Metadata,
    subject_col: str,
    group_col: str,
    baseline_group: str,
    target_group: str
) -> pd.DataFrame:
    """
    Computes paired longitudinal deltas (target minus baseline) of every
    dysbiosis index in `indices` (quality-control columns are skipped), one
    row per subject. When a subject has several samples in a group (e.g.
    several sites), their mean is used. The output index is the subject ID
    taken from `subject_col`.
    """
    metadata_df = metadata.to_dataframe()
    missing = [c for c in (subject_col, group_col)
               if c not in metadata_df.columns]
    if missing:
        raise ValueError("Metadata column(s) not found: %s."
                         % ", ".join(missing))

    # Join only the two needed columns, under private names, to avoid any
    # collision with index or metadata column names.
    design = pd.DataFrame({
        "__subject__": metadata_df[subject_col],
        "__group__": metadata_df[group_col],
    }).dropna()
    design["__subject__"] = design["__subject__"].map(_normalize_label)
    design["__group__"] = design["__group__"].map(_normalize_label)

    baseline_label = str(baseline_group)
    target_label = str(target_group)
    if metadata.columns[group_col].type == "numeric":
        # Let "1" match a numeric column value of 1.0
        try:
            baseline_label = _normalize_label(float(baseline_group))
            target_label = _normalize_label(float(target_group))
        except ValueError:
            raise ValueError(
                "Column '%s' is numeric, so baseline_group and target_group "
                "must be numbers (got %r and %r)."
                % (group_col, baseline_group, target_group)) from None

    metrics = index_columns(indices)
    merged = indices[metrics].join(design, how="inner")
    merged = merged[merged["__group__"].isin([baseline_label, target_label])]

    # Several sites per subject and visit are averaged, so each subject
    # contributes a single observation per visit.
    per_visit = merged.groupby(["__subject__", "__group__"])[metrics].mean()
    groups = per_visit.index.get_level_values("__group__")
    baseline = per_visit[groups == baseline_label].droplevel("__group__")
    target = per_visit[groups == target_label].droplevel("__group__")
    paired = baseline.index.intersection(target.index)
    if len(paired) == 0:
        raise ValueError(
            "No subject in '%s' has samples in both groups '%s' and '%s' of "
            "column '%s'." % (subject_col, baseline_group, target_group,
                              group_col))

    df_res = target.loc[paired] - baseline.loc[paired]
    df_res = df_res.sort_index()
    df_res.index.name = "subject-id"
    return df_res
