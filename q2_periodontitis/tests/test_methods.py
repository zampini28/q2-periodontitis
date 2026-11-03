import math
import unittest
import warnings

import numpy as np
import pandas as pd
import biom
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_periodontitis._methods import (
    COMPLEXES,
    load_complex_reference,
    load_synonyms,
    load_signature,
    PC2014_REFERENCE,
    SMDI_SIMPLIFIED,
    SMDI_SPECIES,
    taxon_key,
    parse_taxonomy_string,
    species_candidates,
    hmt_number,
    build_reference,
    assign_feature_to_complex,
    aggregate_counts,
    clr_microbiome,
    compute_smdi,
    collapse_complexes,
    delta_indices
)


def _taxonomy(taxa, confidence=None):
    df = pd.DataFrame({"Taxon": list(taxa.values())},
                      index=pd.Index(list(taxa.keys()), name="Feature ID"))
    if confidence is not None:
        df["Confidence"] = [confidence[k] for k in taxa]
    return df


class TestParser(TestPluginBase):
    package = "q2_periodontitis.tests"

    def setUp(self):
        super().setUp()
        self.synonyms = load_synonyms()

    def test_silva_binomial(self):
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; "
                "f__Porphyromonadaceae; g__Porphyromonas; s__Porphyromonas_gingivalis"),
            ("Porphyromonas", "Porphyromonas gingivalis"))

    def test_silva_numeric_genus_suffix(self):
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; "
                "f__Prevotellaceae; g__Prevotella_7; s__Prevotella_melaninogenica"),
            ("Prevotella", "Prevotella melaninogenica"))
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Spirochaetota; c__Spirochaetia; o__Spirochaetales; "
                "f__Spirochaetaceae; g__Treponema_2; s__uncultured_bacterium"),
            ("Treponema", None))

    def test_gtdb_letter_suffixes(self):
        # Genus and epithet suffixes (GTDB / Greengenes2)
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Fusobacteriota; c__Fusobacteriia; o__Fusobacteriales; "
                "f__Fusobacteriaceae; g__Fusobacterium_A; s__Fusobacterium_A nucleatum_D"),
            ("Fusobacterium", "Fusobacterium nucleatum"))
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; "
                "f__Bacteroidaceae; g__Prevotella_D; s__Prevotella_D intermedia"),
            ("Prevotella", "Prevotella intermedia"))

    def test_greengenes_epithet_only(self):
        self.assertEqual(
            parse_taxonomy_string(
                "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; "
                "f__Porphyromonadaceae; g__Porphyromonas; s__gingivalis"),
            ("Porphyromonas", "Porphyromonas gingivalis"))

    def test_no_rank_prefixes(self):
        self.assertEqual(
            parse_taxonomy_string(
                "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;"
                "Porphyromonadaceae;Porphyromonas;gingivalis"),
            ("Porphyromonas", "Porphyromonas gingivalis"))

    def test_oral_taxon_number_is_not_a_suffix(self):
        _, species = parse_taxonomy_string(
            "d__Bacteria; g__Fretibacterium; s__Fretibacterium_human_oral_taxon_361")
        self.assertEqual(species, "Fretibacterium human oral taxon 361")

    def test_uninformative_labels(self):
        for label in ["uncultured_bacterium", "metagenome", "unidentified",
                      "human_oral_metagenome", "Ambiguous_taxa", "uncultured"]:
            genus, species = parse_taxonomy_string(
                f"d__Bacteria; g__Porphyromonas; s__{label}")
            self.assertEqual((genus, species), ("Porphyromonas", None), label)
        self.assertEqual(
            parse_taxonomy_string("d__Bacteria; g__uncultured; s__metagenome"),
            (None, None))
        self.assertEqual(parse_taxonomy_string("Unassigned"), (None, None))
        self.assertEqual(parse_taxonomy_string(""), (None, None))

    def test_synonyms(self):
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; g__Tannerella; s__Bacteroides_forsythus", self.synonyms),
            ("Tannerella", "Tannerella forsythia"))
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; g__Wolinella; s__Wolinella_recta", self.synonyms)[1],
            "Campylobacter rectus")
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; g__Peptostreptococcus; s__Peptostreptococcus_micros",
                self.synonyms)[1],
            "Parvimonas micra")


    def test_ehomd_bracketed_genus_and_placeholder_genera(self):
        self.assertEqual(
            parse_taxonomy_string(
                "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;"
                "f__Peptostreptococcaceae;g__Peptostreptococcaceae_[G-6];"
                "s__[Eubacterium]_nodatum"),
            ("Peptostreptococcaceae [G-6]", "Eubacterium nodatum"))
        self.assertEqual(
            parse_taxonomy_string(
                "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;"
                "f__Bacteroidales_[F-2];g__Bacteroidales_[G-2];s__bacterium_HMT_274"),
            ("Bacteroidales [G-2]", "Bacteroidales [G-2] bacterium HMT 274"))

    def test_ehomd_subspecies_and_hmt(self):
        self.assertEqual(
            parse_taxonomy_string(
                "k__Bacteria;g__Fusobacterium;s__nucleatum_subsp._polymorphum")[1],
            "Fusobacterium nucleatum subsp. polymorphum")
        self.assertEqual(
            parse_taxonomy_string("k__Bacteria;g__Fretibacterium;s__sp._HMT_360")[1],
            "Fretibacterium sp. HMT 360")

    def test_multi_species_labels(self):
        _, species = parse_taxonomy_string(
            "d__Bacteria; g__Streptococcus; s__Streptococcus_mitis/oralis")
        self.assertEqual(species_candidates(species),
                         ["Streptococcus mitis", "Streptococcus oralis"])
        self.assertEqual(species_candidates("Streptococcus mitis"),
                         ["Streptococcus mitis"])
        self.assertEqual(species_candidates(None), [])

    def test_ehomd_clades_are_dropped(self):
        self.assertEqual(
            parse_taxonomy_string(
                "k__Bacteria;g__Streptococcus;s__infantis_clade_431")[1],
            "Streptococcus infantis")
        self.assertEqual(
            parse_taxonomy_string(
                "k__Bacteria;g__Streptococcus;s__oralis_subsp._dentisani_clade_058")[1],
            "Streptococcus oralis subsp. dentisani")

    def test_taxon_key(self):
        self.assertEqual(taxon_key("Fretibacterium sp. HMT 360"), "hmt 360")
        self.assertEqual(taxon_key("Fretibacterium sp. oral taxon 360"),
                         "hmt 360")
        self.assertEqual(taxon_key("Bacteroidales [G-2] bacterium HMT 274"),
                         "hmt 274")
        self.assertEqual(taxon_key("Treponema Denticola"), "treponema denticola")

    def test_hmt_number(self):
        self.assertEqual(hmt_number("Fretibacterium sp. HMT 360"), 360)
        self.assertEqual(hmt_number("Desulfobulbus sp. oral taxon 041"), 41)
        self.assertEqual(hmt_number("Fretibacterium human oral taxon 361"), 361)
        self.assertEqual(hmt_number("TM7 sp. HOT-356"), 356)
        self.assertIsNone(hmt_number("Treponema medium"))
        self.assertIsNone(hmt_number(None))


class TestReference(TestPluginBase):
    package = "q2_periodontitis.tests"

    def setUp(self):
        super().setUp()
        self.synonyms = load_synonyms()
        self.species, self.genera = build_reference(load_complex_reference(),
                                                    self.synonyms)

    def assign(self, tax_str):
        return assign_feature_to_complex(tax_str, self.species, self.genera,
                                         self.synonyms)

    def test_reference_is_documented(self):
        ref = load_complex_reference()
        self.assertEqual(list(ref.columns),
                         ["taxon", "complex", "level", "source", "evidence"])
        self.assertTrue(ref["complex"].isin(COMPLEXES).all())
        self.assertNotIn("Other", set(ref["complex"]))
        self.assertFalse(ref["taxon"].duplicated().any())
        self.assertTrue(ref["source"].isin(
            ["socransky1998microbial", "socransky2002dental"]).all())
        self.assertFalse(ref["evidence"].isna().any())
        # The only genus-level entry is the Actinomyces (Blue) group
        genera = ref[ref["level"] == "genus"]
        self.assertEqual(list(genera["taxon"]), ["Actinomyces"])

    def test_red_complex(self):
        self.assertEqual(
            {t for t, c in self.species.items() if c == "Red"},
            {"porphyromonas gingivalis", "tannerella forsythia",
             "treponema denticola"})

    def test_species_only_for_heterogeneous_genera(self):
        for tax in ["g__Prevotella_7; s__Prevotella_melaninogenica",
                    "g__Porphyromonas; s__Porphyromonas_catoniae",
                    "g__Streptococcus; s__Streptococcus_mutans",
                    "g__Fusobacterium; s__nucleatum_subsp._animalis",
                    "g__Treponema; s__Treponema_medium",
                    "g__Treponema",
                    "g__Aggregatibacter; s__actinomycetemcomitans"]:
            self.assertEqual(self.assign(tax), ("Other", "unassigned"), tax)
        self.assertEqual(self.assign("g__Prevotella; s__Prevotella_intermedia"),
                         ("Orange", "species"))

    def test_actinomyces_genus_is_blue(self):
        self.assertEqual(self.assign("g__Actinomyces; s__sp._HMT_169"),
                         ("Blue", "genus"))
        self.assertEqual(self.assign("g__Actinomyces; s__Actinomyces_naeslundii"),
                         ("Blue", "species"))
        # A. odontolyticus (now Schaalia) is Purple, not Blue
        self.assertEqual(self.assign("g__Actinomyces; s__Actinomyces_odontolyticus"),
                         ("Purple", "species"))
        self.assertEqual(self.assign("g__Schaalia; s__odontolytica"),
                         ("Purple", "species"))

    def test_synonyms_resolve_to_complex(self):
        self.assertEqual(self.assign("g__Tannerella; s__Bacteroides_forsythus"),
                         ("Red", "species"))
        self.assertEqual(self.assign("g__Streptococcus; s__sanguis"),
                         ("Yellow", "species"))
        self.assertEqual(
            self.assign("g__Peptostreptococcaceae_[G-6]; s__[Eubacterium]_nodatum"),
            ("Orange", "species"))

    def test_ambiguous_labels(self):
        self.assertEqual(self.assign("g__Streptococcus; s__Streptococcus_mitis/oralis"),
                         ("Yellow", "species"))
        self.assertEqual(self.assign("g__Streptococcus; s__Streptococcus_mitis/mutans"),
                         ("Other", "ambiguous"))
        self.assertEqual(
            self.assign("g__Streptococcus; s__Streptococcus_constellatus/gordonii"),
            ("Other", "ambiguous"))

    def test_invalid_complex(self):
        ref = pd.DataFrame({"taxon": ["Selenomonas noxia"], "complex": ["Pink"]})
        with self.assertRaisesRegex(ValueError, "Invalid complex name.*'Pink'"):
            build_reference(ref, self.synonyms)

    def test_invalid_level(self):
        ref = pd.DataFrame({"taxon": ["Selenomonas"], "complex": ["Red"],
                            "level": ["family"]})
        with self.assertRaisesRegex(ValueError, "Invalid level"):
            build_reference(ref, self.synonyms)

    def test_level_inferred_from_words(self):
        ref = pd.DataFrame({"taxon": ["Selenomonas", "Selenomonas noxia"],
                            "complex": ["Red", "Green"]})
        species, genera = build_reference(ref, self.synonyms)
        self.assertEqual(species, {"selenomonas noxia": "Green"})
        self.assertEqual(genera, {"selenomonas": "Red"})

    def test_builtin_signatures(self):
        pc = pd.read_csv(PC2014_REFERENCE, sep="\t", dtype=str)
        # 17 species/phylotypes with moderate evidence (Table 3)
        self.assertEqual(len(pc), 17)
        self.assertFalse(pc["hmt"].duplicated().any())
        full = load_signature(SMDI_SPECIES)
        # 49 species with MDG >= 0.4: 31 dysbiotic, 18 normobiotic
        self.assertEqual(len(full), 49)
        self.assertEqual((full["direction"] == "dysbiotic").sum(), 31)
        self.assertTrue((full["mean_decrease_gini"].astype(float) >= 0.4).all())
        self.assertFalse(full["taxon"].map(taxon_key).duplicated().any())
        sig = load_signature(SMDI_SIMPLIFIED)
        self.assertEqual(dict(zip(sig["taxon"], sig["direction"])),
                         {"Fretibacterium": "dysbiotic", "Treponema": "dysbiotic",
                          "Actinomyces": "normobiotic"})


class CollapseFixture(TestPluginBase):
    package = "q2_periodontitis.tests"

    # Counts (features x samples); totals S1 = 100, S2 = 1000, S3 = 100, S4 = 0
    COUNTS = np.array([
        [10, 5, 0, 0],      # F0 P. gingivalis -> Red
        [20, 0, 50, 0],     # F1 F. nucleatum -> Orange
        [30, 600, 50, 0],   # F2 S. gordonii -> Yellow
        [20, 395, 0, 0],    # F3 A. naeslundii -> Blue
        [10, 0, 0, 0],      # F4 V. parvula -> Purple
        [10, 0, 0, 0],      # F5 L. casei -> Other
    ])
    TAXA = {
        "F0": "d__Bacteria; g__Porphyromonas; s__Porphyromonas_gingivalis",
        "F1": "d__Bacteria; g__Fusobacterium; s__Fusobacterium_nucleatum",
        "F2": "d__Bacteria; g__Streptococcus; s__Streptococcus_gordonii",
        "F3": "d__Bacteria; g__Actinomyces; s__Actinomyces_naeslundii",
        "F4": "d__Bacteria; g__Veillonella; s__Veillonella_parvula",
        "F5": "d__Bacteria; g__Lactobacillus; s__Lactobacillus_casei",
    }
    # Half the smallest non-zero relative abundance of the table (5/1000)
    PC = 0.0025

    def setUp(self):
        super().setUp()
        self.table = biom.Table(self.COUNTS, list(self.TAXA),
                                ["S1", "S2", "S3", "S4"])
        self.taxonomy = _taxonomy(self.TAXA)


class TestCollapseComplexes(CollapseFixture):

    def test_complexes_hand_computed(self):
        complexes, _, mapping = collapse_complexes(self.table, self.taxonomy)
        self.assertEqual(list(complexes.columns), COMPLEXES)
        # Red Orange Yellow Green Purple Blue Other
        np.testing.assert_allclose(complexes.loc["S1"].values,
                                   [0.1, 0.2, 0.3, 0.0, 0.1, 0.2, 0.1])
        np.testing.assert_allclose(complexes.sum(axis=1).values, [1, 1, 1, 0])
        self.assertEqual(mapping.loc["F5", "Complex"], "Other")
        self.assertEqual(mapping.loc["F5", "Resolution"], "unassigned")
        self.assertEqual(mapping.loc["F0", "Species"], "Porphyromonas gingivalis")

    def test_indices_columns(self):
        _, indices, _ = collapse_complexes(self.table, self.taxonomy)
        self.assertEqual(list(indices.columns),
                         ["SMDI", "SMDI_simplified", "PC2014_moderate",
                          "qc_depth", "qc_assigned_fraction",
                          "qc_species_fraction"])




    def test_smdi_simplified_hand_computed(self):
        _, indices, _ = collapse_complexes(self.table, self.taxonomy)
        # Fretibacterium and Treponema are absent (pseudocount only), so
        # SMDI = ln(pc) - ln(Actinomyces + pc); the CLR means cancel out.
        pc = self.PC
        self.assertAlmostEqual(indices.loc["S1", "SMDI_simplified"],
                               math.log(pc / (0.2 + pc)), places=12)
        self.assertAlmostEqual(indices.loc["S2", "SMDI_simplified"],
                               math.log(pc / (0.395 + pc)), places=12)
        self.assertAlmostEqual(indices.loc["S3", "SMDI_simplified"], 0.0,
                               places=12)
        self.assertTrue(np.isnan(indices.loc["S4", "SMDI_simplified"]))

    def test_custom_smdi_signature(self):
        sig = pd.DataFrame({"taxon": ["Porphyromonas gingivalis",
                                      "Streptococcus gordonii"],
                            "direction": ["dysbiotic", "normobiotic"],
                            "level": ["species", "species"]})
        _, indices, _ = collapse_complexes(self.table, self.taxonomy,
                                           smdi_signature=sig)
        pc = self.PC
        self.assertAlmostEqual(indices.loc["S1", "SMDI_custom"],
                               math.log((0.1 + pc) / (0.3 + pc)), places=12)
        self.assertAlmostEqual(indices.loc["S3", "SMDI_custom"],
                               math.log(pc / (0.5 + pc)), places=12)

    def test_full_smdi_hand_computed(self):
        _, indices, _ = collapse_complexes(self.table, self.taxonomy)
        # Of the 49 signature species only P. gingivalis (dysbiotic),
        # A. naeslundii and V. parvula (normobiotic) are in the table; the
        # other 46 enter at the pseudocount. The CLR mean cancels out.
        pc = self.PC
        ln = math.log
        dys = (ln(0.1 + pc) + 30 * ln(pc)) / 31
        norm = (ln(0.2 + pc) + ln(0.1 + pc) + 16 * ln(pc)) / 18
        self.assertAlmostEqual(indices.loc["S1", "SMDI"], dys - norm, places=12)
        self.assertTrue(np.isnan(indices.loc["S4", "SMDI"]))

    def test_fixed_smdi_pseudocount(self):
        eps = 1e-6
        _, indices, _ = collapse_complexes(self.table, self.taxonomy,
                                           smdi_pseudocount=eps)
        self.assertAlmostEqual(indices.loc["S1", "SMDI_simplified"],
                               math.log(eps / (0.2 + eps)), places=12)
        # S3 has no Actinomyces: all three genera sit at the pseudocount
        self.assertAlmostEqual(indices.loc["S3", "SMDI_simplified"], 0.0,
                               places=12)
        ln = math.log
        dys = (ln(0.1 + eps) + 30 * ln(eps)) / 31
        norm = (ln(0.2 + eps) + ln(0.1 + eps) + 16 * ln(eps)) / 18
        self.assertAlmostEqual(indices.loc["S1", "SMDI"], dys - norm,
                               places=12)
        self.assertTrue(np.isnan(indices.loc["S4", "SMDI"]))

    def test_fixed_pseudocount_is_table_independent(self):
        eps = 1e-5
        _, full, _ = collapse_complexes(self.table, self.taxonomy,
                                        smdi_pseudocount=eps)
        alone = self.table.filter(["S1"], inplace=False)
        _, single, _ = collapse_complexes(alone, self.taxonomy,
                                          smdi_pseudocount=eps)
        self.assertAlmostEqual(full.loc["S1", "SMDI"], single.loc["S1", "SMDI"],
                               places=12)
        # The original pseudocount changes when S2 (1000 reads) is removed
        _, full0, _ = collapse_complexes(self.table, self.taxonomy)
        _, single0, _ = collapse_complexes(alone, self.taxonomy)
        self.assertNotAlmostEqual(full0.loc["S1", "SMDI"],
                                  single0.loc["S1", "SMDI"], places=6)

    def test_invalid_smdi_pseudocount(self):
        for value in (0.0, -1e-6, 1.0):
            with self.assertRaisesRegex(ValueError, "smdi_pseudocount"):
                collapse_complexes(self.table, self.taxonomy,
                                   smdi_pseudocount=value)

    def test_signature_taxa_missing_from_table_warn(self):
        sig = pd.DataFrame({"taxon": ["Filifactor alocis",
                                      "Streptococcus gordonii"],
                            "direction": ["dysbiotic", "normobiotic"],
                            "level": ["species", "species"]})
        with self.assertWarnsRegex(UserWarning, "SMDI_custom.*filifactor alocis"):
            collapse_complexes(self.table, self.taxonomy, smdi_signature=sig)

    def test_invalid_signatures(self):
        base = {"taxon": ["Treponema", "Actinomyces"],
                "direction": ["dysbiotic", "normobiotic"],
                "level": ["genus", "genus"]}
        cases = [
            ({**base, "direction": ["dysbiotic", "up"]}, "Invalid direction"),
            ({**base, "direction": ["dysbiotic", "dysbiotic"]}, "normobiotic"),
            ({**base, "level": ["genus", "species"]}, "one level"),
            ({**base, "taxon": ["Treponema", "treponema"]}, "more than once"),
        ]
        for data, message in cases:
            with self.assertRaisesRegex(ValueError, message):
                compute_smdi(self.table, self.taxonomy, pd.DataFrame(data), {})
        with self.assertRaisesRegex(ValueError, "direction"):
            compute_smdi(self.table, self.taxonomy,
                         pd.DataFrame({"taxon": ["Treponema"]}), {})

    def test_qc_columns(self):
        conf = {"F0": 0.99, "F1": 0.5, "F2": 0.95, "F3": 0.95, "F4": 0.95,
                "F5": 0.95}
        _, indices, _ = collapse_complexes(self.table,
                                           _taxonomy(self.TAXA, conf),
                                           min_confidence=0.7)
        np.testing.assert_allclose(indices["qc_depth"], [100, 1000, 100, 0])
        # S1: F1 (20) is low-confidence, F5 (10) is unassigned
        self.assertAlmostEqual(indices.loc["S1", "qc_assigned_fraction"], 0.7)
        self.assertAlmostEqual(indices.loc["S1", "qc_species_fraction"], 0.8)
        self.assertAlmostEqual(indices.loc["S2", "qc_species_fraction"], 1.0)
        self.assertTrue(np.isnan(indices.loc["S4", "qc_assigned_fraction"]))

    def test_pc2014_fraction(self):
        taxa = {
            "A": "k__Bacteria;g__Fretibacterium;s__sp._HMT_360",
            "B": "k__Bacteria;g__Filifactor;s__alocis",
            "C": "k__Bacteria;g__Bacteroidales_[G-2];s__bacterium_HMT_274",
            "D": "k__Bacteria;g__Treponema;s__denticola",
            "E": "k__Bacteria;g__Fretibacterium;s__sp._HMT_361",
        }
        table = biom.Table(np.array([[10], [20], [5], [40], [25]]), list(taxa),
                           ["S1"])
        _, indices, mapping = collapse_complexes(table, _taxonomy(taxa))
        self.assertAlmostEqual(indices.loc["S1", "PC2014_moderate"], 0.35)
        self.assertEqual(list(mapping["PC2014"]), ["yes", "yes", "yes", "no", "no"])
        self.assertEqual(mapping.loc["D", "Complex"], "Red")

    def test_min_confidence(self):
        conf = {"F0": 0.99, "F1": 0.5, "F2": 0.95, "F3": 0.95, "F4": 0.95,
                "F5": 0.95}
        taxonomy = _taxonomy(self.TAXA, conf)
        complexes, _, mapping = collapse_complexes(self.table, taxonomy,
                                                   min_confidence=0.7)
        self.assertEqual(mapping.loc["F1", "Resolution"], "low-confidence")
        self.assertAlmostEqual(complexes.loc["S1", "Orange"], 0.0)
        self.assertAlmostEqual(complexes.loc["S1", "Other"], 0.3)
        np.testing.assert_allclose(complexes.sum(axis=1).values, [1, 1, 1, 0])

    def test_min_confidence_applies_to_smdi(self):
        conf = {"F0": 0.5, "F1": 0.95, "F2": 0.95, "F3": 0.95, "F4": 0.95,
                "F5": 0.95}
        sig = pd.DataFrame({"taxon": ["Porphyromonas gingivalis",
                                      "Streptococcus gordonii"],
                            "direction": ["dysbiotic", "normobiotic"],
                            "level": ["species", "species"]})
        _, indices, _ = collapse_complexes(self.table,
                                           _taxonomy(self.TAXA, conf),
                                           smdi_signature=sig,
                                           min_confidence=0.7)
        # F0 is low-confidence: P. gingivalis enters at the pseudocount,
        # while its reads stay in the table (the pseudocount is unchanged).
        pc = self.PC
        self.assertAlmostEqual(indices.loc["S1", "SMDI_custom"],
                               math.log(pc / (0.3 + pc)), places=12)

    def test_missing_confidence_is_not_filtered(self):
        taxonomy = _taxonomy(self.TAXA, {k: "-" for k in self.TAXA})
        _, _, mapping = collapse_complexes(self.table, taxonomy,
                                           min_confidence=0.9)
        self.assertEqual(mapping.loc["F0", "Complex"], "Red")

    def test_features_without_taxonomy_raise(self):
        taxonomy = _taxonomy({k: v for k, v in self.TAXA.items() if k != "F5"})
        with self.assertRaisesRegex(ValueError, "no taxonomy.*F5"):
            collapse_complexes(self.table, taxonomy)

    def test_custom_mapping_overrides(self):
        custom = pd.DataFrame({"taxon": ["Lactobacillus casei",
                                         "streptococcus gordonii"],
                               "complex": ["Green", "Green"],
                               "note": ["x", "y"]})
        complexes, _, mapping = collapse_complexes(self.table, self.taxonomy,
                                                   custom_mapping=custom)
        self.assertEqual(mapping.loc["F5", "Complex"], "Green")
        self.assertEqual(mapping.loc["F2", "Complex"], "Green")
        self.assertAlmostEqual(complexes.loc["S1", "Green"], 0.4)

    def test_custom_mapping_invalid(self):
        custom = pd.DataFrame({"taxon": ["Lactobacillus casei"],
                               "complex": ["Extended"]})
        with self.assertRaisesRegex(ValueError, "Invalid complex name"):
            collapse_complexes(self.table, self.taxonomy, custom_mapping=custom)
        custom = pd.DataFrame({"name": ["Lactobacillus casei"],
                               "complex": ["Red"]})
        with self.assertRaisesRegex(ValueError, "taxon"):
            collapse_complexes(self.table, self.taxonomy, custom_mapping=custom)

    def test_action_with_artifacts(self):
        table = qiime2.Artifact.import_data("FeatureTable[Frequency]", self.table)
        taxonomy = qiime2.Artifact.import_data("FeatureData[Taxonomy]",
                                               self.taxonomy)
        reference = qiime2.Artifact.import_data(
            "SocranskyReference", self.get_data_path("reference-valid.tsv"))
        signature = qiime2.Artifact.import_data(
            "DysbiosisSignature", self.get_data_path("signature-valid.tsv"))
        action = self.plugin.actions["collapse_complexes"]
        complexes, indices, mapping = action(
            table=table, taxonomy=taxonomy, custom_mapping=reference,
            smdi_signature=signature)
        self.assertEqual(str(indices.type), "PeriodontalIndices")
        self.assertEqual(mapping.view(pd.DataFrame).loc["F5", "Complex"],
                         "Yellow")
        indices_df = indices.view(pd.DataFrame)
        self.assertIn("SMDI_custom", indices_df.columns)
        # NaN for the empty sample survives the round trip
        self.assertTrue(np.isnan(indices_df.loc["S4", "SMDI"]))
        np.testing.assert_allclose(
            complexes.view(pd.DataFrame).sum(axis=1).values, [1, 1, 1, 0])


class TestCollapseEmpty(TestPluginBase):
    package = "q2_periodontitis.tests"

    INDEX_COLUMNS = [
        "SMDI", "SMDI_simplified", "PC2014_moderate",
        "qc_depth", "qc_assigned_fraction", "qc_species_fraction"
    ]
    MAPPING_COLUMNS = [
        "Taxonomy", "Genus", "Species", "Complex", "Resolution",
        "PC2014", "Confidence"
    ]

    def test_zero_features(self):
        table = biom.Table(np.zeros((0, 2)), [], ["S1", "S2"])
        taxonomy = pd.DataFrame(
            {"Taxon": []},
            index=pd.Index([], name="Feature ID"))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            complexes, indices, mapping = collapse_complexes(table, taxonomy)

        self.assertEqual(list(complexes.index), ["S1", "S2"])
        self.assertEqual(list(complexes.columns), COMPLEXES)
        np.testing.assert_allclose(complexes.to_numpy(), 0.0)

        self.assertEqual(list(indices.index), ["S1", "S2"])
        self.assertEqual(list(indices.columns), self.INDEX_COLUMNS)
        for column in ["SMDI", "SMDI_simplified", "PC2014_moderate",
                       "qc_assigned_fraction", "qc_species_fraction"]:
            self.assertTrue(indices[column].isna().all(), column)
        self.assertTrue((indices["qc_depth"] == 0).all())

        self.assertTrue(mapping.empty)
        self.assertEqual(list(mapping.columns), self.MAPPING_COLUMNS)

    def test_zero_samples(self):
        taxa = {
            "F0": "d__Bacteria; g__Porphyromonas; s__Porphyromonas_gingivalis",
            "F1": "d__Bacteria; g__Fusobacterium; s__Fusobacterium_nucleatum",
        }
        table = biom.Table(np.zeros((2, 0)), list(taxa), [])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            complexes, indices, mapping = collapse_complexes(
                table, _taxonomy(taxa))

        self.assertEqual(complexes.shape, (0, len(COMPLEXES)))
        self.assertEqual(list(complexes.columns), COMPLEXES)
        self.assertEqual(indices.shape[0], 0)
        self.assertEqual(list(indices.columns), self.INDEX_COLUMNS)
        self.assertEqual(list(mapping.index), ["F0", "F1"])
        self.assertEqual(mapping.loc["F0", "Complex"], "Red")
        self.assertEqual(mapping.loc["F1", "Complex"], "Orange")


class TestClr(CollapseFixture):

    def test_matches_microbiome_definition(self):
        counts = aggregate_counts(self.table, self.taxonomy, {}, "genus")
        self.assertEqual(sorted(counts.columns),
                         ["actinomyces", "fusobacterium", "lactobacillus",
                          "porphyromonas", "streptococcus", "veillonella"])
        clr = clr_microbiome(counts)
        # The pseudocount is added to every entry, not only to zeros
        rel = counts.loc["S1"] / 100 + self.PC
        expected = np.log(rel) - np.log(rel).mean()
        np.testing.assert_allclose(clr.loc["S1"], expected)
        np.testing.assert_allclose(clr.loc[["S1", "S2", "S3"]].sum(axis=1), 0,
                                   atol=1e-12)
        self.assertTrue(clr.loc["S4"].isna().all())

    def test_fixed_pseudocount(self):
        counts = aggregate_counts(self.table, self.taxonomy, {}, "genus")
        clr = clr_microbiome(counts, pseudocount=1e-3)
        rel = counts.loc["S2"] / 1000 + 1e-3
        np.testing.assert_allclose(clr.loc["S2"],
                                   np.log(rel) - np.log(rel).mean())
        with self.assertRaisesRegex(ValueError, "greater than 0"):
            clr_microbiome(counts, pseudocount=0)

    def test_unresolved_features_keep_their_reads(self):
        taxa = dict(self.TAXA, F5="d__Bacteria; g__uncultured")
        counts = aggregate_counts(self.table, _taxonomy(taxa), {}, "species")
        np.testing.assert_allclose(counts.sum(axis=1), [100, 1000, 100, 0])
        self.assertIn("unresolved: d__Bacteria; g__uncultured", counts.columns)


class TestDeltaIndices(TestPluginBase):
    package = "q2_periodontitis.tests"

    def setUp(self):
        super().setUp()
        self.indices = pd.DataFrame({
            "SMDI": [1.5, -0.2, 0.5, 0.0],
            "SMDI_simplified": [2.0, 1.0, 0.0, 0.5],
            "qc_depth": [1000, 2000, 1500, 900],
        }, index=["a", "b", "c", "d"])

    def metadata(self, subject, group, extra=None):
        df = pd.DataFrame({"subject": subject, "visit": group},
                          index=pd.Index(["a", "b", "c", "d"][:len(subject)],
                                         name="id"))
        for k, v in (extra or {}).items():
            df[k] = v
        return qiime2.Metadata(df)

    def test_success(self):
        md = self.metadata(["P1", "P1", "P2", "P2"], ["Pre", "Post", "Pre", "Post"],
                           extra={"SMDI": ["x"] * 4})  # colliding column name
        res = delta_indices(self.indices, md, subject_col="subject",
                            group_col="visit", baseline_group="Pre",
                            target_group="Post")
        self.assertEqual(list(res.index), ["P1", "P2"])
        self.assertEqual(res.index.name, "subject-id")
        # QC columns are not differenced
        self.assertEqual(list(res.columns), ["SMDI", "SMDI_simplified"])
        self.assertAlmostEqual(res.loc["P1", "SMDI"], -1.7)
        self.assertAlmostEqual(res.loc["P2", "SMDI_simplified"], 0.5)

    def test_numeric_group_column(self):
        md = self.metadata(["P1", "P1", "P2", "P2"], [0, 1, 0, 1])
        self.assertEqual(md.columns["visit"].type, "numeric")
        res = delta_indices(self.indices, md, "subject", "visit", "0", "1")
        self.assertAlmostEqual(res.loc["P2", "SMDI"], -0.5)

    def test_unpaired_subjects_are_skipped(self):
        md = self.metadata(["P1", "P1", "P2", "P3"], ["Pre", "Post", "Pre", "Post"])
        res = delta_indices(self.indices, md, "subject", "visit", "Pre", "Post")
        self.assertEqual(list(res.index), ["P1"])

    def test_repeated_samples_are_averaged(self):
        # P1 has two baseline sites: their mean (1.5 + 0.5) / 2 = 1.0 is used
        md = self.metadata(["P1", "P1", "P1", "P2"], ["Pre", "Post", "Pre", "Pre"])
        res = delta_indices(self.indices, md, "subject", "visit", "Pre", "Post")
        self.assertEqual(list(res.index), ["P1"])
        self.assertAlmostEqual(res.loc["P1", "SMDI"], -0.2 - 1.0)

    def test_no_pairs_raise(self):
        md = self.metadata(["P1", "P2", "P3", "P4"], ["Pre", "Post", "Pre", "Post"])
        with self.assertRaisesRegex(ValueError, "No subject"):
            delta_indices(self.indices, md, "subject", "visit", "Pre", "Post")

    def test_missing_column_raise(self):
        md = self.metadata(["P1", "P1"], ["Pre", "Post"])
        with self.assertRaisesRegex(ValueError, "not found: patient"):
            delta_indices(self.indices, md, "patient", "visit", "Pre", "Post")

    def test_action_output_type(self):
        indices = qiime2.Artifact.import_data("PeriodontalIndices", self.indices)
        md = self.metadata(["P1", "P1", "P2", "P2"], ["Pre", "Post", "Pre", "Post"])
        delta, = self.plugin.actions["delta_indices"](
            indices=indices, metadata=md, subject_col="subject",
            group_col="visit", baseline_group="Pre", target_group="Post")
        self.assertEqual(str(delta.type), "SampleData[PeriodontalIndicesDelta]")
        self.assertEqual(list(delta.view(pd.DataFrame).index), ["P1", "P2"])


if __name__ == "__main__":
    unittest.main()
