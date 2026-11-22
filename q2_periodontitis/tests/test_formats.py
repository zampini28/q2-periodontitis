import unittest

import pandas as pd
import qiime2
from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_periodontitis._format import (
    SocranskyComplexesFormat,
    SocranskyComplexesDirectoryFormat,
    PeriodontalIndicesFormat,
    PeriodontalIndicesDirectoryFormat,
    PeriodontalIndicesDeltaFormat,
    PeriodontalIndicesDeltaDirectoryFormat,
    SocranskyMappingFormat,
    SocranskyMappingDirectoryFormat,
    SocranskyReferenceFormat,
    SocranskyReferenceDirectoryFormat,
    DysbiosisSignatureFormat,
    DysbiosisSignatureDirectoryFormat
)
from q2_periodontitis._type import (
    SocranskyComplexes,
    PeriodontalIndices,
    PeriodontalIndicesDelta,
    SocranskyMapping,
    SocranskyReference,
    DysbiosisSignature
)
from q2_types.sample_data import SampleData


class TestTypes(TestPluginBase):
    package = "q2_periodontitis.tests"

    def test_types_registered(self):
        for semantic_type in (SocranskyComplexes, PeriodontalIndices,
                              PeriodontalIndicesDelta, SocranskyMapping,
                              SocranskyReference, DysbiosisSignature):
            self.assertRegisteredSemanticType(semantic_type)

    def test_types_to_formats(self):
        self.assertSemanticTypeRegisteredToFormat(
            SocranskyComplexes, SocranskyComplexesDirectoryFormat)
        self.assertSemanticTypeRegisteredToFormat(
            PeriodontalIndices, PeriodontalIndicesDirectoryFormat)
        self.assertSemanticTypeRegisteredToFormat(
            SampleData[PeriodontalIndicesDelta],
            PeriodontalIndicesDeltaDirectoryFormat)
        self.assertSemanticTypeRegisteredToFormat(
            SocranskyMapping, SocranskyMappingDirectoryFormat)
        self.assertSemanticTypeRegisteredToFormat(
            SocranskyReference, SocranskyReferenceDirectoryFormat)
        self.assertSemanticTypeRegisteredToFormat(
            DysbiosisSignature, DysbiosisSignatureDirectoryFormat)


class TestFormatValidation(TestPluginBase):
    package = "q2_periodontitis.tests"

    def test_valid_formats(self):
        for fmt, name in [(SocranskyComplexesFormat, "complexes.tsv"),
                          (PeriodontalIndicesFormat, "indices.tsv"),
                          (PeriodontalIndicesDeltaFormat, "deltas.tsv"),
                          (SocranskyMappingFormat, "mapping.tsv"),
                          (SocranskyReferenceFormat, "reference-valid.tsv"),
                          (SocranskyReferenceFormat,
                           "reference-extra-column.tsv"),
                          (DysbiosisSignatureFormat, "signature-valid.tsv")]:
            fmt(self.get_data_path(name), mode="r").validate()

    def test_indices_non_numeric(self):
        with self.assertRaisesRegex(ValidationError, "non-numeric"):
            PeriodontalIndicesFormat(
                self.get_data_path("indices-non-numeric.tsv"),
                mode="r").validate()

    def test_indices_ragged(self):
        with self.assertRaisesRegex(ValidationError, "fields"):
            PeriodontalIndicesFormat(self.get_data_path("indices-ragged.tsv"),
                                     mode="r").validate()

    def test_indices_missing_values_are_valid(self):
        # indices.tsv has empty fields for an empty sample
        PeriodontalIndicesFormat(self.get_data_path("indices.tsv"),
                                 mode="r").validate()

    def test_signature_invalid(self):
        with self.assertRaisesRegex(ValidationError,
                                    "invalid direction 'bad'"):
            DysbiosisSignatureFormat(
                self.get_data_path("signature-invalid-direction.tsv"),
                mode="r").validate()
        with self.assertRaisesRegex(ValidationError, "invalid level 'phylum'"):
            DysbiosisSignatureFormat(
                self.get_data_path("signature-invalid-level.tsv"),
                mode="r").validate()
        with self.assertRaisesRegex(ValidationError, "taxon.*direction"):
            DysbiosisSignatureFormat(
                self.get_data_path("reference-valid.tsv"), mode="r").validate()

    def test_reference_invalid_complex(self):
        with self.assertRaisesRegex(ValidationError, "invalid complex 'Pink'"):
            SocranskyReferenceFormat(
                self.get_data_path("reference-invalid-complex.tsv"),
                mode="r").validate()

    def test_reference_requires_header(self):
        with self.assertRaisesRegex(ValidationError, "taxon.*complex"):
            SocranskyReferenceFormat(
                self.get_data_path("reference-no-header.tsv"),
                mode="r").validate()

    def test_reference_import_rejects_invalid(self):
        with self.assertRaises(ValidationError):
            qiime2.Artifact.import_data(
                "SocranskyReference",
                self.get_data_path("reference-invalid-complex.tsv"))


class TestTransformers(TestPluginBase):
    package = "q2_periodontitis.tests"

    def roundtrip(self, fmt, filename):
        _, df = self.transform_format(fmt, pd.DataFrame, filename=filename)
        transformer = self.get_transformer(pd.DataFrame, fmt)
        df2 = self.get_transformer(fmt, pd.DataFrame)(transformer(df))
        pd.testing.assert_frame_equal(df, df2)
        return df

    def test_complexes_roundtrip(self):
        df = self.roundtrip(SocranskyComplexesFormat, "complexes.tsv")
        self.assertEqual(list(df.index), ["S1", "S2"])
        self.assertAlmostEqual(df.loc["S1", "Yellow"], 0.3)

    def test_indices_roundtrip(self):
        df = self.roundtrip(PeriodontalIndicesFormat, "indices.tsv")
        self.assertEqual(list(df.columns),
                         ["SMDI", "SMDI_simplified", "PC2014_moderate",
                          "qc_depth"])
        self.assertTrue(df.loc["S2", ["SMDI", "SMDI_simplified"]].isna().all())

    def test_delta_roundtrip_keeps_string_subject_ids(self):
        df = self.roundtrip(PeriodontalIndicesDeltaFormat, "deltas.tsv")
        self.assertEqual(list(df.index), ["001", "P2"])
        self.assertEqual(df.index.name, "subject-id")

    def test_mapping_roundtrip(self):
        df = self.roundtrip(SocranskyMappingFormat, "mapping.tsv")
        self.assertEqual(df.loc["F1", "Complex"], "Red")

    def test_reference_roundtrip(self):
        df = self.roundtrip(SocranskyReferenceFormat, "reference-valid.tsv")
        self.assertEqual(list(df["taxon"]),
                         ["Selenomonas noxia", "Lactobacillus casei"])
        self.assertEqual(list(df["complex"]), ["Green", "Yellow"])

    def test_signature_roundtrip(self):
        df = self.roundtrip(DysbiosisSignatureFormat, "signature-valid.tsv")
        self.assertEqual(list(df["direction"]), ["dysbiotic", "normobiotic"])

    def test_to_metadata(self):
        for fmt, name in [(SocranskyComplexesFormat, "complexes.tsv"),
                          (PeriodontalIndicesFormat, "indices.tsv"),
                          (PeriodontalIndicesDeltaFormat, "deltas.tsv"),
                          (SocranskyMappingFormat, "mapping.tsv")]:
            _, md = self.transform_format(fmt, qiime2.Metadata, filename=name)
            self.assertIsInstance(md, qiime2.Metadata)
        _, md = self.transform_format(PeriodontalIndicesDeltaFormat,
                                      qiime2.Metadata, filename="deltas.tsv")
        self.assertEqual(list(md.ids), ["001", "P2"])


if __name__ == "__main__":
    unittest.main()
