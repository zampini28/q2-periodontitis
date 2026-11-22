from qiime2.plugin import (
    Plugin,
    Citations,
    Str,
    Float,
    Range,
    Metadata
)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy
from q2_types.sample_data import SampleData
from ._version import __version__
from ._type import (
    SocranskyComplexes,
    PeriodontalIndices,
    PeriodontalIndicesDelta,
    SocranskyMapping,
    SocranskyReference,
    DysbiosisSignature
)
from ._format import (
    SocranskyComplexesDirectoryFormat,
    PeriodontalIndicesDirectoryFormat,
    PeriodontalIndicesDeltaDirectoryFormat,
    SocranskyMappingDirectoryFormat,
    SocranskyReferenceDirectoryFormat,
    DysbiosisSignatureDirectoryFormat
)
from ._methods import collapse_complexes, delta_indices
from ._visualizer import summarize

citations = Citations.load('citations.bib', package='q2_periodontitis')

plugin = Plugin(
    name='periodontitis',
    version=__version__,
    website='https://github.com/LaBiOmicS/q2-periodontitis',
    package='q2_periodontitis',
    citations=[citations['socransky1998microbial'],
               citations['socransky2002dental']],
    description='QIIME 2 plugin that collapses subgingival amplicon profiles '
                'into the Socransky microbial complexes and computes the '
                'published Subgingival Microbial Dysbiosis Index (SMDI).',
    short_description='Plugin for Periodontal Biofilm & Dysbiosis '
                      'Characterization.'
)

# Register Directory Formats
plugin.register_formats(
    SocranskyComplexesDirectoryFormat,
    PeriodontalIndicesDirectoryFormat,
    PeriodontalIndicesDeltaDirectoryFormat,
    SocranskyMappingDirectoryFormat,
    SocranskyReferenceDirectoryFormat,
    DysbiosisSignatureDirectoryFormat
)

plugin.register_semantic_types(
    SocranskyComplexes,
    PeriodontalIndices,
    PeriodontalIndicesDelta,
    SocranskyMapping,
    SocranskyReference,
    DysbiosisSignature
)

# Register Semantic Type to Format Mappings
plugin.register_semantic_type_to_format(
    SocranskyComplexes,
    artifact_format=SocranskyComplexesDirectoryFormat
)

plugin.register_semantic_type_to_format(
    PeriodontalIndices,
    artifact_format=PeriodontalIndicesDirectoryFormat
)

plugin.register_semantic_type_to_format(
    SampleData[PeriodontalIndicesDelta],
    artifact_format=PeriodontalIndicesDeltaDirectoryFormat
)

plugin.register_semantic_type_to_format(
    SocranskyMapping,
    artifact_format=SocranskyMappingDirectoryFormat
)

plugin.register_semantic_type_to_format(
    SocranskyReference,
    artifact_format=SocranskyReferenceDirectoryFormat
)

plugin.register_semantic_type_to_format(
    DysbiosisSignature,
    artifact_format=DysbiosisSignatureDirectoryFormat
)

# Register Methods
plugin.methods.register_function(
    function=collapse_complexes,
    inputs={
        'table': FeatureTable[Frequency],
        'taxonomy': FeatureData[Taxonomy],
        'custom_mapping': SocranskyReference,
        'smdi_signature': DysbiosisSignature
    },
    parameters={
        'min_confidence': Float % Range(0, 1, inclusive_end=True),
        'smdi_pseudocount': Float % Range(0, 1, inclusive_start=False)
    },
    outputs=[
        ('complex_table', SocranskyComplexes),
        ('indices', PeriodontalIndices),
        ('mapping', SocranskyMapping)
    ],
    input_descriptions={
        'table': 'Feature table of raw ASV/OTU counts (not rarefied or '
                 'normalized: the indices are computed from read counts).',
        'taxonomy': 'Taxonomic classifications for features, ideally to '
                    'species level against an oral reference such as eHOMD.',
        'custom_mapping': 'Optional taxon-to-complex reference table '
                          '(SocranskyReference: TSV with a header containing '
                          'the columns "taxon" and "complex", and optionally '
                          '"level" = species|genus). Entries are added to, '
                          'and override, the built-in reference.',
        'smdi_signature': 'Optional additional dysbiosis signature '
                          '(DysbiosisSignature: TSV with the columns "taxon", '
                          '"direction" = dysbiotic|normobiotic and "level" = '
                          'species|genus). When given, the indices include an '
                          '"SMDI_custom" column computed from it with the '
                          'SMDI method.'
    },
    parameter_descriptions={
        'min_confidence': 'Features whose taxonomy Confidence is below this '
                          'value are assigned to "Other" (Resolution '
                          '"low-confidence") and never count as SMDI or '
                          'PC2014 taxa. Features without a numeric '
                          'Confidence are never filtered. Default 0.0 (no '
                          'filtering).',
        'smdi_pseudocount': 'Fixed pseudocount added to every relative '
                            'abundance before the CLR of the SMDI indices. '
                            'By default (not set) the original definition '
                            'is used: half the smallest non-zero relative '
                            'abundance of the whole table, so SMDI values '
                            'depend on the other samples of the table. Set '
                            'it (e.g. 1e-6) to make each sample\'s SMDI '
                            'independent of the table it is computed in; '
                            'values are then not identical to the '
                            'published definition.'
    },
    output_descriptions={
        'complex_table': 'Relative abundance of each Socransky complex '
                         '(Red, Orange, Yellow, Green, Purple, Blue) and of '
                         'unassigned reads (Other); rows sum to 1.',
        'indices': 'Per-sample indices: SMDI (49 species), '
                   'SMDI_simplified (3 genera), SMDI_custom (only with '
                   '--i-smdi-signature), PC2014_moderate, and the '
                   'quality-control columns qc_depth, qc_assigned_fraction '
                   'and qc_species_fraction.',
        'mapping': 'Per-feature parsed genus/species, complex, resolution '
                   '(species, genus, ambiguous, unassigned, low-confidence) '
                   'and membership of the Perez-Chaparro et al. (2014) list.'
    },
    citations=[citations['socransky1998microbial'],
               citations['socransky2002dental'],
               citations['chen2022smdi'],
               citations['perezchaparro2014newly']],
    name='Collapse features into Socransky complexes and compute dysbiosis '
         'indices',
    description='Assigns each feature to a Socransky complex using only '
                'species (and the Actinomyces genus) documented by Socransky '
                'et al. (1998) and Socransky & Haffajee (2002), and computes '
                'the SMDI (Chen et al. 2022) and the abundance of the '
                'moderate-evidence taxa of Perez-Chaparro et al. (2014).'
)

plugin.methods.register_function(
    function=delta_indices,
    inputs={
        'indices': PeriodontalIndices
    },
    parameters={
        'metadata': Metadata,
        'subject_col': Str,
        'group_col': Str,
        'baseline_group': Str,
        'target_group': Str
    },
    outputs=[
        ('delta_table', SampleData[PeriodontalIndicesDelta])
    ],
    input_descriptions={
        'indices': 'Calculated dysbiosis indices table.'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata containing subject IDs and treatment '
                    'groups.',
        'subject_col': 'Metadata column identifying the paired unit: a '
                       'subject, or a subject-and-site label for site-level '
                       'pairing.',
        'group_col': 'Metadata column containing treatment/time-point group '
                     'identifiers.',
        'baseline_group': 'Baseline group name (e.g. Pre-treatment / '
                          'Baseline).',
        'target_group': 'Target group name (e.g. Post-treatment / Follow-up).'
    },
    output_descriptions={
        'delta_table': 'Per-unit deltas (target minus baseline) of every '
                       'dysbiosis index; quality-control (qc_) columns are '
                       'not differenced. Indexed by the values of '
                       'subject_col (subject-id), not by sample ID.'
    },
    name='Longitudinal treatment deltas',
    description='Calculates paired deltas of the dysbiosis indices between '
                'a baseline and a follow-up group. Several samples of the '
                'same unit in one group (e.g. several sites) are averaged '
                'first; units without both groups are dropped, and an error '
                'is raised if none is left.'
)

# Register Visualizers
plugin.visualizers.register_function(
    function=summarize,
    inputs={
        'complexes': SocranskyComplexes,
        'indices': PeriodontalIndices,
        'mapping': SocranskyMapping
    },
    parameters={
        'metadata': Metadata,
        'group_column': Str,
        'subject_column': Str
    },
    input_descriptions={
        'complexes': 'Collapsed Socransky complexes table.',
        'indices': 'Calculated dysbiosis indices table.',
        'mapping': 'Feature resolution mapping table.'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata.',
        'group_column': 'Categorical metadata column whose groups are '
                        'compared. Without it no statistical test is run.',
        'subject_column': 'Metadata column with subject IDs. When given, the '
                          'samples of each subject within a group are '
                          'averaged before testing, so that repeated '
                          'samples (sites, visits) are not counted as '
                          'independent observations.'
    },
    citations=[citations['benjamini1995controlling']],
    name='Summarize complexes and dysbiosis indices',
    description='Generates an offline dashboard (.qzv) with complex '
                'profiles, index distributions, taxonomic-resolution '
                'quality control and, for a chosen metadata column, '
                'non-parametric comparisons with effect sizes (AUC) and '
                'Benjamini-Hochberg q-values.'
)

import q2_periodontitis._transformer  # noqa: E402, F401 (registers them)
