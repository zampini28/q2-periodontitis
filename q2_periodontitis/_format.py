import qiime2.plugin.model as model

from ._methods import COMPLEXES, DIRECTIONS, LEVELS


class NumericTableFormat(model.TextFileFormat):
    """Tab-separated table: an ID column followed by numeric columns; empty
    fields are missing values."""

    def _validate_(self, level):
        max_lines = {'min': 100, 'max': None}[level]
        with self.open() as fh:
            header = fh.readline().rstrip('\r\n')
            if not header:
                raise model.ValidationError("File is empty.")
            columns = header.split('\t')
            if len(columns) < 2:
                raise model.ValidationError(
                    "Header line must be tab-separated with an ID column "
                    "followed by at least one value column.")
            for line_no, line in enumerate(fh, start=2):
                if max_lines is not None and line_no > max_lines + 1:
                    break
                line = line.rstrip('\r\n')
                if not line:
                    continue
                fields = line.split('\t')
                if len(fields) != len(columns):
                    raise model.ValidationError(
                        "Line %d has %d fields, expected %d."
                        % (line_no, len(fields), len(columns)))
                for value in fields[1:]:
                    if value == '':
                        # Missing value (e.g. an index undefined for an
                        # empty sample)
                        continue
                    try:
                        float(value)
                    except ValueError:
                        raise model.ValidationError(
                            "Line %d: non-numeric value %r."
                            % (line_no, value))


class SocranskyComplexesFormat(NumericTableFormat):
    pass


SocranskyComplexesDirectoryFormat = model.SingleFileDirectoryFormat(
    'SocranskyComplexesDirectoryFormat',
    'complexes.tsv',
    SocranskyComplexesFormat
)


class PeriodontalIndicesFormat(NumericTableFormat):
    pass


PeriodontalIndicesDirectoryFormat = model.SingleFileDirectoryFormat(
    'PeriodontalIndicesDirectoryFormat',
    'indices.tsv',
    PeriodontalIndicesFormat
)


class PeriodontalIndicesDeltaFormat(NumericTableFormat):
    """Per-subject deltas; the first column (``subject-id``) holds the
    subject IDs."""
    pass


PeriodontalIndicesDeltaDirectoryFormat = model.SingleFileDirectoryFormat(
    'PeriodontalIndicesDeltaDirectoryFormat',
    'deltas.tsv',
    PeriodontalIndicesDeltaFormat
)


class SocranskyMappingFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            header = fh.readline()
            if not header:
                raise model.ValidationError("File is empty.")
            if "\t" not in header:
                raise model.ValidationError(
                    "Header line must be tab-separated.")


SocranskyMappingDirectoryFormat = model.SingleFileDirectoryFormat(
    'SocranskyMappingDirectoryFormat',
    'mapping.tsv',
    SocranskyMappingFormat
)


class SocranskyReferenceFormat(model.TextFileFormat):
    """Tab-separated taxon -> complex table with a mandatory header that
    contains (at least) the columns ``taxon`` and ``complex``."""

    def _validate_(self, level):
        max_lines = {'min': 100, 'max': None}[level]
        with self.open() as fh:
            header = fh.readline().rstrip('\r\n')
            if not header:
                raise model.ValidationError("File is empty.")
            columns = [c.strip() for c in header.split('\t')]
            missing = [c for c in ('taxon', 'complex') if c not in columns]
            if missing:
                raise model.ValidationError(
                    "Header must be tab-separated and contain the columns "
                    "'taxon' and 'complex'; missing: %s." % ', '.join(missing))
            taxon_idx = columns.index('taxon')
            complex_idx = columns.index('complex')

            for line_no, line in enumerate(fh, start=2):
                if max_lines is not None and line_no > max_lines + 1:
                    break
                line = line.rstrip('\r\n')
                if not line.strip():
                    continue
                fields = line.split('\t')
                if len(fields) != len(columns):
                    raise model.ValidationError(
                        "Line %d has %d fields, expected %d."
                        % (line_no, len(fields), len(columns)))
                if not fields[taxon_idx].strip():
                    raise model.ValidationError(
                        "Line %d has an empty taxon." % line_no)
                complex_name = fields[complex_idx].strip()
                if complex_name not in COMPLEXES:
                    raise model.ValidationError(
                        "Line %d: invalid complex %r. Valid complexes: %s."
                        % (line_no, complex_name, ', '.join(COMPLEXES)))


SocranskyReferenceDirectoryFormat = model.SingleFileDirectoryFormat(
    'SocranskyReferenceDirectoryFormat',
    'reference.tsv',
    SocranskyReferenceFormat
)


class DysbiosisSignatureFormat(model.TextFileFormat):
    """Tab-separated list of signature taxa with a mandatory header that
    contains the columns ``taxon``, ``direction`` (dysbiotic or
    normobiotic) and ``level`` (species or genus)."""

    def _validate_(self, level):
        max_lines = {'min': 100, 'max': None}[level]
        with self.open() as fh:
            header = fh.readline().rstrip('\r\n')
            if not header:
                raise model.ValidationError("File is empty.")
            columns = [c.strip() for c in header.split('\t')]
            missing = [c for c in ('taxon', 'direction', 'level')
                       if c not in columns]
            if missing:
                raise model.ValidationError(
                    "Header must be tab-separated and contain the columns "
                    "'taxon', 'direction' and 'level'; missing: %s."
                    % ', '.join(missing))
            idx = {c: columns.index(c)
                   for c in ('taxon', 'direction', 'level')}

            for line_no, line in enumerate(fh, start=2):
                if max_lines is not None and line_no > max_lines + 1:
                    break
                line = line.rstrip('\r\n')
                if not line.strip():
                    continue
                fields = line.split('\t')
                if len(fields) != len(columns):
                    raise model.ValidationError(
                        "Line %d has %d fields, expected %d."
                        % (line_no, len(fields), len(columns)))
                if not fields[idx['taxon']].strip():
                    raise model.ValidationError(
                        "Line %d has an empty taxon." % line_no)
                direction = fields[idx['direction']].strip().lower()
                if direction not in DIRECTIONS:
                    raise model.ValidationError(
                        "Line %d: invalid direction %r. Valid: %s."
                        % (line_no, direction, ', '.join(DIRECTIONS)))
                tax_level = fields[idx['level']].strip().lower()
                if tax_level not in LEVELS:
                    raise model.ValidationError(
                        "Line %d: invalid level %r. Valid: %s."
                        % (line_no, tax_level, ', '.join(LEVELS)))


DysbiosisSignatureDirectoryFormat = model.SingleFileDirectoryFormat(
    'DysbiosisSignatureDirectoryFormat',
    'signature.tsv',
    DysbiosisSignatureFormat
)
