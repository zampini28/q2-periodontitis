import pandas as pd
import qiime2

from .plugin_setup import plugin
from ._format import (
    SocranskyComplexesFormat,
    PeriodontalIndicesFormat,
    PeriodontalIndicesDeltaFormat,
    SocranskyMappingFormat,
    SocranskyReferenceFormat,
    DysbiosisSignatureFormat,
)


@plugin.register_transformer
def _1(format: SocranskyComplexesFormat) -> pd.DataFrame:
    with format.open() as fh:
        return pd.read_csv(fh, sep='\t', index_col=0)


@plugin.register_transformer
def _2(df: pd.DataFrame) -> SocranskyComplexesFormat:
    format = SocranskyComplexesFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index_label='sample-id')
    return format


@plugin.register_transformer
def _3(format: PeriodontalIndicesFormat) -> pd.DataFrame:
    with format.open() as fh:
        return pd.read_csv(fh, sep='\t', index_col=0)


@plugin.register_transformer
def _4(df: pd.DataFrame) -> PeriodontalIndicesFormat:
    format = PeriodontalIndicesFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index_label='sample-id')
    return format


@plugin.register_transformer
def _5(format: SocranskyMappingFormat) -> pd.DataFrame:
    with format.open() as fh:
        return pd.read_csv(fh, sep='\t', index_col=0)


@plugin.register_transformer
def _6(df: pd.DataFrame) -> SocranskyMappingFormat:
    format = SocranskyMappingFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index_label='Feature ID')
    return format


@plugin.register_transformer
def _8(format: PeriodontalIndicesFormat) -> qiime2.Metadata:
    with format.open() as fh:
        return qiime2.Metadata(pd.read_csv(fh, sep='\t', index_col=0))


@plugin.register_transformer
def _9(format: SocranskyComplexesFormat) -> qiime2.Metadata:
    with format.open() as fh:
        return qiime2.Metadata(pd.read_csv(fh, sep='\t', index_col=0))


@plugin.register_transformer
def _10(format: SocranskyMappingFormat) -> qiime2.Metadata:
    with format.open() as fh:
        return qiime2.Metadata(pd.read_csv(fh, sep='\t', index_col=0))


@plugin.register_transformer
def _11(format: SocranskyReferenceFormat) -> pd.DataFrame:
    with format.open() as fh:
        df = pd.read_csv(fh, sep='\t', dtype=str, skip_blank_lines=True)
    df.columns = [c.strip() for c in df.columns]
    for col in ('taxon', 'complex'):
        df[col] = df[col].str.strip()
    return df


@plugin.register_transformer
def _12(df: pd.DataFrame) -> SocranskyReferenceFormat:
    format = SocranskyReferenceFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index=False)
    return format


@plugin.register_transformer
def _13(format: PeriodontalIndicesDeltaFormat) -> pd.DataFrame:
    with format.open() as fh:
        df = pd.read_csv(fh, sep='\t', index_col=0, dtype={0: str})
    df.index = df.index.astype(str)
    df.index.name = 'subject-id'
    return df


@plugin.register_transformer
def _14(df: pd.DataFrame) -> PeriodontalIndicesDeltaFormat:
    format = PeriodontalIndicesDeltaFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index_label='subject-id')
    return format


@plugin.register_transformer
def _15(format: PeriodontalIndicesDeltaFormat) -> qiime2.Metadata:
    with format.open() as fh:
        df = pd.read_csv(fh, sep='\t', index_col=0, dtype={0: str})
    # QIIME 2 Metadata only accepts its own ID header names
    df.index.name = 'id'
    return qiime2.Metadata(df)


@plugin.register_transformer
def _16(format: DysbiosisSignatureFormat) -> pd.DataFrame:
    with format.open() as fh:
        df = pd.read_csv(fh, sep='\t', dtype=str, skip_blank_lines=True)
    df.columns = [c.strip() for c in df.columns]
    for col in ('taxon', 'direction', 'level'):
        df[col] = df[col].str.strip()
    return df


@plugin.register_transformer
def _17(df: pd.DataFrame) -> DysbiosisSignatureFormat:
    format = DysbiosisSignatureFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index=False)
    return format
