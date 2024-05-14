# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from skbio import DNA

from q2_types.feature_data import DNAIterator
from q2_dwq2 import (SingleRecordDNAFASTAFormat,
                     LocalAlignmentSearchResultsFormat)

from .plugin_setup import plugin


# Define and register transformers
@plugin.register_transformer
def _1(ff: SingleRecordDNAFASTAFormat) -> DNA:
    # by default, DNA.read will read the first sequence in the file
    with ff.open() as fh:
        return DNA.read(fh)


@plugin.register_transformer
def _2(seq: DNA) -> SingleRecordDNAFASTAFormat:
    ff = SingleRecordDNAFASTAFormat()
    seq.write(str(ff.path))
    return ff


@plugin.register_transformer
def _3(ff: SingleRecordDNAFASTAFormat) -> DNAIterator:
    with ff.open() as fh:
        return DNAIterator([DNA.read(fh)])


@plugin.register_transformer
def _4(ff: LocalAlignmentSearchResultsFormat) -> pd.DataFrame:
    df = pd.read_csv(str(ff.path), sep='\t')
    df.set_index(['query id', 'reference id'], inplace=True)
    df = df.astype({'percent similarity': 'float',
                    'alignment length': 'int',
                    'score': 'float'})
    return df


@plugin.register_transformer
def _5(df: pd.DataFrame) -> LocalAlignmentSearchResultsFormat:
    ff = LocalAlignmentSearchResultsFormat()
    df.to_csv(str(ff.path), sep='\t')
    return ff
