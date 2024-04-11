# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio

import qiime2


def seq1_factory():
    seq = skbio.DNA("AACCGGTTGGCCAA", metadata={"id": "seq1"})
    return qiime2.Artifact.import_data(
        "SingleDNASequence", seq, view_type=skbio.DNA)


def seq2_factory():
    seq = skbio.DNA("AACCGCTGGCGAA", metadata={"id": "seq2"})
    return qiime2.Artifact.import_data(
        "SingleDNASequence", seq, view_type=skbio.DNA)


def nw_align_example_1(use):
    seq1 = use.init_artifact('seq1', seq1_factory)
    seq2 = use.init_artifact('seq2', seq2_factory)

    msa, = use.action(
        use.UsageAction(plugin_id='dwq2',
                        action_id='nw_align'),
        use.UsageInputs(seq1=seq1, seq2=seq2),
        use.UsageOutputNames(aligned_sequences='msa'),
    )
