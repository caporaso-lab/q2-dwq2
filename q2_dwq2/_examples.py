# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile

import skbio

import qiime2

from q2_dwq2 import SingleRecordDNAFASTAFormat


def seq1_factory():
    seq = skbio.DNA("AACCGGTTGGCCAA", metadata={"id": "seq1"})
    return _create_seq_artifact(seq)


def seq2_factory():
    seq = skbio.DNA("AACCGCTGGCGAA", metadata={"id": "seq2"})
    return _create_seq_artifact(seq)


def _create_seq_artifact(seq: skbio.DNA):
    with tempfile.NamedTemporaryFile() as f:
        # write our skbio.DNA object to file in fasta format
        seq.write(f)
        # reset to the beginning of the file
        f.seek(0)
        # instantiate our file format (ff) object with the fasta file
        ff = SingleRecordDNAFASTAFormat(f.name, mode='r')
        # return the sequence packaged in a "SingleDNASequence" qiime2.Artifact
        return qiime2.Artifact.import_data("SingleDNASequence", ff)


def nw_align_example_1(use):
    seq1 = use.init_artifact('seq1', seq1_factory)
    seq2 = use.init_artifact('seq2', seq2_factory)

    msa, = use.action(
        use.UsageAction(plugin_id='dwq2',
                        action_id='nw_align'),
        use.UsageInputs(seq1=seq1, seq2=seq2),
        use.UsageOutputNames(aligned_sequences='msa'),
    )
