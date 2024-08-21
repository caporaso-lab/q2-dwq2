# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib.resources

import skbio

import qiime2
from q2_types.feature_data import DNAFASTAFormat


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


def align_and_summarize_example_1(use):
    seq1 = use.init_artifact('seq1', seq1_factory)
    seq2 = use.init_artifact('seq2', seq2_factory)

    msa, msa_summary, = use.action(
        use.UsageAction(plugin_id='dwq2',
                        action_id='align_and_summarize'),
        use.UsageInputs(seq1=seq1, seq2=seq2),
        use.UsageOutputNames(aligned_sequences='msa',
                             msa_summary='msa_summary'),
    )


def get_filepath_from_package(path):
    # this utility function may be transferred to the QIIME 2 framework
    # track progress here: https://github.com/qiime2/qiime2/issues/792
    return importlib.resources.files('q2_dwq2') / path


def query1_factory():
    fp = get_filepath_from_package(
        'examples/data/search-and-summarize/query.fasta')
    return qiime2.Artifact.import_data(
        "FeatureData[Sequence]", fp, view_type=DNAFASTAFormat
    )


def reference1_factory():
    fp = get_filepath_from_package(
        'examples/data/search-and-summarize/reference.fasta')
    return qiime2.Artifact.import_data(
        "FeatureData[Sequence]", fp, view_type=DNAFASTAFormat
    )


def search_and_summarize_example_1(use):
    query_seqs = use.init_artifact('query_seqs', query1_factory)
    reference_seqs = use.init_artifact('reference_seqs', reference1_factory)
    use.comment("This is an example of running this Pipeline serially.")
    use.comment("The modification to run this in parallel depends on the "
                "interface you're using (for example, using q2cli you would "
                "append the --parallel flag).")

    hits, hits_table, = use.action(
        use.UsageAction(plugin_id='dwq2',
                        action_id='search_and_summarize'),
        use.UsageInputs(query_seqs=query_seqs,
                        reference_seqs=reference_seqs,
                        split_size=1),
        use.UsageOutputNames(hits='hits', hits_table='hits-table')
    )
