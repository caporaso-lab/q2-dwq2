# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._methods import _nw_align_defaults, _las_defaults, _split_seqs_defaults
from ._visualizers import _tabulate_las_defaults


def align_and_summarize(
        ctx, seq1, seq2,
        gap_open_penalty=_nw_align_defaults['gap_open_penalty'],
        gap_extend_penalty=_nw_align_defaults['gap_extend_penalty'],
        match_score=_nw_align_defaults['match_score'],
        mismatch_score=_nw_align_defaults['mismatch_score']):
    nw_align_action = ctx.get_action('dwq2', 'nw_align')
    summarize_alignment_action = ctx.get_action('dwq2', 'summarize_alignment')

    msa, = nw_align_action(
                    seq1, seq2, gap_open_penalty=gap_open_penalty,
                    gap_extend_penalty=gap_extend_penalty,
                    match_score=match_score, mismatch_score=mismatch_score)
    msa_summary, = summarize_alignment_action(msa)

    return (msa, msa_summary)


def search_and_summarize(
        ctx, query_seqs, reference_seqs,
        reference_metadata=None,
        n=_las_defaults['n'],
        split_size=_split_seqs_defaults['split_size'],
        title=_tabulate_las_defaults['title'],
        gap_open_penalty=_las_defaults['gap_open_penalty'],
        gap_extend_penalty=_las_defaults['gap_extend_penalty'],
        match_score=_las_defaults['match_score'],
        mismatch_score=_las_defaults['mismatch_score']):
    split_action = ctx.get_action('dwq2', 'split_sequences')
    combine_action = ctx.get_action('dwq2', 'combine_las_reports')
    las_action = ctx.get_action('dwq2', 'local_alignment_search')
    tabulate_las_results_action = ctx.get_action('dwq2', 'tabulate_las_results')

    query_splits, = split_action(query_seqs, split_size=split_size)

    las_results = []
    for q in query_splits.values():
        las_result, = las_action(
            q, reference_seqs, n=n, gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty, match_score=match_score,
            mismatch_score=mismatch_score)
        las_results.append(las_result)

    las_results, = combine_action(las_results)
    result_table, = tabulate_las_results_action(
        las_results, title=title, reference_metadata=reference_metadata)

    return (las_results, result_table)
