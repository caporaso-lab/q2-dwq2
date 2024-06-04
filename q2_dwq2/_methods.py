# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools

import pandas as pd

from skbio.alignment import (global_pairwise_align_nucleotide,
                             local_pairwise_align_nucleotide,
                             TabularMSA)
from skbio import DNA

from q2_types.feature_data import DNAIterator

_nw_align_defaults = {
    'gap_open_penalty': 5,
    'gap_extend_penalty': 2,
    'match_score': 1,
    'mismatch_score': -2
}

_las_defaults = {
    'n': 5,
    'gap_open_penalty': 5,
    'gap_extend_penalty': 2,
    'match_score': 2,
    'mismatch_score': -3
}

_split_seqs_defaults = {
    'split_size': 5
}


def nw_align(
        seq1: DNA,
        seq2: DNA,
        gap_open_penalty: float = _nw_align_defaults['gap_open_penalty'],
        gap_extend_penalty: float = _nw_align_defaults['gap_extend_penalty'],
        match_score: float = _nw_align_defaults['match_score'],
        mismatch_score: float = _nw_align_defaults['mismatch_score']) \
        -> TabularMSA:
    msa, _, _ = global_pairwise_align_nucleotide(
        seq1=seq1, seq2=seq2, gap_open_penalty=gap_open_penalty,
        gap_extend_penalty=gap_extend_penalty, match_score=match_score,
        mismatch_score=mismatch_score
    )

    return msa


def local_alignment_search(
        query_seqs: DNAIterator,
        reference_seqs: DNAIterator,
        n: int = _las_defaults['n'],
        gap_open_penalty: float = _las_defaults['gap_open_penalty'],
        gap_extend_penalty: float = _las_defaults['gap_extend_penalty'],
        match_score: float = _las_defaults['match_score'],
        mismatch_score: float = _las_defaults['mismatch_score']) \
        -> pd.DataFrame:
    # compute all of the alignments and their associated percent
    # similarities, lengths, and scores
    results = []
    for q, r in itertools.product(query_seqs, reference_seqs):
        aln, score, _ = local_pairwise_align_nucleotide(
            q, r,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
            match_score=match_score,
            mismatch_score=mismatch_score)
        percent_similarity = (100 * (1. - aln[0].distance(aln[1])))
        aln_length = aln.shape[1]
        results.append([q.metadata['id'],
                        r.metadata['id'],
                        percent_similarity,
                        aln_length,
                        score,
                        str(aln[0]),
                        str(aln[1])])

    if len(results) == 0:
        raise ValueError("At least one query sequence and one reference "
                         "sequence must be provided.")

    # package up the best hits and return them
    columns = ['query id',
               'reference id',
               'percent similarity',
               'alignment length',
               'score',
               'aligned query',
               'aligned reference']
    results = pd.DataFrame(results, columns=columns)

    # Filter results to the top n scores for each query; if n == 0,
    # retain all hits
    if n != 0:
        top_results = results.groupby('query id')['score']
        top_results = top_results.nlargest(n).reset_index(level=1, drop=True)
        top_results = pd.merge(results, top_results, how='inner',
                               on=['query id', 'score'])
    else:
        top_results = results

    # sort results by descending score
    top_results.sort_values(by=['query id', 'score'], ascending=[True, False],
                            inplace=True)
    top_results.set_index(['query id', 'reference id'], inplace=True)
    return top_results


# this is adapted from the itertools.batched documentation. when
# Python 3.12 is supported, this can be replaced with a call to
# itertools.batched
# https://docs.python.org/3/library/itertools.html#itertools.batched
def _batched(iterable, n):
    # _batched('ABCDEFG', 3) → ABC DEF G
    iterator = iter(iterable)
    while batch := tuple(itertools.islice(iterator, n)):
        yield batch


def split_sequences(
        seqs: DNAIterator,
        split_size: int = _split_seqs_defaults['split_size']) \
        -> DNAIterator:
    result = {i : DNAIterator(split)
              for i, split in enumerate(_batched(seqs, split_size))}
    return result


def combine_las_reports(reports: pd.DataFrame) -> pd.DataFrame:
    results = pd.concat(reports.values())
    return results
