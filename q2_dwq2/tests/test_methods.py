# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt

from skbio.alignment import TabularMSA
from skbio.sequence import DNA

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform

from q2_dwq2._methods import (nw_align, local_alignment_search, split_sequences,
                              combine_las_reports)
from q2_dwq2._types_and_formats import SingleRecordDNAFASTAFormat


class NWAlignTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_simple1(self):
        # test alignment of a pair of sequences
        sequence1 = DNA('AAAAAAAAGGTGGCCTTTTTTTT')
        sequence2 = DNA('AAAAAAAAGGGGCCTTTTTTTT')
        observed = nw_align(sequence1, sequence2)

        aligned_sequence1 = DNA('AAAAAAAAGGTGGCCTTTTTTTT')
        aligned_sequence2 = DNA('AAAAAAAAGG-GGCCTTTTTTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)

    def test_simple2(self):
        # test alignment of a different pair of sequences
        # loaded from file this time, for demonstration purposes
        sequence1 = transform(
            self.get_data_path('seq-1.fasta'),
            from_type=SingleRecordDNAFASTAFormat,
            to_type=DNA)
        sequence2 = transform(
            self.get_data_path('seq-2.fasta'),
            from_type=SingleRecordDNAFASTAFormat,
            to_type=DNA)
        observed = nw_align(sequence1, sequence2)

        aligned_sequence1 = DNA('ACCGGTGGAACCGG-TAACACCCAC')
        aligned_sequence2 = DNA('ACCGGT--AACCGGTTAACACCCAC')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertNotEqual(observed, expected)

    def test_alt_match_score(self):
        sequence1 = DNA('AAAATTT')
        sequence2 = DNA('AAAAGGTTT')
        # call with default value for match score
        observed = nw_align(sequence1, sequence2)

        aligned_sequence1 = DNA('--AAAATTT')
        aligned_sequence2 = DNA('AAAAGGTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)

        # call with non-default value for match_score
        observed = nw_align(sequence1, sequence2, match_score=10)

        # the following expected outcome was determined by calling
        # skbio.alignment.global_pairwise_align_nucleotide directly. the
        # goal isn't to test that the underlying library code (i.e.,
        # skbio.alignment.global_pairwise_align_nucleotide) is working, b/c
        # I trust that that is already tested (or I wouldn't use it). rather,
        # the goal is to test that my wrapper of it is working. in this case,
        # specifically, i'm testing that passing an alternative value for
        # match_score changes the output alignment
        aligned_sequence1 = DNA('AAAA--TTT')
        aligned_sequence2 = DNA('AAAAGGTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)

    def test_alt_gap_open_penalty(self):
        sequence1 = DNA('AAAATTT')
        sequence2 = DNA('AAAAGGTTT')
        observed = nw_align(sequence1, sequence2, gap_open_penalty=0.01)

        aligned_sequence1 = DNA('AAAA-T-TT-')
        aligned_sequence2 = DNA('AAAAG-GTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)

        observed = nw_align(sequence1, sequence2)

        aligned_sequence1 = DNA('--AAAATTT')
        aligned_sequence2 = DNA('AAAAGGTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)

    def test_alt_gap_extend_penalty(self):
        sequence1 = DNA('AAAATTT')
        sequence2 = DNA('AAAAGGTTT')
        observed = nw_align(sequence1, sequence2, gap_open_penalty=0.01)

        aligned_sequence1 = DNA('AAAA-T-TT-')
        aligned_sequence2 = DNA('AAAAG-GTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)

        observed = nw_align(sequence1, sequence2, gap_open_penalty=0.01,
                            gap_extend_penalty=0.001)

        aligned_sequence1 = DNA('AAAA--TTT')
        aligned_sequence2 = DNA('AAAAGGTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)

    def test_alt_mismatch_score(self):
        sequence1 = DNA('AAAATTT')
        sequence2 = DNA('AAAAGGTTT')
        observed = nw_align(sequence1, sequence2, gap_open_penalty=0.01)

        aligned_sequence1 = DNA('AAAA-T-TT-')
        aligned_sequence2 = DNA('AAAAG-GTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)

        observed = nw_align(sequence1, sequence2, gap_open_penalty=0.1,
                            mismatch_score=-0.1)

        aligned_sequence1 = DNA('-AAA-ATTT')
        aligned_sequence2 = DNA('AAAAGGTTT')
        expected = TabularMSA([aligned_sequence1, aligned_sequence2])

        self.assertEqual(observed, expected)


class LocalAlignmentSearchTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_simple_single_query(self):
        query_sequence = [DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q1'})]
        reference_sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),  # 90% match
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),  # 100% match
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),  # 95% match
        ]
        observed = local_alignment_search(query_sequence, reference_sequences)
        expected = pd.DataFrame([
          ['q1', 'r2', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
           'ACACTCTCCACCCATTTGCT'],
          ['q1', 'r3', 95., 20, 35., 'ACACTCTCCACCCATTTGCT',
           'ACACTCTCCAGCCATTTGCT'],
          ['q1', 'r1', 90., 20, 30., 'ACACTCTCCACCCATTTGCT',
           'ACACTCACCACCCAATTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        expected.set_index(['query id', 'reference id'], inplace=True)
        pdt.assert_frame_equal(observed, expected)

    def test_simple_multi_query(self):
        query_sequences = [DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q1'}),
                           DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'q2'})]
        reference_sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),  # == q2
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),  # == q1
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),  # 95% match
        ]
        observed = local_alignment_search(query_sequences, reference_sequences)
        expected = pd.DataFrame([
          ['q1', 'r2', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
                                      'ACACTCTCCACCCATTTGCT'],
          ['q1', 'r3', 95., 20, 35., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCTCCAGCCATTTGCT'],
          ['q1', 'r1', 90., 20, 30., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCACCACCCAATTGCT'],
          ['q2', 'r1', 100., 20, 40., 'ACACTCACCACCCAATTGCT',
                                      'ACACTCACCACCCAATTGCT'],
          ['q2', 'r2', 90., 20, 30., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCACCCATTTGCT'],
          ['q2', 'r3', 85., 20, 25., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCAGCCATTTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        expected.set_index(['query id', 'reference id'], inplace=True)
        pdt.assert_frame_equal(observed, expected)

    def test_alt_n(self):
        query_sequences = [DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q1'})]
        reference_sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),  # 90% match
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),  # 100% match
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),  # 95% match
        ]
        observed = local_alignment_search(
            query_sequences, reference_sequences, n=1)
        expected = pd.DataFrame([
          ['q1', 'r2', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
           'ACACTCTCCACCCATTTGCT']
         ],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        expected.set_index(['query id', 'reference id'], inplace=True)
        pdt.assert_frame_equal(observed, expected)

        n = 3
        observed = local_alignment_search(
            query_sequences, reference_sequences, n=n)
        self.assertEqual(len(observed), n * len(query_sequences))

        query_sequences = [DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q1'}),
                           DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q2'}),
                           DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q3'})]
        observed = local_alignment_search(
            query_sequences, reference_sequences, n=n)
        self.assertEqual(len(observed), n * len(query_sequences))

    def test_retain_all_hits(self):
        query_sequences = [DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q1'})]
        reference_sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),  # 90% match
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),  # 100% match
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),  # 95% match
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r4'}),  # 95% match
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r5'}),  # 95% match
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r6'}),  # 95% match
        ]
        observed = local_alignment_search(
            query_sequences, reference_sequences, n=0)
        self.assertEqual(len(observed),
                         len(query_sequences) * len(reference_sequences))

    def test_empty_input_edge_cases(self):
        s = [DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q1'})]
        self.assertRaisesRegex(ValueError, "At least one",
                               local_alignment_search, s, [])

    def test_alt_alignment_parameters(self):
        query_sequence = [DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'q1'})]
        reference_sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),  # 90% match
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),  # 100% match
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),  # 95% match
            DNA('ACACTCTCCCCATTTGCT', metadata={'id': 'r4'}),  # two deletions
        ]
        observed_w_defaults = \
            local_alignment_search(query_sequence, reference_sequences)

        observed_w_alt_match_score = \
            local_alignment_search(query_sequence, reference_sequences,
                                   match_score=42)
        self.assertNotEqual(observed_w_defaults['score'][('q1', 'r1')],
                            observed_w_alt_match_score['score'][('q1', 'r1')])

        observed_w_alt_mismatch_score = \
            local_alignment_search(query_sequence, reference_sequences,
                                   mismatch_score=-42)
        self.assertNotEqual(
            observed_w_defaults['score'][('q1', 'r1')],
            observed_w_alt_mismatch_score['score'][('q1', 'r1')])

        observed_w_alt_gap_open_penalty = \
            local_alignment_search(query_sequence, reference_sequences,
                                   gap_open_penalty=42)
        self.assertNotEqual(
            observed_w_defaults['score'][('q1', 'r4')],
            observed_w_alt_gap_open_penalty['score'][('q1', 'r4')])

        observed_w_alt_gap_extend_penalty = \
            local_alignment_search(query_sequence, reference_sequences,
                                   gap_extend_penalty=42)
        self.assertNotEqual(
            observed_w_defaults['score'][('q1', 'r4')],
            observed_w_alt_gap_extend_penalty['score'][('q1', 'r4')])


class SplitSequencesTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def _flatten_splits(self, splits):
        result = []
        for k, dna_iter in splits.items():
            result.extend(list(dna_iter))
        return result

    def test_one_seq_per_split(self):
        sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r4'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r5'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r6'}),
        ]

        # confirm the expected number of splits are returned
        observed = split_sequences(sequences, split_size=1)
        self.assertEqual(len(observed), 6)
        # confirm that all input sequences are contained in the
        # split sequences
        self.assertEqual(self._flatten_splits(observed), sequences)

        sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),
        ]

        observed = split_sequences(sequences, split_size=1)
        self.assertEqual(len(observed), 2)
        self.assertEqual(self._flatten_splits(observed), sequences)

    def test_one_split(self):
        sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r4'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r5'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r6'}),
        ]

        observed = split_sequences(sequences, split_size=6)
        self.assertEqual(len(observed), 1)
        self.assertEqual(self._flatten_splits(observed), sequences)

        sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),
        ]

        observed = split_sequences(sequences, split_size=2)
        self.assertEqual(len(observed), 1)
        self.assertEqual(self._flatten_splits(observed), sequences)

    def test_multiple_seqs_multiple_splits(self):
        sequences = [
            DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),
            DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r4'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r5'}),
            DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r6'}),
        ]

        observed = split_sequences(sequences, split_size=2)
        self.assertEqual(len(observed), 3)
        self.assertEqual(self._flatten_splits(observed), sequences)

        observed = split_sequences(sequences, split_size=3)
        self.assertEqual(len(observed), 2)
        self.assertEqual(self._flatten_splits(observed), sequences)

        observed = split_sequences(sequences, split_size=4)
        self.assertEqual(len(observed), 2)
        self.assertEqual(self._flatten_splits(observed), sequences)

        observed = split_sequences(sequences, split_size=5)
        self.assertEqual(len(observed), 2)
        self.assertEqual(self._flatten_splits(observed), sequences)


class CombineLasReportsTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_combine_two_las_reports(self):
        expected = pd.DataFrame([
          ['q1', 'r2', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
                                      'ACACTCTCCACCCATTTGCT'],
          ['q1', 'r3', 95., 20, 35., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCTCCAGCCATTTGCT'],
          ['q1', 'r1', 90., 20, 30., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCACCACCCAATTGCT'],
          ['q2', 'r1', 100., 20, 40., 'ACACTCACCACCCAATTGCT',
                                      'ACACTCACCACCCAATTGCT'],
          ['q2', 'r2', 90., 20, 30., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCACCCATTTGCT'],
          ['q2', 'r3', 85., 20, 25., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCAGCCATTTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        expected.set_index(['query id', 'reference id'], inplace=True)

        part1 = pd.DataFrame([
          ['q1', 'r2', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
                                      'ACACTCTCCACCCATTTGCT'],
          ['q1', 'r3', 95., 20, 35., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCTCCAGCCATTTGCT'],
          ['q1', 'r1', 90., 20, 30., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCACCACCCAATTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        part1.set_index(['query id', 'reference id'], inplace=True)

        part2 = pd.DataFrame([
          ['q2', 'r1', 100., 20, 40., 'ACACTCACCACCCAATTGCT',
                                      'ACACTCACCACCCAATTGCT'],
          ['q2', 'r2', 90., 20, 30., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCACCCATTTGCT'],
          ['q2', 'r3', 85., 20, 25., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCAGCCATTTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        part2.set_index(['query id', 'reference id'], inplace=True)

        observed = combine_las_reports({0: part1, 1: part2})

        pdt.assert_frame_equal(observed, expected)

    def test_combine_three_las_reports(self):
        expected = pd.DataFrame([
          ['q1', 'r1', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
                                      'ACACTCTCCACCCATTTGCT'],
          ['q1', 'r2', 95., 20, 35., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCTCCAGCCATTTGCT'],
          ['q2', 'r3', 90., 20, 30., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCACCACCCAATTGCT'],
          ['q2', 'r4', 100., 20, 40., 'ACACTCACCACCCAATTGCT',
                                      'ACACTCACCACCCAATTGCT'],
          ['q3', 'r2', 90., 20, 30., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCACCCATTTGCT'],
          ['q3', 'r3', 85., 20, 25., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCAGCCATTTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        expected.set_index(['query id', 'reference id'], inplace=True)

        part1 = pd.DataFrame([
          ['q1', 'r1', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
                                      'ACACTCTCCACCCATTTGCT'],
          ['q1', 'r2', 95., 20, 35., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCTCCAGCCATTTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        part1.set_index(['query id', 'reference id'], inplace=True)

        part2 = pd.DataFrame([
          ['q2', 'r3', 90., 20, 30., 'ACACTCTCCACCCATTTGCT',
                                     'ACACTCACCACCCAATTGCT'],
          ['q2', 'r4', 100., 20, 40., 'ACACTCACCACCCAATTGCT',
                                      'ACACTCACCACCCAATTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        part2.set_index(['query id', 'reference id'], inplace=True)

        part3 = pd.DataFrame([
          ['q3', 'r2', 90., 20, 30., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCACCCATTTGCT'],
          ['q3', 'r3', 85., 20, 25., 'ACACTCACCACCCAATTGCT',
                                     'ACACTCTCCAGCCATTTGCT']],
         columns=['query id', 'reference id', 'percent similarity',
                  'alignment length', 'score', 'aligned query',
                  'aligned reference'])
        part3.set_index(['query id', 'reference id'], inplace=True)

        observed = combine_las_reports({0: part1, 1: part2, 2: part3})

        pdt.assert_frame_equal(observed, expected)
