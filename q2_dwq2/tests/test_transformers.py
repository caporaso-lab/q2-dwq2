# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt
from skbio import DNA

from qiime2.plugin.testing import TestPluginBase

from q2_dwq2 import (SingleRecordDNAFASTAFormat,
                     LocalAlignmentSearchResultsFormat)


class SingleDNASequenceTransformerTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_single_record_fasta_to_DNA_simple1(self):
        _, observed = self.transform_format(
            SingleRecordDNAFASTAFormat, DNA, filename='seq-1.fasta')

        expected = DNA('ACCGGTGGAACCGGTAACACCCAC',
                       metadata={'id': 'example-sequence-1', 'description': ''})

        self.assertEqual(observed, expected)

    def test_single_record_fasta_to_DNA_simple2(self):
        _, observed = self.transform_format(
            SingleRecordDNAFASTAFormat, DNA, filename='seq-2.fasta')

        expected = DNA('ACCGGTAACCGGTTAACACCCAC',
                       metadata={'id': 'example-sequence-2', 'description': ''})

        self.assertEqual(observed, expected)

    def test_DNA_to_single_record_fasta_simple1(self):
        in_ = DNA('ACCGGTGGAACCGGTAACACCCAC',
                  metadata={'id': 'example-sequence-1', 'description': ''})
        tx = self.get_transformer(DNA, SingleRecordDNAFASTAFormat)

        observed = tx(in_)
        # confirm "round-trip" of DNA -> SingleRecordDNAFASTAFormat -> DNA
        # results in an observed sequence that is the same as the starting
        # sequence
        self.assertEqual(observed.view(DNA), in_)

    def test_DNA_to_single_record_fasta_simple2(self):
        in_ = DNA('ACCGGTAACCGGTTAACACCCAC',
                  metadata={'id': 'example-sequence-2', 'description': ''})
        tx = self.get_transformer(DNA, SingleRecordDNAFASTAFormat)

        observed = tx(in_)
        self.assertEqual(observed.view(DNA), in_)


class LocalAlignmentSearchResultsTransformerTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_las_results_format_to_df_simple1(self):
        _, observed = self.transform_format(
            LocalAlignmentSearchResultsFormat, pd.DataFrame,
            filename='las-results-1.tsv')
        expected = pd.DataFrame([
          ['que1', 'ref1', 100., 5, 20., 'AAAAA', 'AAAAA'],
          ['que1', 'ref2', 90., 10, 35., 'AAAAATTTTT', 'AAAAACTTTT'],
          ['que2', 'ref1', 95., 20, 70., 'AAAAATTTTTAAAAATTTTT',
           'AAAAATTTT-AAAAATTTTT']],
          columns=['query id', 'reference id', 'percent similarity',
                   'alignment length', 'score', 'aligned query',
                   'aligned reference'])
        expected.set_index(['query id', 'reference id'], inplace=True)

        pdt.assert_frame_equal(observed, expected)

    def test_las_results_format_to_df_simple2(self):
        _, observed = self.transform_format(
            LocalAlignmentSearchResultsFormat, pd.DataFrame,
            filename='las-results-2.tsv')
        expected = pd.DataFrame([
          ['que99', 'ref42', 100., 5, 20., 'AAAAA', 'AAAAA']],
          columns=['query id', 'reference id', 'percent similarity',
                   'alignment length', 'score', 'aligned query',
                   'aligned reference'])
        expected.set_index(['query id', 'reference id'], inplace=True)

        pdt.assert_frame_equal(observed, expected)

    def test_df_to_las_results_format_simple1(self):
        in_ = pd.DataFrame([
          ['q1', 'r2', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
           'ACACTCTCCACCCATTTGCT'],
          ['q1', 'r3', 95., 20, 35., 'ACACTCTCCACCCATTTGCT',
           'ACACTCTCCAGCCATTTGCT'],
          ['q1', 'r1', 90., 20, 30., 'ACACTCTCCACCCATTTGCT',
           'ACACTCACCACCCAATTGCT']],
          columns=['query id', 'reference id', 'percent similarity',
                   'alignment length', 'score', 'aligned query',
                   'aligned reference'])
        in_.set_index(['query id', 'reference id'], inplace=True)

        tx = self.get_transformer(pd.DataFrame,
                                  LocalAlignmentSearchResultsFormat)
        observed = tx(in_)
        # confirm "round-trip" of DataFrame ->
        # LocalAlignmentSearchResultsFormat -> DataFrame results in equal
        # input and output DataFrames
        pdt.assert_frame_equal(observed.view(pd.DataFrame), in_)

    def test_df_to_las_results_format_simple2(self):
        in_ = pd.DataFrame([
          ['q99', 'r42', 100., 20, 40., 'ACACTCTCCACCCATTTGCT',
           'ACACTCTCCACCCATTTGCT'],
          ['q1', 'r3', 95., 20, 35., 'ACACTCTCCACCCATTTGCT',
           'ACACTCTCCAGCCATTTGCT']],
          columns=['query id', 'reference id', 'percent similarity',
                   'alignment length', 'score', 'aligned query',
                   'aligned reference'])
        in_.set_index(['query id', 'reference id'], inplace=True)

        tx = self.get_transformer(pd.DataFrame,
                                  LocalAlignmentSearchResultsFormat)
        observed = tx(in_)
        pdt.assert_frame_equal(observed.view(pd.DataFrame), in_)
