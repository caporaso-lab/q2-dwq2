# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt

import skbio

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAIterator
from qiime2.sdk.parallel_config import ParallelConfig


class AlignAndSummarizeTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_simple1(self):
        # access the pipeline as QIIME 2 sees it,
        # for correct assignment of `ctx` variable
        align_and_summarize_pipeline = \
            self.plugin.pipelines['align_and_summarize']

        sequence1 = skbio.DNA('AAAAAAAAGGTGGCCTTTTTTTT',
                              metadata={'id': 's1', 'description': ''})
        sequence2 = skbio.DNA('AAAAAAAAGGGGCCTTTTTTTT',
                              metadata={'id': 's2', 'description': ''})
        sequence1_art = qiime2.Artifact.import_data(
            "SingleDNASequence", sequence1, view_type=skbio.DNA)
        sequence2_art = qiime2.Artifact.import_data(
            "SingleDNASequence", sequence2, view_type=skbio.DNA)
        observed_msa, observed_viz = align_and_summarize_pipeline(
            sequence1_art, sequence2_art)

        aligned_sequence1 = skbio.DNA('AAAAAAAAGGTGGCCTTTTTTTT',
                                      metadata={'id': 's1', 'description': ''})
        aligned_sequence2 = skbio.DNA('AAAAAAAAGG-GGCCTTTTTTTT',
                                      metadata={'id': 's2', 'description': ''})
        expected_msa = skbio.TabularMSA([aligned_sequence1, aligned_sequence2])

        # observed_msa output is a qiime2.Artifact, so view it as a
        # skbio.TabularMSA for comparison to expected_msa.
        self.assertEqual(observed_msa.view(skbio.TabularMSA), expected_msa)

        # observed_viz is a qiime2.Visualization.
        # access its index.html file for testing.
        index_fp = observed_viz.get_index_paths(relative=False)['html']
        with open(index_fp, 'r') as fh:
            observed_index = fh.read()
            self.assertIn(str(aligned_sequence1), observed_index)
            self.assertIn(str(aligned_sequence2), observed_index)


class SearchAndSummarizeTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def setUp(self):
        super().setUp()

        self.search_and_summarize_pipeline = \
            self.plugin.pipelines['search_and_summarize']
        query_sequences = [skbio.DNA('ACACTCTCCACCCATTTGCT',
                                     metadata={'id': 'q1'}),
                           skbio.DNA('ACACTCACCACCCAATTGCT',
                                     metadata={'id': 'q2'})]
        query_sequences = DNAIterator(query_sequences)
        self.query_sequences_art = qiime2.Artifact.import_data(
            "FeatureData[Sequence]", query_sequences, view_type=DNAIterator
        )
        reference_sequences = [
            skbio.DNA('ACACTCACCACCCAATTGCT', metadata={'id': 'r1'}),  # == q2
            skbio.DNA('ACACTCTCCACCCATTTGCT', metadata={'id': 'r2'}),  # == q1
            skbio.DNA('ACACTCTCCAGCCATTTGCT', metadata={'id': 'r3'}),
        ]
        reference_sequences = DNAIterator(reference_sequences)
        self.reference_sequences_art = qiime2.Artifact.import_data(
            "FeatureData[Sequence]", reference_sequences, view_type=DNAIterator
        )
        metadata_df = pd.DataFrame([
            ['r1', 'abc;def', 10101],
            ['r2', 'ghi;jkl;mno', 101],
            ['r3', 'pqr; stu', 110011],
            ['r4', 'vwxyz', 424242]],
            columns=['sample-id', 'some column', 'another'])
        metadata_df.set_index('sample-id', inplace=True)
        self.metadata = qiime2.Metadata(metadata_df)

    def _test_simple1_helper(self, observed_hits, observed_viz):
        expected_hits = pd.DataFrame([
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
        expected_hits.set_index(['query id', 'reference id'], inplace=True)

        pdt.assert_frame_equal(observed_hits.view(pd.DataFrame), expected_hits)

        # observed_viz is a qiime2.Visualization.
        # access its index.html file for testing.
        index_fp = observed_viz.get_index_paths(relative=False)['html']
        with open(index_fp, 'r') as fh:
            observed_index = fh.read()
            self.assertIn('q1', observed_index)
            self.assertIn('q2', observed_index)
            self.assertIn('r1', observed_index)
            self.assertIn('r2', observed_index)
            self.assertIn('ACACTCACCACCCAATTGCT', observed_index)
            self.assertIn('ACACTCTCCACCCATTTGCT', observed_index)

    def test_simple1_serial(self):

        observed_hits, observed_viz = self.search_and_summarize_pipeline(
            self.query_sequences_art, self.reference_sequences_art)
        self._test_simple1_helper(observed_hits, observed_viz)

    def test_metadata1_serial(self):

        observed_hits, observed_viz = self.search_and_summarize_pipeline(
            self.query_sequences_art, self.reference_sequences_art,
            reference_metadata=self.metadata)
        self._test_simple1_helper(observed_hits, observed_viz)

        # check that metadata is included in output Visualization
        index_fp = observed_viz.get_index_paths(relative=False)['html']
        with open(index_fp, 'r') as fh:
            observed_index = fh.read()
            self.assertIn('ghi;jkl;mno', observed_index)
            self.assertIn('some column', observed_index)
            self.assertIn('110011', observed_index)
            # metadata associated with ids not referenced in hits does
            # not show up in observed
            self.assertNotIn('424242', observed_index)

    def test_simple1_parallel(self):
        with ParallelConfig():
            observed_hits, observed_viz = \
                self.search_and_summarize_pipeline.parallel(
                    self.query_sequences_art, self.reference_sequences_art)
            self._test_simple1_helper(observed_hits, observed_viz)
