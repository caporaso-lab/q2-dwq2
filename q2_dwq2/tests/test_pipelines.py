# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio

import qiime2
from qiime2.plugin.testing import TestPluginBase


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
