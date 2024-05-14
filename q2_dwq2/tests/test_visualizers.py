# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path

import pandas as pd
from skbio.alignment import TabularMSA
from skbio.sequence import DNA

from qiime2.plugin.testing import TestPluginBase

from q2_dwq2._visualizers import summarize_alignment, tabulate_las_results


class SummarizeAlignmentTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_simple1(self):
        aligned_sequence1 = DNA('AAAAAAAAGGTGGCCTTTTTTTT')
        aligned_sequence2 = DNA('AAAAAAAAGG-GGCCTTTTTTTT')
        msa = TabularMSA([aligned_sequence1, aligned_sequence2])

        with self.temp_dir as output_dir:
            summarize_alignment(output_dir, msa)

            with open(os.path.join(output_dir, 'index.html'), 'r') as fh:
                observed = fh.read()

            self.assertIn('AAAAAAAAGGTGGCCTTTTTTTT', observed)
            self.assertIn('AAAAAAAAGG-GGCCTTTTTTTT', observed)


class TabulateLasResultsTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_simple1(self):
        hits = pd.DataFrame([
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

        with self.temp_dir as output_dir:
            observed = tabulate_las_results(output_dir, hits)

            with open(os.path.join(output_dir, 'index.html'), 'r') as fh:
                observed = fh.read()

                self.assertIn('q1', observed)
                self.assertIn('q2', observed)
                self.assertIn('r1', observed)
                self.assertIn('r2', observed)
                self.assertIn('ACACTCACCACCCAATTGCT', observed)
                self.assertIn('ACACTCTCCACCCATTTGCT', observed)
