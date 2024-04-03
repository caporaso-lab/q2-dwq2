# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.alignment import TabularMSA
from skbio.sequence import DNA

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform

from q2_dwq2._methods import nw_align
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
