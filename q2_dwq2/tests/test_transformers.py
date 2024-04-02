# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from skbio import DNA

from qiime2.plugin.testing import TestPluginBase

from q2_dwq2 import SingleRecordDNAFASTAFormat


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
