# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from skbio import DNA

from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_dwq2._types_formats_transformers import (
    SingleDNASequence, SingleRecordDNAFASTAFormat
)


class SingleDNASequenceTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_semantic_type_registration(self):
        self.assertRegisteredSemanticType(SingleDNASequence)


class SingleRecordDNAFASTAFormatTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_simple1(self):
        filenames = ['seq-1.fasta', 'seq-2.fasta', 't-thermophilis-rrna.fasta']
        filepaths = [self.get_data_path(fn) for fn in filenames]

        for fp in filepaths:
            format = SingleRecordDNAFASTAFormat(fp, mode='r')
            format.validate()

    def test_invalid_default_validation(self):
        fp = self.get_data_path('bad-sequence-1.fasta')
        format = SingleRecordDNAFASTAFormat(fp, mode='r')
        self.assertRaisesRegex(ValidationError,
                               "4 non-ACGT characters.*171 positions.",
                               format.validate)

    def test_invalid_max_validation(self):
        fp = self.get_data_path('bad-sequence-1.fasta')
        format = SingleRecordDNAFASTAFormat(fp, mode='r')
        self.assertRaisesRegex(ValidationError,
                               "4 non-ACGT characters.*171 positions.",
                               format.validate,
                               level='max')

    def test_invalid_min_validation(self):
        fp = self.get_data_path('bad-sequence-1.fasta')
        format = SingleRecordDNAFASTAFormat(fp, mode='r')
        # min validation is successful
        format.validate(level='min')
        # but max validation raises an error
        self.assertRaisesRegex(ValidationError,
                               "4 non-ACGT characters.*171 positions.",
                               format.validate,
                               level='max')

        fp = self.get_data_path('bad-sequence-2.fasta')
        format = SingleRecordDNAFASTAFormat(fp, mode='r')
        self.assertRaisesRegex(ValidationError,
                               "4 non-ACGT characters.*50 positions.",
                               format.validate,
                               level='min')


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
