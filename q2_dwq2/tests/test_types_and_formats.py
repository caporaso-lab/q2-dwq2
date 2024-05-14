# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase

from q2_dwq2 import (
    SingleDNASequence, SingleRecordDNAFASTAFormat, LocalAlignmentSearchResults,
    LocalAlignmentSearchResultsFormat
)


class SingleDNASequenceTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_semantic_type_registration(self):
        self.assertRegisteredSemanticType(SingleDNASequence)


class LocalAlignmentSearchResultsTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_semantic_type_registration(self):
        self.assertRegisteredSemanticType(LocalAlignmentSearchResults)


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


class LocalAlignmentSearchResultsFormatTests(TestPluginBase):
    package = 'q2_dwq2.tests'

    def test_simple1(self):
        filenames = ['las-results-1.tsv']
        filepaths = [self.get_data_path(fn) for fn in filenames]

        for fp in filepaths:
            format = LocalAlignmentSearchResultsFormat(fp, mode='r')
            format.validate()

    def test_invalid_headers(self):
        fp = self.get_data_path('bad-las-results-header.tsv')
        format = LocalAlignmentSearchResultsFormat(fp, mode='r')
        self.assertRaisesRegex(ValidationError,
                               " id .* query id\\.",
                               format.validate)

    def test_missing_fields(self):
        fp = self.get_data_path('bad-las-results-missing-field.tsv')
        format = LocalAlignmentSearchResultsFormat(fp, mode='r')
        self.assertRaisesRegex(ValidationError,
                               "Line 4: Found 6 .* 7\\.",
                               format.validate)

    def test_validation_level(self):
        fp = self.get_data_path('las-results-validation-level-test.tsv')

        # min validation is successful
        format = LocalAlignmentSearchResultsFormat(fp, mode='r')
        format.validate(level='min')

        # but max validation raises an error
        format = LocalAlignmentSearchResultsFormat(fp, mode='r')
        self.assertRaisesRegex(ValidationError,
                               "Line 7: Found 6 .* 7\\.",
                               format.validate,
                               level='max')
