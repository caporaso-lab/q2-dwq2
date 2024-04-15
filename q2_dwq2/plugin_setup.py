# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import Citations, Plugin, Float, Range, Visualization
from q2_types.feature_data import FeatureData, AlignedSequence
from q2_dwq2 import __version__
from q2_dwq2._methods import nw_align
from q2_dwq2._visualizers import summarize_alignment
from q2_dwq2._pipelines import align_and_summarize
from q2_dwq2._examples import nw_align_example_1, align_and_summarize_example_1
from q2_dwq2 import (
    SingleDNASequence, SingleRecordDNAFASTAFormat,
    SingleRecordDNAFASTADirectoryFormat)

citations = Citations.load("citations.bib", package="q2_dwq2")

plugin = Plugin(
    name="dwq2",
    version=__version__,
    website="https://cap-lab.bio",
    package="q2_dwq2",
    description="Plugin for use with *Developing with QIIME 2* (DWQ2).",
    short_description="Educational plugin.",
    citations=[citations['Caporaso-Bolyen-2024']]
)

# Register semantic types
plugin.register_semantic_types(SingleDNASequence)

# Register formats
plugin.register_formats(SingleRecordDNAFASTAFormat,
                        SingleRecordDNAFASTADirectoryFormat)

# Define and register new ArtifactClass
plugin.register_artifact_class(SingleDNASequence,
                               SingleRecordDNAFASTADirectoryFormat,
                               description="A single DNA sequence.")


# Register methods
_nw_align_inputs = {'seq1': SingleDNASequence,
                    'seq2': SingleDNASequence}
_nw_align_parameters = {
        'gap_open_penalty': Float % Range(0, None, inclusive_start=False),
        'gap_extend_penalty': Float % Range(0, None, inclusive_start=False),
        'match_score': Float % Range(0, None, inclusive_start=False),
        'mismatch_score': Float % Range(None, 0, inclusive_end=True)}
_nw_align_outputs = {'aligned_sequences': FeatureData[AlignedSequence]}
_nw_align_input_descriptions = {'seq1': 'The first sequence to align.',
                                'seq2': 'The second sequence to align.'}
_nw_align_parameter_descriptions = {
        'gap_open_penalty': ('The penalty incurred for opening a new gap. By '
                             'convention this is a positive number.'),
        'gap_extend_penalty': ('The penalty incurred for extending an existing '
                               'gap. By convention this is a positive number.'),
        'match_score': ('The score for matching characters at an alignment '
                        'position. By convention, this is a positive number.'),
        'mismatch_score': ('The score for mismatching characters at an '
                           'alignment position. By convention, this is a '
                           'negative number.')}
_nw_align_output_descriptions = {
        'aligned_sequences': 'The pairwise aligned sequences.'
    }

plugin.methods.register_function(
    function=nw_align,
    inputs=_nw_align_inputs,
    parameters=_nw_align_parameters,
    outputs=_nw_align_outputs,
    input_descriptions=_nw_align_input_descriptions,
    parameter_descriptions=_nw_align_parameter_descriptions,
    output_descriptions=_nw_align_output_descriptions,
    name='Pairwise global sequence alignment.',
    description=("Align two DNA sequences using Needleman-Wunsch (NW). "
                 "This is a Python implementation of NW, so it is very slow! "
                 "This action is for demonstration purposes only. 🐌"),
    citations=[citations['Needleman1970']],
    examples={'Align two DNA sequences.': nw_align_example_1}
)

# Register visualizers
plugin.visualizers.register_function(
    function=summarize_alignment,
    inputs={'msa': FeatureData[AlignedSequence]},
    input_descriptions={'msa': 'The multiple sequence alignment to summarize.'},
    parameters={},
    parameter_descriptions={},
    name='Summarize an alignment.',
    description='Summarize a multiple sequence alignment.',
    citations=[],
)

# Register pipelines
# Order is important in the outputs dict unlike for inputs and parameters,
# because these outputs are mapped onto the return values of the registered
# function.
_align_and_summarize_outputs = {}
_align_and_summarize_outputs.update(_nw_align_outputs)
_align_and_summarize_outputs['msa_summary'] = Visualization

_align_and_summarize_output_descriptions = {}
_align_and_summarize_output_descriptions.update(_nw_align_output_descriptions)
_align_and_summarize_output_descriptions['msa_summary'] = \
    "Visual summary of the pairwise alignment."

plugin.pipelines.register_function(
    function=align_and_summarize,
    inputs=_nw_align_inputs,
    parameters=_nw_align_parameters,
    outputs=_align_and_summarize_outputs,
    input_descriptions=_nw_align_input_descriptions,
    parameter_descriptions=_nw_align_parameter_descriptions,
    output_descriptions=_align_and_summarize_output_descriptions,
    name="Pairwise global alignment and summarization.",
    description=("Perform global pairwise sequence alignment using a slow "
                 "Needleman-Wunsch (NW) implementation, and generate a "
                 "visual summary of the alignment."),
    # Only citations new to this Pipeline need to be defined. Citations for
    # the Actions called from the Pipeline are automatically included.
    citations=[],
    examples={'Align two sequences and summarize the alignment.':
              align_and_summarize_example_1}
)

importlib.import_module('q2_dwq2._transformers')
