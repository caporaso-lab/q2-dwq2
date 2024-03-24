# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Citations, Plugin
from q2_types.feature_table import FeatureTable, Frequency
from q2_dwq2 import __version__
from q2_dwq2._methods import duplicate_table

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

plugin.methods.register_function(
    function=duplicate_table,
    inputs={'table': FeatureTable[Frequency]},
    parameters={},
    outputs=[('new_table', FeatureTable[Frequency])],
    input_descriptions={'table': 'The feature table to be duplicated.'},
    parameter_descriptions={},
    output_descriptions={'new_table': 'The duplicated feature table.'},
    name='Duplicate table',
    description=("Create a copy of a feature table with a new uuid. "
                 "This is for demonstration purposes only. 🧐"),
    citations=[]
)
