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


def summarize_alignment(output_dir: str, msa: TabularMSA) -> None:
    title = "Alignment summary"
    with open(os.path.join(output_dir, "index.html"), "w") as fh:
        fh.write(_html_template % (title, repr(msa)))


def tabulate_las_results(output_dir: str,
                         hits: pd.DataFrame,
                         title: str = "Local Alignment Search Results") \
                         -> None:
    with open(os.path.join(output_dir, "index.html"), "w") as fh:
        fh.write(_html_template % (title, hits.to_html()))


_html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>%s</title>
    <style>
        body {
            padding: 20px;
        }
        p.alignment {
            font-family: 'Courier New', Courier, monospace;
        }
    </style>
</head>
<body>
    <pre>
%s
    </pre>
</body>
</html>
"""
