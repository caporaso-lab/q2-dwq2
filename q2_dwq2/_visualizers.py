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

import qiime2


_tabulate_las_defaults = {'title': 'Local Alignment Search Results'}


def summarize_alignment(output_dir: str, msa: TabularMSA) -> None:
    title = "Alignment summary"
    with open(os.path.join(output_dir, "index.html"), "w") as fh:
        fh.write(_html_template % (title, repr(msa)))


def tabulate_las_results(output_dir: str,
                         hits: pd.DataFrame,
                         title: str = _tabulate_las_defaults['title'],
                         reference_metadata: qiime2.Metadata = None) \
                         -> None:
    if reference_metadata is not None:
        reference_metadata = reference_metadata.to_dataframe()

        hits.reset_index(inplace=True)

        metadata_index = reference_metadata.index.name
        metadata_columns = reference_metadata.columns.to_list()
        reference_metadata.reset_index(inplace=True)

        missing_ids = \
            set(hits['reference id']) - set(reference_metadata[metadata_index])
        if len(missing_ids) > 0:
            raise KeyError(
                f"The following {len(missing_ids)} IDs are missing from "
                f"reference metadata: {missing_ids}.")

        hits = pd.merge(hits, reference_metadata,
                        left_on='reference id',
                        right_on=metadata_index,
                        how='inner')

        hits.set_index(['query id', 'reference id'], inplace=True)
        column_order = \
            ['percent similarity', 'alignment length', 'score'] + \
            metadata_columns + ['aligned query', 'aligned reference']
        hits = hits.reindex(columns=column_order)
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
