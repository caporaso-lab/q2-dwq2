# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

description = ("A template QIIME 2 plugin.")

setup(
    name="q2-dwq2",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Greg Caporaso",
    author_email="greg.caporaso@nau.edu",
    description=description,
    url="https://cap-lab.bio",
    entry_points={
        "qiime2.plugins": ["q2-dwq2=q2_dwq2.plugin_setup:plugin"]
    },
    package_data={
        "q2_dwq2": ["citations.bib"],
        "q2_dwq2.tests": ["data/*", "data/search-and-summarize/*"],
    },
    zip_safe=False,
)
