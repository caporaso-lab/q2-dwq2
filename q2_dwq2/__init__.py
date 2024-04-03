# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2024, Greg Caporaso.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions

from ._types_and_formats import (
    SingleDNASequence, SingleRecordDNAFASTAFormat,
    SingleRecordDNAFASTADirectoryFormat)

__version__ = get_versions()["version"]
del get_versions

from . import _version
__version__ = _version.get_versions()['version']

__all__ = [
    "SingleDNASequence", "SingleRecordDNAFASTAFormat",
    "SingleRecordDNAFASTADirectoryFormat"]