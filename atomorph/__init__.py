#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Atomorph: A tool for converting atomic structure files.

Atomorph provides functionalities for converting atomic structure files between
different formats, with support for constraints, sorting, and multi-frame processing.
"""

__version__ = "1.0.2"

from .converter import StructureConverter
from .utils import setup_logging, get_version, list_available_formats

__all__ = [
    "StructureConverter",
    "setup_logging",
    "get_version",
    "list_available_formats",
    "__version__"
]