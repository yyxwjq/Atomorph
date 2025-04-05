#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Atomorph: A tool for converting atomic structure files.

Atomorph provides functionalities for converting atomic structure files between
different formats, with support for constraints, sorting, and multi-frame processing.
It also supports batch conversion of multiple files in a directory structure.
"""

__version__ = "1.0.3"

from .converter import StructureConverter
from .utils import setup_logging, get_version, list_available_formats

# Shortcut functions
def convert_file(*args, **kwargs):
    """Shorthand for creating a StructureConverter and calling convert()"""
    converter = StructureConverter()
    return converter.convert(*args, **kwargs)

def batch_convert(*args, **kwargs):
    """
    Shorthand for creating a StructureConverter and calling batch_convert()
    
    This function provides a convenient way to batch convert structure files
    without explicitly creating a StructureConverter instance.
    
    Args:
        input_dir: Input directory path
        output_dir: Output directory path
        pattern: File pattern to match (default: "POSCAR")
        recursive: Whether to search recursively (default: True)
        output_format: Output file format (default: "vasp")
        merge_output: Whether to merge all structures into a single multi-frame file (default: False)
        output_merged_file: Output path for merged file (used only if merge_output is True)
        **kwargs: Additional conversion parameters
        
    Returns:
        Dictionary with statistics about the conversion process
        
    Example:
        >>> from atomorph import batch_convert
        >>> stats = batch_convert(
        ...     input_dir="./calculations",
        ...     output_dir="./converted",
        ...     pattern="POSCAR",
        ...     output_format="xyz"
        ... )
    """
    converter = StructureConverter()
    return converter.batch_convert(*args, **kwargs)

__all__ = [
    "StructureConverter",
    "setup_logging",
    "get_version",
    "list_available_formats",
    "convert_file",
    "batch_convert",
    "__version__"
]