#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Utility functions for the Atomorph package.
"""

import os
import re
import logging
from typing import List, Dict, Union, Optional
from pathlib import Path

# Configure logging
logger = logging.getLogger("atomorph")

def setup_logging(level: int = logging.INFO, log_file: Optional[str] = None) -> None:
    """
    Set up logging configuration.
    
    Args:
        level: Logging level
        log_file: Path to log file
    """
    logger.setLevel(level)
    
    # Create handlers
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    
    # Create formatters
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    console_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(console_handler)
    
    # Add file handler if log_file is specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

def ensure_dir(path: Union[str, Path]) -> Path:
    """
    Ensure directory exists, creating it if necessary.
    
    Args:
        path: Directory path
        
    Returns:
        Path object of the directory
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path

def get_version() -> str:
    """
    Get package version from __init__.py.
    
    Returns:
        Version string
    """
    package_init = Path(__file__).parent / "__init__.py"
    with open(package_init, 'r') as f:
        content = f.read()
    
    version_match = re.search(r"__version__ = ['\"]([^'\"]+)['\"]", content)
    if version_match:
        return version_match.group(1)
    return "unknown"

def list_available_formats() -> Dict[str, List[str]]:
    """
    List available file formats for conversion.
    
    Returns:
        Dictionary of input and output formats
    """
    from ase.io.formats import all_formats
    
    input_formats = []
    output_formats = []
    
    for name, desc in all_formats.items():
        if desc.can_read:
            input_formats.append(name)
        if desc.can_write:
            output_formats.append(name)
    
    # Add our custom mappings
    custom_formats = {
        "xyz": "extxyz",
        "vasp": "vasp",
        "poscar": "vasp"
    }
    
    for name in custom_formats:
        if name not in input_formats:
            input_formats.append(name)
        if name not in output_formats:
            output_formats.append(name)
    
    return {
        "input_formats": sorted(input_formats),
        "output_formats": sorted(output_formats)
    }

def validate_constraints(constraints: Union[str, List, None]) -> bool:
    """
    Validate constraint format.
    
    Args:
        constraints: Constraints specification
        
    Returns:
        True if constraints are valid
    """
    if constraints is None:
        return True
        
    if constraints == "fixed":
        return True
        
    if not isinstance(constraints, list):
        return False
        
    if len(constraints) < 2:
        return False
        
    constraint_type = constraints[0]
    if constraint_type not in ["elements", "layers", "indices"]:
        return False
        
    return True

def format_file_size(size_bytes: int) -> str:
    """
    Format file size in human-readable format.
    
    Args:
        size_bytes: File size in bytes
        
    Returns:
        Formatted file size string
    """
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024 * 1024:
        return f"{size_bytes / 1024:.2f} KB"
    elif size_bytes < 1024 * 1024 * 1024:
        return f"{size_bytes / (1024 * 1024):.2f} MB"
    else:
        return f"{size_bytes / (1024 * 1024 * 1024):.2f} GB" 