#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Command-line interface for Atomorph structure converter.
"""

import os
import sys
import argparse
from typing import List, Optional, Union, Any
from pathlib import Path

from .converter import StructureConverter

def parse_constraints(constraint_str: str) -> List[str]:
    """
    Parse constraint string into a list of constraints.
    
    Args:
        constraint_str: Constraint string (e.g., "fixed", "elements:C,H", "layers:0.5,0.8")
        
    Returns:
        List of constraints
    """
    if not constraint_str:
        return None
        
    if constraint_str == "fixed":
        return "fixed"
        
    parts = constraint_str.split(":")
    if len(parts) != 2:
        raise ValueError(f"Invalid constraint format: {constraint_str}")
        
    constraint_type = parts[0].strip().lower()
    constraint_values = parts[1].strip()
    
    if constraint_type == "elements":
        elements = [e.strip() for e in constraint_values.split(",")]
        return ["elements"] + elements
    elif constraint_type == "layers":
        return ["layers", constraint_values]
    elif constraint_type == "indices":
        indices = []
        for idx_part in constraint_values.split(","):
            idx_part = idx_part.strip()
            if "-" in idx_part:
                start, end = map(int, idx_part.split("-"))
                indices.append(f"{start}-{end}")
            else:
                indices.append(idx_part)
        return ["indices"] + indices
    else:
        raise ValueError(f"Unknown constraint type: {constraint_type}")

def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Args:
        args: Command-line arguments
        
    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Convert atomic structure files between different formats.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/output
    parser.add_argument("input", help="Input file path")
    parser.add_argument("output", help="Output file path")
    parser.add_argument("-if", "--input-format", help="Input file format (auto-detect if not specified)")
    parser.add_argument("-of", "--output-format", help="Output file format (auto-detect if not specified)")
    
    # Mode options
    parser.add_argument("--mode", choices=["single", "multi"], 
                       help="Conversion mode (single frame or multi-frame)")
    
    # Frame selection
    parser.add_argument("--frame", help="Frame selection (e.g., '1 2 3-5')")
    
    # Element order
    parser.add_argument("--element-order", help="Custom element order (comma-separated)")
    
    # Constraints
    parser.add_argument("--constraints", help="""Atomic constraints:
                         fixed (fix all atoms),
                         elements:C,H (fix specific elements),
                         layers:0.3,0.5 (fix atoms in specific fractional z-range),
                         indices:1,2,5-10 (fix specific atoms by index)""")
    
    # Sorting
    parser.add_argument("--sort", choices=["ascending", "descending"], 
                       help="Sort atoms by chemical symbol")
    
    # Parallel processing
    parser.add_argument("--parallel", action="store_true",
                       help="Use parallel processing for multi-frame conversion")
    
    # Multiple frames
    parser.add_argument("--multi-frame", action="store_true",
                       help="Process multiple frames")
    
    # Separate directories
    parser.add_argument("--separate-dirs", action="store_true",
                       help="Save frames in separate directories")
    
    # Use fractional coordinates for layer constraints
    parser.add_argument("--use-fractional", action="store_true",
                       help="Use fractional coordinates for layer constraints")
    
    # Debug mode
    parser.add_argument("--debug", action="store_true",
                       help="Enable debug mode")
    
    return parser.parse_args(args)

def main(args: Optional[List[str]] = None) -> int:
    """
    Main function.
    
    Args:
        args: Command-line arguments
        
    Returns:
        Exit code
    """
    try:
        # Parse arguments
        args = parse_args(args)
        
        # Create converter
        converter = StructureConverter()
        
        # Parse element order
        element_mapping = None
        if args.element_order:
            element_mapping = [e.strip() for e in args.element_order.split(",")]
        
        # Parse constraints
        constraints = None
        if args.constraints:
            constraints = parse_constraints(args.constraints)
        
        # Print debug info
        if args.debug:
            print(f"Input file: {args.input}")
            print(f"Output file: {args.output}")
            print(f"Input format: {args.input_format}")
            print(f"Output format: {args.output_format}")
            print(f"Mode: {args.mode}")
            print(f"Frame selection: {args.frame}")
            print(f"Element order: {element_mapping}")
            print(f"Constraints: {constraints}")
            print(f"Sort type: {args.sort}")
            print(f"Use fractional: {args.use_fractional}")
        
        # Convert
        converter.convert(
            input_path=args.input,
            output_path=args.output,
            input_format=args.input_format,
            output_format=args.output_format,
            mode=args.mode,
            frame=args.frame,
            element_mapping=element_mapping,
            constraints=constraints,
            sort_type=args.sort,
            parallel=args.parallel,
            multi_frame=args.multi_frame,
            separate_dirs=args.separate_dirs,
            use_fractional=args.use_fractional
        )
        
        return 0
        
    except ValueError as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        return 1
    except Exception as e:
        if args and args.debug:
            import traceback
            traceback.print_exc()
        print(f"Unexpected error: {str(e)}", file=sys.stderr)
        return 2

if __name__ == "__main__":
    sys.exit(main()) 