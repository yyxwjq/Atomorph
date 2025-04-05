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
    
    # Debug mode
    parser.add_argument("--debug", action="store_true",
                     help="Enable debug mode")
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # Single file conversion command
    convert_parser = subparsers.add_parser("convert", aliases=["conv", "c"], help="Convert a single file (short alias: c)")
    
    # Input/output
    convert_parser.add_argument("input", help="Input file path")
    convert_parser.add_argument("output", help="Output file path")
    convert_parser.add_argument("-if", "--input-format", help="Input file format (auto-detect if not specified)")
    convert_parser.add_argument("-of", "--output-format", help="Output file format (auto-detect if not specified)")
    
    # Mode options
    convert_parser.add_argument("-m", "--mode", choices=["single", "multi"], 
                       help="Conversion mode (single frame or multi-frame)")
    
    # Frame selection
    convert_parser.add_argument("-f", "--frame", help="Frame selection (e.g., '1 2 3-5')")
    
    # Element order
    convert_parser.add_argument("-e", "--element-order", help="Custom element order (comma-separated)")
    
    # Constraints
    convert_parser.add_argument("-c", "--constraints", help="""Atomic constraints:
                         fixed (fix all atoms),
                         elements:C,H (fix specific elements),
                         layers:0.3,0.5 (fix atoms in specific fractional z-range),
                         indices:1,2,5-10 (fix specific atoms by index)""")
    
    # Sorting
    convert_parser.add_argument("-s", "--sort", choices=["ascending", "descending", "asc", "desc"], 
                       help="Sort atoms by chemical symbol (shorthand: asc, desc)")
    
    # Parallel processing
    convert_parser.add_argument("-p", "--parallel", action="store_true",
                       help="Use parallel processing for multi-frame conversion")
    
    # Multiple frames
    convert_parser.add_argument("--multi-frame", action="store_true",
                       help="Process multiple frames")
    
    # Separate directories
    convert_parser.add_argument("--separate-dirs", action="store_true",
                       help="Save frames in separate directories")
    
    # Use fractional coordinates for layer constraints
    convert_parser.add_argument("--use-fractional", action="store_true",
                       help="Use fractional coordinates for layer constraints")
    
    # Batch conversion command
    batch_parser = subparsers.add_parser("batch", aliases=["b"], help="Batch convert files in a directory (short alias: b)")
    
    # Input/output directories
    batch_parser.add_argument("input_dir", help="Input directory path")
    batch_parser.add_argument("output_dir", help="Output directory path")
    
    # File pattern
    batch_parser.add_argument("-p", "--pattern", default="POSCAR", 
                        help="File pattern to match (default: 'POSCAR', or use '*.ext' for extensions)")
    
    # Recursive search
    batch_parser.add_argument("-nr", "--no-recursive", action="store_true",
                        help="Do not search recursively (search only in top directory)")
    
    # Output format
    batch_parser.add_argument("-of", "--output-format", default="vasp",
                        help="Output file format (default: vasp)")
    
    # Merge output
    batch_parser.add_argument("-m", "--merge", action="store_true",
                        help="Merge all structures into a single multi-frame file")
    
    # Merged output file
    batch_parser.add_argument("-o", "--merge-output", 
                        help="Output path for merged file (only used with --merge)")
    
    # Element order
    batch_parser.add_argument("-e", "--element-order", help="Custom element order (comma-separated)")
    
    # Constraints
    batch_parser.add_argument("-c", "--constraints", help="""Atomic constraints:
                         fixed (fix all atoms),
                         elements:C,H (fix specific elements),
                         layers:0.3,0.5 (fix atoms in specific fractional z-range),
                         indices:1,2,5-10 (fix specific atoms by index)""")
    
    # Sorting
    batch_parser.add_argument("-s", "--sort", choices=["ascending", "descending", "asc", "desc"], 
                       help="Sort atoms by chemical symbol (shorthand: asc, desc)")
    
    args = parser.parse_args(args)
    
    # If no command is specified, default to convert
    if args.command is None:
        # For backward compatibility, check if required convert args are present
        if hasattr(args, 'input') and hasattr(args, 'output'):
            args.command = "convert"
        else:
            parser.print_help()
            sys.exit(1)
    
    return args

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
        
        # 处理sort参数的简写
        if hasattr(args, 'sort') and args.sort:
            if args.sort == 'asc':
                args.sort = 'ascending'
            elif args.sort == 'desc':
                args.sort = 'descending'
        
        # Handle different commands
        if args.command == "convert" or args.command == "conv" or args.command == "c":
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
                print(f"Command: {args.command}")
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
            
            # Convert single file
            print(f"Converting {args.input} to {args.output}...")
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
            print(f"Conversion completed! Output file: {args.output}")
        
        elif args.command == "batch" or args.command == "b":
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
                print(f"Command: {args.command}")
                print(f"Input directory: {args.input_dir}")
                print(f"Output directory: {args.output_dir}")
                print(f"Pattern: {args.pattern}")
                print(f"Recursive: {not args.no_recursive}")
                print(f"Output format: {args.output_format}")
                print(f"Merge output: {args.merge}")
                print(f"Merge output file: {args.merge_output}")
                print(f"Element order: {element_mapping}")
                print(f"Constraints: {constraints}")
                print(f"Sort type: {args.sort}")
            
            # Batch convert files
            print(f"Starting batch conversion from {args.input_dir} to {args.output_dir}...")
            converter.batch_convert(
                input_dir=args.input_dir,
                output_dir=args.output_dir,
                pattern=args.pattern,
                recursive=not args.no_recursive,
                output_format=args.output_format,
                merge_output=args.merge,
                output_merged_file=args.merge_output,
                # Additional conversion parameters
                element_mapping=element_mapping,
                constraints=constraints,
                sort_type=args.sort
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