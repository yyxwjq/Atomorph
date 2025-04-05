#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Core conversion functionality for atomic structure files.
"""

import os
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
import ase.io
from ase import Atoms
import sys
from ase.io import read, write
from ase.constraints import FixAtoms
import numpy as np
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
import multiprocessing

# Filter out ASE spacegroup warnings
warnings.filterwarnings("ignore", category=UserWarning, module="ase.spacegroup.spacegroup")

class StructureConverter:
    """A class for converting atomic structure files between different formats."""
    
    # Default element order (periodic table order)
    DEFAULT_ELEMENT_ORDER = [
        "H", "He",
        "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
        "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
        "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]

    def __init__(self):
        """Initialize the StructureConverter."""
        self.single_frame_only_formats = ["vasp"]
        self.format_mapping = {
            "xyz": "extxyz",  # Map xyz to extxyz
            "extxyz": "extxyz"
        }
        self._element_order = None
        self._constraints = None
        self._layer_constraints = None
        self._use_fractional = False  # Default to using cartesian coordinates
        self.sort_order = "ascending"  # Default ascending sort
        self.MAX_FILE_SIZE = 100 * 1024 * 1024 * 1024  # 100GB
        self.MAX_WORKERS = os.cpu_count() or 4  # Use CPU cores or default to 4 threads
    
    def _get_element_sort(self, element: str) -> int:
        """
        Get the position of an element in the periodic table.
        
        Args:
            element: Element symbol
            
        Returns:
            Position in the periodic table
        """
        try:
            # Return element position in periodic table (0-based)
            return self.DEFAULT_ELEMENT_ORDER.index(element)
        except ValueError:
            # If element not in periodic table, return a large number to ensure it's last
            return 999
    
    def _check_file_size(self, file_path):
        """Check if file size exceeds limit"""
        file_size = os.path.getsize(file_path)
        if file_size > self.MAX_FILE_SIZE:
            raise ValueError(f"""File size ({file_size/1024/1024/1024:.2f}GB) 
                             , exceeds limit ({self.MAX_FILE_SIZE/1024/1024/1024:.2f}GB)""")

    def _show_progress(self, total, desc="Processing"):
        """Show progress bar"""
        return tqdm(total=total, desc=desc, unit="frame")

    def _process_frame(self, frame, output_path, constraints):
        """Process a single structure frame"""
        try:
            self._write_vasp(frame, output_path, constraints)
            return True, output_path
        except Exception as e:
            return False, str(e)

    def convert(
        self,
        input_path: Union[str, Path],
        output_path: Union[str, Path],
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        mode: Optional[str] = None,
        frame: Optional[str] = None,
        element_mapping: Optional[List[str]] = None,
        constraints: Optional[Union[str, List[str]]] = None,
        layer_constraints: Optional[List[Dict[str, float]]] = None,
        sort_type: Optional[str] = None,
        parallel: bool = True,
        multi_frame: bool = False,
        separate_dirs: bool = False,
        use_fractional: bool = False
    ) -> None:
        """
        Convert atomic structure file.
        
        Args:
            input_path: Input file path
            output_path: Output file path
            input_format: Input file format
            output_format: Output file format
            mode: Conversion mode ('single' or 'multi')
            frame: Frame selection string
            element_mapping: Element mapping list
            constraints: Atomic constraints
            layer_constraints: Layer-based constraints
            sort_type: Sort type for atoms
            parallel: Whether to use parallel processing
            multi_frame: Whether to process multiple frames
            separate_dirs: Whether to save frames in separate directories
            use_fractional: Whether to use fractional coordinates for layer constraints
        """
        try:
            # Check file size
            self._check_file_size(input_path)
            
            # Set options
            self._element_order = element_mapping
            self._constraints = constraints
            self._layer_constraints = layer_constraints
            self._use_fractional = use_fractional  # Store the flag
            
            # Convert paths to Path objects
            input_path = Path(input_path)
            output_path = Path(output_path)
            
            # Auto-detect file format
            if input_format is None:
                input_format = input_path.suffix[1:] if input_path.suffix else "vasp"
            if output_format is None:
                output_format = output_path.suffix[1:] if output_path.suffix else "vasp"
            
            # Map file formats
            input_format = self.format_mapping.get(input_format, input_format)
            output_format = self.format_mapping.get(output_format, output_format)
            
            # Auto-detect mode
            if mode is None:
                mode = "multi" if multi_frame else "single"
            
            # Validate mode
            self._validate_mode(mode, output_format)
            
            # Parse frame selection
            frame_indices = self._parse_frame_selection(frame) if frame else None
            
            # Read structures
            structures = self._read_structures(input_path, input_format, frame)
            
            # Apply transformations
            if sort_type:
                structures = self._sort_atoms(structures, sort_type)
            if self._element_order:
                structures = self._apply_element_order(structures)
            
            # Handle output
            if output_format == "vasp":
                # For VASP format, use custom write function
                if mode == "single" and len(structures) == 1:
                    # Single frame, single output
                    self._write_vasp(structures[0], output_path, self._constraints)
                else:
                    # Multiple frames or multi-frame mode
                    # Create output directory
                    output_path.mkdir(parents=True, exist_ok=True)
                    if parallel:
                        # Parallel processing for multi-frame
                        print(f"Using {self.MAX_WORKERS} threads for parallel processing...")
                        with ThreadPoolExecutor(max_workers=self.MAX_WORKERS) as executor:
                            # Create task list
                            futures = []
                            for i, structure in enumerate(structures):
                                if separate_dirs:
                                    # Use original frame index for directory name
                                    frame_idx = frame_indices[i] if frame_indices else i
                                    frame_dir = output_path / f"frame_{frame_idx+1}"
                                    frame_dir.mkdir(exist_ok=True)
                                    frame_output = frame_dir / "POSCAR"
                                else:
                                    frame_output = output_path / f"frame_{i+1}.vasp"
                                future = executor.submit(self._process_frame, structure, frame_output, self._constraints)
                                futures.append(future)
                            
                            # Show progress with progress bar
                            with self._show_progress(len(futures), "Processing frames") as pbar:
                                for future in as_completed(futures):
                                    success, result = future.result()
                                    if not success:
                                        print(f"Warning: Processing failed - {result}")
                                    pbar.update(1)
                    else:
                        # Serial processing for multi-frame
                        with self._show_progress(len(structures), "Processing frames") as pbar:
                            for i, structure in enumerate(structures):
                                if separate_dirs:
                                    # Use original frame index for directory name
                                    frame_idx = frame_indices[i] if frame_indices else i
                                    frame_dir = output_path / f"frame_{frame_idx+1}"
                                    frame_dir.mkdir(exist_ok=True)
                                    frame_output = frame_dir / "POSCAR"
                                else:
                                    frame_output = output_path / f"frame_{i+1}.vasp"
                                self._write_vasp(structure, frame_output, self._constraints)
                                pbar.update(1)
            else:
                # For other formats, use ASE's write function
                if mode == "single" and len(structures) == 1:
                    # Single frame, single output
                    ase.io.write(output_path, structures[0], format=output_format)
                else:
                    # Multiple frames or multi-frame mode
                    if output_format in self.single_frame_only_formats or mode == "multi":
                        # Create output directory
                        output_path.mkdir(parents=True, exist_ok=True)
                        for i, structure in enumerate(structures):
                            if separate_dirs:
                                # Use original frame index for directory name
                                frame_idx = frame_indices[i] if frame_indices else i
                                frame_dir = output_path / f"frame_{frame_idx+1}"
                                frame_dir.mkdir(exist_ok=True)
                                frame_output = frame_dir / f"structure.{output_format}"
                            else:
                                frame_output = output_path / f"frame_{i+1}.{output_format}"
                            ase.io.write(frame_output, structure, format=output_format)
                    else:
                        # Write all frames to a single file
                        ase.io.write(output_path, structures, format=output_format)
            
            print(f"Conversion completed! Output file: {output_path}")
            
        except Exception as e:
            raise ValueError(f"Conversion failed: {str(e)}")
    
    def _detect_mode(self, input_path: Path, input_format: str) -> str:
        """Detect conversion mode"""
        try:
            structures = ase.io.read(input_path, format=input_format, index=":")
            return "multi" if len(structures) > 1 else "single"
        except Exception:
            return "single"
    
    def _validate_mode(self, mode: str, output_format: str) -> None:
        """Validate conversion mode"""
        if mode == "multi" and output_format in self.single_frame_only_formats:
            if isinstance(output_format, str) and Path(output_format).suffix:
                raise ValueError(f"Output format {output_format} does not support multi-frame structures. Please choose another output format or use single-frame mode.")

    def _parse_frame_selection(self, frame_str: str) -> List[int]:
        """
        Parse frame selection string into a list of frame indices.
        
        Args:
            frame_str: Frame selection string (e.g., "1 2 3-6")
            
        Returns:
            List of frame indices (0-based)
        """
        if not frame_str:
            return []
            
        indices = set()
        parts = frame_str.split()
        
        for part in parts:
            if '-' in part:
                start, end = map(int, part.split('-'))
                indices.update(range(start-1, end))
            else:
                try:
                    indices.add(int(part)-1)
                except ValueError:
                    raise ValueError(f"Invalid frame number: {part}")
                
        return sorted(list(indices))

    def _read_structures(
        self,
        input_path: Path,
        input_format: str,
        frame: Optional[str] = None,
    ) -> List[Atoms]:
        """Read structures"""
        try:
            # Read all frames first
            all_structures = ase.io.read(input_path, format=input_format, index=":")
            if not isinstance(all_structures, list):
                all_structures = [all_structures]
            
            if frame is not None:
                # Parse frame selection
                frame_indices = self._parse_frame_selection(frame)
                if not frame_indices:
                    raise ValueError("Invalid frame selection format")
                
                # Select specified frames
                structures = []
                for idx in frame_indices:
                    if 0 <= idx < len(all_structures):
                        structures.append(all_structures[idx])
                    else:
                        print(f"Warning: Frame {idx+1} out of range, skipping")
                
                if not structures:
                    raise ValueError("No valid frames selected")
            else:
                structures = all_structures
            
            # Check if structures have valid cells - only issue warnings now
            for idx, structure in enumerate(structures):
                if structure.cell.rank < 3:
                    # XYZ files typically don't preserve lattice information
                    if str(input_path).lower().endswith('.xyz'):
                        raise ValueError(f"XYZ files do not contain lattice information. Please provide a file format with lattice data (e.g., VASP POSCAR) or set appropriate lattice parameters.")
                    else:
                        raise ValueError(f'Image {idx} lose the Lattice informations,')
            
            return structures
        except Exception as e:
            raise ValueError(f"Failed to read file: {str(e)}")
    
    def _sort_atoms(self, structures: List[Atoms], sort_type: str) -> List[Atoms]:
        """Sort atoms by element type"""
        for structure in structures:
            # Get current element order and positions
            elements = structure.get_chemical_symbols()
            positions = structure.positions.copy()
            numbers = structure.numbers.copy()
            
            # Create list of (element, index) tuples for sorting
            element_indices = list(enumerate(elements))
            
            # Sort based on element order
            if sort_type == "ascending":
                # Sort by alphabetical order
                element_indices.sort(key=lambda x: x[1])
            else:  # descending
                # Sort by alphabetical order
                element_indices.sort(key=lambda x: x[1], reverse=True)
            
            # Get sorted indices
            indices = [i for i, _ in element_indices]
            
            # Update structure
            structure.positions = positions[indices]
            structure.numbers = numbers[indices]
            structure.set_chemical_symbols([elements[i] for i in indices])
            
            # Update additional properties if they exist
            if structure.has('initial_magmoms'):
                magmoms = structure.get_initial_magmoms()
                structure.set_initial_magmoms(magmoms[indices])
            if structure.has('initial_charges'):
                charges = structure.get_initial_charges()
                structure.set_initial_charges(charges[indices])
            
            # Update atomic type order
            unique_elements = list(dict.fromkeys(structure.get_chemical_symbols()))
            structure.new_element_order = unique_elements
        
        return structures
    
    def _apply_element_order(self, structures: List[Atoms]) -> List[Atoms]:
        """Apply custom element order"""
        for structure in structures:
            # Get current element order and positions
            elements = structure.get_chemical_symbols()
            positions = structure.positions.copy()
            numbers = structure.numbers.copy()
            
            # Create mapping from element to desired position
            element_positions = {elem: i for i, elem in enumerate(self._element_order)}
            
            # Create list of (element, index) tuples for sorting
            element_indices = list(enumerate(elements))
            
            # Sort based on custom element order
            element_indices.sort(key=lambda x: element_positions.get(x[1], len(self._element_order)))
            
            # Get sorted indices
            indices = [i for i, _ in element_indices]
            
            # Update structure
            structure.positions = positions[indices]
            structure.numbers = numbers[indices]
            structure.set_chemical_symbols([elements[i] for i in indices])
            
            # Update additional properties if they exist
            if structure.has('initial_magmoms'):
                magmoms = structure.get_initial_magmoms()
                structure.set_initial_magmoms(magmoms[indices])
            if structure.has('initial_charges'):
                charges = structure.get_initial_charges()
                structure.set_initial_charges(charges[indices])
            
            # Update atomic type order
            structure.new_element_order = self._element_order
        
        return structures
    
    def _apply_constraints(self, structures: List[Atoms]) -> List[Atoms]:
        """Apply atomic constraints"""
        for structure in structures:
            if self._constraints == "fixed":
                # Fix all atoms
                structure.set_constraint(FixAtoms(indices=list(range(len(structure)))))
            elif isinstance(self._constraints, list):
                if all(isinstance(x, str) for x in self._constraints):
                    # Fix specific elements
                    indices = []
                    for i, symbol in enumerate(structure.get_chemical_symbols()):
                        if symbol in self._constraints:
                            indices.append(i)
                    if indices:
                        structure.set_constraint(FixAtoms(indices=indices))
                elif all(isinstance(x, int) for x in self._constraints):
                    # Fix specific atom indices
                    structure.set_constraint(FixAtoms(indices=self._constraints))
        
        return structures
    
    def _apply_layer_constraints(self, structures: List[Atoms]) -> List[Atoms]:
        """Apply layer-based constraints"""
        for structure in structures:
            indices = []
            for layer in self._layer_constraints:
                # Get atoms in layer
                z_positions = structure.positions[:, 2]
                layer_indices = np.where((z_positions >= layer['start']) & (z_positions < layer['end']))[0]
                indices.extend(layer_indices)
            if indices:
                structure.set_constraint(FixAtoms(indices=indices))
        
        return structures
    
    def _write_vasp(self, structure: Atoms, output_path: Union[str, Path], constraints: Optional[Union[str, List[str]]] = None) -> None:
        """Write structure in VASP format."""
        # Create output directory if needed
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Get chemical symbols and counts
        elements = structure.get_chemical_symbols()
        if hasattr(structure, 'new_element_order'):
            # Use custom element order
            unique_elements = [e for e in structure.new_element_order if e in elements]
        else:
            # Default to periodic table order
            unique_elements = sorted(set(elements), key=self._get_element_sort)
        
        # Statistics variables
        fixed_count = 0
        total_atoms = len(structure)
        
        # Calculate fixed atoms
        fixed_atoms = []
        if constraints:
            if constraints == 'fixed':
                # Fix all atoms
                fixed_atoms = list(range(total_atoms))
                fixed_count = total_atoms
            elif isinstance(constraints, list) and len(constraints) >= 2:
                constraint_type = constraints[0]
                
                if constraint_type == 'elements':
                    # Fix specific elements
                    fixed_elements = constraints[1:]
                    for i, symbol in enumerate(elements):
                        if symbol in fixed_elements:
                            fixed_atoms.append(i)
                    fixed_count = len(fixed_atoms)
                    
                elif constraint_type == 'layers':
                    # Fix specific layers
                    use_fractional = getattr(self, '_use_fractional', True)
                    
                    # Default to only process the first layer range parameter
                    if len(constraints) > 1:
                        layer_range = constraints[1]
                        try:
                            # Ensure parameter format is correct
                            if ',' in layer_range:
                                start, end = map(float, layer_range.split(','))
                                
                                if use_fractional:
                                    # Use fractional coordinates (recommended way)
                                    cell_inv = np.linalg.inv(structure.cell)
                                    positions = structure.positions
                                    frac_positions = np.dot(positions, cell_inv)
                                    
                                    for i, frac_pos in enumerate(frac_positions):
                                        frac_z = frac_pos[2]
                                        if start <= frac_z <= end:
                                            fixed_atoms.append(i)
                                else:
                                    # Use Cartesian coordinates
                                    for i, pos in enumerate(structure.positions):
                                        z = pos[2]
                                        if start <= z <= end:
                                            fixed_atoms.append(i)
                            else:
                                print(f"Warning: Layer range parameter format incorrect '{layer_range}', should be 'start,end'")
                                
                        except ValueError as e:
                            print(f"Warning: Unable to parse layer range '{layer_range}': {str(e)}")
                    else:
                        print("Warning: No layer range parameter provided")
                    
                    fixed_count = len(fixed_atoms)
                    
                elif constraint_type == 'indices':
                    # Fix specific atom indices
                    for idx_str in constraints[1:]:
                        try:
                            if '-' in str(idx_str):
                                start, end = map(int, idx_str.split('-'))
                                for idx in range(start-1, end):  # Convert to 0-based index
                                    if 0 <= idx < total_atoms:
                                        fixed_atoms.append(idx)
                            else:
                                idx = int(idx_str) - 1  # Convert to 0-based index
                                if 0 <= idx < total_atoms:
                                    fixed_atoms.append(idx)
                        except ValueError as e:
                            print(f"Warning: Unable to parse index '{idx_str}': {str(e)}")
                    
                    fixed_count = len(fixed_atoms)
        
        # Remove duplicates
        fixed_atoms = list(set(fixed_atoms))
        fixed_count = len(fixed_atoms)
        
        # Write VASP file
        with open(output_path, 'w') as f:
            # Write elements
            f.write(' '.join(unique_elements) + '\n')
            
            # Write scale factor
            f.write('1.0\n')
            
            # Write cell
            cell = structure.get_cell()
            for i in range(3):
                f.write(f'{cell[i,0]:20.16f} {cell[i,1]:20.16f} {cell[i,2]:20.16f}\n')
            
            # Write elements and counts
            f.write(' '.join(unique_elements) + '\n')
            counts = [elements.count(e) for e in unique_elements]
            f.write(' '.join(map(str, counts)) + '\n')
            
            # Write selective dynamics and coordinates
            if constraints:
                f.write('Selective dynamics\n')
            f.write('Cartesian\n')
            
            # Write positions with constraints
            # Write atomic coordinates in element order
            for element in unique_elements:
                element_indices = [i for i, symbol in enumerate(elements) if symbol == element]
                for i in element_indices:
                    pos = structure.positions[i]
                    is_fixed = i in fixed_atoms
                    
                    if constraints:
                        if is_fixed:
                            f.write(f'{pos[0]:20.16f} {pos[1]:20.16f} {pos[2]:20.16f} F F F\n')
                        else:
                            f.write(f'{pos[0]:20.16f} {pos[1]:20.16f} {pos[2]:20.16f} T T T\n')
                    else:
                        f.write(f'{pos[0]:20.16f} {pos[1]:20.16f} {pos[2]:20.16f}\n')
        
        print(f"VASP file written to: {output_path}")
        print(f"Total atoms: {total_atoms}, Fixed atoms: {fixed_count}")

    def batch_convert(
        self,
        input_dir: Union[str, Path],
        output_dir: Union[str, Path],
        pattern: str = "POSCAR",
        recursive: bool = True,
        output_format: str = "vasp",
        merge_output: bool = False,
        output_merged_file: Optional[Union[str, Path]] = None,
        **convert_kwargs
    ) -> Dict[str, Any]:
        """
        Batch convert all files matching pattern in a directory and its subdirectories.
        
        This function traverses the given directory, finds all files matching the specified
        pattern, and converts them to the target format. It can operate recursively through
        subdirectories and optionally merge all structures into a single multi-frame file.
        
        Args:
            input_dir: Input directory path
            output_dir: Output directory path
            pattern: File pattern to match (default: "POSCAR"). Can be a specific filename
                     or a pattern with wildcard (e.g., "*.vasp")
            recursive: Whether to search recursively through subdirectories (default: True)
            output_format: Output file format (default: "vasp")
            merge_output: Whether to merge all structures into a single multi-frame file (default: False)
            output_merged_file: Output path for merged file (used only if merge_output is True)
            **convert_kwargs: Additional conversion parameters to pass to convert() method,
                              such as element_mapping, constraints, sort_type, etc.
            
        Returns:
            Dictionary with statistics about the conversion process, including:
            - total_files: Total number of files processed
            - successful: Number of files successfully converted
            - failed: Number of files that failed to convert
            - failed_files: List of failed files with error messages
            - merged_file: Path to the merged output file (if merge_output=True)
            - merged_count: Number of structures in the merged file
            
        Example:
            >>> converter = StructureConverter()
            >>> # Convert all POSCAR files to XYZ
            >>> stats = converter.batch_convert(
            ...     input_dir="./calculations",
            ...     output_dir="./converted",
            ...     pattern="POSCAR",
            ...     output_format="xyz"
            ... )
            >>> print(f"Converted {stats['successful']} files")
            >>> 
            >>> # Merge all VASP files into one XYZ file
            >>> stats = converter.batch_convert(
            ...     input_dir="./calculations",
            ...     output_dir="./merged",
            ...     pattern="*.vasp",
            ...     output_format="xyz",
            ...     merge_output=True
            ... )
            >>> print(f"Merged {stats['merged_count']} structures")
        """
        input_dir = Path(input_dir)
        output_dir = Path(output_dir)
        
        if not input_dir.exists() or not input_dir.is_dir():
            raise ValueError(f"Input directory '{input_dir}' does not exist or is not a directory")
        
        # Create output directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find all matching files
        matching_files = []
        
        if recursive:
            for root, _, files in os.walk(input_dir):
                for file in files:
                    if file == pattern or (pattern.startswith('*.') and file.endswith(pattern[1:])):
                        file_path = Path(root) / file
                        matching_files.append(file_path)
        else:
            for file in input_dir.iterdir():
                if file.is_file() and (file.name == pattern or 
                                    (pattern.startswith('*.') and file.name.endswith(pattern[1:]))):
                    matching_files.append(file)
        
        print(f"Found {len(matching_files)} files matching pattern '{pattern}'")
        
        if not matching_files:
            return {"status": "warning", "message": f"No files matching pattern '{pattern}' found"}
        
        # Statistics
        stats = {
            "total_files": len(matching_files),
            "successful": 0,
            "failed": 0,
            "skipped": 0,
            "failed_files": []
        }
        
        # Structures for merging
        all_structures = []
        
        # Process each file
        progress_bar = tqdm(matching_files, desc="Converting files", unit="file")
        for input_file in progress_bar:
            try:
                # Determine output path
                rel_path = input_file.relative_to(input_dir)
                parent_dirs = str(rel_path.parent)
                
                if merge_output:
                    # Read structures and add to merge list
                    try:
                        structures = self._read_structures(input_file, None)
                        
                        # Add source file info to each structure
                        for structure in structures:
                            structure.info['source_file'] = str(rel_path)
                        
                        all_structures.extend(structures)
                        stats["successful"] += 1
                    except Exception as e:
                        print(f"Error reading {input_file}: {str(e)}")
                        stats["failed"] += 1
                        stats["failed_files"].append({"file": str(input_file), "error": str(e)})
                else:
                    # Create output subdirectory
                    if parent_dirs and parent_dirs != '.':
                        output_subdir = output_dir / parent_dirs
                        output_subdir.mkdir(parents=True, exist_ok=True)
                        output_file = output_subdir / f"{input_file.stem}.{output_format}"
                    else:
                        output_file = output_dir / f"{input_file.stem}.{output_format}"
                    
                    # Convert file
                    try:
                        self.convert(
                            input_path=input_file,
                            output_path=output_file,
                            output_format=output_format,
                            **convert_kwargs
                        )
                        stats["successful"] += 1
                    except Exception as e:
                        print(f"Error converting {input_file}: {str(e)}")
                        stats["failed"] += 1
                        stats["failed_files"].append({"file": str(input_file), "error": str(e)})
            
            except Exception as e:
                print(f"Error processing {input_file}: {str(e)}")
                stats["failed"] += 1
                stats["failed_files"].append({"file": str(input_file), "error": str(e)})
        
        # If merging, write all structures to a single file
        if merge_output and all_structures:
            if not output_merged_file:
                output_merged_file = output_dir / f"merged.{output_format}"
            else:
                output_merged_file = Path(output_merged_file)
            
            try:
                print(f"Writing {len(all_structures)} structures to {output_merged_file}")
                
                # Check if directory exists
                output_merged_file.parent.mkdir(parents=True, exist_ok=True)
                
                # Write using ASE
                ase.io.write(str(output_merged_file), all_structures, format=output_format)
                
                stats["merged_file"] = str(output_merged_file)
                stats["merged_count"] = len(all_structures)
            except Exception as e:
                print(f"Error writing merged file: {str(e)}")
                stats["merge_error"] = str(e)
        
        # Print summary
        print(f"\nConversion summary:")
        print(f"  Total files:   {stats['total_files']}")
        print(f"  Successful:    {stats['successful']}")
        print(f"  Failed:        {stats['failed']}")
        
        if merge_output and "merged_file" in stats:
            print(f"  Merged to:     {stats['merged_file']} ({stats['merged_count']} structures)")
        
        return stats