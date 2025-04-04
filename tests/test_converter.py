#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the StructureConverter class.
"""

import os
import sys
import unittest
from pathlib import Path
import tempfile
import shutil
import numpy as np

# Add package to path for testing
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from atomorph import StructureConverter
from ase import Atoms


class TestStructureConverter(unittest.TestCase):
    """Tests for the StructureConverter class."""
    
    def setUp(self):
        """Set up the test environment."""
        self.converter = StructureConverter()
        self.test_dir = Path(tempfile.mkdtemp())
        
        # Create a simple structure for testing
        self.structure = Atoms(
            symbols=['C', 'H', 'H', 'H', 'H'],
            positions=[
                [0.0, 0.0, 0.0],
                [0.6, 0.6, 0.6],
                [-0.6, -0.6, 0.6],
                [0.6, -0.6, -0.6],
                [-0.6, 0.6, -0.6]
            ],
            cell=[[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]
        )
        
        # Save the structure as XYZ
        self.xyz_file = self.test_dir / 'test.xyz'
        self.structure.write(self.xyz_file, format='xyz')
        
    def tearDown(self):
        """Clean up after the test."""
        shutil.rmtree(self.test_dir)
    
    def test_element_sorting(self):
        """Test element sorting functionality."""
        # Create a structure with mixed elements
        structure = Atoms(
            symbols=['H', 'C', 'O', 'H', 'C'],
            positions=np.random.random((5, 3)) * 3,
            cell=[[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]]
        )
        
        # Save to temporary file
        xyz_file = self.test_dir / 'mixed.xyz'
        vasp_file = self.test_dir / 'mixed.vasp'
        structure.write(xyz_file, format='xyz')
        
        # Convert with ascending sorting
        self.converter.convert(
            input_path=xyz_file,
            output_path=vasp_file,
            sort_type="ascending"
        )
        
        # Check if the VASP file was created
        self.assertTrue(vasp_file.exists())
        
        # Read back the structure
        from ase.io import read
        sorted_structure = read(vasp_file, format='vasp')
        
        # Check if elements are sorted
        symbols = sorted_structure.get_chemical_symbols()
        self.assertEqual(symbols[0], 'C')  # C comes before H and O alphabetically
        
    def test_constraints(self):
        """Test constraints functionality."""
        # Test fixed elements constraint
        vasp_fixed_file = self.test_dir / 'fixed_elements.vasp'
        
        self.converter.convert(
            input_path=self.xyz_file,
            output_path=vasp_fixed_file,
            constraints=["elements", "H"]
        )
        
        # Check if the file exists
        self.assertTrue(vasp_fixed_file.exists())
        
        # Check the content for selective dynamics
        with open(vasp_fixed_file, 'r') as f:
            content = f.read()
            self.assertIn('Selective dynamics', content)
            self.assertIn('F F F', content)  # Fixed atoms
            self.assertIn('T T T', content)  # Movable atoms
    
    def test_multi_frame(self):
        """Test multi-frame processing."""
        # Create a multi-frame XYZ file
        multi_xyz = self.test_dir / 'multi.xyz'
        with open(multi_xyz, 'w') as f:
            # Write the same structure twice
            for _ in range(2):
                f.write('5\n')
                f.write('Atoms\n')
                f.write('C 0.0 0.0 0.0\n')
                f.write('H 0.6 0.6 0.6\n')
                f.write('H -0.6 -0.6 0.6\n')
                f.write('H 0.6 -0.6 -0.6\n')
                f.write('H -0.6 0.6 -0.6\n')
        
        # Output directory
        output_dir = self.test_dir / 'multi_output'
        
        # Convert with multi-frame
        self.converter.convert(
            input_path=multi_xyz,
            output_path=output_dir,
            multi_frame=True
        )
        
        # Check if directory was created with frame files
        self.assertTrue(output_dir.is_dir())
        self.assertTrue((output_dir / 'frame_1.vasp').exists())
        self.assertTrue((output_dir / 'frame_2.vasp').exists())
    
    def test_parse_frame_selection(self):
        """Test frame selection parsing."""
        # Test single frame
        indices = self.converter._parse_frame_selection("2")
        self.assertEqual(indices, [1])  # 0-based index
        
        # Test multiple frames
        indices = self.converter._parse_frame_selection("1 3 5")
        self.assertEqual(indices, [0, 2, 4])
        
        # Test range
        indices = self.converter._parse_frame_selection("2-4")
        self.assertEqual(indices, [1, 2, 3])
        
        # Test mixed
        indices = self.converter._parse_frame_selection("1 3-5 7")
        self.assertEqual(indices, [0, 2, 3, 4, 6])


if __name__ == '__main__':
    unittest.main() 