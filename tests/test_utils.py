#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for the utils module.
"""

import os
import sys
import unittest
from pathlib import Path
import tempfile
import shutil

# Add package to path for testing
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from atomorph.utils import (
    ensure_dir, get_version, validate_constraints, format_file_size
)


class TestUtils(unittest.TestCase):
    """Tests for the utils module."""
    
    def setUp(self):
        """Set up the test environment."""
        self.test_dir = Path(tempfile.mkdtemp())
    
    def tearDown(self):
        """Clean up after the test."""
        shutil.rmtree(self.test_dir)
    
    def test_ensure_dir(self):
        """Test ensure_dir function."""
        # Test creating a new directory
        new_dir = self.test_dir / 'new_dir'
        result = ensure_dir(new_dir)
        
        self.assertTrue(new_dir.exists())
        self.assertTrue(new_dir.is_dir())
        self.assertEqual(result, new_dir)
        
        # Test with existing directory
        result = ensure_dir(new_dir)
        self.assertEqual(result, new_dir)
        
        # Test nested directory
        nested_dir = self.test_dir / 'parent' / 'child' / 'grandchild'
        result = ensure_dir(nested_dir)
        
        self.assertTrue(nested_dir.exists())
        self.assertTrue(nested_dir.is_dir())
        self.assertEqual(result, nested_dir)
    
    def test_get_version(self):
        """Test get_version function."""
        # Version should be a string
        version = get_version()
        self.assertIsInstance(version, str)
        
        # Version should match expected format (e.g., "1.0.0")
        import re
        version_pattern = r'^\d+\.\d+\.\d+$'
        self.assertTrue(re.match(version_pattern, version))
    
    def test_validate_constraints(self):
        """Test validate_constraints function."""
        # Test None constraints
        self.assertTrue(validate_constraints(None))
        
        # Test string constraints
        self.assertTrue(validate_constraints("fixed"))
        self.assertFalse(validate_constraints("invalid"))
        
        # Test list constraints
        self.assertTrue(validate_constraints(["elements", "H"]))
        self.assertTrue(validate_constraints(["layers", "0.3,0.5"]))
        self.assertTrue(validate_constraints(["indices", "1", "2", "5-10"]))
        
        # Test invalid list constraints
        self.assertFalse(validate_constraints([]))
        self.assertFalse(validate_constraints(["elements"]))
        self.assertFalse(validate_constraints(["invalid", "value"]))
    
    def test_format_file_size(self):
        """Test format_file_size function."""
        # Test bytes
        self.assertEqual(format_file_size(512), "512 B")
        
        # Test kilobytes
        self.assertEqual(format_file_size(1024), "1.00 KB")
        self.assertEqual(format_file_size(2048), "2.00 KB")
        
        # Test megabytes
        mb_size = 1024 * 1024 * 5.25
        self.assertEqual(format_file_size(mb_size), "5.25 MB")
        
        # Test gigabytes
        gb_size = 1024 * 1024 * 1024 * 2.5
        self.assertEqual(format_file_size(gb_size), "2.50 GB")


if __name__ == '__main__':
    unittest.main() 