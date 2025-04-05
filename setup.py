#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Setup script for Atomorph package.
"""

from setuptools import setup, find_packages
import os
import re

# Read the version from __init__.py
with open(os.path.join('atomorph', '__init__.py'), 'r') as f:
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M)
    if version_match:
        version = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string in atomorph/__init__.py")

# Read README.md for long_description
with open('README.md', 'r') as f:
    long_description = f.read()

# Installation requirements
install_requires = [
    'numpy>=1.20.0',
    'ase>=3.22.1',
    'tqdm>=4.62.0',
]

# Development requirements
dev_requires = [
    'pytest>=6.0.0',
    'pytest-cov>=2.10.0',
    'black>=21.5b2',
    'flake8>=3.9.0',
    'mypy>=0.812',
    'sphinx>=4.0.0',
    'sphinx-rtd-theme>=0.5.0',
]

setup(
    name="atomorph",
    version=version,
    description="A tool for converting atomic structure files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="The Atomorph Team",
    author_email="contact@atomorph.org",
    url="https://github.com/atomorph/atomorph",
    packages=find_packages(),
    install_requires=install_requires,
    extras_require={
        "dev": dev_requires,
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Natural Language :: Chinese (Simplified)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "atomorph=atomorph.cli:main",
        ],
    },
    keywords='atomic structure converter VASP XYZ ASE',
    project_urls={
        "Documentation": "https://atomorph.readthedocs.io/",
        "Source": "https://github.com/atomorph/atomorph",
        "Bug Reports": "https://github.com/atomorph/atomorph/issues",
    },
    include_package_data=True,
    license='MIT',
) 