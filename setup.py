#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Setup script for Atomorph package.
"""

from setuptools import setup, find_packages
import os
import re

# Get version number
def get_version():
    init_py = open('atomorph/__init__.py').read()
    return re.search(r"__version__ = ['\"]([^'\"]+)['\"]", init_py).group(1)

# Read long description
def get_long_description():
    with open('README.md', encoding='utf-8') as f:
        return f.read()

# Installation requirements
install_requires = [
    'numpy>=1.20.0',
    'ase>=3.22.1',
    'tqdm>=4.62.0',
]

# Development requirements
dev_requires = [
    'pytest>=6.0.0',
    'pytest-cov>=2.12.0',
    'black>=21.5b2',
    'isort>=5.9.0',
    'flake8>=3.9.0',
    'mypy>=0.812',
    'sphinx>=4.0.0',
    'sphinx-rtd-theme>=0.5.2',
    'twine>=3.4.0',
    'build>=0.7.0',
]

setup(
    name="atomorph",
    version=get_version(),
    description="A tool for converting atomic structure files",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    author="Atomorph Team",
    author_email="atomorph@example.com",
    url="https://github.com/atomorph-team/atomorph",
    packages=find_packages(),
    install_requires=install_requires,
    extras_require={
        "dev": dev_requires,
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "atomorph=atomorph.cli:main",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/atomorph-team/atomorph/issues",
        "Source": "https://github.com/atomorph-team/atomorph",
        "Documentation": "https://github.com/atomorph-team/atomorph",
    },
) 