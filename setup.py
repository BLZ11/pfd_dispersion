#!/usr/bin/env python3
"""
Setup script for the PFD Dispersion Calculator package.

Installation:
    pip install .

Development installation:
    pip install -e .
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pfd-dispersion",
    version="1.0.0",
    author="Barbaro Zulueta",
    author_email="baz14@pitt.edu",
    description="Petersson-Frisch Dispersion (PFD) calculator for two-body and three-body dispersion energies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BLZ11/pfd_dispersion",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
    ],
    entry_points={
        "console_scripts": [
            "pfd-calc=pfd_calculator:main",
        ],
    },
    keywords=[
        "dispersion",
        "DFT",
        "density functional theory",
        "quantum chemistry",
        "computational chemistry",
        "PFD",
        "Axilrod-Teller-Muto",
    ],
    project_urls={
        "Bug Reports": "https://github.com/BLZ11/pfd_dispersion/issues",
        "Source": "https://github.com/BLZ11/pfd_dispersion",
    },
)
