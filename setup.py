#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='biomabuilder',
    version='0.1.0',
    description='BioMaBuilder is a program for modelling the macro-complex structure of biomolecules formed by proteins and DNA/RNA, starting with each pairing interaction of a complex: protein-protein, protein-DNA/RNA',
    author ='Carolina Hernandez-Oliver, Ricard Lambea-Jane and Jose Vicente Roig-Genoves',
    author_email='carolina.hernandez01@estudiant.upf.edu, ricard.lambea01@estudiant.upf.edu, josevicente.roig01@estudiant.upf.edu',
    url='https://github.com/ricardlambea/SBI-PYT_project',
    packages=setuptools.find_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
