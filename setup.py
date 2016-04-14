#!/usr/bin/env python

from setuptools import setup
setup(
    name='PhyloGenClust',
    version='0.1dev',
    author='Mehmood Alam Khan',
    author_email='malagori@kth.se',
    url='https://github.com/malagori/PhyloGenClust',
    packages=['phylogenclust',],
    scripts = ['phyloGenClust', 'scripts/convertSpeciesTreeLabels.py', 'scripts/cutGeneTree' ,
               'scripts/prepareSeqNames.py', 
               'scripts/cutReconciledGeneTree.py',
               'scripts/greedyApproach.py', 'scripts/prepareSeqNames.py' ],
    license='GPLv3',
    long_description=open('README.md').read(),
    install_requires= ['Biopython >=1.62'], 
)
