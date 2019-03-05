#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Before run any test, all the executables of programs that you are going to
use need to be set to their actual paths on your system accordingly, otherwise
tests that using these programs will be skipped.

If you want the test script do cleanup after testing, set True for CLEANUP,
then all the files and directories generated during testing will be removed
upon all tests finished.
"""

MUSCLE = r"C:\Users\tianz\bin\muscle.exe"
MAFFT = None
CLUSTAL = None

# Path to the executables of ancestral states reconstruction programs
CODEML = r"C:\Users\tianz\bin\codeml.exe"
RAXML = None

# Path to the executables of ML tree inference programs
FASTTREE = None
IQTREE = r"C:\Users\tianz\bin\iqtree.exe"
PHYML = None

# Path to the executables of ML tree inference programs
EVOLVER = r"C:\Users\tianz\bin\evolver.exe"
SEQGEN = None

CLEANUP = True

# The following executables are used by myself, I have all these executables
# installed into $HOME/bin directory and the directory is a part of $PATH.

# Path to the executables of multiple sequence alignment programs
# MUSCLE = 'muscle'  # MUSCLE v3.8.31
# MAFFT = 'mafft'  # MAFFT v7.407 (2018/Jul/23)
# CLUSTAL = 'clustal'  # Clustal Omega - 1.2.4
#
# # Path to the executables of ancestral states reconstruction programs
# CODEML = 'codeml'  # PAML 4.9e
# RAXML = 'raxml'  # RAxML version 8.2.12
#
# # Path to the executables of ML tree inference programs
# FASTTREE = 'FastTree'  # FastTree 2.1.10
# IQTREE = 'iqtree'  # IQ-TREE multicore version 1.5.4
# PHYML = 'PhyML'  # PhyML version 20120412
#
# # Path to the executables of ML tree inference programs
# EVOLVER = 'evolver'  # PAML 4.9e
# SEQGEN = 'seq-gen'  # Version 1.3.4


if __name__ == '__main__':
    pass
