#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ProtParCon import msa, asr, imc, sim, detect

CLEANUP = True

sequence = 'lysozyme.fa'
tree = 'lysozyme.newick'

# Replace the path with real path on your system
muscle = 'muscle'
codeml = 'codeml'
evolver = 'evolver'

result = imc(sequence, tree, aligner=muscle, ancestor=codeml, simulator=evolver,
             verbose=True)
comparisons = detect()

