#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ProtParCon import msa, asr, imc, sim, detect

sequence = 'lysozyme.fa'
tree = 'lysozyme.newick'

# Replace the path with real path on your system
muscle = 'muscle'
codeml = 'codeml'
evolver = 'evolver'

alignment = msa(muscle, sequence, trimming=True, verbose=True)
ancestors = asr(codeml, alignment, tree, 'JTT', verbose=True)
obs = imc(ancestors, verbose=True)
simulations = sim(evolver, tree=tree, sequence=alignment, verbose=True)
exp = imc(simulations, verbose=True)

# If you choose to identify changes step by step, after get obs and exp, you
# need to handel these results by your self. Each of them consists
# of five objects: pars (dict), cons (dict), divs (dict), details (list),
# and length (int). The three dicts have the same keys, and their values are
# lists. See function imc() for details, refer function detect() to see an
# example for handling pairwise comparison.
