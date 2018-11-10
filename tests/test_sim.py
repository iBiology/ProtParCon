#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import logging
import unittest

PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, PATH)
from imc.sim import _guess, sim

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP

EVOLVER = settings.EVOLVER
SEQGEN = settings.SEQGEN

ASR = {'evolver': EVOLVER, 'seq-gen': SEQGEN}


class TestSim(unittest.TestCase):
    def setUp(self):
        self.msa = os.path.join(PATH, 'tests', 'data', 'sim', 'msa.fa')
        self.tree = os.path.join(PATH, 'tests', 'data', 'sim', 'tree.newick')
        self.rm = ''
    
    def tearDown(self):
        if CLEANUP and self.rm and os.path.exists(self.rm):
            try:
                if os.path.isfile(self.rm):
                    os.remove(self.rm)
                else:
                    shutil.rmtree(self.rm)
            except OSError as err:
                logging.warning('Failed to cleanup testing data {} due to:'
                                '\n\t{}'.format(self.rm, err))
            
    @unittest.skipIf(not all(ASR.values()), 'No ASR program provided.')
    def test__guess(self):
        aligners = set(ASR.keys())
        self.assertSetEqual(aligners, set(_guess(a)[0] for a in ASR.values()))
    
    @unittest.skipIf(EVOLVER is None, 'No EVOLVER executable provided.')
    def test_evolver_default(self):
        out = sim(EVOLVER, self.tree)
        self.assertTrue(os.path.isfile(out))
        self.rm = out
        
    @unittest.skipIf(EVOLVER is None, 'No EVOLVER executable provided.')
    def test_evolver_msa(self):
        out = sim(EVOLVER, self.tree, msa=self.msa)
        self.assertTrue(os.path.isfile(out))
        self.rm = out
        
    @unittest.skipIf(EVOLVER is None, 'No EVOLVER executable provided.')
    def test_evolver_msa_outfile(self):
        outfile = os.path.join(PATH, 'tests', 'data', 'sim',
                               'evolver.output.simulation.tsv')
        out = sim(EVOLVER, self.tree, msa=self.msa, outfile=outfile, n=10)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(outfile, out)
        self.rm = outfile

    @unittest.skipIf(SEQGEN is None, 'No Seq-Gen executable provided.')
    def test_seqgen_default(self):
        out = sim(SEQGEN, self.tree)
        self.assertTrue(os.path.isfile(out))
        self.rm = out

    @unittest.skipIf(SEQGEN is None, 'No Seq-Gen executable provided.')
    def test_seqgen_msa(self):
        out = sim(SEQGEN, self.tree, msa=self.msa)
        self.assertTrue(os.path.isfile(out))
        self.rm = out

    @unittest.skipIf(SEQGEN is None, 'No Seq-Gen executable provided.')
    def test_seqgen_msa_wd(self):
        outfile = os.path.join(PATH, 'tests', 'data', 'sim',
                               'seqgen.output.simulation.tsv')
        out = sim(SEQGEN, self.tree, msa=self.msa, outfile=outfile, n=10)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(outfile, out)
        self.rm = outfile
        
  
if __name__ == '__main__':
    unittest.main()
