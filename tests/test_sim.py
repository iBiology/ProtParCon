#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import logging
import unittest

PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, PATH)
from ProtParCon.sim import _guess, sim

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP

EVOLVER = settings.EVOLVER
SEQGEN = settings.SEQGEN

ASR = {'evolver': EVOLVER, 'seq-gen': SEQGEN}


class TestSim(unittest.TestCase):
    def setUp(self):
        self.msa = os.path.join(PATH, 'tests', 'data', 'sim', 'msa.fa')
        self.anc = os.path.join(PATH, 'tests', 'data', 'sim', 'asr.tsv')
        self.tree = os.path.join(PATH, 'tests', 'data', 'sim', 'tree.newick')
        self.out = os.path.join(PATH, 'tests', 'data', 'sim', 'sim.out')
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
            
    @unittest.skipIf(not all(ASR.values()), 'No simulate program provided.')
    def test__guess(self):
        aligners = set(ASR.keys())
        self.assertSetEqual(aligners, set(_guess(a)[0] for a in ASR.values()))
    
    @unittest.skipIf(EVOLVER is None, 'No EVOLVER executable provided.')
    def test_evolver_default(self):
        out = sim(EVOLVER, self.tree, length=100)
        self.assertTrue(os.path.isfile(out))
        self.rm = out

    @unittest.skipIf(EVOLVER is None, 'No EVOLVER executable provided.')
    def test_evolver_out(self):
        out = sim(EVOLVER, self.tree, length=100, outfile=self.out)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(out, self.out)
        self.rm = out

    @unittest.skipIf(EVOLVER is None, 'No EVOLVER executable provided.')
    def test_evolver_msa(self):
        out = sim(EVOLVER, self.tree, sequence=self.msa, outfile=self.out)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(out, self.out)
        self.rm = out

    @unittest.skipIf(EVOLVER is None, 'No EVOLVER executable provided.')
    def test_evolver_anc(self):
        out = sim(EVOLVER, self.tree, sequence=self.anc, outfile=self.out)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(out, self.out)
        self.rm = out

    @unittest.skipIf(SEQGEN is None, 'No Seq-Gen executable provided.')
    def test_seqgen_default(self):
        out = sim(SEQGEN, self.tree, length=100, outfile=self.out)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(out, self.out)
        self.rm = out

    @unittest.skipIf(SEQGEN is None, 'No Seq-Gen executable provided.')
    def test_seqgen_out(self):
        out = sim(SEQGEN, self.tree, msa=self.msa, outfile=self.out)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(out, self.out)
        self.rm = out

    @unittest.skipIf(SEQGEN is None, 'No Seq-Gen executable provided.')
    def test_seqgen_msa(self):
        out = sim(SEQGEN, self.tree, msa=self.msa, outfile=self.out, n=10)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(out, self.out)
        self.rm = out

    @unittest.skipIf(SEQGEN is None, 'No Seq-Gen executable provided.')
    def test_seqgen_anc(self):
        outfile = os.path.join(PATH, 'tests', 'data', 'sim',
                               'seqgen.output.simulation.tsv')
        out = sim(SEQGEN, self.tree, msa=self.anc, outfile=outfile, n=10)
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(outfile, out)
        self.rm = outfile
        
  
if __name__ == '__main__':
    unittest.main()
