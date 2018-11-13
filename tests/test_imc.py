#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import logging
import unittest

PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, PATH)
from ProtParCon.imc import imc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP

MUSCLE = settings.MUSCLE
CODEML = settings.CODEML
IQTREE = settings.IQTREE
EVOLVER = settings.EVOLVER
SEQGEN = settings.SEQGEN


class TestAUT(unittest.TestCase):
    def setUp(self):
        self.seq = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'sequence.fa')
        self.msa = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'msa.fa')
        self.tree = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'tree.newick')
        self.asr = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'asr.tsv')
        self.sim = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'sim.tsv')
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

    @unittest.skipIf((MUSCLE is None) or (CODEML is None),
                     'No valid executable provided.')
    def test_imc_seq_muscle_codeml(self):
        wd = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'simple')
        outs = ['ancestors.tsv', 'ProtParCon.counts.tsv', 'ProtParCon.details.tsv',
                'msa.fa', 'trimmed.msa.fa']
        pars, cons, details, aup = imc(self.seq, self.tree, aligner=MUSCLE,
                                       ancestor=CODEML, wd=wd)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(details)
        self.assertIsNone(aup)
        self.assertListEqual(outs, sorted(os.listdir(wd)))
        self.rm = wd
    
    @unittest.skipIf((MUSCLE is None) or (CODEML is None) or (EVOLVER is None),
                     'No valid executable provided.')
    def test_imc_seq_muscle_codeml_evolver_iqtree(self):
        wd = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'seq')
        outs = ['ancestors.tsv', 'ProtParCon.aut.txt', 'ProtParCon.counts.tsv',
                'ProtParCon.details.tsv', 'msa.fa', 'simulations.tsv',
                'trimmed.msa.fa']
        pars, cons, details, aup = imc(self.seq, tree=self.tree, aligner=MUSCLE,
                                       ancestor=CODEML,  simulator=EVOLVER,
                                       toper=IQTREE, wd=wd)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(details)
        self.assertLessEqual(0.0, aup)
        self.assertListEqual(outs, sorted(os.listdir(wd)))
        self.rm = wd

    @unittest.skipIf((CODEML is None) or (EVOLVER is None),
                     'No valid executable provided.')
    def test_imc_msa_codeml_evolver_iqtree(self):
        wd = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'msa')
        outs = ['ancestors.tsv', 'ProtParCon.aut.txt', 'ProtParCon.counts.tsv',
                'ProtParCon.details.tsv', 'simulations.tsv', 'trimmed.msa.fa']
        pars, cons, details, aup = imc(self.msa, tree=self.tree,
                                       ancestor=CODEML,  simulator=EVOLVER,
                                       toper=IQTREE, wd=wd)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(details)
        self.assertLessEqual(0.0, aup)
        self.assertListEqual(outs, sorted(os.listdir(wd)))
        self.rm = wd

    @unittest.skipIf(EVOLVER is None, 'No valid executable provided.')
    def test_imc_asr_evolver_iqtree(self):
        wd = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'asr')
        outs = ['ProtParCon.aut.txt', 'ProtParCon.counts.tsv',
                'ProtParCon.details.tsv', 'simulations.tsv']
        pars, cons, details, aup = imc(self.asr, simulator=EVOLVER,
                                       toper=IQTREE, wd=wd)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(details)
        self.assertLessEqual(0.0, aup)
        self.assertListEqual(outs, sorted(os.listdir(wd)))
        self.rm = wd
    
    def test_imc_obs(self):
        wd = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'obs')
        outs = ['ProtParCon.counts.tsv', 'ProtParCon.details.tsv']
        pars, cons, details, aup = imc(self.asr, wd=wd)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(details)
        self.assertIsNone(aup)
        self.assertListEqual(outs, sorted(os.listdir(wd)))
        self.rm = wd
        
    def test_imc_sim(self):
        wd = os.path.join(PATH, 'tests', 'data', 'ProtParCon', 'sim')
        outs = ['ProtParCon.counts.tsv', 'ProtParCon.details.tsv']
        pars, cons, details, aup = imc(self.sim, wd=wd)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(details)
        self.assertIsNone(aup)
        self.assertListEqual(outs, sorted(os.listdir(wd)))
        self.rm = wd

  
if __name__ == '__main__':
    unittest.main()
