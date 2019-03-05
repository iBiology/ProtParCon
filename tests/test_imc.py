#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import logging
import unittest

PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
wd = os.path.join(PATH, 'tests', 'data', 'imc')

sys.path.insert(0, PATH)
from ProtParCon.imc import imc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP
# CLEANUP = True

MUSCLE = settings.MUSCLE
CODEML = settings.CODEML
IQTREE = settings.IQTREE
EVOLVER = settings.EVOLVER
SEQGEN = settings.SEQGEN


class TestAUT(unittest.TestCase):
    def setUp(self):
        self.seq = os.path.join(wd, 'sequence.fa')
        self.msa = os.path.join(wd, 'msa.fa')
        self.tree = os.path.join(wd, 'tree.newick')
        self.asr = os.path.join(wd, 'asr.tsv')
        self.sim = os.path.join(wd, 'sim.tsv')
        self.rms = []

    def tearDown(self):
        if CLEANUP and self.rms:
            for rm in self.rms:
                try:
                    os.remove(os.path.join(PATH, 'tests', 'data', 'imc', rm))
                except OSError as err:
                    logging.warning('Failed to cleanup testing data {} due to:'
                                    '\n\t{}'.format(rm, err))

    @unittest.skipIf((MUSCLE is None) or (CODEML is None),
                     'No valid executable provided.')
    def test_imc_seq_muscle_codeml(self):
        outs = {'sequence.muscle.codeml.tsv',
                'sequence.muscle.codeml.counts.tsv',
                'sequence.muscle.codeml.details.tsv',
                'sequence.muscle.fa', 'sequence.muscle.trimmed.fa'}
        pars, cons, divs, details, _ = imc(self.seq, self.tree, aligner=MUSCLE,
                                           ancestor=CODEML, save=True)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(divs)
        self.assertIsNotNone(details)
        self.assertTrue(outs.issubset(set(os.listdir(wd))))
        self.rms = outs

    @unittest.skipIf((MUSCLE is None) or (CODEML is None) or (EVOLVER is None),
                     'No valid executable provided.')
    def test_imc_seq_muscle_codeml_evolver_iqtree(self):
        outs = {'sequence.evolver.tsv',
                'sequence.muscle.codeml.tsv',
                'sequence.muscle.codeml.counts.tsv',
                'sequence.muscle.codeml.details.tsv',
                'sequence.muscle.fa',
                'sequence.muscle.trimmed.fa'}
        pars, cons, divs, details, _ = imc(self.seq, tree=self.tree,
                                           aligner=MUSCLE, ancestor=CODEML,
                                           simulator=EVOLVER,
                                           save=True, verbose=True)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(divs)
        self.assertIsNotNone(details)
        self.assertTrue(outs.issubset(set(os.listdir(wd))))
        self.rms = outs

    @unittest.skipIf((CODEML is None) or (EVOLVER is None),
                     'No valid executable provided.')
    def test_imc_msa_codeml_evolver_iqtree(self):
        outs = {'msa.evolver.tsv',
                'msa.codeml.tsv',
                'msa.codeml.counts.tsv',
                'msa.codeml.details.tsv',
                'msa.trimmed.fa'}
        pars, cons, divs, details, _ = imc(self.msa, tree=self.tree,
                                           ancestor=CODEML,  simulator=EVOLVER,
                                           save=True, verbose=True)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(divs)
        self.assertIsNotNone(details)
        self.assertTrue(outs.issubset(set(os.listdir(wd))))
        self.rms = outs

    @unittest.skipIf(EVOLVER is None, 'No valid executable provided.')
    def test_imc_asr_evolver_iqtree(self):
        outs = {'asr.counts.tsv', 'asr.details.tsv'}
        pars, cons, divs, details, _ = imc(self.asr, save=True, verbose=True)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(divs)
        self.assertIsNotNone(details)
        self.assertTrue(outs.issubset(set(os.listdir(wd))))
        self.rms = outs

    def test_imc_threshold(self):
        outs = {'sequence.muscle.codeml.tsv',
                'sequence.muscle.codeml.counts.tsv',
                'sequence.muscle.codeml.details.tsv', 'sequence.muscle.fa',
                'sequence.muscle.trimmed.fa'}
        pars, cons, divs, details, _ = imc(self.seq, self.tree, aligner=MUSCLE,
                                           ancestor=CODEML,  threshold=0.5,
                                           save=True, verbose=True)
        self.assertIsNotNone(pars)
        self.assertIsNotNone(cons)
        self.assertIsNotNone(divs)
        self.assertIsNotNone(details)
        self.assertTrue(outs.issubset(set(os.listdir(wd))))
        self.rms = outs

  
if __name__ == '__main__':
    unittest.main()
