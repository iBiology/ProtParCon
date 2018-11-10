#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import logging
import unittest

PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, PATH)
from imc.aut import aut

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import settings

CLEANUP = settings.CLEANUP

IQTREE = settings.IQTREE

class TestAUT(unittest.TestCase):
    def setUp(self):
        self.msa = os.path.join(PATH, 'tests', 'data', 'aut', 'msa.fa')
        self.tree = os.path.join(PATH, 'tests', 'data', 'aut', 'tree.newick')
        self.trees = os.path.join(PATH, 'tests', 'data', 'aut', 'trees.newick')
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

    @unittest.skipIf(IQTREE is None, 'No IQ-TREE executable was provided.')
    def test_aut_default(self):
        aup, _, _ = aut(IQTREE, self.msa, self.tree)
        print('AU TEST P-value: {} (test default)'.format(aup))
        self.assertLessEqual(0.0, aup)

    @unittest.skipIf(IQTREE is None, 'No IQ-TREE executable was provided.')
    def test_aut_model(self):
        aup, _, _ = aut(IQTREE, self.msa, self.tree, model='JTT')
        print('AU TEST P-value: {} (test JTT model)'.format(aup))
        self.assertLessEqual(0.0, aup)

    @unittest.skipIf(IQTREE is None, 'No IQ-TREE executable was provided.')
    def test_aut_trees(self):
        aup, _, _ = aut(IQTREE, self.msa, self.trees, model='JTT')
        print('AU TEST P-value: {} (test multiple trees)'.format(aup))
        self.assertLessEqual(0.0, aup)

    @unittest.skipIf(IQTREE is None, 'No IQ-TREE executable was provided.')
    def test_aut_outfile(self):
        out = os.path.join(PATH, 'tests', 'data', 'aut', 'iMC.autest.txt')
        aup, _, _ = aut(IQTREE, self.msa, self.tree, outfile=out)
        print('AU TEST P-value: {} (test outfile)'.format(aup))
        self.assertLessEqual(0.0, aup)
        self.assertTrue(os.path.isfile(out))
        self.rm = out

  
if __name__ == '__main__':
    unittest.main()
